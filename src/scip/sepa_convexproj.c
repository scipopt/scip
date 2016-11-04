/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_convexproj.c
 * @brief  convexproj separator
 * @author Felipe Serrano
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_convexproj.h"
#include "scip/nlp.h"
#include "scip/prop_nlobbt.h"
#include "nlpi/exprinterpret.h"
#include "nlpi/nlpi.h"


#define SEPA_NAME              "convexproj"
#define SEPA_DESC              "separate at projection of point onto convex region"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE      /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                 TRUE      /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXDEPTH             -1      /* maximum depth at which the separator is applied; -1 means no limit */
#define DEFAULT_NLPTIMELIMIT        0.0      /**< default time limit of NLP solver; 0.0 for no limit */
#define DEFAULT_NLPITERLIM          250      /**< default NLP iteration limit */

#define MIN_VIOLATION              1e-4      /* minimum violation of point to be separated */

#define NLPFEASTOL                 1e-6      /**< NLP feasibility tolerance */
#define NLPVERBOSITY                  0      /**< NLP solver verbosity */

/*
 * Data structures
 */

/** side that makes an nlrow convex */
enum ConvexSide
{
   LHS = 0,                                  /**< left hand side */
   RHS = 1,                                  /**< right hand side */
   BOTH = 2                                  /**< both sides, ie, nlrow is linear */
};
typedef enum ConvexSide CONVEXSIDE;

/** separator data */
struct SCIP_SepaData
{
   SCIP_NLPI*            nlpi;               /**< nlpi used to create the nlpi problem */
   SCIP_NLPIPROBLEM*     nlpiprob;           /**< nlpi problem representing the convex NLP relaxation */
   SCIP_VAR**            nlpivars;           /**< array containing all variables of the nlpi */
   SCIP_HASHMAP*         var2nlpiidx;        /**< mapping between variables and nlpi indices */
   int                   nlpinvars;          /**< total number of nlpi variables */

   SCIP_Bool             cutoff;             /**< whether the separator detected a cutoff */
   SCIP_Bool             skipsepa;           /**< should separator be skipped? */

   SCIP_NLROW**          nlrows;             /**< convex nlrows */
   CONVEXSIDE*           convexsides;        /**< which sides make the nlrows convex */
   SCIP_Real*            constraintviolation;/**< array storing the violation of constraint by current solution; 0.0 if it is not violated */

   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */

   /* parameter */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0.0 for no limit */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for no limit */
   int                   maxdepth;           /**< maximal depth at which the separator is applied */

   /* statistic variables */
   int                   nconvexnlrows;      /**< total number of convex nonlinear nlrows */
   int                   nlinearnlrows;      /**< total number of linear nlrows */
   int                   nnlrows;            /**< total number of nlrows */
   int                   ncuts;              /**< number of cuts generated */
   int                   ncutsadded;         /**< number of cuts added in current call */
};


/*
 * Local methods
 */

/** clears the sepadata data */
static
SCIP_RETCODE sepadataClear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   assert(sepadata != NULL);

   if( sepadata->nlpiprob != NULL )
   {
      assert(sepadata->nlpi != NULL);

      SCIPfreeBlockMemoryArray(scip, &sepadata->nlpivars, sepadata->nlpinvars);
      SCIPfreeBlockMemoryArray(scip, &sepadata->nlrows, sepadata->nnlrows);
      SCIPfreeBlockMemoryArray(scip, &sepadata->convexsides, sepadata->nnlrows);
      SCIPfreeBlockMemoryArray(scip, &sepadata->constraintviolation, sepadata->nnlrows);
      SCIPhashmapFree(&sepadata->var2nlpiidx);
      SCIPnlpiFreeProblem(sepadata->nlpi, &sepadata->nlpiprob);
      SCIP_CALL( SCIPexprintFree(&sepadata->exprinterpreter) );

      sepadata->nlpinvars = 0;
   }
   assert(sepadata->nlpinvars == 0);

   sepadata->skipsepa = FALSE;

   return SCIP_OKAY;
}

/** computes gradient of exprtree at projection */
static
SCIP_RETCODE computeGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_SOL*             projection,         /**< point where we compute gradient */
   SCIP_EXPRTREE*        exprtree,           /**< exprtree to which we compute the gradient */
   SCIP_Real*            grad                /**< buffer to store the gradient; order ?? */
   )
{
   SCIP_Real* x;
   SCIP_Real val;
   int nvars;
   int i;

   /*
   SCIPdebugMsg(scip, "Computing gradient of: ");
   SCIPexprtreePrint(exprtree, SCIPgetMessagehdlr(scip), NULL, NULL ,NULL);
   SCIPdebugMsg(scip, "\n");
   */

   nvars = SCIPexprtreeGetNVars(exprtree);
   SCIP_CALL( SCIPallocBufferArray(scip, &x, nvars) );

   /* compile expression exprtree, if not done before */
   if( SCIPexprtreeGetInterpreterData(exprtree) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );
   }

   for( i = 0; i < nvars; ++i )
   {
      x[i] = SCIPgetSolVal(scip, projection, SCIPexprtreeGetVars(exprtree)[i]);
   }

   SCIP_CALL( SCIPexprintGrad(exprint, exprtree, x, TRUE, &val, grad) );

   /*SCIPdebug( for( i = 0; i < nvars; ++i ) printf("%e [%s]\n", grad[i], SCIPvarGetName(SCIPexprtreeGetVars(exprtree)[i])) );*/


   SCIPfreeBufferArray(scip, &x);

   return SCIP_OKAY;
}

/** computes gradient cut (linearization) of nlrow at projection */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_SOL*             projection,         /**< point where we compute gradient cut */
   SCIP_NLROW*           nlrow,              /**< constraint for which we generate gradient cut */
   SCIP_Real             activity,           /**< activity of constraint at projection */
   SCIP_ROW**            row                 /**< storage for cut */
   )
{
   char rowname[SCIP_MAXSTRLEN];
   SCIP_SEPADATA* sepadata;
   SCIP_Real gradx0; /* <grad f(x_0), x_0> */
   int i;

   assert(scip != NULL);
   assert(sepa != NULL);
   assert(exprint != NULL);
   assert(nlrow != NULL);
   assert(row != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   gradx0 = 0.0;
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "proj_cut_%s_%u", SCIPnlrowGetName(nlrow), ++(sepadata->ncuts));

   /* gradient cuts are globally valid whenever the constraint from which they are deduced is globally valid
    * since we build the convex relaxation using only globally valid constraints, the cuts are globally valid
    */
   SCIP_CALL( SCIPcreateEmptyRowSepa(scip, row, sepa, rowname, SCIPnlrowGetLhs(nlrow), SCIPnlrowGetRhs(nlrow),
           TRUE, FALSE , TRUE) );

   /* an nlrow has a linear part, quadratic part and expression tree; ideally I would just build the gradient but I do not
    * know if some variable of the linear part belong to the quadratic part or worse, to the expression tree, so I can't
    * just build the gradient; the solution for this is to create the row right away and compute the gradients of each
    * part and add them to the row; this way, the row will take care to add coefficients corresponding to the same variable
    * that belong to different parts of the nlrow
    */

   /* linear part */
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, SCIPnlrowGetNLinearVars(nlrow), SCIPnlrowGetLinearVars(nlrow), SCIPnlrowGetLinearCoefs(nlrow)) );
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
      gradx0 += SCIPgetSolVal(scip, projection, SCIPnlrowGetLinearVars(nlrow)[i]) * SCIPnlrowGetLinearCoefs(nlrow)[i];

   /* quadratic part; it seems (and probably is the only thing that makes sense) that the variables belonging to a quadelem
    * are quadvars[quadelem.idx1] and quadvars[quadelem.idx2]
    */
   SCIP_CALL( SCIPcacheRowExtensions(scip, *row) );
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real grad1;
      SCIP_Real grad2;

      var1  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];
      grad1 = SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, projection, var2);
      grad2 = SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, projection, var1);

      SCIP_CALL( SCIPaddVarToRow(scip, *row, var1, grad1) );
      SCIP_CALL( SCIPaddVarToRow(scip, *row, var2, grad2) );

      gradx0 += grad1 * SCIPgetSolVal(scip, projection, var1) + grad2 * SCIPgetSolVal(scip, projection, var2);
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, *row) );

   /* expression tree part */
   {
      SCIP_Real* grad;
      SCIP_EXPRTREE* tree;

      tree = SCIPnlrowGetExprtree(nlrow);

      if( tree != NULL )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &grad, SCIPexprtreeGetNVars(tree)) );

         SCIP_CALL( computeGradient(scip, sepadata->exprinterpreter, projection, tree, grad) );
         SCIP_CALL( SCIPaddVarsToRow(scip, *row, SCIPexprtreeGetNVars(tree), SCIPexprtreeGetVars(tree), grad) );

         for( i = 0; i < SCIPexprtreeGetNVars(tree); i++ )
            gradx0 +=  grad[i] * SCIPgetSolVal(scip, projection, SCIPexprtreeGetVars(tree)[i]);

         SCIPfreeBufferArray(scip, &grad);
      }
   }

   SCIPdebugPrintf("gradient: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );
   SCIPdebugPrintf("gradient dot x_0: %g\n", gradx0);

   /* gradient cut is lhs <= f(x_0) - <grad f(x_0), x_0> + <grad f(x_0), x> <= rhs */
   if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)) )
   {
      assert(SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));
      SCIP_CALL( SCIPchgRowRhs(scip, *row, SCIPnlrowGetRhs(nlrow) - activity + gradx0) );
   }
   else
   {
      assert(!SCIPisInfinity(scip, -SCIPnlrowGetLhs(nlrow)));
      SCIP_CALL( SCIPchgRowLhs(scip, *row, SCIPnlrowGetLhs(nlrow) - activity + gradx0) );
   }

   SCIPdebugPrintf("gradient cut: ");
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, *row, NULL) ) );

   return SCIP_OKAY;
}

/** set quadratic part of objective function: \sum_i x_i^2; the objective function is ||x - x_0||^2, where x_0 is the
 * point to separate; the only part that changes is the term -2 x_0 \cdot x which is linear. The linear part is set
 * every time we want to separate a point, see separateCuts
 */
static
SCIP_RETCODE setQuadraticObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< the cut separator data */
   )
{
   SCIP_QUADELEM* quadelems;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, sepadata->nlpinvars) );
   for( i = 0; i < sepadata->nlpinvars; i++ )
   {
      SCIP_VAR* var;

      var = sepadata->nlpivars[i];
      assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

      quadelems[i].idx1 = (int)(size_t)SCIPhashmapGetImage(sepadata->var2nlpiidx, (void*)var);
      quadelems[i].idx2 = quadelems[i].idx1;
      quadelems[i].coef = 1.0;
   }

   /* set quadratic part of objective function */
   SCIP_CALL( SCIPnlpiSetObjective(sepadata->nlpi, sepadata->nlpiprob,
            0, NULL, NULL, sepadata->nlpinvars, quadelems, NULL, NULL, 0.0) );

   /* free memory */
   SCIPfreeBufferArray(scip, &quadelems);

   return SCIP_OKAY;
}

/** projects sol onto convex relaxation (stored in sepadata) and tries to generate gradient cuts at the projection
 * it generates cuts only for the constraints that were violated by the LP solution and are now active or still
 * violated (in case we don't solve to optimality).
 * @note: this method modifies the sepadata, for instance sepadata->ncutsadded stores how many cuts where added.
 * @todo: store a feasible solution if one is found to use as warmstart
 */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< the cut separator itself */
   SCIP_SOL*             sol                 /**< solution that should be separated */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_NLROW*    nlrow;
   SCIP_SOL*      projection;
   SCIP_Real*     linvals;
   SCIP_Real*     nlpisol;
   SCIP_Real      timelimit;
   int            nlpinvars;
   int            i;
   int            iterlimit;
   int*           lininds;
   SCIP_Bool      nlpunstable;

   nlpunstable = FALSE;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);
   assert(sepadata->constraintviolation != NULL);

   nlpinvars = sepadata->nlpinvars;
   /* set linear part of objective function: \norm(x - x^0)^2 = \norm(x)^2 - \sum 2 * x_i * x^0_i + const
    * we ignore the constant; x0 is `sol`
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &linvals, nlpinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lininds, nlpinvars) );
   for( i = 0; i < nlpinvars; i++ )
   {
      SCIP_VAR* var;

      var = sepadata->nlpivars[i];
      assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

      lininds[i] = (int)(size_t)SCIPhashmapGetImage(sepadata->var2nlpiidx, (void*)var);
      linvals[i] = - 2.0 * SCIPgetSolVal(scip, sol, var);

      /* if coefficient is too big, don't separate */
      if( SCIPisInfinity(scip, REALABS(linvals[i])) )
      {
         SCIPfreeBufferArray(scip, &linvals);
         SCIPfreeBufferArray(scip, &lininds);
         SCIPdebugMsg(scip, "Don't separate points too close to infinity\n");

         return SCIP_OKAY;
      }
   }

   /* set linear part of objective function */
   SCIP_CALL( SCIPnlpiChgLinearCoefs(sepadata->nlpi, sepadata->nlpiprob,
            -1, nlpinvars, lininds, linvals) );

   /* set parameters in nlpi; time and iterations limit, tolerance, verbosity; for time limit, get time limit of scip;
    * if scip doesn't have much time left, don't run separator. otherwise, timelimit is the minimum between whats left
    * for scip and the timelimit setting
    */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 1.0 )
      {
         SCIPdebugMsg(scip, "skip NLP solve; no time left\n");
         return SCIP_OKAY;
      }
   }
   if( sepadata->nlptimelimit > 0.0 )
      timelimit = MIN(sepadata->nlptimelimit, timelimit);
   SCIP_CALL( SCIPnlpiSetRealPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );

   iterlimit = sepadata->nlpiterlimit > 0 ? sepadata->nlpiterlimit : INT_MAX;
   SCIP_CALL( SCIPnlpiSetIntPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_ITLIM, iterlimit) );
   SCIP_CALL( SCIPnlpiSetRealPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_FEASTOL, NLPFEASTOL) );
   SCIP_CALL( SCIPnlpiSetIntPar(sepadata->nlpi, sepadata->nlpiprob, SCIP_NLPPAR_VERBLEVEL, NLPVERBOSITY) );

   /* compute the projection onto the convex NLP relaxation */
   SCIP_CALL( SCIPnlpiSolve(sepadata->nlpi, sepadata->nlpiprob) );
   SCIPdebugMsg(scip, "NLP solstat = %d\n", SCIPnlpiGetSolstat(sepadata->nlpi, sepadata->nlpiprob));

   /* if solution is feasible, add a cuts */
   switch( SCIPnlpiGetSolstat(sepadata->nlpi, sepadata->nlpiprob) )
   {
      case SCIP_NLPSOLSTAT_GLOBOPT:
      case SCIP_NLPSOLSTAT_LOCOPT:
         /* @todo: if solution is optimal, we might as well add the cut <x - P(x_0), x_0 - P(x_0)> <= 0
          * even though this cut is implied by all the gradient cuts of the rows active at the projection,
          * we do not add them all (only the gradient cuts of constraints that violated the LP solution */
      case SCIP_NLPSOLSTAT_FEASIBLE:

         /* get solution: build SCIP_SOL out of nlpi sol */
         SCIP_CALL( SCIPnlpiGetSolution(sepadata->nlpi, sepadata->nlpiprob, &nlpisol, NULL, NULL, NULL) );
         SCIP_CALL( SCIPcreateSol(scip, &projection, NULL) );
         for( i = 0; i < nlpinvars; i++ )
         {
            SCIP_VAR* var;

            var = sepadata->nlpivars[i];
            assert(SCIPhashmapExists(sepadata->var2nlpiidx, (void*)var) );

            SCIP_CALL( SCIPsetSolVal(scip, projection, var,
                     nlpisol[(int)(size_t)SCIPhashmapGetImage(sepadata->var2nlpiidx, (void *)var)]) );
         }
         SCIPprintSol(scip, projection, NULL, TRUE);

         /* generate cuts for constraints that violate the solution that we want to separate (sol)
          * NOTE: such a constraint could be not active at the projection; we only add cuts for
          * constraints that are also active at the projection. because of numerics, it could
          * happen that no constraint is active at the projection; therefore, we consider constraints
          * that are: violated at sol, active or violate at the projection (nlpsol);
          * if the projection is on the interior of the region, we do nothing
          */
         for( i = 0; i < sepadata->nnlrows; ++i )
         {
            SCIP_Real activity;

            /* ignore constraints that are not violated by `sol` */
            if( SCIPisZero(scip, sepadata->constraintviolation[i]) )
               continue;

            nlrow = sepadata->nlrows[i];
            assert(nlrow != NULL);

            /* ignore linear constraints, since there shouldn't be any violated linear constraint by `sol` */
            if( sepadata->convexsides[i] == BOTH )
               continue;

            /* check for currently active constraints at projected point */
            SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, projection, &activity) );

            SCIPdebugMsg(scip, "NlRow activity at nlpi solution: %g <= %g <= %g\n",
                  SCIPnlrowGetLhs(nlrow), activity, SCIPnlrowGetRhs(nlrow) );
            /* in case the projection is computed optimally, we should assert this
            assert(SCIPisFeasLE(scip, SCIPnlrowGetLhs(nlrow), activity));
            assert(SCIPisFeasLE(scip, activity, SCIPnlrowGetRhs(nlrow)));
            */

            /* if nlrow is active or violates the projection, build gradient cut at projection */
            if( SCIPisFeasGE(scip, activity, SCIPnlrowGetRhs(nlrow)) || SCIPisFeasLE(scip, activity, SCIPnlrowGetLhs(nlrow)) )
            {
               SCIP_ROW* row;

               SCIPdebugMsg(scip, "active nlrow: (%e) ", sepadata->constraintviolation[i]);
               SCIP_CALL( SCIPprintNlRow(scip, nlrow, NULL) );

               SCIP_CALL( generateCut(scip, sepa, sepadata->exprinterpreter, projection, nlrow, activity, &row) );
               SCIPprintRow(scip, row, NULL);

               /* add cut if it is efficacious for the point we want to separate (sol) */
               SCIPdebugMsg(scip, "cut with efficacy %g generated\n", SCIPgetCutEfficacy(scip, sol, row));
               if( SCIPisCutEfficacious(scip, sol, row) )
               {
                  SCIP_Bool infeasible;

                  SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE, &infeasible) );

                  if( infeasible )
                     sepadata->cutoff = TRUE;
                  else
                     sepadata->ncutsadded++;
               }

               /* release the row */
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }

         }

#ifdef SCIP_DEBUG
         {
            SCIP_Real distance;

            /* compute distance between LP sol and its projection (only makes sense when it is optimal) */
            distance = 0.0;
            for( i = 0; i < SCIPgetNNLPVars(scip); ++i )
            {
               SCIP_VAR* var;

               var = SCIPgetNLPVars(scip)[i];
               assert(var != NULL);

               /* assert NLP solution is within the bounds of the variable (only make sense when sol is optimal) */
               if( !SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
                  assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(var), SCIPvarGetNLPSol(var)));
               if( !SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
                  assert(SCIPisFeasLE(scip, SCIPvarGetNLPSol(var), SCIPvarGetUbLocal(var)));

               /*SCIPdebugMsg(scip, "NLP sol (LP sol): %s = %f (%g)\n", SCIPvarGetName(var),
                *     SCIPvarGetNLPSol(var), SCIPgetSolVal(scip, sol, var));
                */

               distance += SQR( SCIPvarGetNLPSol(var) - SCIPgetSolVal(scip, sol, var) );
            }

            SCIPdebugMsg(scip, "NLP objval: %e, distance: %e\n", SCIPgetNLPObjval(scip), distance);
         }
#endif

         /* free solution */
         SCIP_CALL( SCIPfreeSol(scip, &projection) );
         break;


      case SCIP_NLPSOLSTAT_LOCINFEASIBLE:
         /* one would like to do something smart here, e.g. cut off the node;
          * but most likely we are here because of numerics, so fallthrough
          */
      case SCIP_NLPSOLSTAT_UNKNOWN:
         /* unknown... assume numerical issues */
         nlpunstable = TRUE;
         break;

      case SCIP_NLPSOLSTAT_UNBOUNDED:
      default:
         SCIPerrorMessage("Projection NLP is not unbounded by construction, should not get here!\n");
         nlpunstable = TRUE;
         SCIPABORT();
   }


   /* if nlp is detected to be unstable, don't try to separate again */
   if( nlpunstable )
   {
      /* TODO: maybe change objective function to \sum [(x_i - x_i^*)/max(|x_i^*|, 1)]^2
       * or some other scaling when unstable and try again.
       *       maybe free it here */
      sepadata->skipsepa = TRUE;
   }

   /* reset objective
    * TODO: BMS clear buffer memory?
    */
   for( i = 0; i < nlpinvars; i++ )
      linvals[i] = 0.0;
   SCIP_CALL( SCIPnlpiChgLinearCoefs(sepadata->nlpi, sepadata->nlpiprob,
            -1, nlpinvars, lininds, linvals) );
   SCIPfreeBufferArray(scip, &linvals);
   SCIPfreeBufferArray(scip, &lininds);


   return SCIP_OKAY;
}

/** computes the violation and maximum violation of the convex nlrows stored in sepadata wrt sol */
static
SCIP_RETCODE computeMaxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< solution that should be separated */
   SCIP_Real*            maxviolation        /**< buffer to store maximum violation */
   )
{
   SCIP_NLROW*    nlrow;
   SCIP_Real      sum;
   int            i;

   assert(sepadata != NULL);
   assert(sepadata->constraintviolation != NULL);

   *maxviolation = 0.0;
   sum = 0.0;
   for( i = 0; i < sepadata->nnlrows; i++ )
   {
      SCIP_Real activity;

      /* skip linear constraints */
      if( sepadata->convexsides[i] == BOTH )
         continue;

      nlrow = sepadata->nlrows[i];

      /* get activity of nlrow */
      SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrow, sol, &activity) );

      /* violation = max{activity - rhs, 0.0} when convex and max{lhs - activity, 0.0} when concave */
      if( sepadata->convexsides[i] == RHS )
      {
         assert(SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX);
         sepadata->constraintviolation[i] = MAX(activity - SCIPnlrowGetRhs(nlrow), 0.0);
      }
      if( sepadata->convexsides[i] == LHS )
      {
         assert(SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE);
         sepadata->constraintviolation[i] = MAX(SCIPnlrowGetLhs(nlrow) - activity, 0.0);
      }

      sum += sepadata->constraintviolation[i];

      /* compute maximum */
      if( *maxviolation < sepadata->constraintviolation[i] )
         *maxviolation = sepadata->constraintviolation[i];
   }

   SCIPdebugMsg(scip, "Maximum violation %g\n", *maxviolation);
   SCIPdebugMsg(scip, "Average violation %g\n", sum/sepadata->nnlrows);

   return SCIP_OKAY;
}


/** stores, from the constraints represented by nlrows, the convex ones in sepadata
 * counts the number of linear and nonlinear constraints as well */
static
SCIP_RETCODE storeConvexNlrows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_NLROW**          nlrows,             /**< nlrows from which to store convex ones */
   int                   nnlrows             /**< number of nlrows */
   )
{
   int i;

   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIPdebugMsg(scip, "storing convex nlrows\n");

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->nlrows), nnlrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->convexsides), nnlrows) );

   /* count the number of convex, linear and non-convex nonlinear rows; store the convex ones */
   sepadata->nconvexnlrows = 0;
   sepadata->nlinearnlrows = 0;
   sepadata->nnlrows = 0;
   for( i = 0; i < nnlrows; ++i )
   {
      SCIP_NLROW* nlrow;

      nlrow = nlrows[i];
      assert(nlrow != NULL);

      /* linear case */
      /* TODO: this should be equivalent to the curvature being convex */
      if( SCIPnlrowGetNQuadElems(nlrow) == 0 && SCIPnlrowGetExprtree(nlrow) == NULL )
      {
         ++(sepadata->nlinearnlrows);
         sepadata->convexsides[sepadata->nnlrows] = BOTH;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
      /* nonlinear case */
      else if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrow)) && SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONVEX )
      {
         ++(sepadata->nconvexnlrows);
         sepadata->convexsides[sepadata->nnlrows] = RHS;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
      else if( !SCIPisInfinity(scip, SCIPnlrowGetLhs(nlrow)) && SCIPnlrowGetCurvature(nlrow) == SCIP_EXPRCURV_CONCAVE )
      {
         ++(sepadata->nconvexnlrows);
         sepadata->convexsides[sepadata->nnlrows] = LHS;
         sepadata->nlrows[sepadata->nnlrows] = nlrow;
         ++(sepadata->nnlrows);
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( sepadataClear(scip, sepadata) );

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}



/** solving process initialization method of separator (called when branch and bound process is about to begin) */
static
SCIP_DECL_SEPAINITSOL(sepaInitsolConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   sepadata->nlinearnlrows = 0;
   sepadata->nconvexnlrows = 0;

   return SCIP_OKAY;
}


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolConvexproj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);

   assert(sepadata != NULL);

   SCIP_CALL( sepadataClear(scip, sepadata) );

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpConvexproj)
{  /*lint --e{715}*/

   SCIP_Real maxviolation;
   SCIP_SOL* lpsol;
   SCIP_SEPADATA* sepadata;

   *result = SCIP_DIDNOTRUN;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* do not run if there is no interesting convex relaxation (with at least one nonlinear convex constraint),
    * or if we have found it to be numerically unstable
    * @todo: should it be with at least 2 nonlinear convex constraints?
    */
   if( sepadata->skipsepa )
   {
      SCIPdebugMsg(scip, "not running because convex relaxation is uninteresting or numerically unstable\n");
      return SCIP_OKAY;
   }

   /* only call separator up to a maximum depth */
   if( sepadata->maxdepth >= 0 && SCIPgetDepth(scip) > sepadata->maxdepth )
      return SCIP_OKAY;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* do not run if SCIP does not have constructed an NLP */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMsg(scip, "NLP not constructed, skipping convex projection separator\n");
      return SCIP_OKAY;
   }

   /* recompute convex NLP relaxation if the variable set changed and we are still at the root node */
   if( sepadata->nlpiprob != NULL && SCIPgetNVars(scip) != sepadata->nlpinvars  && SCIPgetDepth(scip) == 0 )
   {
      SCIP_CALL( sepadataClear(scip, sepadata) );
      assert(sepadata->nlpiprob == NULL);
   }

   /* create or update convex NLP relaxation */
   if( sepadata->nlpiprob == NULL )
   {
      /* store convex nonlinear constraints */
      SCIP_CALL( storeConvexNlrows(scip, sepadata, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip)) );

      /* check that convex NLP relaxation is interesting (more than one nonlinear constraint) */
      if( sepadata->nconvexnlrows < 1 )
      {
         SCIPdebugMsg(scip, "convex relaxation uninteresting, don't run\n");
         sepadata->skipsepa = TRUE;
         return SCIP_OKAY;
      }

      /* initialize some of the sepadata */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->constraintviolation), sepadata->nnlrows) );
      SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &sepadata->exprinterpreter) );

      sepadata->nlpinvars = SCIPgetNVars(scip);
      sepadata->nlpi = SCIPgetNlpis(scip)[0];
      assert(sepadata->nlpi != NULL);

      SCIP_CALL( SCIPnlpiCreateProblem(sepadata->nlpi, &sepadata->nlpiprob, "convexproj-nlp") );
      SCIP_CALL( SCIPhashmapCreate(&sepadata->var2nlpiidx, SCIPblkmem(scip),
            SCIPcalcHashtableSize(sepadata->nlpinvars)) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &sepadata->nlpivars, SCIPgetVars(scip), sepadata->nlpinvars) );

      /* I shouldn't care about the cutoff, just assert that the lp solution satisfies the cutoff bound */
      SCIP_CALL( SCIPcreateConvexNlpNlobbt(scip, sepadata->nlpi, SCIPgetNLPNlRows(scip), SCIPgetNNLPNlRows(scip),
            sepadata->nlpiprob, sepadata->var2nlpiidx, NULL, SCIPgetCutoffbound(scip)) );

      /* add rows of the LP */
      if( SCIPgetDepth(scip) == 0 )
      {
         SCIP_CALL( SCIPaddConvexNlpRowsNlobbt(scip, sepadata->nlpi, sepadata->nlpiprob, sepadata->var2nlpiidx, SCIPgetLPRows(scip),
               SCIPgetNLPRows(scip)) );
      }

      /* set quadratic part of objective function */
      SCIP_CALL( setQuadraticObj(scip, sepadata) );
   }
   else
   {
      SCIP_CALL( SCIPupdateConvexNlpNlobbt(scip, sepadata->nlpi, sepadata->nlpiprob, sepadata->var2nlpiidx,
            sepadata->nlpivars, sepadata->nlpinvars) );
   }

   /* get current sol: LP or pseudo solution if LP sol is not available */
   SCIP_CALL( SCIPcreateCurrentSol(scip, &lpsol, NULL) );

   /* do not run if current solution's violation is small */
   SCIP_CALL( computeMaxViolation(scip, sepadata, lpsol, &maxviolation) );
   if( maxviolation < MIN_VIOLATION )
   {
      SCIPdebugMsg(scip, "solution doesn't violate constraints enough, do not separate\n");
      SCIP_CALL( SCIPfreeSol(scip, &lpsol) );
      return SCIP_OKAY;
   }

   /* run the separator */
   *result = SCIP_DIDNOTFIND;
   sepadata->cutoff = FALSE;
   sepadata->ncutsadded = 0;

   /* separateCuts computes the projection and then gradient cuts on each constraint that was originally violated;
    * when adding the cuts determines whether there was a cutoff, but the result is set here
    */
   SCIP_CALL( separateCuts(scip, sepa, lpsol) );

   if( sepadata->ncutsadded > 0 )
      *result = SCIP_SEPARATED;
   if( sepadata->cutoff )
      *result = SCIP_CUTOFF;

   /* free memory */
   SCIP_CALL( SCIPfreeSol(scip, &lpsol) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the convexproj separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaConvexproj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create convexproj separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpConvexproj, NULL,
         sepadata) );
   sepadata->ncuts = 0;
   sepadata->ncutsadded = 0;
   sepadata->skipsepa = FALSE;
   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeConvexproj) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolConvexproj) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolConvexproj) );

   /* add convexproj separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxdepth",
         "maximal depth at which the separator is applied (-1: unlimited)",
         &sepadata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/"SEPA_NAME"/nlpiterlimit",
         "iteration limit of NLP solver; 0 for no limit",
         &sepadata->nlpiterlimit, TRUE, DEFAULT_NLPITERLIM, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/"SEPA_NAME"/nlptimelimit",
         "time limit of NLP solver; 0.0 for no limit",
         &sepadata->nlptimelimit, TRUE, DEFAULT_NLPTIMELIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
