#define SCIP_DEBUG
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

/**@file   heur_multistart.c
 * @brief  multistart heuristic for convex and nonconvex MINLPs
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_multistart.h"
#include "scip/heur_subnlp.h"

#include "scip/nlp.h"
#include "nlpi/exprinterpret.h"


#define HEUR_NAME             "multistart"
#define HEUR_DESC             "multistart heuristic for convex and nonconvex MINLPs"
#define HEUR_DISPCHAR         'm'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_RANDOMSEED    0              /**< default random seed */
#define DEFAULT_NRNDPOINTS  100              /**< default number of generated random points per call */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */
   int                   nrndpoints;         /**< number of random points generated per execution call */
   unsigned int          randseed;           /**< seed value for random number generator */
};


/*
 * Local methods
 */

static
SCIP_RETCODE computeGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_SOL*             sol,                /**< solution to compute the gradient for */
   SCIP_Real*            lingrad,            /**< array to store the gradient belonging to the linear part */
   SCIP_Real*            quadgrad,           /**< array to store the gradient belonging to the quadratic part */
   SCIP_Real*            treegrad,           /**< array to store the gradient belonging to the tree */
   SCIP_Real*            norm                /**< buffer to store the square of the norm of the gradient (assume that variables are different in the linear, quadratic, and tree part) */
   )
{
   SCIP_EXPRTREE* tree;
   int i;

   assert(scip != NULL);
   assert(nlrow != NULL);
   assert(exprint != NULL);
   assert(sol != NULL);
   assert(lingrad != NULL);
   assert(quadgrad != NULL);
   assert(treegrad != NULL);
   assert(norm != NULL);

   *norm = 0.0;

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      lingrad[i] = SCIPnlrowGetLinearCoefs(nlrow)[i];
      *norm += SQR(SCIPnlrowGetLinearCoefs(nlrow)[i]);
   }

   /* quadratic part */
   BMSclearMemoryArray(quadgrad, SCIPnlrowGetNQuadVars(nlrow));
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      var1  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];
      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx1 < SCIPnlrowGetNQuadVars(nlrow));
      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx2 < SCIPnlrowGetNQuadVars(nlrow));

      quadgrad[SCIPnlrowGetQuadElems(nlrow)[i].idx1] += SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, sol, var2);
      quadgrad[SCIPnlrowGetQuadElems(nlrow)[i].idx2] += SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, sol, var1);
   }

   for( i = 0; i < SCIPnlrowGetNQuadVars(nlrow); ++i )
      *norm += SQR(quadgrad[i]);

   /* tree part */
   tree = SCIPnlrowGetExprtree(nlrow);
   if( tree != NULL )
   {
      SCIP_Real* x;
      SCIP_Real val;

      SCIP_CALL( SCIPallocBufferArray(scip, &x, SCIPexprtreeGetNVars(tree)) );
      assert(SCIPexprtreeGetNVars(tree) <= SCIPgetNVars(scip));

      /* compile expression tree, if not done before */
      if( SCIPexprtreeGetInterpreterData(tree) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(exprint, tree) );
      }

      /* sets the solution value */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
         x[i] = SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(tree)[i]);

      SCIP_CALL( SCIPexprintGrad(exprint, tree, x, TRUE, &val, treegrad) );
      SCIPfreeBufferArray(scip, &x);

      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
         *norm += SQR(treegrad[i]);
   }

   return SCIP_OKAY;
}

/** computes the maximum violation of a given point; a negative value means that there is a violation */
static
SCIP_RETCODE getMaxViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_Real*            maxviol             /**< buffer to store the maximum violation */
   )
{
   SCIP_NLROW** nlrows;
   SCIP_Real tmp;
   int i;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(maxviol != NULL);

   *maxviol = SCIPinfinity(scip);

   nlrows = SCIPgetNLPNlRows(scip);
   assert(nlrows != NULL);

   for( i = 0; i < SCIPgetNNLPNlRows(scip); ++i )
   {
      assert(nlrows[i] != NULL);

      SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrows[i], sol, &tmp) );
      *maxviol = MIN(*maxviol, tmp);
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE addSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             val                 /**< value to add */
   )
{
   SCIP_Real newval;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(var != NULL);

   newval = SCIPgetSolVal(scip, sol, var) + val;

   /* adjust value */
   newval = MAX(newval, SCIPvarGetLbLocal(var));
   newval = MIN(newval, SCIPvarGetUbLocal(var));

   /* set new solution value */
   SCIPdebugMessage("  change value of %s from %e -> %e\n", SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), newval);
   SCIP_CALL( SCIPsetSolVal(scip, sol, var, newval) );

   return SCIP_OKAY;
}

/** moves the point in the direction of the feasible set (for the continuous relaxation) */
static
SCIP_RETCODE movePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_SOL*             startsol,           /**< starting point */
   SCIP_Real*            maxviol            /**< buffer to store the resulting maximum violation */
   )
{
   SCIP_NLROW** nlrows;
   SCIP_Real* lingrad;
   SCIP_Real* quadgrad;
   SCIP_Real* treegrad;
   int nnlrows;
   int nrounds;

   assert(scip != NULL);
   assert(heurdata != NULL);
   assert(startsol != NULL);
   assert(maxviol != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &lingrad, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &quadgrad, SCIPgetNVars(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &treegrad, SCIPgetNVars(scip)) );
   nrounds = 0;

   nnlrows = SCIPgetNNLPNlRows(scip);
   nlrows = SCIPgetNLPNlRows(scip);
   assert(nlrows != NULL);

   if( heurdata->exprinterpreter == NULL )
   {
      SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &heurdata->exprinterpreter) );
   }

   /* compute start violation */
   SCIP_CALL( getMaxViolation(scip, startsol, maxviol) );

   SCIPdebugMessage("start maxviol = %e\n", *maxviol);

   /* main loop */
   while( !SCIPisFeasGE(scip, *maxviol, 0.0) && nrounds < 100 )
   {
      int i;

      for( i = 0; i < nnlrows; ++i )
      {
         SCIP_NLROW* nlrow;
         SCIP_EXPRTREE* tree;
         SCIP_VAR* var;
         SCIP_Real feasibility;
         SCIP_Real scale;
         SCIP_Real norm;
         int j;

         nlrow = nlrows[i];
         assert(nlrow != NULL);
         SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrow, startsol, &feasibility) );

         /* do not consider non-violated constraints */
         if( SCIPisFeasGE(scip, feasibility, 0.0) )
            continue;

         SCIP_CALL( computeGradient(scip, nlrow, heurdata->exprinterpreter, startsol, lingrad, quadgrad, treegrad, &norm) );
         SCIPdebugMessage("  nlrow %d ||gradient|| = %e\n", i, norm);

         /* skip nonlinear rows with a too small gradient */
         if( SCIPisEQ(scip, norm, 0.0) )
            continue;

         scale = feasibility / norm;
         SCIPdebugMessage("  nlrow %d scale = %e\n", i, scale);

         /* skip nonliner row of the scaler is too small or too large */
         if( SCIPisEQ(scip, scale, 0.0) || SCIPisInfinity(scip, REALABS(scale)) )
            continue;

         /* update point */
         for( j = 0; j < SCIPnlrowGetNLinearVars(nlrow); ++j )
         {
            var = SCIPnlrowGetLinearVars(nlrow)[j];
            SCIP_CALL( addSolVal(scip, startsol, var, scale * lingrad[j]) );
         }

         for( j = 0; j < SCIPnlrowGetNQuadVars(nlrow); ++j )
         {
            var = SCIPnlrowGetQuadVars(nlrow)[j];
            SCIP_CALL( addSolVal(scip, startsol, var, scale * quadgrad[j]) );
         }

         tree = SCIPnlrowGetExprtree(nlrow);
         for( j = 0; tree != NULL && j < SCIPexprtreeGetNVars(tree); ++j )
         {
            var = SCIPexprtreeGetVars(tree)[j];
            SCIP_CALL( addSolVal(scip, startsol, var, scale * treegrad[j]) );
         }
      }

      /* update maximum violation */
      SCIP_CALL( getMaxViolation(scip, startsol, maxviol) );
      SCIPdebugMessage("violation after round %d = %e\n", nrounds, *maxviol);

      ++nrounds;
   }

   SCIPfreeBufferArray(scip, &treegrad);
   SCIPfreeBufferArray(scip, &quadgrad);
   SCIPfreeBufferArray(scip, &lingrad);

   return SCIP_OKAY;
}

static
SCIP_RETCODE applyHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* startsol;
   SCIP_Real maxviol;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPcreateSol(scip, &startsol, heur) );

   SCIP_CALL( movePoint(scip, heurdata, startsol, &maxviol) );

   SCIP_CALL( SCIPfreeSol(scip, &startsol) );

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyMultistart)
{  /*lint --e{715}*/
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurMultistart(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeMultistart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);

   if( heurdata->exprinterpreter != NULL )
   {
      SCIP_CALL( SCIPexprintFree(&heurdata->exprinterpreter) );
   }

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitMultistart)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitMultistart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of multistart primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitMultistart NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolMultistart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of multistart primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolMultistart NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolMultistart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of multistart primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolMultistart NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecMultistart)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   /* check cases for which the heuristic is not applicable */
   if( !SCIPisNLPConstructed(scip) || SCIPfindHeur(scip, "subnlp") == NULL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( applyHeur(scip, heur) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the multistart primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurMultistart(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create multistart primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* initialize random seed */
   heurdata->randseed = DEFAULT_RANDOMSEED;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecMultistart, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyMultistart) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeMultistart) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitMultistart) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitMultistart) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolMultistart) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolMultistart) );

   /* add multistart primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/nrndpoints",
         "number of random points generated per execution call",
         &heurdata->nrndpoints, FALSE, DEFAULT_NRNDPOINTS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
