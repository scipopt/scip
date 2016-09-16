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

#include "nlpi/exprinterpret.h"


#define HEUR_NAME             "multistart"
#define HEUR_DESC             "multistart heuristic for convex and nonconvex MINLPs"
#define HEUR_DISPCHAR         'm'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
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


/** returns an unique index of a variable in the range of 0,..,SCIPgetNVars(scip)-1 */
#ifndef NDEBUG
static
int getVarIndex(
   SCIP_HASHMAP*         varindex,           /**< maps variables to indicies between 0,..,SCIPgetNVars(scip)-1 */
   SCIP_VAR*             var                 /**< variable */
   )
{
   assert(varindex != NULL);
   assert(var != NULL);
   assert(SCIPhashmapExists(varindex, (void*)var));

   return (int)(size_t)SCIPhashmapGetImage(varindex, (void*)var);
}
#else
#define getVarIndex(varindex,var) ((int)(size_t)SCIPhashmapGetImage((varindex), (void*)(var)))
#endif

/** samples and stores random points */
static
SCIP_RETCODE sampleRandomPoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            rndpoints,          /**< array to store all random points */
   int                   nrndpoints,         /**< total number of random points to compute */
   unsigned int*         rndseed,            /**< random seed */
   SCIP_Real             intervalsize        /**< interval size for unbounded variables */
   )
{
   SCIP_VAR** vars;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   int nvars;
   int i;
   int k;

   assert(scip != NULL);
   assert(rndpoints != NULL);
   assert(nrndpoints > 0);
   assert(rndseed != NULL);
   assert(intervalsize > 0.0);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   for( k = 0; k < nrndpoints; ++k )
   {
      SCIP_CALL( SCIPcreateSol(scip, &rndpoints[k], NULL) );

      for( i = 0; i < nvars; ++i )
      {
         lb = MIN(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]));
         ub = MAX(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]));

         /* use a smaller domain for unbounded variables */
         if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub) )
            val = SCIPgetRandomReal(lb, ub, rndseed);
         else if( !SCIPisInfinity(scip, -lb) )
            val = SCIPgetRandomReal(lb, lb + intervalsize, rndseed);
         else if( !SCIPisInfinity(scip, ub) )
            val = SCIPgetRandomReal(ub - intervalsize, ub, rndseed);
         else
         {
            assert(SCIPisInfinity(scip, -lb) && SCIPisInfinity(scip, ub));
            val = SCIPgetRandomReal( -0.5*intervalsize, 0.5*intervalsize, rndseed);
         }
         assert(val >= lb && val <= ub);

         /* set solution value */
         SCIP_CALL( SCIPsetSolVal(scip, rndpoints[k], vars[i], val) );
      }

      assert(rndpoints[k] != NULL);
   }

   return SCIP_OKAY;
}

/** computes the maximum violation of a given point; a negative value means that there is a violation */
static
SCIP_RETCODE getMaxViol(
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

/** computes the gradient for a given point and nonlinear row */
static
SCIP_RETCODE computeGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< nonlinear row */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_SOL*             sol,                /**< solution to compute the gradient for */
   SCIP_HASHMAP*         varindex,           /**< maps variables to indicies between 0,..,SCIPgetNVars(scip)-1 uniquely */
   SCIP_Real*            grad,               /**< buffer to store the gradient; grad[varindex(i)] corresponds to SCIPgetVars(scip)[i] */
   SCIP_Real*            norm                /**< buffer to store ||grad||^2  */
   )
{
   SCIP_EXPRTREE* tree;
   SCIP_VAR* var;
   int i;

   assert(scip != NULL);
   assert(nlrow != NULL);
   assert(varindex != NULL);
   assert(exprint != NULL);
   assert(sol != NULL);
   assert(norm != NULL);

   BMSclearMemoryArray(grad, SCIPgetNVars(scip));
   *norm = 0.0;

   /* linear part */
   for( i = 0; i < SCIPnlrowGetNLinearVars(nlrow); i++ )
   {
      var = SCIPnlrowGetLinearVars(nlrow)[i];
      assert(var != NULL);
      assert(getVarIndex(varindex, var) >= 0 && getVarIndex(varindex, var) < SCIPgetNVars(scip));

      grad[getVarIndex(varindex, var)] += SCIPnlrowGetLinearCoefs(nlrow)[i];
   }

   /* quadratic part */
   for( i = 0; i < SCIPnlrowGetNQuadElems(nlrow); i++ )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;

      var1  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx1];
      var2  = SCIPnlrowGetQuadVars(nlrow)[SCIPnlrowGetQuadElems(nlrow)[i].idx2];

      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx1 < SCIPnlrowGetNQuadVars(nlrow));
      assert(SCIPnlrowGetQuadElems(nlrow)[i].idx2 < SCIPnlrowGetNQuadVars(nlrow));
      assert(getVarIndex(varindex, var1) >= 0 && getVarIndex(varindex, var1) < SCIPgetNVars(scip));
      assert(getVarIndex(varindex, var2) >= 0 && getVarIndex(varindex, var2) < SCIPgetNVars(scip));

      grad[getVarIndex(varindex, var1)] += SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, sol, var2);
      grad[getVarIndex(varindex, var2)] += SCIPnlrowGetQuadElems(nlrow)[i].coef * SCIPgetSolVal(scip, sol, var1);
   }

   /* tree part */
   tree = SCIPnlrowGetExprtree(nlrow);
   if( tree != NULL )
   {
      SCIP_Real* treegrad;
      SCIP_Real* x;
      SCIP_Real val;

      assert(SCIPexprtreeGetNVars(tree) <= SCIPgetNVars(scip));

      SCIP_CALL( SCIPallocBufferArray(scip, &x, SCIPexprtreeGetNVars(tree)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &treegrad, SCIPexprtreeGetNVars(tree)) );

      /* compile expression tree, if not done before */
      if( SCIPexprtreeGetInterpreterData(tree) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(exprint, tree) );
      }

      /* sets the solution value */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
         x[i] = SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(tree)[i]);

      SCIP_CALL( SCIPexprintGrad(exprint, tree, x, TRUE, &val, treegrad) );

      /* update corresponding gradient entry */
      for( i = 0; i < SCIPexprtreeGetNVars(tree); ++i )
      {
         var = SCIPexprtreeGetVars(tree)[i];
         assert(var != NULL);
         assert(getVarIndex(varindex, var) >= 0 && getVarIndex(varindex, var) < SCIPgetNVars(scip));

         grad[getVarIndex(varindex, var)] += treegrad[i];
      }

      SCIPfreeBufferArray(scip, &treegrad);
      SCIPfreeBufferArray(scip, &x);
   }

   /* compute ||grad||^2 */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
      *norm += SQR(grad[i]);

   return SCIP_OKAY;
}

/** use consensus vectors to improve feasibility for a given starting point */
static
SCIP_RETCODE improvePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HASHMAP*         varindex,           /**< maps variables to indicies between 0,..,SCIPgetNVars(scip)-1 */
   SCIP_EXPRINT*         exprinterpreter,    /**< expression interpreter */
   SCIP_SOL*             point,                /**< random generated point */
   int                   nmaxiter,           /**< maximum number of iterations */
   SCIP_Real             minimprfac,         /**< minimum required improving factor to proceed */
   SCIP_Real*            maxviol             /**< pointer to store the maximum violation */
   )
{
   SCIP_NLROW** nlrows;
   SCIP_Real* grad;
   SCIP_Real lastmaxviol;
   int nnlrows;
   int r;
   int i;

   assert(varindex != NULL);
   assert(exprinterpreter != NULL);
   assert(point != NULL);
   assert(nmaxiter > 0);
   assert(minimprfac >= 0.0);
   assert(maxviol != NULL);

   nnlrows = SCIPgetNNLPNlRows(scip);
   nlrows = SCIPgetNLPNlRows(scip);
   assert(nlrows != NULL);

   SCIP_CALL( getMaxViol(scip, point, maxviol) );
   SCIPdebugMessage("maxviol = %e\n", *maxviol);

   /* stop since start point is feasible */
   if( !SCIPisFeasLT(scip, *maxviol, 0.0) )
   {
      SCIPdebugMessage("start point is feasible");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &grad, SCIPgetNVars(scip)) );

   /* main loop */
   do
   {
      SCIP_Real feasibility;
      SCIP_Real activity;
      SCIP_Real nlrownorm;
      SCIP_Real scale;

      lastmaxviol = *maxviol;

      /* move point; compute acitivity and gradient first */
      for( i = 0; i < nnlrows; ++i )
      {
         int j;

         SCIP_CALL( SCIPgetNlRowSolFeasibility(scip, nlrows[i], point, &feasibility) );

         /* do not consider non-violated constraints */
         if( SCIPisFeasGE(scip, feasibility, 0.0) )
            continue;

         SCIP_CALL( SCIPgetNlRowSolActivity(scip, nlrows[i], point, &activity) );
         SCIP_CALL( computeGradient(scip, nlrows[i], exprinterpreter, point, varindex, grad, &nlrownorm) );

         /* compute -g(x_k) / ||grad(g)(x_k)||^2 for a constraint g(x_k) <= 0 */
         scale = -feasibility / nlrownorm;
         if( !SCIPisInfinity(scip, SCIPnlrowGetRhs(nlrows[i])) && SCIPisGT(scip, activity, SCIPnlrowGetRhs(nlrows[i])) )
            scale *= -1.0;

         /* skip nonliner row of the scaler is too small or too large */
         if( SCIPisEQ(scip, scale, 0.0) || SCIPisHugeValue(scip, REALABS(scale)) )
            continue;

         /* update point */
         for( j = 0; j < SCIPgetNVars(scip); ++j )
         {
            SCIP_Real val;

            val = SCIPgetSolVal(scip, point, SCIPgetVars(scip)[j]) + scale * grad[j];
            SCIP_CALL( SCIPsetSolVal(scip, point, SCIPgetVars(scip)[j], val) );
         }
      }

      /* update violations */
      SCIP_CALL( getMaxViol(scip, point, maxviol) );
      SCIPdebugMessage("maxviol = %e\n", *maxviol);
   }
   while( ++r < nmaxiter
      && SCIPisFeasLT(scip, *maxviol, 0.0)
      && (*maxviol - lastmaxviol) / MAX(REALABS(*maxviol), REALABS(lastmaxviol)) >= minimprfac );

   SCIPfreeBufferArray(scip, &grad);

   return SCIP_OKAY;
}

/** sort points w.r.t their violations; filters out points with too large violation */
static
SCIP_RETCODE filterPoints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            points,             /**< array containing improved points */
   SCIP_Real*            violations,         /**< array containing violations (sorted) */
   int                   npoints,            /**< total number of points */
   int*                  nusefulpoints       /**< pointer to store the total number of useful points */
   )
{
   SCIP_Real maxviolation;
   int i;

   assert(points != NULL);
   assert(violations != NULL);
   assert(npoints > 0);
   assert(nusefulpoints != NULL);

   /* sort points w.r.t their violations; non-negative violations correspond to feasible points for the NLP */
   SCIPsortDownRealPtr(violations, (void**)points, npoints);

   maxviolation = violations[npoints - 1];
   *nusefulpoints = 0;

   for( i = 0; i < npoints; ++i )
   {
      /* consider feasible points always as useful; otherwise we take the a point if the violation is not too large,
       * compared to the maximum violation; the latter criterion becomes more strict with the number of selected points
       */
      if( SCIPisFeasGE(scip, violations[i], 0.0)
         || (1.0 - *nusefulpoints / (2.0 * npoints)) * maxviolation < violations[i] )
         ++(*nusefulpoints);
   }
   assert(*nusefulpoints >= 0);

   return SCIP_OKAY;
}

/** returns the relative distance between two points */
static
SCIP_Real getRelDistance(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             x,                  /**< first point */
   SCIP_SOL*             y                   /**< second point */
   )
{
   SCIP_VAR** vars;
   SCIP_Real distance;
   int i;

   vars = SCIPgetVars(scip);
   distance = 0.0;

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      distance += (SCIPgetSolVal(scip, x, vars[i]) - SCIPgetSolVal(scip, y, vars[i]))
         / (MAX(1.0, SCIPvarGetUbLocal(vars[i]) - SCIPvarGetLbLocal(vars[i])));
   }

   return distance;
}

/** cluster (useful) points with a greedy algorithm */
static
SCIP_RETCODE clusterPointsGreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            points,             /**< array containing improved points */
   int                   npoints,            /**< total number of points */
   int*                  clusteridx,         /**< array to store for each point the index of the cluster */
   SCIP_Real             maxdist             /**< maximum distance to consider points in the same cluster */
   )
{
   int ncluster;
   int i;

   assert(points != NULL);
   assert(npoints > 0);
   assert(clusteridx != NULL);
   assert(maxdist >= 0.0);

   BMSclearMemoryArray(clusteridx, npoints);

   ncluster = 0;

   for( i = 0; i < npoints; ++i )
   {
      int j;

      /* point is already in a cluster */
      if( clusteridx[i] != 0 )
         continue;

      /* create a new cluster for i */
      clusteridx[i] = ++ncluster;

      for( j = i + 1; j < npoints; ++j )
      {
         if( clusteridx[j] == 0 && getRelDistance(scip, points[i], points[j]) <= maxdist )
            clusteridx[j] = ncluster;
      }
   }

#ifndef NDEBUG
   for( i = 0; i < npoints; ++i )
   {
      assert(clusteridx[i] > 0);
      assert(clusteridx[i] <= ncluster);
   }
#endif

   return SCIP_OKAY;
}

/** calls the sub-NLP heuristic for a given cluster */
static
SCIP_RETCODE solveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< multi-start heuristic */
   SCIP_HEUR*            nlpheur,            /**< pointer to NLP local search heuristics */
   SCIP_SOL**            points,             /**< array containing improved points */
   int                   npoints,            /**< total number of points */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver */
   SCIP_Real             timelimit,          /**< time limit for NLP solver */
   SCIP_Real             minimprove,         /**< desired minimal relative improvement in objective function value */
   SCIP_Bool*            success             /**< pointer to store if we could find a solution */
   )
{
   SCIP_VAR** vars;
   SCIP_SOL* refpoint;
   SCIP_RESULT nlpresult;
   int i;

   assert(points != NULL);
   assert(npoints > 0);

   vars = SCIPgetVars(scip);
   *success = FALSE;

   SCIP_CALL( SCIPcreateSol(scip, &refpoint, heur) );

   /* compute reference point */
   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      SCIP_Real val;
      int p;

      val = 0.0;

      for( p = 0; p < npoints; ++p )
      {
         assert(points[p] != NULL);
         val += SCIPgetSolVal(scip, points[p], vars[i]);
      }

      SCIP_CALL( SCIPsetSolVal(scip, refpoint, vars[i], val / npoints) );
   }

   /* call sub-NLP heuristic */
   SCIP_CALL( SCIPapplyHeurSubNlp(scip, nlpheur, &nlpresult, refpoint, itercontingent, timelimit, minimprove, NULL) );

   if( nlpresult == SCIP_FOUNDSOL )
      *success = TRUE;

   /* free reference point */
   SCIP_CALL( SCIPfreeSol(scip, &refpoint) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE execHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< heuristic */
   )
{
   SCIP_HASHMAP* varindex;
   int i;

   assert(scip != NULL);
   assert(heur != NULL);

   SCIP_CALL( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(varindex, (void*)SCIPgetVars(scip)[i], (void*)(size_t)i) );
   }

   SCIPhashmapFree(&varindex);

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
