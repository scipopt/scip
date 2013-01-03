/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_pscostdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the pseudo cost values
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_pscostdiving.h"


#define HEUR_NAME             "pscostdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the pseudo cost values"
#define HEUR_DISPCHAR         'p'
#define HEUR_PRIORITY         -1002000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          2
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */

#define MINLPITER                 10000 /**< minimal number of LP iterations allowed in each LP solving call */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             minreldepth;        /**< minimal relative depth to start diving */
   SCIP_Real             maxreldepth;        /**< maximal relative depth to start diving */
   SCIP_Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int                   maxlpiterofs;       /**< additional number of allowed LP iterations */
   SCIP_Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed (0.0: no limit) */
   SCIP_Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool             backtrack;          /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Longint          nlpiterations;      /**< LP iterations used in this heuristic */
   int                   nsuccess;           /**< number of runs that produced at least one feasible solution */
};


/*
 * local methods
 */

/** calculates the pseudocost score for a given variable w.r.t. a given solution value and a given rounding direction */
static
void calcPscostQuot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             primsol,            /**< primal solution of variable */
   SCIP_Real             frac,               /**< fractionality of variable */
   int                   rounddir,           /**< -1: round down, +1: round up, 0: select due to pseudo cost values */
   SCIP_Real*            pscostquot,         /**< pointer to store pseudo cost quotient */
   SCIP_Bool*            roundup             /**< pointer to store whether the variable should be rounded up */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   assert(pscostquot != NULL);
   assert(roundup != NULL);
   assert(SCIPisEQ(scip, frac, primsol - SCIPfeasFloor(scip, primsol)));

   /* bound fractions to not prefer variables that are nearly integral */
   frac = MAX(frac, 0.1);
   frac = MIN(frac, 0.9);

   /* get pseudo cost quotient */
   pscostdown = SCIPgetVarPseudocostVal(scip, var, 0.0-frac);
   pscostup = SCIPgetVarPseudocostVal(scip, var, 1.0-frac);
   assert(pscostdown >= 0.0 && pscostup >= 0.0);

   /* choose rounding direction */
   if( rounddir == -1 )
      *roundup = FALSE;
   else if( rounddir == +1 )
      *roundup = TRUE;
   else if( primsol < SCIPvarGetRootSol(var) - 0.4 )
      *roundup = FALSE;
   else if( primsol > SCIPvarGetRootSol(var) + 0.4 )
      *roundup = TRUE;
   else if( frac < 0.3 )
      *roundup = FALSE;
   else if( frac > 0.7 )
      *roundup = TRUE;
   else if( pscostdown < pscostup )
      *roundup = FALSE;
   else
      *roundup = TRUE;

   /* calculate pseudo cost quotient */
   if( *roundup )
      *pscostquot = sqrt(frac) * (1.0+pscostdown) / (1.0+pscostup);
   else
      *pscostquot = sqrt(1.0-frac) * (1.0+pscostup) / (1.0+pscostdown);

   /* prefer decisions on binary variables */
   if( SCIPvarIsBinary(var) )
      (*pscostquot) *= 1000.0;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyPscostdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurPscostdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreePscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;
   heurdata->nsuccess = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_VAR* var;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real primsol;
   SCIP_Real frac;
   SCIP_Real pscostquot;
   SCIP_Real bestpscostquot;
   SCIP_Bool bestcandmayrounddown;
   SCIP_Bool bestcandmayroundup;
   SCIP_Bool bestcandroundup;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_Bool backtracked;
   SCIP_Longint ncalls;
   SCIP_Longint nsolsfound;
   SCIP_Longint nlpiterations;
   SCIP_Longint maxnlpiterations;
   int nlpcands;
   int startnlpcands;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP solution is basic (which allows fast resolve in diving) */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) && SCIPgetDepth(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only try to dive, if we are in the correct part of the tree, given by minreldepth and maxreldepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->minreldepth*maxdepth || depth > heurdata->maxreldepth*maxdepth )
      return SCIP_OKAY;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   nlpiterations = SCIPgetNNodeLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = 10*SCIPheurGetNBestSolsFound(heur) + heurdata->nsuccess;
   maxnlpiterations = (SCIP_Longint)((1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * nlpiterations);
   maxnlpiterations += heurdata->maxlpiterofs;

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + MINLPITER);

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   /* don't try to dive, if there are no fractional variables */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      if( heurdata->maxdiveubquotnosol > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquotnosol > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   else
   {
      if( heurdata->maxdiveubquot > 0.0 )
         searchubbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      else
         searchubbound = SCIPinfinity(scip);
      if( heurdata->maxdiveavgquot > 0.0 )
         searchavgbound = SCIPgetLowerbound(scip)
            + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
      else
         searchavgbound = SCIPinfinity(scip);
   }
   searchbound = MIN(searchubbound, searchavgbound);
   if( SCIPisObjIntegral(scip) )
      searchbound = SCIPceil(scip, searchbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;


   *result = SCIP_DIDNOTFIND;

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* enables collection of variable statistics during probing */
   SCIPenableVarHistory(scip);

   /* get LP objective value */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing pscostdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   cutoff = FALSE;
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startnlpcands = nlpcands;
   while( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound))
      && !SCIPisStopped(scip) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round variable with largest rel. difference of pseudo cost values in corresponding
       *     direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable in the objective value direction
       */
      bestcand = -1;
      bestpscostquot = -1.0;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         primsol = lpcandssol[c];
         frac = lpcandsfrac[c];
         if( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
            if( bestcandmayrounddown || bestcandmayroundup )
            {
               /* choose rounding direction:
                * - if variable may be rounded in both directions, round corresponding to the pseudo cost values
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               roundup = FALSE;
               if( mayrounddown && mayroundup )
                  calcPscostQuot(scip, var, primsol, frac, 0, &pscostquot, &roundup);
               else if( mayrounddown )
                  calcPscostQuot(scip, var, primsol, frac, +1, &pscostquot, &roundup);
               else
                  calcPscostQuot(scip, var, primsol, frac, -1, &pscostquot, &roundup);

               /* check, if candidate is new best candidate */
               if( pscostquot > bestpscostquot )
               {
                  bestcand = c;
                  bestpscostquot = pscostquot;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded: calculate pseudo cost quotient and preferred direction */
            calcPscostQuot(scip, var, primsol, frac, 0, &pscostquot, &roundup);

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if( bestcandmayrounddown || bestcandmayroundup || pscostquot > bestpscostquot )
            {
               bestcand = c;
               bestpscostquot = pscostquot;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
         }
      }
      assert(bestcand != -1);

      /* if all candidates are roundable, try to round the solution */
      if( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            SCIPdebugMessage("pscostdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = lpcands[bestcand];

      backtracked = FALSE;
      do
      {
         /* if the variable is already fixed or if the solution value is outside the domain, numerical troubles may have
          * occured or variable was fixed by propagation while backtracking => Abort diving!
          */
         if( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
         {
            SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), lpcandssol[bestcand]);
            cutoff = TRUE;
            break;
         }
         if( SCIPisFeasLT(scip, lpcandssol[bestcand], SCIPvarGetLbLocal(var)) || SCIPisFeasGT(scip, lpcandssol[bestcand], SCIPvarGetUbLocal(var)) )
         {
            SCIPdebugMessage("selected variable's <%s> solution value is outside the domain [%g,%g] (solval: %.9f), diving aborted\n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), lpcandssol[bestcand]);
            assert(backtracked);
            break;
         }

         /* apply rounding of best candidate */
         if( bestcandroundup == !backtracked )
         {
            /* round variable up */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPfeasCeil(scip, lpcandssol[bestcand]), SCIPvarGetUbLocal(var));
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );
         }
         else
         {
            /* round variable down */
            SCIPdebugMessage("  dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               lpcandssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPvarGetLbLocal(var), SCIPfeasFloor(scip, lpcandssol[bestcand]));
            SCIP_CALL( SCIPchgVarUbProbing(scip, lpcands[bestcand], SCIPfeasFloor(scip, lpcandssol[bestcand])) );
         }

         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if( !cutoff )
         {
            /* resolve the diving LP */
            /* Errors in the LP solver should not kill the overall solving process, if the LP is just needed for a heuristic.
             * Hence in optimized mode, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
             */
#ifdef NDEBUG
            SCIP_RETCODE retstat;
            nlpiterations = SCIPgetNLPIterations(scip);
            retstat = SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror);
            if( retstat != SCIP_OKAY )
            {
               SCIPwarningMessage(scip, "Error while solving LP in Pscostdiving heuristic; LP solve terminated with code <%d>\n",retstat);
            }
#else
            nlpiterations = SCIPgetNLPIterations(scip);
            SCIP_CALL( SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror) );
#endif

            if( lperror )
               break;

            /* update iteration count */
            heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

            /* get LP solution status, objective value, and fractional variables, that should be integral */
            lpsolstat = SCIPgetLPSolstat(scip);
            cutoff = (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
         }

         /* perform backtracking if a cutoff was detected */
         if( cutoff && !backtracked && heurdata->backtrack )
         {
            SCIPdebugMessage("  *** cutoff detected at level %d - backtracking\n", SCIPgetProbingDepth(scip));
            SCIP_CALL( SCIPbacktrackProbing(scip, SCIPgetProbingDepth(scip)-1) );
            SCIP_CALL( SCIPnewProbingNode(scip) );
            backtracked = TRUE;
         }
         else
            backtracked = FALSE;
      }
      while( backtracked );

      if( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update pseudo cost values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, lpcands[bestcand], 1.0-lpcandsfrac[bestcand],
                     objval - oldobjval, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, lpcands[bestcand], 0.0-lpcandsfrac[bestcand],
                     objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );
      }
      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d, lpiterations=%"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT"\n",
         lpsolstat, objval, searchbound, nlpcands, heurdata->nlpiterations, maxnlpiterations);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMessage("pscostdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );

   if( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMessage("pscostdiving heuristic finished\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the pscostdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurPscostdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Pscostdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPscostdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyPscostdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreePscostdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitPscostdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitPscostdiving) );

   /* pscostdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/pscostdiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/pscostdiving/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );

   return SCIP_OKAY;
}

