/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_indcoefdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the matrix coefficients, indicator variables are used as well
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_indcoefdiving.h"
#include <scip/cons_indicator.h>

#define HEUR_NAME             "indcoefdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the matrix coefficients, indicator variables are used as well"
#define HEUR_DISPCHAR         '@'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE



/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH            0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH            1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT         0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS          1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT          0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                            *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT         0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                            *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL     0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL    0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK             TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_ONLYINDICATORS       FALSE /**< use only indicator variables as branching candidates */

#define MINLPITER                    10000 /**< minimal number of LP iterations allowed in each LP solving call */



/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*        sol;                 /**< working solution */
   SCIP_Real        minreldepth;         /**< minimal relative depth to start diving */
   SCIP_Real        maxreldepth;         /**< maximal relative depth to start diving */
   SCIP_Real        maxlpiterquot;       /**< maximal fraction of diving LP iterations compared to node LP iterations */
   int              maxlpiterofs;        /**< additional number of allowed LP iterations */
   SCIP_Real        maxdiveubquot;       /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                          *   where diving is performed (0.0: no limit) */
   SCIP_Real        maxdiveavgquot;      /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                          *   where diving is performed (0.0: no limit) */
   SCIP_Real        maxdiveubquotnosol;  /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Real        maxdiveavgquotnosol; /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
   SCIP_Bool        backtrack;           /**< use one level of backtracking if infeasibility is encountered? */
   SCIP_Bool        onlyindicators;      /**< use only indicator variables as branching candidates */
   SCIP_Longint     nlpiterations;       /**< LP iterations used in this heuristic */
   int              nsuccess;            /**< number of runs that produced at least one feasible solution */
};




/*
 * local methods
 */

/** get candidate variables */
static
SCIP_RETCODE getCandVars(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_VAR**            cands,               /**< candidate variables */
   SCIP_Real*            candssol,            /**< solution values of candidates */
   SCIP_Real*            candsfrac,           /**< fractionalities of candidates */
   int*                  ncands,              /**< number of candidates */
   int*                  nlpcands,            /**< number of fractional candidates */
   SCIP_Bool             onlyindicators       /**< use only indicator variables */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   int c;
#ifndef NDEBUG
   int nMaxVars;

   nMaxVars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
#endif

   *nlpcands = 0;
   if ( !onlyindicators )
   {
      SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, nlpcands, NULL, NULL) );
      assert( *nlpcands <= nMaxVars );

      /* copy variables and enter them into a hash set */
      for (c = 0; c < *nlpcands; ++c)
      {
         cands[c] = lpcands[c];
         candsfrac[c] = lpcandsfrac[c];
         candssol[c] = lpcandssol[c];
      }
      *ncands = *nlpcands;
   }
   else
      *ncands = 0;

   /* now check indicator constraints */
   conshdlr = SCIPfindConshdlr(scip, "indicator");
   assert( conshdlr != NULL );

   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   for (c = 0; c < nconss; ++c)
   {
      /* check whether constraint is violated */
      if ( SCIPisViolatedIndicator(scip, conss[c], NULL) )
      {
         SCIP_VAR* binvar;
         SCIP_Real val;

         binvar = SCIPgetBinaryVarIndicator(conss[c]);
         val = SCIPgetSolVal(scip, NULL, binvar);
         if ( SCIPisFeasIntegral(scip, val) )
         {
            cands[*ncands] = binvar;
            candssol[*ncands] = val;
            candsfrac[*ncands] = SCIPfrac(scip, val);
            ++(*ncands);
            assert( *ncands <= nMaxVars );
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIndcoefdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitIndcoefdiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXIT(heurExitIndcoefdiving) /*lint --e{715}*/
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


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolIndcoefdiving NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolIndcoefdiving NULL

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIndcoefdiving)
{
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIndcoefdiving(scip) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIndcoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_LPSOLSTAT lpsolstat;
   SCIP_VAR* var;
   SCIP_VAR** cands;
   SCIP_Real* candssol;
   SCIP_Real* candsfrac;
   SCIP_Real searchubbound;
   SCIP_Real searchavgbound;
   SCIP_Real searchbound;
   SCIP_Real objval;
   SCIP_Real oldobjval;
   SCIP_Real frac;
   SCIP_Real bestfrac;
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
   int ncands;
   int nlpcands;
   int startncands;
   int nMaxVars;
   int depth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int nlocksdown;
   int nlocksup;
   int nviolrows;
   int bestnviolrows;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if ( ! SCIPhasCurrentNodeLP(scip) || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
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
   nMaxVars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = nMaxVars;
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   *result = SCIP_DIDNOTFIND;

   /* get storage for candidate variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &cands, nMaxVars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candsfrac, nMaxVars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candssol, nMaxVars) );

   /* get candidate variables - *before* probing is started */
   SCIP_CALL( getCandVars(scip, cands, candssol, candsfrac, &ncands, &nlpcands, heurdata->onlyindicators) );

   /* start diving */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* get LP objective value, and fractional variables, that should be integral */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") executing indcoefdiving heuristic: depth=%d, %d fractionals, %d candidates, dualbound=%g, avgbound=%g, cutoffbound=%g, searchbound=%g\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, ncands, SCIPgetDualbound(scip), SCIPgetAvgDualbound(scip),
      SCIPretransformObj(scip, SCIPgetCutoffbound(scip)), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   cutoff = FALSE;
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startncands = ncands;
   while( !SCIPisStopped(scip) && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && ncands > 0
      && (divedepth < 10
         || ncands <= startncands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound)) )
   {
      SCIP_CALL( SCIPnewProbingNode(scip) );
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round variable with least number of locks in corresponding direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable with least number of locks in opposite of its feasible rounding direction
       */
      bestcand = -1;
      bestnviolrows = INT_MAX;
      bestfrac = SCIP_INVALID;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      for (c = 0; c < ncands; ++c)
      {
	 var = cands[c];
	 frac = candsfrac[c];

         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);

         if ( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
            if ( bestcandmayrounddown || bestcandmayroundup )
            {
               /* choose rounding direction:
                * - if variable may be rounded in both directions, round corresponding to the fractionality
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               if ( mayrounddown && mayroundup )
                  roundup = (frac > 0.5);
               else
                  roundup = mayrounddown;

               if ( roundup )
               {
                  frac = 1.0 - frac;
                  nviolrows = SCIPvarGetNLocksUp(var);
               }
               else
                  nviolrows = SCIPvarGetNLocksDown(var);

               /* penalize too small fractions */
               if ( frac < 0.01 )
                  nviolrows *= 100;

               /* prefer decisions on binary variables */
               if ( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
                  nviolrows *= 100;

               /* check, if candidate is new best candidate */
               if ( nviolrows + frac < bestnviolrows + bestfrac )
               {
                  bestcand = c;
                  bestnviolrows = nviolrows;
                  bestfrac = frac;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded */
            nlocksdown = SCIPvarGetNLocksDown(var);
            nlocksup = SCIPvarGetNLocksUp(var);
            roundup = (nlocksdown > nlocksup || (nlocksdown == nlocksup && frac > 0.5));
            if ( roundup )
            {
               nviolrows = nlocksup;
               frac = 1.0 - frac;
            }
            else
               nviolrows = nlocksdown;

            /* penalize too small fractions */
            if ( frac < 0.01 )
               nviolrows *= 100;

            /* prefer decisions on binary variables */
            if ( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
               nviolrows *= 100;

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if ( bestcandmayrounddown || bestcandmayroundup || nviolrows + frac < bestnviolrows + bestfrac )
            {
               bestcand = c;
               bestnviolrows = nviolrows;
               bestfrac = frac;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
            assert(bestfrac < SCIP_INVALID);
         }
      }
      assert( bestcand >= 0 );

      /* if all candidates are roundable, try to round the solution */
      if ( bestcandmayrounddown || bestcandmayroundup )
      {
         SCIP_Bool success;

         /* create solution from diving LP and try to round it */
         SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
         SCIP_CALL( SCIProundSol(scip, heurdata->sol, &success) );

         if ( success )
         {
            SCIPdebugMessage("indcoefdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

            /* try to add solution to SCIP */
            SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );

            /* check, if solution was feasible and good enough */
            if ( success )
            {
               SCIPdebugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = cands[bestcand];

      backtracked = FALSE;
      do
      {
         /* if the variable is already fixed, numerical troubles may have occured or
          * variable was fixed by propagation while backtracking => Abort diving!
          */
         if ( SCIPvarGetLbLocal(var) >= SCIPvarGetUbLocal(var) - 0.5 )
         {
            SCIPdebugMessage("Selected variable <%s> already fixed to [%g,%g] (solval: %.9f), diving aborted \n",
               SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), candssol[bestcand]);
            cutoff = TRUE;
            break;
         }

         /* apply rounding of best candidate */
         if ( bestcandroundup == !backtracked )
         {
	    SCIP_Real value = SCIPfeasCeil(scip, candssol[bestcand]);

	    if ( SCIPisFeasIntegral(scip, candssol[bestcand]) )
	    {
	       assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
	       value = 1.0;
	    }

            /* round variable up */
            SCIPdebugMessage("  up: dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               candssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               value, SCIPvarGetUbLocal(var));

	    SCIP_CALL( SCIPchgVarLbProbing(scip, var, value) );
         }
         else
         {
	    SCIP_Real value = SCIPfeasFloor(scip, candssol[bestcand]);
	    if ( SCIPisFeasIntegral(scip, candssol[bestcand]) )
	    {
	       assert( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY );
	       value = 0.0;
	    }

            /* round variable down */
            SCIPdebugMessage("  down: dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT": var <%s>, round=%u/%u, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
               divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
               SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
               candssol[bestcand], SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var),
               SCIPvarGetLbLocal(var), value);

	    SCIP_CALL( SCIPchgVarUbProbing(scip, cands[bestcand], value) );
         }

         /* apply domain propagation */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if ( !cutoff )
         {
            /* resolve the diving LP */
            nlpiterations = SCIPgetNLPIterations(scip);
            SCIP_CALL( SCIPsolveProbingLP(scip, MAX((int)(maxnlpiterations - heurdata->nlpiterations), MINLPITER), &lperror, &cutoff) );
            if( lperror )
               break;

            /* update iteration count */
            heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

            /* get LP solution status, objective value, and fractional variables, that should be integral */
            lpsolstat = SCIPgetLPSolstat(scip);
            cutoff = (lpsolstat == SCIP_LPSOLSTAT_OBJLIMIT || lpsolstat == SCIP_LPSOLSTAT_INFEASIBLE);
         }

         /* perform backtracking if a cutoff was detected */
         if ( cutoff && !backtracked && heurdata->backtrack )
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

      if ( !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update pseudo cost values */
         if ( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, cands[bestcand], 1.0-candsfrac[bestcand],
                     objval - oldobjval, 1.0) );
            }
            else
            {
               SCIP_CALL( SCIPupdateVarPseudocost(scip, cands[bestcand], 0.0-candsfrac[bestcand],
                     objval - oldobjval, 1.0) );
            }
         }

	 /* get new candidate variables */
	 SCIP_CALL( getCandVars(scip, cands, candssol, candsfrac, &ncands, &nlpcands, heurdata->onlyindicators) );
      }
      SCIPdebugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d, ncand=%d\n", lpsolstat, objval, searchbound, nlpcands, ncands);
   }

   /* check if a solution has been found */
   if ( ncands == 0 && !lperror && !cutoff && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SCIP_Bool success;

      /* create solution from diving LP */
      SCIP_CALL( SCIPlinkLPSol(scip, heurdata->sol) );
      SCIPdebugMessage("indcoefdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, heurdata->onlyindicators, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if ( success )
      {
         SCIPdebugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   SCIPfreeBufferArray(scip, &cands);
   SCIPfreeBufferArray(scip, &candsfrac);
   SCIPfreeBufferArray(scip, &candssol);

   /* end diving */
   SCIP_CALL( SCIPendProbing(scip) );

   if ( *result == SCIP_FOUNDSOL )
      heurdata->nsuccess++;

   SCIPdebugMessage("(node %"SCIP_LONGINT_FORMAT") finished indcoefdiving heuristic: %d fractionals, %d candidates, dive %d/%d, LP iter %"SCIP_LONGINT_FORMAT"/%"SCIP_LONGINT_FORMAT", objval=%g/%g, lpsolstat=%d, cutoff=%u\n",
      SCIPgetNNodes(scip), nlpcands, ncands, divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
      SCIPretransformObj(scip, objval), SCIPretransformObj(scip, searchbound), lpsolstat, cutoff);

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the indcoefdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIndcoefdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyIndcoefdiving, heurFreeIndcoefdiving, heurInitIndcoefdiving, heurExitIndcoefdiving,
         heurInitsolIndcoefdiving, heurExitsolIndcoefdiving, heurExecIndcoefdiving,
         heurdata) );

   /* indcoefdiving heuristic parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/minreldepth",
         "minimal relative depth to start diving",
         &heurdata->minreldepth, TRUE, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/maxreldepth",
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, TRUE, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/maxlpiterquot",
         "maximal fraction of diving LP iterations compared to node LP iterations",
         &heurdata->maxlpiterquot, FALSE, DEFAULT_MAXLPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "heuristics/indcoefdiving/maxlpiterofs",
         "additional number of allowed LP iterations",
         &heurdata->maxlpiterofs, FALSE, DEFAULT_MAXLPITEROFS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveubquot, TRUE, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/maxdiveavgquot",
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed (0.0: no limit)",
         &heurdata->maxdiveavgquot, TRUE, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/maxdiveubquotnosol",
         "maximal UBQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveubquotnosol, TRUE, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "heuristics/indcoefdiving/maxdiveavgquotnosol",
         "maximal AVGQUOT when no solution was found yet (0.0: no limit)",
         &heurdata->maxdiveavgquotnosol, TRUE, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/indcoefdiving/backtrack",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->backtrack, FALSE, DEFAULT_BACKTRACK, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/indcoefdiving/onlyindicators",
         "use one level of backtracking if infeasibility is encountered?",
         &heurdata->onlyindicators, FALSE, DEFAULT_ONLYINDICATORS, NULL, NULL) );

   return SCIP_OKAY;
}
