/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_pscostdiving.c,v 1.18 2005/02/03 16:57:45 bzfpfend Exp $"

/**@file   heur_pscostdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the pseudo cost values
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_pscostdiving.h"


#define HEUR_NAME             "pscostdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the pseudo cost values"
#define HEUR_DISPCHAR         'p'
#define HEUR_PRIORITY         -1002000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          2
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      FALSE     /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   FALSE     /* call heuristic during plunging? (should be FALSE for diving heuristics!) */



/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0  /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0  /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT       0.02 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_MAXDIVEUBQUOT       0.8  /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                          *   where diving is performed */
#define DEFAULT_MAXDIVEAVGQUOT      4.0  /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                          *   where diving is performed */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1  /**< maximal UBQUOT when no solution was found yet */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 8.0  /**< maximal AVGQUOT when no solution was found yet */



/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   Real             minreldepth;        /**< minimal relative depth to start diving */
   Real             maxreldepth;        /**< maximal relative depth to start diving */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   Real             maxdiveubquot;      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveavgquot;     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet */
   Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet */
   Longint          nlpiterations;      /**< LP iterations used in this heuristic */
};




/*
 * local methods
 */

static
void calcPscostQuot(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             primsol,            /**< primal solution of variable */
   Real             frac,               /**< fractionality of variable */
   int              rounddir,           /**< -1: round down, +1: round up, 0: select due to pseudo cost values */
   Real*            pscostquot,         /**< pointer to store pseudo cost quotient */
   Bool*            roundup             /**< pointer to store whether the variable should be rounded up */
   )
{
   Real pscostdown;
   Real pscostup;

   assert(pscostquot != NULL);
   assert(roundup != NULL);
   assert(SCIPisEQ(scip, frac, primsol - SCIPfeasFloor(scip, primsol)));

   /* bound fractions to not prefer variables that are nearly integral */
   frac = MAX(frac, 0.1);
   frac = MIN(frac, 0.9);
   
   /* get pseudo cost quotient */
   pscostdown = SCIPgetVarPseudocost(scip, var, 0.0-frac);
   pscostup = SCIPgetVarPseudocost(scip, var, 1.0-frac);
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
   if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      (*pscostquot) *= 1000.0;
}




/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
DECL_HEURFREE(heurFreePscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

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
DECL_HEURINIT(heurInitPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   CHECK_OKAY( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data */
   heurdata->nlpiterations = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
DECL_HEUREXIT(heurExitPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecPscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;
   LPSOLSTAT lpsolstat;
   VAR* var;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real searchubbound;
   Real searchavgbound;
   Real searchbound;
   Real objval;
   Real oldobjval;
   Real primsol;
   Real frac;
   Real pscostquot;
   Real bestpscostquot;
   Bool bestcandmayrounddown;
   Bool bestcandmayroundup;
   Bool bestcandroundup;
   Bool mayrounddown;
   Bool mayroundup;
   Bool roundup;
   Bool lperror;
   Longint ncalls;
   Longint nsolsfound;
   Longint nlpiterations;
   Longint maxnlpiterations;
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

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) )
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
   nlpiterations = SCIPgetNLPIterations(scip);
   ncalls = SCIPheurGetNCalls(heur);
   nsolsfound = SCIPheurGetNSolsFound(heur);
   maxnlpiterations = (1.0 + 10.0*(nsolsfound+1.0)/(ncalls+1.0)) * heurdata->maxlpiterquot * (nlpiterations + 10000);

   /* don't try to dive, if we took too many LP iterations during diving */
   if( heurdata->nlpiterations >= maxnlpiterations )
      return SCIP_OKAY;

   /* allow at least a certain number of LP iterations in this dive */
   maxnlpiterations = MAX(maxnlpiterations, heurdata->nlpiterations + 10000);

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      searchubbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveubquotnosol * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      searchavgbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   }
   else
   {
      searchubbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveubquot * (SCIPgetCutoffbound(scip) - SCIPgetLowerbound(scip));
      searchavgbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   }
   searchbound = MIN(searchubbound, searchavgbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;


   *result = SCIP_DIDNOTFIND;

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );

   /* get LP objective value, and fractional variables, that should be integral */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   debugMessage("(node %lld) executing pscostdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n", 
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective, depth and iteration limits and fractional variables exist, but
    * - if the last rounding was in a direction, that never destroys feasibility, we continue in any case
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   lperror = FALSE;
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startnlpcands = nlpcands;
   while( !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (bestcandmayrounddown || bestcandmayroundup
         || divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations && objval < searchbound)) )
   {
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
         Bool success;
         
         /* create solution from diving LP and try to round it */
         CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
         CHECK_OKAY( SCIProundSol(scip, heurdata->sol, &success) );

         if( success )
         {
            debugMessage("pscostdiving found roundable primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));
         
            /* try to add solution to SCIP */
            CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, &success) );
            
            /* check, if solution was feasible and good enough */
            if( success )
            {
               debugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = lpcands[bestcand];

      if( SCIPgetVarLbDive(scip, var) >= SCIPgetVarUbDive(scip, var) - 0.5 )
      {
         /* the variable is already fixed -> numerical troubles -> abort diving */
         debugMessage("numerical troubles: selected variable <%s> already fixed to [%g,%g] (solval: %.9f)\n",
            SCIPvarGetName(var), SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var), lpcandssol[bestcand]);
         break;
      }

      /* apply rounding of best candidate */
      if( bestcandroundup )
      {
         /* round variable up */
         debugMessage("  dive %d/%d, LP iter %lld/%lld: var <%s>, round=%d/%d, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var),
            SCIPfeasCeil(scip, lpcandssol[bestcand]), SCIPgetVarUbDive(scip, var));
         CHECK_OKAY( SCIPchgVarLbDive(scip, var, SCIPfeasCeil(scip, lpcandssol[bestcand])) );
      }
      else
      {
         /* round variable down */
         debugMessage("  dive %d/%d, LP iter %lld/%lld: var <%s>, round=%d/%d, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var),
            SCIPgetVarLbDive(scip, var), SCIPfeasFloor(scip, lpcandssol[bestcand]));
         CHECK_OKAY( SCIPchgVarUbDive(scip, lpcands[bestcand], SCIPfeasFloor(scip, lpcandssol[bestcand])) );
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip, maxnlpiterations, &lperror) );
      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
      nlpiterations = SCIPgetNLPIterations(scip);

      /* get LP solution status, objective value, and fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update pseudo cost values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[bestcand], 1.0-lpcandsfrac[bestcand], 
                     objval - oldobjval, 1.0) );
            }
            else
            {
               CHECK_OKAY( SCIPupdateVarPseudocost(scip, lpcands[bestcand], 0.0-lpcandsfrac[bestcand], 
                     objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );
      }
      debugMessage("   -> lpsolstat=%d, objval=%g/%g, nfrac=%d, lpiterations=%lld/%lld\n", 
         lpsolstat, objval, searchbound, nlpcands, heurdata->nlpiterations, maxnlpiterations);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      Bool success;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      debugMessage("pscostdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         debugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   CHECK_OKAY( SCIPendDive(scip) );

   debugMessage("pscostdiving heuristic finished\n");

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the pscostdiving heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurPscostdiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING,
         heurFreePscostdiving, heurInitPscostdiving, heurExitPscostdiving, heurExecPscostdiving,
         heurdata) );

   /* pscostdiving heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/minreldepth", 
         "minimal relative depth to start diving",
         &heurdata->minreldepth, DEFAULT_MINRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxreldepth", 
         "maximal relative depth to start diving",
         &heurdata->maxreldepth, DEFAULT_MAXRELDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxlpiterquot", 
         "maximal fraction of diving LP iterations compared to total iteration number",
         &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveubquot",
         "maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound) where diving is performed",
         &heurdata->maxdiveubquot, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveavgquot", 
         "maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed",
         &heurdata->maxdiveavgquot, DEFAULT_MAXDIVEAVGQUOT, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveubquotnosol", 
         "maximal UBQUOT when no solution was found yet",
         &heurdata->maxdiveubquotnosol, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
         "heuristics/pscostdiving/maxdiveavgquotnosol", 
         "maximal AVGQUOT when no solution was found yet",
         &heurdata->maxdiveavgquotnosol, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

