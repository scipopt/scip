/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_objpscostdiving.c,v 1.4 2004/06/02 07:39:07 bzfpfend Exp $"

/**@file   heur_objpscostdiving.c
 * @brief  LP diving heuristic that changes variable's objective value instead of bounds, using pseudo cost values as guide
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_objpscostdiving.h"


#define HEUR_NAME         "objpscostdiving"
#define HEUR_DESC         "LP diving heuristic that changes variable's objective values instead of bounds, using pseudo costs as guide"
#define HEUR_DISPCHAR     'o'
#define HEUR_PRIORITY     -1003000
#define HEUR_FREQ         10
#define HEUR_FREQOFS      6
#define HEUR_MAXDEPTH     -1
#define HEUR_PSEUDONODES  FALSE         /** call heuristic at nodes where only a pseudo solution exist? */



/*
 * Default parameter settings
 */

#define DEFAULT_DIVESTARTDEPTH      0.5 /**< minimal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT       0.1 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_DEPTHFAC            0.5 /**< maximal diving depth: number of binary/integer variables times depthfac */
#define DEFAULT_DEPTHFACNOSOL       2.0 /**< maximal diving depth factor if no feasible solution was found yet */


/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   Real             divestartdepth;     /**< minimal relative depth to start diving */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   Real             depthfac;           /**< maximal diving depth: number of binary/integer variables times depthfac */
   Real             depthfacnosol;      /**< maximal diving depth factor if no feasible solution was found yet */
   Longint          nlpiterations;      /**< LP iterations used in this heuristic */
};




/*
 * local methods
 */

static
void calcPscostQuot(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
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
DECL_HEURFREE(heurFreeObjpscostdiving) /*lint --e{715}*/
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
DECL_HEURINIT(heurInitObjpscostdiving) /*lint --e{715}*/
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
DECL_HEUREXIT(heurExitObjpscostdiving) /*lint --e{715}*/
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
DECL_HEUREXEC(heurExecObjpscostdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;
   LPSOLSTAT lpsolstat;
   VAR* var;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real frac;
   Real pscostquot;
   Real bestpscostquot;
   Real oldobj;
   Real newobj;
   Real objscale;
   Bool bestcandmayrounddown;
   Bool bestcandmayroundup;
   Bool bestcandroundup;
   Bool mayrounddown;
   Bool mayroundup;
   Bool roundup;
   Bool lperror;
   Longint nlpiterations;
   Longint maxnlpiterations;
   int* roundings;
   int nvars;
   int varidx;
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
   assert(SCIPhasActNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNNodes(scip) )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't try to dive, if we are in the higher fraction of the tree, given by divestartdepth */
   depth = SCIPgetDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   maxdepth = MAX(maxdepth, 30);
   if( depth < heurdata->divestartdepth*maxdepth )
      return SCIP_OKAY;

   /* don't try to dive, if we took too many LP iterations during diving */
   nlpiterations = SCIPgetNLPIterations(scip);
   if( heurdata->nlpiterations > heurdata->maxlpiterquot*nlpiterations + 1000 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* calculate the maximal number of LP iterations until heuristic is aborted */
   maxnlpiterations = heurdata->maxlpiterquot*nlpiterations - heurdata->nlpiterations;
   maxnlpiterations = MAX(maxnlpiterations, 10000);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   if( SCIPgetNSolsFound(scip) == 0 )
      maxdivedepth = heurdata->depthfacnosol * nvars;
   else
      maxdivedepth = heurdata->depthfac * nvars;

   /* get temporary memory for remembering the current soft roundings */
   CHECK_OKAY( SCIPallocBufferArray(scip, &roundings, nvars) );
   clearMemoryArray(roundings, nvars);

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );

   /* get LP objective value, and fractional variables, that should be integral */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );

   debugMessage("(node %lld) executing objpscostdiving heuristic: depth=%d, %d fractionals, dualbound=%g, maxnlpiterations=%lld, maxdivedepth=%d\n", 
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nlpcands, SCIPgetDualbound(scip), maxnlpiterations, maxdivedepth);

   /* dive as long we are in the given diving depth and iteration limits and fractional variables exist, but
    * - if the last objective change was in a direction, that corresponds to a feasibile rounding, we continue in any case
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
         || (divedepth < maxdivedepth && heurdata->nlpiterations < maxnlpiterations)) )
   {
      divedepth++;

      /* choose variable for objective change:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, change objective value of variable with largest rel. difference of pseudo cost values
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - change objective value of variable with largest rel. difference of pseudo cost values
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
                  calcPscostQuot(scip, var, frac, 0, &pscostquot, &roundup);
               else if( mayrounddown )
                  calcPscostQuot(scip, var, frac, +1, &pscostquot, &roundup);
               else
                  calcPscostQuot(scip, var, frac, -1, &pscostquot, &roundup);

               /* prefer variables, that have already been softrounded but failed to get integral */
               varidx = SCIPvarGetProbindex(var);
               assert(0 <= varidx && varidx < nvars);
               if( roundings[varidx] != 0 )
                  pscostquot *= 1000.0;

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
            calcPscostQuot(scip, var, frac, 0, &pscostquot, &roundup);

            /* prefer variables, that have already been softrounded but failed to get integral */
            varidx = SCIPvarGetProbindex(var);
            assert(0 <= varidx && varidx < nvars);
            if( roundings[varidx] != 0 )
               pscostquot *= 1000.0;

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
            debugMessage("objpscostdiving found roundable primal solution: obj=%g\n", 
               SCIPgetSolOrigObj(scip, heurdata->sol));
         
            /* try to add solution to SCIP */
            CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, &success) );
            
            /* check, if solution was feasible and good enough */
            if( success )
            {
               debugMessage(" -> solution was feasible and good enough\n");
               *result = SCIP_FOUNDSOL;
            }
         }
      }

      var = lpcands[bestcand];

      /* check, if the best candidate was already subject to soft rounding */
      varidx = SCIPvarGetProbindex(var);
      assert(0 <= varidx && varidx < nvars);
      if( roundings[varidx] == +1 )
      {
         /* variable was already soft rounded upwards: hard round it downwards */
         CHECK_OKAY( SCIPchgVarUbDive(scip, var, SCIPfloor(scip, lpcandssol[bestcand])) );
         debugMessage("  dive %d/%d: var <%s>, round=%d/%d, sol=%g, was already soft rounded upwards -> bounds=[%g,%g]\n",
            divedepth, maxdivedepth, SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var));
      }
      else if( roundings[varidx] == -1 )
      {
         /* variable was already soft rounded downwards: hard round it upwards */
         CHECK_OKAY( SCIPchgVarLbDive(scip, var, SCIPceil(scip, lpcandssol[bestcand])) );
         debugMessage("  dive %d/%d: var <%s>, round=%d/%d, sol=%g, was already soft rounded downwards -> bounds=[%g,%g]\n",
            divedepth, maxdivedepth, SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var));
      }
      else
      {
         assert(roundings[varidx] == 0);

         /* apply soft rounding of best candidate via a change in the objective value */
         objscale = divedepth * 1000.0;
         oldobj = SCIPgetVarObjDive(scip, var);
         if( bestcandroundup )
         {
            /* soft round variable up: make objective value (more) negative */
            if( oldobj < 0.0 )
               newobj = objscale * oldobj;
            else
               newobj = -objscale * oldobj;
            newobj = MIN(newobj, -objscale);
            
            /* remember, that this variable was soft rounded upwards */
            roundings[varidx] = +1;
         }
         else
         {
            /* soft round variable down: make objective value (more) positive */
            if( oldobj > 0.0 )
               newobj = objscale * oldobj;
            else
               newobj = -objscale * oldobj;
            newobj = MAX(newobj, objscale);
            
            /* remember, that this variable was soft rounded downwards */
            roundings[varidx] = -1;
         }
         CHECK_OKAY( SCIPchgVarObjDive(scip, var, newobj) );
         debugMessage("  dive %d/%d, LP iter %lld/%lld: var <%s>, round=%d/%d, sol=%g, bounds=[%g,%g], obj=%g, newobj=%g\n",
            divedepth, maxdivedepth, heurdata->nlpiterations, maxnlpiterations,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var), oldobj, newobj);
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip, maxnlpiterations, &lperror) );
      if( lperror )
         break;

      /* update iteration count */
      heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;
      nlpiterations = SCIPgetNLPIterations(scip);

      /* get LP solution status  and fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new fractional variables */
         CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL) );
      }
      debugMessage("   -> lpsolstat=%d, nfrac=%d\n", lpsolstat, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && !lperror && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      Bool success;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      debugMessage("objpscostdiving found primal solution: obj=%g\n", SCIPgetSolOrigObj(scip, heurdata->sol));

      /* try to add solution to SCIP */
      CHECK_OKAY( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, &success) );

      /* check, if solution was feasible and good enough */
      if( success )
      {
         debugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   CHECK_OKAY( SCIPendDive(scip) );

   /* free temporary memory for remembering the current soft roundings */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &roundings) );

   debugMessage("objpscostdiving heuristic finished\n");

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the objpscostdiving heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurObjpscostdiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                  HEUR_MAXDEPTH, HEUR_PSEUDONODES,
                  heurFreeObjpscostdiving, heurInitObjpscostdiving, heurExitObjpscostdiving, heurExecObjpscostdiving,
                  heurdata) );

   /* objpscostdiving heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objpscostdiving/divestartdepth", 
                  "minimal relative depth to start diving",
                  &heurdata->divestartdepth, DEFAULT_DIVESTARTDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objpscostdiving/maxlpiterquot", 
                  "maximal fraction of diving LP iterations compared to total iteration number",
                  &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objpscostdiving/depthfac",
                  "maximal diving depth: number of binary/integer variables times depthfac",
                  &heurdata->depthfac, DEFAULT_DEPTHFAC, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objpscostdiving/depthfacnosol",
                  "maximal diving depth factor if no feasible solution was found yet",
                  &heurdata->depthfacnosol, DEFAULT_DEPTHFACNOSOL, 0.0, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

