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
#pragma ident "@(#) $Id: heur_objhistdiving.c,v 1.2 2004/02/04 17:27:25 bzfpfend Exp $"

/**@file   heur_objhistdiving.c
 * @brief  LP diving heuristic that changes variable's objective value instead of bounds, using history values as guide
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_objhistdiving.h"


#define HEUR_NAME         "objhistdiving"
#define HEUR_DESC         "LP diving heuristic that changes variable's objective values instead of bounds, using history"
#define HEUR_DISPCHAR     'o'
#define HEUR_PRIORITY     -1003000
#define HEUR_FREQ         10
#define HEUR_FREQOFS      6
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
   Longint          nlpiterations;      /**< LP iterations used in history diving */
};




/*
 * local methods
 */

static
void calcHistQuot(
   SCIP*            scip,               /**< SCIP data structure */
   VAR*             var,                /**< problem variable */
   Real             frac,               /**< fractionality of variable */
   int              rounddir,           /**< -1: round down, +1: round up, 0: select due to history values */
   Real*            histquot,           /**< pointer to store history quotient */
   Bool*            roundup             /**< pointer to store whether the variable should be rounded up */
   )
{
   Real histdown;
   Real histup;

   assert(histquot != NULL);
   assert(roundup != NULL);

   /* bound fractions to not prefer variables that are nearly integral */
   frac = MAX(frac, 0.1);
   frac = MIN(frac, 0.9);
   
   /* get history quotient */
   histdown = SCIPgetVarLPHistory(scip, var, 0.0-frac);
   histup = SCIPgetVarLPHistory(scip, var, 1.0-frac);
   assert(histdown >= 0.0 && histup >= 0.0);
   
   /* choose rounding direction */
   if( rounddir == -1 )
      *roundup = FALSE;
   else if( rounddir == +1 )
      *roundup = TRUE;
   else if( frac < 0.3 )
      *roundup = FALSE;
   else if( frac > 0.7 )
      *roundup = TRUE;
   else if( histdown < histup )
      *roundup = FALSE;
   else
      *roundup = TRUE;
   
   /* calculate history quotient */
   if( *roundup )
      *histquot = sqrt(frac) * (1.0+histdown) / (1.0+histup);
   else
      *histquot = sqrt(1.0-frac) * (1.0+histup) / (1.0+histdown);
   
   /* prefer decisions on binary variables */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      (*histquot) *= 1000.0;
}




/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
DECL_HEURFREE(heurFreeObjhistdiving) /*lint --e{715}*/
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


/** initialization method of primal heuristic (called when problem solving starts) */
static
DECL_HEURINIT(heurInitObjhistdiving) /*lint --e{715}*/
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


/** deinitialization method of primal heuristic (called when problem solving exits) */
static
DECL_HEUREXIT(heurExitObjhistdiving) /*lint --e{715}*/
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
DECL_HEUREXEC(heurExecObjhistdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;
   LPSOLSTAT lpsolstat;
   VAR* var;
   VAR** lpcands;
   Real* lpcandssol;
   Real* lpcandsfrac;
   Real frac;
   Real histquot;
   Real besthistquot;
   Real oldobj;
   Real newobj;
   Real objscale;
   Bool bestcandmayrounddown;
   Bool bestcandmayroundup;
   Bool bestcandroundup;
   Bool mayrounddown;
   Bool mayroundup;
   Bool roundup;
   Longint nlpiterations;
   Longint ndivinglpiterations;
   int* roundings;
   int nvars;
   int varidx;
   int nlpcands;
   int startnlpcands;
   int actdepth;
   int maxdepth;
   int maxdivedepth;
   int divedepth;
   int bestcand;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasActnodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* don't dive two times at the same node */
   if( SCIPgetLastDivenode(scip) == SCIPgetNodenum(scip) )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't try to dive, if we are in the higher fraction of the tree, given by divestartdepth */
   actdepth = SCIPgetActDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   if( actdepth < heurdata->divestartdepth*maxdepth )
      return SCIP_OKAY;

   /* don't try to dive, if we took too many LP iterations during diving */
   nlpiterations = SCIPgetNLPIterations(scip);
   if( heurdata->nlpiterations > heurdata->maxlpiterquot*nlpiterations )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

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
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );

   debugMessage("(node %lld) executing objhistdiving heuristic: depth=%d, %d fractionals, dualbound=%g, maxdivedepth=%d\n", 
      SCIPgetNodenum(scip), SCIPgetActDepth(scip), nlpcands, SCIPgetDualbound(scip), maxdivedepth);

   /* dive as long we are in the given diving depth limits and fractional variables exist, but
    * - if the last objective change was in a direction, that corresponds to a feasibile rounding, we continue in any case
    * - if possible, we dive at least with the depth 10
    * - if the number of fractional variables decreased at least with 1 variable per 2 dive depths, we continue diving
    */
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   startnlpcands = nlpcands;
   while( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (bestcandmayrounddown || bestcandmayroundup
         || divedepth < 10
         || nlpcands <= startnlpcands - divedepth/2
         || divedepth < maxdivedepth) )
   {
      divedepth++;

      /* choose variable for objective change:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, change objective value of variable with largest rel. difference of history values
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - change objective value of variable with largest rel. difference of history values
       */
      bestcand = -1;
      besthistquot = -1.0;
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
                * - if variable may be rounded in both directions, round corresponding to the history values
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               roundup = FALSE;
               if( mayrounddown && mayroundup )
                  calcHistQuot(scip, var, frac, 0, &histquot, &roundup);
               else if( mayrounddown )
                  calcHistQuot(scip, var, frac, +1, &histquot, &roundup);
               else
                  calcHistQuot(scip, var, frac, -1, &histquot, &roundup);

               /* prefer variables, that have already been softrounded but failed to get integral */
               varidx = SCIPvarGetProbindex(var);
               assert(0 <= varidx && varidx < nvars);
               if( roundings[varidx] != 0 )
                  histquot *= 1000.0;

               /* check, if candidate is new best candidate */
               if( histquot > besthistquot )
               {
                  bestcand = c;
                  besthistquot = histquot;
                  bestcandmayrounddown = mayrounddown;
                  bestcandmayroundup = mayroundup;
                  bestcandroundup = roundup;
               }
            }
         }
         else
         {
            /* the candidate may not be rounded: calculate history quotient and preferred direction */
            calcHistQuot(scip, var, frac, 0, &histquot, &roundup);

            /* prefer variables, that have already been softrounded but failed to get integral */
            varidx = SCIPvarGetProbindex(var);
            assert(0 <= varidx && varidx < nvars);
            if( roundings[varidx] != 0 )
               histquot *= 1000.0;

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if( bestcandmayrounddown || bestcandmayroundup || histquot > besthistquot )
            {
               bestcand = c;
               besthistquot = histquot;
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
            debugMessage("objhistdiving found roundable primal solution: obj=%g\n", SCIPgetSolObj(scip, heurdata->sol));
         
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
         debugMessage("  dive %d/%d: var <%s>, round=%d/%d, sol=%g, bounds=[%g,%g], obj=%g, newobj=%g\n",
            divedepth, maxdivedepth,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var), oldobj, newobj);
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip) );

      /* get LP solution status  and fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new fractional variables */
         CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );
      }
      debugMessage("   -> lpsolstat=%d, nfrac=%d\n", lpsolstat, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      Bool success;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      debugMessage("objhistdiving found primal solution: obj=%g\n", SCIPgetSolObj(scip, heurdata->sol));

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

   /* update iteration count */
   heurdata->nlpiterations += SCIPgetNLPIterations(scip) - nlpiterations;

   /* free temporary memory for remembering the current soft roundings */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &roundings) );

   debugMessage("objhistdiving heuristic finished\n");

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the objhistdiving heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurObjhistdiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                  HEUR_PSEUDONODES,
                  heurFreeObjhistdiving, heurInitObjhistdiving, heurExitObjhistdiving, heurExecObjhistdiving,
                  heurdata) );

   /* objhistdiving heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objhistdiving/divestartdepth", 
                  "minimal relative depth to start diving",
                  &heurdata->divestartdepth, DEFAULT_DIVESTARTDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objhistdiving/maxlpiterquot", 
                  "maximal fraction of diving LP iterations compared to total iteration number",
                  &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objhistdiving/depthfac",
                  "maximal diving depth: number of binary/integer variables times depthfac",
                  &heurdata->depthfac, DEFAULT_DEPTHFAC, 0.0, REAL_MAX, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/objhistdiving/depthfacnosol",
                  "maximal diving depth factor if no feasible solution was found yet",
                  &heurdata->depthfacnosol, DEFAULT_DEPTHFACNOSOL, 0.0, REAL_MAX, NULL, NULL) );
   
   return SCIP_OKAY;
}

