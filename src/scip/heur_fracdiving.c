/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_fracdiving.c,v 1.1 2004/01/24 17:21:10 bzfpfend Exp $"

/**@file   heur_fracdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_fracdiving.h"


#define HEUR_NAME         "fracdiving"
#define HEUR_DESC         "LP diving heuristic that chooses fixings w.r.t. the fractionalities"
#define HEUR_DISPCHAR     'f'
#define HEUR_PRIORITY     -1000000
#define HEUR_FREQ         10
#define HEUR_PSEUDONODES  FALSE         /** call heuristic at nodes where only a pseudo solution exist? */



/*
 * Default parameter settings
 */

#define DEFAULT_DIVESTARTDEPTH      0.5 /**< minimal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT       0.1 /**< maximal fraction of diving LP iterations compared to total iteration number */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (actlowerbound - lowerbound)/(upperbound - lowerbound)
                                         *   where diving is performed */
#define DEFAULT_MAXDIVEAVGQUOT      4.0 /**< maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 8.0 /**< maximal AVGQUOT when no solution was found yet */



/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
   Real             divestartdepth;     /**< minimal relative depth to start diving */
   Real             maxlpiterquot;      /**< maximal fraction of diving LP iterations compared to total iteration number */
   Real             maxdiveubquot;      /**< maximal quotient (actlowerbound - lowerbound)/(upperbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveavgquot;     /**< maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet */
   Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet */
   Longint          nlpiterations;      /**< LP iterations used in fractional diving */
};




/*
 * local methods
 */





/*
 * Callback methods
 */

static
DECL_HEURFREE(SCIPheurFreeFracdiving) /*lint --e{715}*/
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

static
DECL_HEURINIT(SCIPheurInitFracdiving) /*lint --e{715}*/
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

static
DECL_HEUREXIT(SCIPheurExitFracdiving) /*lint --e{715}*/
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

static
DECL_HEUREXEC(SCIPheurExecFracdiving) /*lint --e{715}*/
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
   Real obj;
   Real objgain;
   Real bestobjgain;
   Real frac;
   Real bestfrac;
   Bool bestcandmayrounddown;
   Bool bestcandmayroundup;
   Bool bestcandroundup;
   Bool mayrounddown;
   Bool mayroundup;
   Bool roundup;
   Longint nlpiterations;
   Longint ndivinglpiterations;
   Longint nsolsfound;
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

   /* calculate the objective search bound */
   nsolsfound = SCIPgetNSolsFound(scip);
   if( nsolsfound == 0 )
   {
      searchubbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveubquotnosol * (SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip));
      searchavgbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveavgquotnosol * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   }
   else
   {
      searchubbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveubquot * (SCIPgetUpperbound(scip) - SCIPgetLowerbound(scip));
      searchavgbound = SCIPgetLowerbound(scip)
         + heurdata->maxdiveavgquot * (SCIPgetAvgLowerbound(scip) - SCIPgetLowerbound(scip));
   }
   searchbound = MIN(searchubbound, searchavgbound);

   /* calculate the maximal diving depth: 10 * min{number of integer variables, max depth} */
   maxdivedepth = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
   maxdivedepth = MIN(maxdivedepth, maxdepth);
   maxdivedepth *= 10;

   /* start diving */
   CHECK_OKAY( SCIPstartDive(scip) );

   /* get LP objective value, and fractional variables, that should be integral */
   lpsolstat = SCIP_LPSOLSTAT_OPTIMAL;
   objval = SCIPgetLPObjval(scip);
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );

   debugMessage("(node %lld) executing fracdiving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n", 
      SCIPgetNodenum(scip), SCIPgetActDepth(scip), nlpcands, SCIPgetDualbound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective limits and fractional variables exist, but
    * - if the last rounding was in a direction, that never destroys feasibility, we continue in any case
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
         || (divedepth < maxdivedepth && objval < searchbound)) )
   {
      divedepth++;

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round least fractional variable in corresponding direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable with least increasing objective value
       */
      bestcand = -1;
      bestobjgain = SCIPinfinity(scip);
      bestfrac = SCIP_INVALID;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestcandroundup = FALSE;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         frac = lpcandsfrac[c];
         obj = SCIPvarGetObj(var);
         if( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded: choose this candidate only, if the best candidate may also be rounded */
            if( bestcandmayrounddown || bestcandmayroundup )
            {
               /* choose rounding direction:
                * - if variable may be rounded in both directions, round corresponding to the fractionality
                * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
                *   the current fractional solution
                */
               if( mayrounddown && mayroundup )
                  roundup = (frac > 0.5);
               else
                  roundup = mayrounddown;

               if( roundup )
               {
                  frac = 1.0 - frac;
                  objgain = frac*obj;
               }
               else
                  objgain = -frac*obj;

               /* penalize too small fractions */
               if( frac < 0.01 )
                  objgain *= 1000.0;
            
               /* prefer decisions on binary variables */
               if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
                  objgain *= 1000.0;

               /* check, if candidate is new best candidate */
               if( SCIPisLT(scip, objgain, bestobjgain) || (SCIPisEQ(scip, objgain, bestobjgain) && frac < bestfrac) )
               {
                  bestcand = c;
                  bestobjgain = objgain;
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
            if( frac > 0.5 )
            {
               roundup = TRUE;
               frac = 1.0 - frac;
            }
            else
               roundup = FALSE;

            /* penalize too small fractions */
            if( frac < 0.01 )
               frac += 10.0;
            
            /* prefer decisions on binary variables */
            if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
               frac *= 1000.0;

            /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
            if( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
            {
               bestcand = c;
               bestfrac = frac;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
            assert(bestfrac < SCIP_INVALID);
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
            debugMessage("fracdiving found roundable primal solution: obj=%g\n", SCIPgetSolObj(scip, heurdata->sol));
         
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

      if( SCIPgetVarLbDive(scip, var) >= SCIPgetVarUbDive(scip, var) - 0.5 )
      {
         /* the variable is already fixed -> numerical troubles -> abort diving */
         break;
      }

      /* apply rounding of best candidate */
      if( bestcandroundup )
      {
         /* round variable up */
         debugMessage("  dive %d/%d: var <%s>, round=%d/%d, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
            divedepth, maxdivedepth,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var),
            SCIPceil(scip, lpcandssol[bestcand]), SCIPgetVarUbDive(scip, var));
         CHECK_OKAY( SCIPchgVarLbDive(scip, var, SCIPceil(scip, lpcandssol[bestcand])) );
      }
      else
      {
         /* round variable down */
         debugMessage("  dive %d/%d: var <%s>, round=%d/%d, sol=%g, oldbounds=[%g,%g], newbounds=[%g,%g]\n",
            divedepth, maxdivedepth,
            SCIPvarGetName(var), bestcandmayrounddown, bestcandmayroundup,
            lpcandssol[bestcand], SCIPgetVarLbDive(scip, var), SCIPgetVarUbDive(scip, var),
            SCIPgetVarLbDive(scip, var), SCIPfloor(scip, lpcandssol[bestcand]));
         CHECK_OKAY( SCIPchgVarUbDive(scip, lpcands[bestcand], SCIPfloor(scip, lpcandssol[bestcand])) );
      }

      /* resolve the diving LP */
      CHECK_OKAY( SCIPsolveDiveLP(scip) );

      /* get LP solution status, objective value, and fractional variables, that should be integral */
      lpsolstat = SCIPgetLPSolstat(scip);
      if( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
      {
         /* get new objective value */
         oldobjval = objval;
         objval = SCIPgetLPObjval(scip);

         /* update history values */
         if( SCIPisGT(scip, objval, oldobjval) )
         {
            if( bestcandroundup )
            {
               CHECK_OKAY( SCIPupdateVarLPHistory(scip, lpcands[bestcand], 1.0-lpcandsfrac[bestcand], 
                              objval - oldobjval, 1.0) );
            }
            else
            {
               CHECK_OKAY( SCIPupdateVarLPHistory(scip, lpcands[bestcand], 0.0-lpcandsfrac[bestcand], 
                              objval - oldobjval, 1.0) );
            }
         }

         /* get new fractional variables */
         CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );
      }
      debugMessage("   -> lpsolstat=%d, objval=%g, nfrac=%d\n", lpsolstat, objval, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      Bool success;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPlinkLPSol(scip, heurdata->sol) );
      debugMessage("fracdiving found primal solution: obj=%g\n", SCIPgetSolObj(scip, heurdata->sol));

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

   debugMessage("fracdiving heuristic finished\n");

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the fracdiving heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurFracdiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_PSEUDONODES,
                  SCIPheurFreeFracdiving, SCIPheurInitFracdiving, SCIPheurExitFracdiving, SCIPheurExecFracdiving,
                  heurdata) );

   /* fracdiving heuristic parameters */
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/fracdiving/divestartdepth", 
                  "minimal relative depth to start diving",
                  &heurdata->divestartdepth, DEFAULT_DIVESTARTDEPTH, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/fracdiving/maxlpiterquot", 
                  "maximal fraction of diving LP iterations compared to total iteration number",
                  &heurdata->maxlpiterquot, DEFAULT_MAXLPITERQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/fracdiving/maxdiveubquot",
                  "maximal quotient (actlowerbound - lowerbound)/(upperbound - lowerbound) where diving is performed",
                  &heurdata->maxdiveubquot, DEFAULT_MAXDIVEUBQUOT, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/fracdiving/maxdiveavgquot", 
                  "maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound) where diving is performed",
                  &heurdata->maxdiveavgquot, DEFAULT_MAXDIVEAVGQUOT, 0.0, SCIP_INVALID, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/fracdiving/maxdiveubquotnosol", 
                  "maximal UBQUOT when no solution was found yet",
                  &heurdata->maxdiveubquotnosol, DEFAULT_MAXDIVEUBQUOTNOSOL, 0.0, 1.0, NULL, NULL) );
   CHECK_OKAY( SCIPaddRealParam(scip,
                  "heuristics/fracdiving/maxdiveavgquotnosol", 
                  "maximal AVGQUOT when no solution was found yet",
                  &heurdata->maxdiveavgquotnosol, DEFAULT_MAXDIVEAVGQUOTNOSOL, 0.0, SCIP_INVALID, NULL, NULL) );
   
   return SCIP_OKAY;
}

