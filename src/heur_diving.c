/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_diving.c
 * @brief  LP diving heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_diving.h"


#define HEUR_NAME         "diving"
#define HEUR_DESC         "LP diving heuristic"
#define HEUR_DISPCHAR     'd'
#define HEUR_PRIORITY     -1000000
#define HEUR_FREQ         10
#define HEUR_PSEUDONODES  FALSE         /** call heuristic at nodes where only a pseudo solution exist? */



/*
 * Default parameter settings
 */

#define SCIP_DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (actlowerbound - lowerbound)/(upperbound - lowerbound)
                                              *   where diving is performed */
#define SCIP_DEFAULT_MAXDIVEAVGQUOT      4.0 /**< maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                              *   where diving is performed */
#define SCIP_DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet */
#define SCIP_DEFAULT_MAXDIVEAVGQUOTNOSOL 8.0 /**< maximal AVGQUOT when no solution was found yet */



/* locally defined heuristic data */
struct HeurData
{
   Real             maxdiveubquot;      /**< maximal quotient (actlowerbound - lowerbound)/(upperbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveavgquot;     /**< maximal quotient (actlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed */
   Real             maxdiveubquotnosol; /**< maximal UBQUOT when no solution was found yet */
   Real             maxdiveavgquotnosol;/**< maximal AVGQUOT when no solution was found yet */
};




/*
 * local methods
 */





/*
 * Callback methods
 */

static
DECL_HEURFREE(SCIPheurFreeDiving)
{
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
DECL_HEUREXEC(SCIPheurExecDiving)
{
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
   Real bestfrac;
   Real frac;
   Bool bestcandmayrounddown;
   Bool bestcandmayroundup;
   Bool bestcandroundup;
   Bool mayrounddown;
   Bool mayroundup;
   Bool roundup;
   int nlpcands;
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

   /* only try to dive, if we are in the lower 20% of the tree */
   actdepth = SCIPgetActDepth(scip);
   maxdepth = SCIPgetMaxDepth(scip);
   if( actdepth < 0.8*maxdepth )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* calculate the objective search bound */
   if( SCIPgetNSolsFound(scip) == 0 )
   {
      searchubbound = SCIPgetTransLowerBound(scip)
         + heurdata->maxdiveubquotnosol * (SCIPgetTransUpperBound(scip) - SCIPgetTransLowerBound(scip));
      searchavgbound = SCIPgetTransLowerBound(scip)
         + heurdata->maxdiveavgquotnosol * (SCIPgetAvgTransLowerBound(scip) - SCIPgetTransLowerBound(scip));
   }
   else
   {
      searchubbound = SCIPgetTransLowerBound(scip)
         + heurdata->maxdiveubquot * (SCIPgetTransUpperBound(scip) - SCIPgetTransLowerBound(scip));
      searchavgbound = SCIPgetTransLowerBound(scip)
         + heurdata->maxdiveavgquot * (SCIPgetAvgTransLowerBound(scip) - SCIPgetTransLowerBound(scip));
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

   debugMessage("executing diving heuristic: depth=%d, %d fractionals, dualbound=%g, searchbound=%g\n", 
      SCIPgetActDepth(scip), nlpcands, SCIPgetDualBound(scip), SCIPretransformObj(scip, searchbound));

   /* dive as long we are in the given objective limits and fractional variables exist
    * if the last rounding was in a direction, that never destroys feasibility, we continue in any case
    */
   divedepth = 0;
   bestcandmayrounddown = FALSE;
   bestcandmayroundup = FALSE;
   while( lpsolstat == SCIP_LPSOLSTAT_OPTIMAL && nlpcands > 0
      && (bestcandmayrounddown || bestcandmayroundup || (divedepth < maxdivedepth && objval < searchbound)) )
   {
      divedepth++;

      todoMessage("use a better variable selection/rounding criteria in diving (e.g. history depending)");

      /* choose variable fixing:
       * - prefer variables that may not be rounded without destroying LP feasibility:
       *   - of these variables, round least fractional variable in corresponding direction
       * - if all remaining fractional variables may be rounded without destroying LP feasibility:
       *   - round variable in the feasible direction, that is closest to the rounded value
       */
      bestcand = -1;
      bestcandmayrounddown = TRUE;
      bestcandmayroundup = TRUE;
      bestfrac = 1.0;
      for( c = 0; c < nlpcands; ++c )
      {
         var = lpcands[c];
         mayrounddown = SCIPvarMayRoundDown(var);
         mayroundup = SCIPvarMayRoundUp(var);
         if( mayrounddown || mayroundup )
         {
            /* the candidate may be rounded */
            if( bestcandmayrounddown || bestcandmayroundup )
            {
               frac = lpcandsfrac[c];
               roundup = FALSE;
               if( (mayrounddown && mayroundup) || divedepth < 10 )
               {
                  if( frac > 0.5 )
                  {
                     frac = MIN(frac, 1.0-frac);
                     roundup = TRUE;
                  }
               }
               else if( mayroundup )
               {
                  frac = 1.0 - frac;
                  roundup = TRUE;
               }
               if( frac < bestfrac )
               {
                  bestcand = c;
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
            frac = lpcandsfrac[c];
            roundup = FALSE;
            if( frac > 0.5 )
            {
               frac = 1.0-frac;
               roundup = TRUE;
            }
            if( bestcandmayrounddown || bestcandmayroundup || frac < bestfrac )
            {
               bestcand = c;
               bestfrac = frac;
               bestcandmayrounddown = FALSE;
               bestcandmayroundup = FALSE;
               bestcandroundup = roundup;
            }
         }
      }
      assert(bestcand != -1);

      var = lpcands[bestcand];

      if( SCIPvarGetLb(var) >= SCIPvarGetUb(var) - 0.5 )
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
         objval = SCIPgetLPObjval(scip);
         CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands) );
      }
      debugMessage("   -> lpsolstat=%d, objval=%g, nfrac=%d\n", lpsolstat, objval, nlpcands);
   }

   /* check if a solution has been found */
   if( nlpcands == 0 && lpsolstat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      SOL* sol;
      Bool stored;

      /* create solution from diving LP */
      CHECK_OKAY( SCIPcreateLPSol(scip, &sol, heur) );
      debugMessage("diving found primal solution: obj=%g\n", SCIPgetSolObj(scip, sol));

      /* try to add solution to SCIP */
      CHECK_OKAY( SCIPtrySolFree(scip, &sol, FALSE, FALSE, &stored) );

      /* check, if solution was feasible and good enough */
      if( stored )
      {
         debugMessage(" -> solution was feasible and good enough\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   /* end diving */
   CHECK_OKAY( SCIPendDive(scip) );

   debugMessage("diving heuristic finished\n");

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the diving heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurDiving(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* allocate and initialise heuristic data; this has to be freed in the destructor */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );
   heurdata->maxdiveubquot = SCIP_DEFAULT_MAXDIVEUBQUOT;
   heurdata->maxdiveavgquot = SCIP_DEFAULT_MAXDIVEAVGQUOT;
   heurdata->maxdiveubquotnosol = SCIP_DEFAULT_MAXDIVEUBQUOTNOSOL;
   heurdata->maxdiveavgquotnosol = SCIP_DEFAULT_MAXDIVEAVGQUOTNOSOL;

   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_PSEUDONODES,
                  SCIPheurFreeDiving, NULL, NULL, SCIPheurExecDiving,
                  heurdata) );

   return SCIP_OKAY;
}

