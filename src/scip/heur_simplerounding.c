/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_simplerounding.c
 * @ingroup PRIMALHEURISTICS
 * @brief  simple and fast LP rounding heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_simplerounding.h"


#define HEUR_NAME             "simplerounding"
#define HEUR_DESC             "simple and fast LP rounding heuristic"
#define HEUR_DISPCHAR         'r'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Longint          lastlp;             /**< last LP number where the heuristic was applied */
   int                   nroundablevars;     /**< number of variables that can be rounded (-1 if not yet calculated) */
};




/*
 * Callback methods
 */
/** copy method for primal heuristic plugins (called when SCIP copies plugins) */

static
SCIP_DECL_HEURCOPY(heurCopySimplerounding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSimplerounding(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#define heurFreeSimplerounding NULL


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(SCIPheurGetData(heur) == NULL);

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->lastlp = -1;
   heurdata->nroundablevars = -1;
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolSimplerounding)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->lastlp = -1;

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolSimplerounding NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_SOL* sol;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Longint nlps;
   int nlpcands;
   int c;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* on our first call or after each pricing round, calculate the number of roundable variables */
   if( heurdata->nroundablevars == -1  || heurtiming == SCIP_HEURTIMING_DURINGPRICINGLOOP )
   {
      SCIP_VAR** vars;
      int nvars;
      int nroundablevars;
      int i;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);
      nroundablevars = 0;
      for( i = 0; i < nvars; ++i )
      {
         if( SCIPvarMayRoundDown(vars[i]) || SCIPvarMayRoundUp(vars[i]) )
            nroundablevars++;
      }
      heurdata->nroundablevars = nroundablevars;
   }

   /* don't call heuristic if there are no roundable variables; except we are called during pricing, in this case we
    * want to detect a (mixed) integer (LP) solution which is primal feasible */
   if( heurdata->nroundablevars == 0 && heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP )
      return SCIP_OKAY;

   /* don't call heuristic, if we have already processed the current LP solution */
   nlps = SCIPgetNLPs(scip);
   if( nlps == heurdata->lastlp )
      return SCIP_OKAY;
   heurdata->lastlp = nlps;

   /* get fractional variables, that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL) );

   /* only call heuristic, if LP solution is fractional; except we are called during pricing, in this case we
    * want to detect a (mixed) integer (LP) solution which is primal feasible */
   if( nlpcands == 0  && heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP )
      return SCIP_OKAY;

   /* don't call heuristic, if there are more fractional variables than roundable ones */
   if( nlpcands > heurdata->nroundablevars )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("executing simple rounding heuristic: %d fractionals\n", nlpcands);

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert(sol != NULL);

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* round all roundable fractional columns in the corresponding direction as long as no unroundable column was found */
   for( c = 0; c < nlpcands; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real oldsolval;
      SCIP_Real newsolval;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;

      oldsolval = lpcandssol[c];
      assert(!SCIPisFeasIntegral(scip, oldsolval));
      var = lpcands[c];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      SCIPdebugMessage("simple rounding heuristic: var <%s>, val=%g, rounddown=%u, roundup=%u\n",
         SCIPvarGetName(var), oldsolval, mayrounddown, mayroundup);

      /* choose rounding direction */
      if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if( SCIPvarGetObj(var) >= 0.0 )
            newsolval = SCIPfeasFloor(scip, oldsolval);
         else
            newsolval = SCIPfeasCeil(scip, oldsolval);
      }
      else if( mayrounddown )
         newsolval = SCIPfeasFloor(scip, oldsolval);
      else if( mayroundup )
         newsolval = SCIPfeasCeil(scip, oldsolval);
      else
         break;

      /* store new solution value */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
   }

   /* check, if rounding was successful */
   if( c == nlpcands )
   {
      SCIP_Bool stored;

      /* check solution for feasibility, and add it to solution store if possible
       * neither integrality nor feasibility of LP rows has to be checked, because all fractional
       * variables were already moved in feasible direction to the next integer
       */
      SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, FALSE, &stored) );

      if( stored )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("found feasible rounded solution:\n");
         SCIPprintSol(scip, sol, NULL, FALSE);
#endif
         *result = SCIP_FOUNDSOL;
      }
   }

   return SCIP_OKAY;
}




/*
 * heuristic specific interface methods
 */

/** creates the simple rounding heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSimplerounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopySimplerounding,
         heurFreeSimplerounding, heurInitSimplerounding, heurExitSimplerounding,
         heurInitsolSimplerounding, heurExitsolSimplerounding, heurExecSimplerounding,
         NULL) );

   return SCIP_OKAY;
}

