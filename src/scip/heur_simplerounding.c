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
#pragma ident "@(#) $Id: heur_simplerounding.c,v 1.13 2005/02/02 19:34:12 bzfpfend Exp $"

/**@file   heur_simplerounding.c
 * @brief  simple and fast LP rounding heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_simplerounding.h"


#define HEUR_NAME             "simplerounding"
#define HEUR_DESC             "simple and fast LP rounding heuristic"
#define HEUR_DISPCHAR         'r'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      FALSE     /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   TRUE      /* call heuristic during plunging? (should be FALSE for diving heuristics!) */


/* locally defined heuristic data */
struct HeurData
{
   SOL*             sol;                /**< working solution */
};




/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#define heurFreeSimplerounding NULL


/** initialization method of primal heuristic (called after problem was transformed) */
static
DECL_HEURINIT(heurInitSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(SCIPheurGetData(heur) == NULL);
   assert(scip != NULL);

   /* create heuristic data */
   CHECK_OKAY( SCIPallocMemory(scip, &heurdata) );
   CHECK_OKAY( SCIPcreateSol(scip, &heurdata->sol, heur) );
   SCIPheurSetData(heur, heurdata);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
DECL_HEUREXIT(heurExitSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   CHECK_OKAY( SCIPfreeSol(scip, &heurdata->sol) );
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecSimplerounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   HEURDATA* heurdata;
   SOL* sol;
   VAR** lpcands;
   Real* lpcandssol;
   int nlpcands;
   int c;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   CHECK_OKAY( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL) );

   /* only call heuristic, if LP solution is fractional */
   if( nlpcands == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   debugMessage("executing simple rounding heuristic: %d fractionals\n", nlpcands);

   /* get the working solution from heuristic's local data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   sol = heurdata->sol;
   assert(sol != NULL);

   /* copy the current LP solution to the working solution */
   CHECK_OKAY( SCIPlinkLPSol(scip, sol) );

   /* round all roundable fractional columns in the corresponding direction as long as no unroundable column was found */
   for( c = 0; c < nlpcands; ++c )
   {
      VAR* var;
      Real oldsolval;
      Real newsolval;
      Bool mayrounddown;
      Bool mayroundup;

      oldsolval = lpcandssol[c];
      assert(!SCIPisIntegral(scip, oldsolval));
      var = lpcands[c];
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);
      debugMessage("simple rounding heuristic: var <%s>, val=%g, rounddown=%d, roundup=%d\n",
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
      CHECK_OKAY( SCIPsetSolVal(scip, sol, var, newsolval) );
   }

   /* check, if rounding was successful */
   if( c == nlpcands )
   {
      Bool stored;

      /* check solution for feasibility, and add it to solution store if possible
       * neither integrality nor feasibility of LP rows has to be checked, because all fractional
       * variables were already moved in feasible direction to the next integer
       */
      CHECK_OKAY( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, &stored) );

      if( stored )
      {
#ifdef DEBUG
         printf("found feasible rounded solution:\n");
         SCIPprintSol(scip, sol, NULL);
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
RETCODE SCIPincludeHeurSimplerounding(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* include heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
                  HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING,
                  heurFreeSimplerounding, heurInitSimplerounding, heurExitSimplerounding, heurExecSimplerounding,
                  NULL) );

   return SCIP_OKAY;
}

