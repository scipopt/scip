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

/**@file   heur_indrounding.c
 * @brief  rounding heuristic that takes care of indicator constraints
 * @author Marc Pfetsch
 *
 * This is a simple rounding algorithm to produce feasible solutions: we round the current variables
 * and set the indicator variables such that we get a feasible solution (if possible).
 *
 * This is a simple idea implementing a 2-approximation for MaxFS.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include <scip/cons_indicator.h>

#include "heur_indrounding.h"


#define HEUR_NAME             "indrounding"
#define HEUR_DESC             "indicator rounding heuristic"
#define HEUR_DISPCHAR         '5'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      FALSE




/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
};



/*
 * Callback methods
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIndrounding)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitIndrounding)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitIndrounding)
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
#define heurInitsolIndrounding NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolIndrounding NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIndrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_CONSHDLR* indconshdlr;
   SCIP_HEURDATA* heurdata;
   SCIP_Bool* isIndicatorVar;
   SCIP_VAR** vars;
   SCIP_SOL* sol;
   SCIP_CONS** conss;
   SCIP_Bool stored;
   SCIP_Bool changed;
   int nconss;
   int nvars;
   int c;
   int v;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;
   assert( SCIPhasCurrentNodeLP(scip) );

   SCIPdebugMessage("Running indicator rounding heuristic ...\n");
   *result = SCIP_DIDNOTFIND;

   /* get variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   assert( vars != NULL );

   /* get indicator constraints */
   indconshdlr = SCIPfindConshdlr(scip, "indicator");
   assert( indconshdlr != NULL );

   conss = SCIPconshdlrGetConss(indconshdlr);
   nconss = SCIPconshdlrGetNConss(indconshdlr);

   /* find indicator variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &isIndicatorVar, nvars) );
   for (v = 0; v < nvars; ++v)
      isIndicatorVar[v] = FALSE;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_VAR* var;
      int idx;

      var = SCIPgetBinaryVarIndicator(conss[c]);
      idx = SCIPvarGetProbindex(var);

      if ( idx >= 0 )
      {
	 assert( idx < nvars );
	 isIndicatorVar[idx] = TRUE;
      }

      var = SCIPgetSlackVarIndicator(conss[c]);
      idx = SCIPvarGetProbindex(var);

      if ( idx >= 0 )
      {
	 assert( idx < nvars );
	 isIndicatorVar[idx] = TRUE;
      }
   }

   /* now round variables other than indicator variables */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   sol = heurdata->sol;
   assert( sol != NULL );
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );
   assert( sol != NULL );

   /* round all variable except indicator variables */
   for (v = 0; v < nvars; ++v)
   {
      assert( SCIPvarGetProbindex(vars[v]) == v );
      if ( ! isIndicatorVar[v] )
      {
	 SCIP_VAR* var;
	 SCIP_Real val;

	 var = vars[v];
	 val = SCIPgetSolVal(scip, sol, var);
	 assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );

	 if ( (SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER)
            && (! SCIPisFeasIntegral(scip, val)) )
	 {
	    SCIP_Real newsolval;
	    SCIP_Real obj;

	    /* if objective coefficient is 0 round to nearest integer */
	    obj = SCIPvarGetObj(var);

	    if ( SCIPisZero(scip, obj) )
	    {
               newsolval = SCIPround(scip, val);
	    }
	    else
	    {
	       if ( obj > 0.0 )
		  newsolval = SCIPfeasFloor(scip, val);
	       else
		  newsolval = SCIPfeasCeil(scip, val);
	    }

	    /* SCIPdebugMessage("var <%s>, val=%g, newval=%g\n", SCIPvarGetName(var), val, newsolval); */

	    /* store new solution value */
	    SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
	 }
      }
   }

   /* set indicator variables */
   SCIP_CALL( SCIPmakeIndicatorsFeasible(scip, indconshdlr, sol, &changed) );

   /* check solution for feasibility, and add it to solution store if possible
    * neither integrality nor feasibility of LP rows has to be checked, because all fractional
    * variables were already moved in feasible direction to the next integer
    */
   SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, TRUE, TRUE, TRUE, &stored) );

   /* SCIP_CALL( SCIPcheckSolOrig(scip, sol, &stored, TRUE, FALSE) ); */

   if ( stored )
   {
      SCIPdebugMessage("found feasible rounded solution of value %f.\n", SCIPgetSolOrigObj(scip, sol));
#ifdef SCIP_DEBUG
      /* SCIPprintSol(scip, sol, NULL, FALSE); */
#endif
      *result = SCIP_FOUNDSOL;
   }

   SCIPfreeBufferArray(scip, &isIndicatorVar);

   return SCIP_OKAY;
}


/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIndrounding)
{
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIndrounding(scip) );

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the indrounding heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIndrounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyIndrounding, heurFreeIndrounding, heurInitIndrounding, heurExitIndrounding,
         heurInitsolIndrounding, heurExitsolIndrounding, heurExecIndrounding,
         heurdata) );

   return SCIP_OKAY;
}
