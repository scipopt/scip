/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_indicator.c
 * @brief  handle partial solutions for linear problems with indicators and otherwise continuous variables
 * @author Marc Pfetsch
 *
 * For linear problems with indicators and otherwise continuous variables, the indicator constraint handler can produce
 * partial solutions, i.e., values for the indicator variables. This partial solution can be passed to this heuristic,
 * which then fixes these values and solves an LP. Additionally a local search for a better solution is added.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/scip.h"
#include "scip/heur_indicator.h"
#include "scip/cons_indicator.h"

#define HEUR_NAME             "indicator"
#define HEUR_DESC             "indicator heuristic to create feasible solutions from values for indicator variables"
#define HEUR_DISPCHAR         'A'
#define HEUR_PRIORITY         -20200
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_ONEOPT        TRUE           /**< whether the one-opt heuristic should be started */


/** primal heuristic data */
struct SCIP_HeurData
{
   int                   nindconss;          /**< number of indicator constraints */
   SCIP_CONS**           indconss;           /**< indicator constraints */
   SCIP_Bool*            solcand;            /**< bitset of indicator variables ind solution candidate */
   SCIP_Bool             oneopt;             /**< whether the one-opt heuristic should be started */
   SCIP_SOL*             lastsol;            /**< last solution considered for improvement */
};

/*
 * Local methods
 */

/** try one-opt on given solution */
static
SCIP_RETCODE tryOneOpt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< indicator heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   SCIP_Bool*            solcand,            /**< values for indicator variables in partial solution */
   SCIP_Bool*            foundsol            /**< whether a solution has been found */
   )
{
   SCIP_Bool cutoff;
   SCIP_Bool lperror;
   SCIP_Bool stored;
   SCIP_SOL* sol;
   int cnt = 0;
   int i;
   int c;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( heurdata != NULL );
   assert( nindconss == 0 || indconss != NULL );
   assert( solcand != NULL );
   assert( foundsol != NULL );

   SCIPdebugMessage("Performing one-opt ...\n");
   *foundsol = FALSE;

   SCIP_CALL( SCIPstartProbing(scip) );

   for (i = 0; i < nindconss; ++i)
   {
      SCIP_VAR* binvar;

      assert( SCIPconsIsActive(indconss[i]) );
      binvar = SCIPgetBinaryVarIndicator(indconss[i]);
      assert( binvar != NULL );

      /* skip constraints with fixed variables */
      if ( SCIPvarGetLbLocal(binvar) < 0.5 && SCIPvarGetUbLocal(binvar) > 0.5 )
         continue;

      /* get rid of all bound changes */
      SCIP_CALL( SCIPbacktrackProbing(scip, 0) );
      SCIP_CALL( SCIPnewProbingNode(scip) );
      ++cnt;

      /* fix variables */
      for (c = 0; c < nindconss; ++c)
      {
         SCIP_Bool s;

         assert( SCIPconsIsActive(indconss[c]) );
         binvar = SCIPgetBinaryVarIndicator(indconss[c]);
         assert( binvar != NULL );

         /* fix variables according to solution candidate, except constraint i */
         if ( c == i )
            s = ! solcand[c];
         else
            s = solcand[c];

         if ( ! s )
         {
            if ( SCIPvarGetLbLocal(binvar) < 0.5 && SCIPvarGetUbLocal(binvar) > 0.5 )
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, binvar, 1.0) );
            }
         }
         else
         {
            if ( SCIPvarGetUbLocal(binvar) > 0.5 && SCIPvarGetLbLocal(binvar) < 0.5 )
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, binvar, 0.0) );
            }
         }
      }

      /* propagate variables */
      SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );
      if ( cutoff )
      {
         SCIPdebugMessage("Hit cutoff in propagation.\n");
         continue;
      }

      /* solve LP to move continuous variables */
      SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );

      /* the LP often reaches the objective limit - we currently do not use such solutions */
      if ( lperror || cutoff || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
#ifdef SCIP_DEBUG
         if ( lperror )
            SCIPdebugMessage("An LP error occured.\n");
         else
            SCIPdebugMessage("Hit cutoff in LP solving.\n");
#endif
         continue;
      }

      /* create solution */
      SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

      /* copy the current LP solution to the working solution */
      SCIP_CALL( SCIPlinkLPSol(scip, sol) );

      /* check solution for feasibility */
      SCIPdebugMessage("Found solution candidate with value %g.\n", SCIPgetSolTransObj(scip, sol));

      /* only check integrality, because we solved an LP */
      SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, FALSE, &stored) );
      if ( stored )
         *foundsol = TRUE;
   }
   SCIP_CALL( SCIPendProbing(scip) );

   SCIPdebugMessage("Finished one-opt (tried variables: %d).\n", cnt);

   return SCIP_OKAY;
}


/** try given solution */
static
SCIP_RETCODE trySolCandidate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< indicator heuristic */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   SCIP_Bool*            solcand,            /**< values for indicator variables in partial solution */
   SCIP_Bool*            foundsol            /**< whether a solution has been found */
   )
{
   SCIP_Bool cutoff;
   SCIP_Bool lperror;
   SCIP_Bool stored;
   SCIP_SOL* sol;
   int c;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( heurdata != NULL );
   assert( nindconss == 0 || indconss != NULL );
   assert( solcand != NULL );
   assert( foundsol != NULL );

   SCIPdebugMessage("Trying to generate feasible solution with indicators from solution candidate ...\n");
   *foundsol = FALSE;

   SCIP_CALL( SCIPstartProbing(scip) );
   SCIP_CALL( SCIPnewProbingNode(scip) );

   /* fix variables */
   for (c = 0; c < nindconss; ++c)
   {
      SCIP_VAR* binvar;

      assert( SCIPconsIsActive(indconss[c]) );
      binvar = SCIPgetBinaryVarIndicator(indconss[c]);
      assert( binvar != NULL );

      /* Fix binary variables not in cover to 1 and corresponding slack variables to 0. The other binary variables are fixed to 0. */
      if ( ! solcand[c] )
      {
         /* to be sure, check for non-fixed variables */
         if ( SCIPvarGetLbLocal(binvar) < 0.5 && SCIPvarGetUbLocal(binvar) > 0.5 )
         {
            SCIP_CALL( SCIPchgVarLbProbing(scip, binvar, 1.0) );
         }
      }
      else
      {
         if ( SCIPvarGetUbLocal(binvar) > 0.5 && SCIPvarGetLbLocal(binvar) < 0.5 )
         {
            SCIP_CALL( SCIPchgVarUbProbing(scip, binvar, 0.0) );
         }
      }
   }

   /* propagate variables */
   SCIP_CALL( SCIPpropagateProbing(scip, -1, &cutoff, NULL) );
   if ( cutoff )
   {
      SCIPdebugMessage("Hit cutoff in propagation.\n");
      SCIP_CALL( SCIPendProbing(scip) );
      return SCIP_OKAY;
   }

   /* solve LP to move continuous variables */
   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );

   /* the LP often reaches the objective limit - we currently do not use such solutions */
   if ( lperror || cutoff || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
   {
#ifdef SCIP_DEBUG
      if ( lperror )
         SCIPdebugMessage("An LP error occured.\n");
      else
         SCIPdebugMessage("Hit cutoff in LP solving.\n");
#endif
      SCIP_CALL( SCIPendProbing(scip) );
      return SCIP_OKAY;
   }

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, heur) );

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* check solution for feasibility */
#ifdef SCIP_DEBUG
   SCIPdebugMessage("Found solution candidate with value %g.\n", SCIPgetSolTransObj(scip, sol));
#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
   SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, TRUE, TRUE, &stored) );
   if ( stored )
   {
      *foundsol = TRUE;
      SCIPdebugMessage("Solution is feasible and stored.\n");
   }
   else
      SCIPdebugMessage("Solution was not stored.\n");
#else
   /* only check integrality, because we solved an LP */
   SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, FALSE, TRUE, FALSE, &stored) );
   if ( stored )
      *foundsol = TRUE;
#endif
   SCIP_CALL( SCIPendProbing(scip) );

   /* possibly perform one-opt */
   if ( stored && heurdata->oneopt )
   {
      SCIP_Bool found = FALSE;
      assert( *foundsol );
      SCIP_CALL( tryOneOpt(scip, heur, heurdata, nindconss, indconss, solcand, &found) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIndicator)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIndicator(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIndicator)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   SCIPfreeBlockMemoryArrayNull(scip, &(heurdata->indconss), heurdata->nindconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(heurdata->solcand), heurdata->nindconss);

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIndicator)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* call heuristic, if solution candidate is available */
   if ( heurdata->solcand != NULL )
   {
      SCIP_Bool foundsol = FALSE;

      assert( heurdata->nindconss > 0 );
      assert( heurdata->indconss != NULL );

      /* The heuristic will only be successful if there are no integral variables and no binary variables except the
       * indicator variables. */
      if ( SCIPgetNIntVars(scip) > 0 || heurdata->nindconss < SCIPgetNBinVars(scip) )
         return SCIP_OKAY;

      SCIP_CALL( trySolCandidate(scip, heur, heurdata, heurdata->nindconss, heurdata->indconss, heurdata->solcand, &foundsol) );

      if ( foundsol )
         *result = SCIP_FOUNDSOL;
      else
         *result = SCIP_DIDNOTFIND;

      /* free memory */
      SCIPfreeBlockMemoryArray(scip, &(heurdata->solcand), heurdata->nindconss);
      SCIPfreeBlockMemoryArray(scip, &(heurdata->indconss), heurdata->nindconss);
   }

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the indicator primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIndicator(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Indicator primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->nindconss = 0;
   heurdata->indconss = NULL;
   heurdata->solcand = NULL;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecIndicator, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyIndicator) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeIndicator) );

   /* add parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heur/"HEUR_NAME"/oneopt",
         "whether the one-opt heuristic should be started",
         &heurdata->oneopt, TRUE, DEFAULT_ONEOPT, NULL, NULL) );

   return SCIP_OKAY;
}


/** pass partial solution for indicator variables to heuristic */
SCIP_RETCODE SCIPheurPassIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< indicator heuristic */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   SCIP_Bool*            solcand             /**< values for indicator variables in partial solution */
   )
{
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( nindconss > 0 );
   assert( indconss != NULL );
   assert( solcand != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* copy indicator information */
   if ( heurdata->indconss != NULL )
      SCIPfreeBlockMemoryArray(scip, &(heurdata->indconss), heurdata->nindconss);

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(heurdata->indconss), indconss, nindconss) );
   heurdata->nindconss = nindconss;

   /* copy partial solution */
   if ( heurdata->solcand != NULL )
      BMScopyMemoryArray(heurdata->solcand, solcand, nindconss);
   else
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(heurdata->solcand), solcand, nindconss) );

   return SCIP_OKAY;
}
