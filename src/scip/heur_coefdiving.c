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

/**@file   heur_coefdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the matrix coefficients
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 *
 * Indicator constraints are taken into account if present.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_coefdiving.h"
#include "scip/cons_indicator.h"


#define HEUR_NAME             "coefdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the matrix coefficients"
#define HEUR_DISPCHAR         'c'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          1
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_MINRELDEPTH         0.0 /**< minimal relative depth to start diving */
#define DEFAULT_MAXRELDEPTH         1.0 /**< maximal relative depth to start diving */
#define DEFAULT_MAXLPITERQUOT      0.05 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_CONSHDLR*        indconshdlr;        /**< indicator constraint handler (or NULL) */
   SCIP_DIVESET*         diveset;            /**< diving settings */
};


/*
 * local methods
 */

/** get indicator candidate variables */
static
SCIP_Real getVarBranchScore(
   SCIP*                 scip,
   SCIP_VAR*             cand,
   SCIP_Real             candfrac
   )
{
   SCIP_Bool roundup;
   SCIP_Real score;
   SCIP_Bool mayrounddown = SCIPvarMayRoundDown(cand);
   SCIP_Bool mayroundup = SCIPvarMayRoundUp(cand);

   if( mayrounddown || mayroundup )
   {
      /* choose rounding direction:
       * - if variable may be rounded in both directions, round corresponding to the fractionality
       * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
       *   the current fractional solution
       */
      if( mayrounddown && mayroundup )
         roundup = (candfrac > 0.5);
      else
         roundup = mayrounddown;

      if( roundup )
      {
         candfrac = 1.0 - candfrac;
         score = SCIPvarGetNLocksUp(cand);
      }
      else
         score = SCIPvarGetNLocksDown(cand);

      /* penalize too small fractions */
      if( candfrac < 0.01 )
         score *= 100;

      /* prefer decisions on binary variables */
      if( !SCIPvarIsBinary(cand) )
         score *= 100;

      return score + candfrac + SCIPgetNLPRows(scip);
   }
   else
   {
      /* the candidate may not be rounded */
      int nlocksdown = SCIPvarGetNLocksDown(cand);
      int nlocksup = SCIPvarGetNLocksUp(cand);
      roundup = (nlocksdown > nlocksup || (nlocksdown == nlocksup && candfrac > 0.5));
      if( roundup )
      {
         score = nlocksup;
         candfrac = 1.0 - candfrac;
      }
      else
         score = nlocksdown;

      /* penalize too small fractions */
      if( candfrac < 0.01 )
         score *= 100;

      /* prefer decisions on binary variables */
      if( !SCIPvarIsBinary(cand) )
         score *= 100;

      /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
      assert( (0.0 < candfrac && candfrac < 1.0) || SCIPvarIsBinary(cand) );

      return score + candfrac;
   }
}

static
SCIP_RETCODE getIndCandVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           indconss,           /**< indicator constraints */
   int                   nindconss,          /**< number of indicator constraints */
   SCIP_VAR**            indcands,           /**< indicator candidate variables */
   SCIP_Real*            indcandssol,        /**< solution values of candidates */
   SCIP_Real*            indcandfrac,        /**< fractionalities of candidates */
   int*                  nindcands           /**< number of candidates */
   )
{
   SCIP_VAR* binvar;
   SCIP_Real val;
   int c;

   assert( scip != NULL );
   assert( indconss != NULL );
   assert( indcands != NULL );
   assert( nindcands != NULL );
   assert( indcandssol != NULL );
   assert( indcandfrac != NULL );

   *nindcands = 0;
   for (c = 0; c < nindconss; ++c)
   {
      /* check whether constraint is violated */
      if ( SCIPisViolatedIndicator(scip, indconss[c], NULL) )
      {
         binvar = SCIPgetBinaryVarIndicator(indconss[c]);
         val = SCIPgetSolVal(scip, NULL, binvar);

         /* fractional indicator variables are treated by lpcands */
         if ( SCIPisFeasIntegral(scip, val) )
         {
            indcands[*nindcands] = binvar;
            indcandssol[*nindcands] = val;
            indcandfrac[*nindcands] = SCIPfrac(scip, val);
            ++(*nindcands);
         }
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyCoefdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeHeurCoefdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeCoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPdivesetFree(&heurdata->diveset);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitCoefdiving) /*lint --e{715}*/
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
   SCIPdivesetReset(heurdata->diveset);

   /* get indicator constraint handler */
   heurdata->indconshdlr = SCIPfindConshdlr(scip, "indicator");

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitCoefdiving) /*lint --e{715}*/
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


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecCoefdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   heurdata = SCIPheurGetData(heur);

   assert(heurdata != NULL);
   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, heurdata->diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

static
SCIP_DECL_DIVESETGETCANDS(divesetGetCandsCoefdiving)
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR** indcands;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Real* lpcandsfrac;
   SCIP_CONS** indconss;
   SCIP_Real* indcandssol;
   SCIP_Real* indcandsfrac;
   int nlpcands;
   int nindcands;
   int nindconss;


   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, &nlpcands, NULL, NULL) );

   heurdata = SCIPheurGetData(SCIPdivesetGetHeur(diveset));
   nindcands = 0;
   if( heurdata->indconshdlr != NULL )
   {
      indconss = SCIPconshdlrGetConss(heurdata->indconshdlr);
      nindconss = SCIPconshdlrGetNConss(heurdata->indconshdlr);

      if ( nindconss > 0 )
      {
         /* get storage for candidate variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &indcands, nindconss) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indcandssol, nindconss) );
         SCIP_CALL( SCIPallocBufferArray(scip, &indcandsfrac, nindconss) );
         /* get indicator candidates */
         SCIP_CALL( getIndCandVars(scip, indconss, nindconss, indcands, indcandssol, indcandsfrac, &nindcands) );
      }
   }
   *nbranchcands = nlpcands + nindcands;

   if( *nbranchcands > 0 )
   {
      SCIP_Real* scores;
      int c;

      SCIP_CALL( SCIPallocBufferArray(scip, branchcands, *nbranchcands) );
      SCIP_CALL( SCIPallocBufferArray(scip, branchcandssol, *nbranchcands) );
      SCIP_CALL( SCIPallocBufferArray(scip, branchcandsfrac, *nbranchcands) );
      if( nlpcands > 0 )
      {
         BMScopyMemoryArray(*branchcands, lpcands, nlpcands);
         BMScopyMemoryArray(*branchcandssol, lpcandssol, nlpcands);
         BMScopyMemoryArray(*branchcandsfrac, lpcandsfrac, nlpcands);
      }
      if( nindcands > 0 )
      {
         BMScopyMemoryArray(&((*branchcands)[nlpcands]), indcands, nindcands);
         BMScopyMemoryArray(&((*branchcandssol)[nlpcands]), indcandssol, nindcands);
         BMScopyMemoryArray(&((*branchcandsfrac)[nlpcands]), indcandsfrac, nindcands);
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &scores, *nbranchcands) );

      for( c = 0; c < *nbranchcands; ++c )
         scores[c] = getVarBranchScore(scip, (*branchcands)[c], (*branchcandsfrac)[c]);

      SCIPsortRealRealRealPtr(scores, *branchcandsfrac, *branchcandssol, (void **)(*branchcands), *nbranchcands);

      SCIPfreeBufferArray(scip, &scores);
   }
   if( nindconss > 0 )
   {
      SCIPfreeBufferArray(scip, &indcandsfrac);
      SCIPfreeBufferArray(scip, &indcandssol);
      SCIPfreeBufferArray(scip, &indcands);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_DIVESETFREECANDS(divesetFreecandsCoefdiving)
{
   if( nbranchcands > 0 )
   {
      SCIPfreeBufferArray(scip, branchcandsfrac);
      SCIPfreeBufferArray(scip, branchcandssol);
      SCIPfreeBufferArray(scip, branchcands);
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_DIVESETCANDBRANCHDIR(divesetCandbranchdirCoefdiving)
{
   SCIP_Bool roundup;
   SCIP_Bool mayrounddown = SCIPvarMayRoundDown(cand);
   SCIP_Bool mayroundup = SCIPvarMayRoundUp(cand);

   if( mayrounddown && mayroundup )
      roundup = (candsfrac > 0.5);
   else if( mayrounddown || mayroundup )
      roundup = mayrounddown;
   else
   {
      /* the candidate may not be rounded */
      int nlocksdown = SCIPvarGetNLocksDown(cand);
      int nlocksup = SCIPvarGetNLocksUp(cand);
      roundup = (nlocksdown > nlocksup || (nlocksdown == nlocksup && candsfrac > 0.5));
   }

   /* return corresponding branching direction */
   if( roundup )
      return SCIP_BRANCHDIR_UPWARDS;
   else
      return SCIP_BRANCHDIR_DOWNWARDS;
}


/*
 * heuristic specific interface methods
 */

/** creates the coefdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCoefdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create coefdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecCoefdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyCoefdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeCoefdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitCoefdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitCoefdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, &heurdata->diveset, heur, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_MAXLPITEROFS,
         DEFAULT_BACKTRACK, divesetGetCandsCoefdiving, divesetFreecandsCoefdiving, divesetCandbranchdirCoefdiving) );
   return SCIP_OKAY;
}

