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

/**@file   heur_fracdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_fracdiving.h"


#define HEUR_NAME             "fracdiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings w.r.t. the fractionalities"
#define HEUR_DISPCHAR         'f'
#define HEUR_PRIORITY         -1003000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          3
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
   SCIP_DIVESET*         diveset;            /**< diving settings */
};


/*
 * local methods
 */

/** score candidate variable
 *
 *  if candidate cannot be trivially rounded,
 *  score is the difference between frac and [frac] (rounding to the nearest integer)
 *
 *  if candidate can be trivially rounded in at least one direction, the objective gain is used to score the variables
 *
 *  candidate which cannot be rounded trivially always have a lower score than roundable candidates
 */
static
SCIP_Real getVarScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< diving candidate for score */
   SCIP_Real             frac                /**< fractionality of candidate in (last) LP solution */
   )
{
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;
   SCIP_Real obj;
   SCIP_Real objnorm;
   SCIP_Real objgain;

   mayrounddown = SCIPvarMayRoundDown(cand);
   mayroundup = SCIPvarMayRoundUp(cand);
   obj = SCIPvarGetObj(cand);
   objnorm = SCIPgetObjNorm(scip);

   /* divide by objective norm to normalize obj into [-1,1] */
   if( SCIPisPositive(scip, objnorm) )
      obj /= objnorm;


   /* choose rounding direction:
    * - if variable may be rounded in both directions, round corresponding to the fractionality
    * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
    *   the current fractional solution
    */
   if( mayrounddown != mayroundup )
      roundup = mayrounddown;
   else
      roundup = (frac > 0.5);

   if( roundup )
   {
      frac = 1.0 - frac;
      objgain = obj * frac;
   }
   else
      objgain = -obj * frac;

   assert(objgain >= -1.0 && objgain <= 1.0);

   /* penalize too small fractions */
   if( frac < 0.01 )
      frac += 10.0;

   /* prefer decisions on binary variables */
   if( !SCIPvarIsBinary(cand) )
      frac *= 1000.0;

   /* prefer variables which cannot be rounded by scoring their fractionality */
   if( !(mayrounddown || mayroundup) )
      return frac;
   else
      return 2.0 + objgain;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyFracdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurFracdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeFracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free diving settings */
   SCIP_CALL( SCIPdivesetFree(&heurdata->diveset) );

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitFracdiving) /*lint --e{715}*/
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
   assert(heurdata->diveset != NULL);
   SCIPdivesetReset(heurdata->diveset);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitFracdiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecFracdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   diveset = heurdata->diveset;

   assert(diveset != NULL);

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

/** returns the preferred branching direction of candidate */
static
SCIP_DECL_DIVESETCANDBRANCHDIR(divesetCandbranchdirFracdiving)
{
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   SCIP_Bool roundup;

   mayrounddown = SCIPvarMayRoundDown(cand);
   mayroundup = SCIPvarMayRoundUp(cand);

   /* choose rounding direction:
    * - if variable may be rounded in both directions, round corresponding to the fractionality
    * - otherwise, round in the infeasible direction, because feasible direction is tried by rounding
    *   the current fractional solution
    */
   if( mayrounddown != mayroundup )
      roundup = mayrounddown;
   else
      roundup = (candsfrac > 0.5);

   return roundup ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS;
}

/** returns a score for the given candidate -- the best candidate minimizes the diving score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreFracdiving)
{
   return getVarScore(scip, cand, candsfrac);
}

/*
 * heuristic specific interface methods
 */

/** creates the fracdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFracdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Fracdiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFracdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyFracdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeFracdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitFracdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitFracdiving) );

   heurdata->diveset = NULL;
   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, &heurdata->diveset, heur, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_MAXLPITEROFS,
         DEFAULT_BACKTRACK, divesetGetScoreFracdiving, divesetCandbranchdirFracdiving, NULL, NULL) );
   return SCIP_OKAY;
}

