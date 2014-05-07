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

/**@file   heur_veclendiving.c
 * @brief  LP diving heuristic that rounds variables with long column vectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_veclendiving.h"


#define HEUR_NAME             "veclendiving"
#define HEUR_DESC             "LP diving heuristic that rounds variables with long column vectors"
#define HEUR_DISPCHAR         'v'
#define HEUR_PRIORITY         -1003100
#define HEUR_FREQ             10
#define HEUR_FREQOFS          4
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
   SCIP_DIVESET*         diveset;            /**< diving settings for diving control */
};


/*
 * local methods
 */

static
SCIP_Real getVarScore(
   SCIP*                 scip,
   SCIP_VAR*             cand,
   SCIP_Real             frac
   )
{
   SCIP_Real obj;
   SCIP_Real objdelta;
   SCIP_Real score;
   SCIP_Bool roundup;
   int colveclen;

   obj = SCIPvarGetObj(cand);
   roundup = (obj >= 0.0);
   objdelta = (roundup ? (1.0-frac)*obj : -frac * obj);
   assert(objdelta >= 0.0);

   colveclen = (SCIPvarGetStatus(cand) == SCIP_VARSTATUS_COLUMN ? SCIPcolGetNNonz(SCIPvarGetCol(cand)) : 0);

   /* smaller score is better */
   score = (objdelta + SCIPsumepsilon(scip))/((SCIP_Real)colveclen+1.0);

   /* prefer decisions on binary variables */
   if( SCIPvarGetType(cand) != SCIP_VARTYPE_BINARY )
      score *= 1000.0;

   return score;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyVeclendiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurVeclendiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeVeclendiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIP_CALL( SCIPdivesetFree(&heurdata->diveset) );

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitVeclendiving) /*lint --e{715}*/
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

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitVeclendiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecVeclendiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   diveset = heurdata->diveset;

   assert(diveset != NULL);

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

static
SCIP_DECL_DIVESETGETCANDS(divesetGetCandsVeclendiving)
{
   *nbranchcands = SCIPgetNLPBranchCands(scip);
   if( *nbranchcands > 0 )
   {
      SCIP_VAR** lpcands;
      SCIP_Real* lpcandssol;
      SCIP_Real* lpcandsfrac;
      SCIP_Real* scores;
      int c;

      SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, &lpcandsfrac, nbranchcands, NULL, NULL) );

      SCIP_CALL( SCIPduplicateBufferArray(scip, branchcands, lpcands, *nbranchcands) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, branchcandssol, lpcandssol, *nbranchcands) );
      SCIP_CALL( SCIPduplicateBufferArray(scip, branchcandsfrac, lpcandsfrac, *nbranchcands) );


      SCIP_CALL( SCIPallocBufferArray(scip, &scores, *nbranchcands) );
      for( c = 0; c < *nbranchcands; ++c )
         scores[c] = getVarScore(scip, lpcands[c], lpcandsfrac[c]);

      SCIPsortRealRealRealPtr(scores, (*branchcandssol), (*branchcandsfrac), (void **)(*branchcands), *nbranchcands);

      SCIPfreeBufferArray(scip, &scores);
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_DIVESETFREECANDS(divesetFreecandsVeclendiving)
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
SCIP_DECL_DIVESETCANDBRANCHDIR(divesetCandbranchdirVeclendiving)
{
   return SCIPvarGetObj(cand) >= 0 ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS;
}

/*
 * heuristic specific interface methods
 */

/** creates the veclendiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurVeclendiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Veclendiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecVeclendiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyVeclendiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeVeclendiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitVeclendiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitVeclendiving) );

   /* veclendiving heuristic parameters */
   heurdata->diveset = NULL;
   /* create a diveset (this will automatically install some additional parameters for the heuristic) */
   SCIP_CALL( SCIPcreateDiveset(scip, &heurdata->diveset, heur, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_MAXLPITEROFS,
         DEFAULT_BACKTRACK, divesetGetCandsVeclendiving, divesetFreecandsVeclendiving, divesetCandbranchdirVeclendiving) );

   return SCIP_OKAY;
}

