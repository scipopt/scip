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

/**@file   heur_guideddiving.c
 * @brief  LP diving heuristic that chooses fixings in direction of incumbent solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_guideddiving.h"
#include "scip/heur_fracdiving.h"

#define HEUR_NAME             "guideddiving"
#define HEUR_DESC             "LP diving heuristic that chooses fixings in direction of incumbent solutions"
#define HEUR_DISPCHAR         'g'
#define HEUR_PRIORITY         -1007000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          7
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

/** returns whether variable should be rounded up */
static
SCIP_Bool varRoundUp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< candidate variable */
   SCIP_Real             candsol             /**< LP solution value of candidate */
)
{
   SCIP_SOL* bestsol;
   SCIP_Real bestsolval;
   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);
   assert(!SCIPsolIsOriginal(bestsol));

   bestsolval = SCIPgetSolVal(scip, bestsol, cand);

   /* variable should be rounded (guided) into the direction of its incumbent solution value */
   if( candsol < bestsolval )
      return TRUE;
   else
      return FALSE;
}

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyGuideddiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurGuideddiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeGuideddiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->diveset != NULL);

   SCIPdivesetFree(&heurdata->diveset);

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitGuideddiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   /* initialize data by resetting the diving settings */
   SCIPdivesetReset(heurdata->diveset);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitGuideddiving) /*lint --e{715}*/
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
SCIP_DECL_HEUREXEC(heurExecGuideddiving) /*lint --e{715}*/
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

  /* don't dive, if no feasible solutions exist */
   if( SCIPgetNSols(scip) == 0 )
      return SCIP_OKAY;

   /* get best solution that should guide the search; if this solution lives in the original variable space,
    * we cannot use it since it might violate the global bounds of the current problem
    */
   if( SCIPsolIsOriginal(SCIPgetBestSol(scip)) )
      return SCIP_OKAY;

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   diveset = heurdata->diveset;

   /* call generic diving algorithm */
   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );


   return SCIP_OKAY;
}

/* callbacks for diving */

/** returns the preferred branching direction of candidate */
static
SCIP_DECL_DIVESETCANDBRANCHDIR(divesetCandbranchdirGuideddiving)
{
   if( varRoundUp(scip, cand, candsol) )
      return SCIP_BRANCHDIR_UPWARDS;
   else
      return SCIP_BRANCHDIR_DOWNWARDS;
}

/** returns a score for the given candidate -- the best candidate minimizes the diving score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreGuideddiving)
{
   SCIP_Bool roundup;

   /* the rounding direction is the direction towards the incumbent solution */
   roundup = varRoundUp(scip, cand, candsol);

   /* use fractionality score as in frac diving */
   return SCIPgetFracdivingVarScore(scip, cand, candsfrac, roundup);
}

/*
 * heuristic specific interface methods
 */

/** creates the guideddiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurGuideddiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Guideddiving primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecGuideddiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyGuideddiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeGuideddiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitGuideddiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitGuideddiving) );

   heurdata->diveset = NULL;
   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, &heurdata->diveset, heur, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, 1.0, 1.0, DEFAULT_MAXLPITEROFS,
         DEFAULT_BACKTRACK, divesetGetScoreGuideddiving, divesetCandbranchdirGuideddiving, NULL, NULL) );

   return SCIP_OKAY;
}

