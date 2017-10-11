/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_conflictdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. soft locks
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_conflictdiving.h"

#define HEUR_NAME                    "conflictdiving"
#define HEUR_DESC                    "LP diving heuristic that chooses fixings w.r.t. soft locks"
#define HEUR_DISPCHAR                '~'
#define HEUR_PRIORITY                -1000250
#define HEUR_FREQ                    -1
#define HEUR_FREQOFS                 0
#define HEUR_MAXDEPTH                -1
#define HEUR_TIMING                  SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP             FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES            SCIP_DIVETYPE_INTEGRALITY | SCIP_DIVETYPE_SOS1VARIABLE /**< bit mask that represents all supported dive types */
#define DEFAULT_RANDSEED             151 /**< default random seed */

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
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15/**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */

#define DEFAULT_ADDSOLUTION        TRUE /**< should the solution be added to the solution pool */
#define DEFAULT_EXPREGION           'l' /**< which region of the search space should be explore? ('l'ocked, 'u'nlocked)*/

#define DEFAULT_MINNUMSOFTLOCKS      0

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   int*                  inds;
   int*                  nsoftlocks;
   int                   ninds;

   SCIP_Bool             addsol;
   int                   minnumsoftlocks;
   char                  expregion;
};

/*
 * local methods
 */

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyConflictdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeHeurConflictdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeConflictdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitConflictdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   heurdata->inds = NULL;
   heurdata->nsoftlocks = NULL;
   heurdata->ninds = -1;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitConflictdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   if( heurdata->ninds > -1 )
   {
      SCIPfreeBlockMemoryArray(scip, &heurdata->inds, SCIPgetNVars(scip));
      SCIPfreeBlockMemoryArray(scip, &heurdata->nsoftlocks, SCIPgetNVars(scip));
   }

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecConflictdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;
   int nvars;
   int i;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DIDNOTRUN;

   if( heurtiming == SCIP_HEURTIMING_DURINGLPLOOP && SCIPgetDepth(scip) != 0 )
      return SCIP_OKAY;

   assert((heurtiming & SCIP_HEURTIMING_AFTERLPPLUNGE) == 0);

   if( heurtiming == SCIP_HEURTIMING_AFTERLPPLUNGE && SCIPgetDepth(scip) != 0 )
      return SCIP_OKAY;

   *result = SCIP_DELAYED;

   vars = SCIPgetVars(scip);
   assert(vars != NULL);

   nvars = SCIPgetNVars(scip);

   if( heurdata->ninds == -1 )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->inds, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->nsoftlocks, nvars) );

      heurdata->ninds = 0;

      for( i = 0; i < nvars; i++ )
      {
         heurdata->nsoftlocks[i] = SCIPvarGetNLocksSoftUp(vars[i]) + SCIPvarGetNLocksSoftDown(vars[i]);

         if( heurdata->nsoftlocks[i] < heurdata->minnumsoftlocks )
            heurdata->inds[heurdata->ninds++] = i;
      }
   }
   else
   {
      /* check number of softlocks */
      for( i = 0; i < heurdata->ninds; )
      {
         int v = heurdata->inds[i];

         heurdata->nsoftlocks[v] = SCIPvarGetNLocksSoftUp(vars[v]) + SCIPvarGetNLocksSoftDown(vars[v]);

         if( heurdata->nsoftlocks[v] >= heurdata->minnumsoftlocks )
         {
            --heurdata->ninds;
            heurdata->inds[i] = heurdata->inds[heurdata->ninds];
         }
         else
            i++;
      }
   }
   assert(heurdata->ninds >= 0);

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

#define MIN_RAND 1e-06
#define MAX_RAND 1e-05
#define LOCKFRAC 1e-04

/** returns a score for the given candidate -- the best candidate maximizes the diving score */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreConflictdiving)
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* rng;
   SCIP_Real softlocksum;
   SCIP_Real locksum;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;

   rng = SCIPdivesetGetRandnumgen(diveset);
   assert(rng != NULL);

   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   softlocksum = SCIPvarGetNLocksSoftDown(cand) + SCIPvarGetNLocksSoftUp(cand);
   locksum = SCIPvarGetNLocksDown(cand) + SCIPvarGetNLocksUp(cand);

   mayrounddown = (SCIPvarGetNLocksSoftDown(cand) == 0);
   mayroundup = (SCIPvarGetNLocksSoftUp(cand) == 0);

   /* variable can be rounded in exactly one direction */
   if( mayrounddown != mayroundup )
   {
      *roundup = mayroundup;
   }
   else
   {
      *roundup = (SCIPvarGetNLocksSoftUp(cand) >= SCIPvarGetNLocksSoftDown(cand));
   }

   if( *roundup )
   {
      switch( divetype )
      {
         case SCIP_DIVETYPE_INTEGRALITY:
            candsfrac = 1.0 - candsfrac;
            break;
         case SCIP_DIVETYPE_SOS1VARIABLE:
            if ( SCIPisFeasPositive(scip, candsol) )
               candsfrac = 1.0 - candsfrac;
            break;
         default:
            SCIPerrorMessage("Error: Unsupported diving type\n");
            SCIPABORT();
            return SCIP_INVALIDDATA; /*lint !e527*/
      } /*lint !e788*/


      if( heurdata->expregion == 'l' )
         *score = SCIPvarGetNLocksSoftUp(cand)/MAX(1.0, softlocksum)
               + (LOCKFRAC + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND)) * SCIPvarGetNLocksUp(cand)/MAX(1.0, locksum);
      else
         *score = SCIPvarGetNLocksSoftDown(cand)/MAX(1.0, softlocksum)
               + (LOCKFRAC + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND)) * SCIPvarGetNLocksDown(cand)/MAX(1.0, locksum);
   }
   else
   {
      if ( divetype == SCIP_DIVETYPE_SOS1VARIABLE && SCIPisFeasNegative(scip, candsol) )
         candsfrac = 1.0 - candsfrac;

      if( heurdata->expregion == 'l' )
         *score = SCIPvarGetNLocksSoftDown(cand)/MAX(1.0, softlocksum)
               + (LOCKFRAC + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND)) * SCIPvarGetNLocksDown(cand)/MAX(1.0, locksum);
      else
         *score = SCIPvarGetNLocksSoftUp(cand)/MAX(1.0, softlocksum)
               + (LOCKFRAC + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND)) * SCIPvarGetNLocksUp(cand)/MAX(1.0, locksum);
   }

   /* penalize too less softlocks */
   if( heurdata->nsoftlocks[SCIPvarGetProbindex(cand)] < heurdata->minnumsoftlocks )
      (*score) *= 0.1;

   /* penalize too small fractions */
   if( candsfrac < 0.01 )
      (*score) *= 0.1;

   /* prefer decisions on binary variables */
   if( !SCIPvarIsBinary(cand) )
      (*score) *= 0.1;

   /* penalize the variable if it may be rounded. */
   if( mayrounddown || mayroundup )
      (*score) -= SCIPgetNLPRows(scip);

   /* check, if candidate is new best candidate: prefer unroundable candidates in any case */
   assert( (0.0 < candsfrac && candsfrac < 1.0) || SCIPvarIsBinary(cand) || divetype == SCIP_DIVETYPE_SOS1VARIABLE );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** creates the conflictdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurConflictdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create conflictdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ,
         HEUR_FREQOFS, HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecConflictdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyConflictdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeConflictdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitConflictdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitConflictdiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS, DIVESET_DIVETYPES, divesetGetScoreConflictdiving) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/addsol", "should solution be added to the solution pool",
         &heurdata->addsol, TRUE, DEFAULT_ADDSOLUTION, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minnumsoftlocks", "minimal number of softlocks per variable",
         &heurdata->minnumsoftlocks, TRUE, DEFAULT_MINNUMSOFTLOCKS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/" HEUR_NAME "/expregion", "which region should be explored",
         &heurdata->expregion, TRUE, DEFAULT_EXPREGION, "lu", NULL, NULL) );

   return SCIP_OKAY;
}

