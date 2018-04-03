/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_conflictdiving.c
 * @brief  LP diving heuristic that chooses fixings w.r.t. conflict locks
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_conflictdiving.h"

#define HEUR_NAME                    "conflictdiving"
#define HEUR_DESC                    "LP diving heuristic that chooses fixings w.r.t. conflict locks"
#define HEUR_DISPCHAR                '~'
#define HEUR_PRIORITY                1
#define HEUR_FREQ                    10
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
#define DEFAULT_MAXLPITERQUOT      0.15 /**< maximal fraction of diving LP iterations compared to node LP iterations */
#define DEFAULT_MAXLPITEROFS       1000 /**< additional number of allowed LP iterations */
#define DEFAULT_MAXDIVEUBQUOT       0.8 /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOT      0.0 /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                         *   where diving is performed (0.0: no limit) */
#define DEFAULT_MAXDIVEUBQUOTNOSOL  0.1 /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_MAXDIVEAVGQUOTNOSOL 0.0 /**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
#define DEFAULT_BACKTRACK          TRUE /**< use one level of backtracking if infeasibility is encountered? */
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.45/**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_MAXVIOL            TRUE /**< prefer rounding direction with most violation */
#define DEFAULT_MAXNNZOBJFAC      0.025 /**< maximal portion of nonzero objective coeffcients */
#define DEFAULT_MAXVARSFAC         -1.0 /**< maximal fraction of variables involved in a conflict constraint (< 0: auto) */
#define DEFAULT_MINMAXVARS           -1 /**< minimal absolute maximum of variables involved in a conflict constraint (-1: auto) */
#define DEFAULT_MINCONFLICTLOCKS      0 /**< threshold for penalizing the score */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */

   SCIP_Bool             maxviol;            /**< rounding into potentially infeasible or feasible direction */
   SCIP_Real             maxvarsfac;         /**< maximal fraction of variables involved in a conflict constraint */
   SCIP_Real             maxnnzobjfac;       /**< maximal portion of nonzero objective coeffcients */
   int                   minmaxvars;         /**< minimal absolute maximum of variables involved in a conflict constraint */
   int                   minconflictlocks;   /**< threshold for penalizing the score */
};

static
SCIP_Bool shouldRun(
   SCIP*                 scip,
   SCIP_HEURDATA*        heurdata,
   SCIP_RESULT*          result
   )
{
   SCIP_Real dualboundroot;
   SCIP_Real firstdualboundroot;

   assert(heurdata != NULL);

   *result = SCIP_DELAYED;

   /* don't run if no conflict constraints where found */
   if( SCIPgetNConflictConssFound(scip) == 0 )
      return FALSE;

   *result = SCIP_DIDNOTRUN;

   /* don't run if to many nonzero objective coefficients are present */
   if( SCIPgetNObjVars(scip) < SCIPceil(scip, heurdata->maxnnzobjfac * SCIPgetNVars(scip)) )
      return TRUE;

   dualboundroot = SCIPgetDualboundRoot(scip);
   firstdualboundroot = SCIPgetFirstLPDualboundRoot(scip);

   /* check whether the dual bound has not changed too much at the root node */
   if( !SCIPisInfinity(scip, REALABS(firstdualboundroot)) )
   {
      if( SCIPisPositive(scip, dualboundroot) == SCIPisPositive(scip, firstdualboundroot) )
      {
         if( SCIPisZero(scip, dualboundroot) )
         {
            if( SCIPisZero(scip, firstdualboundroot) )
               return TRUE;
         }
         else if( firstdualboundroot / dualboundroot < 1.005 && firstdualboundroot / dualboundroot > 0.995 )
            return TRUE;
      }
   }

   return FALSE;
}

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

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecConflictdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;
   SCIP_Real maxvarsfac;
   int minmaxvars;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   if( !shouldRun(scip, heurdata, result) )
      return SCIP_OKAY;

   maxvarsfac = SCIP_INVALID;
   minmaxvars = INT_MAX;

   if( heurtiming == SCIP_HEURTIMING_DURINGLPLOOP && SCIPgetDepth(scip) != 0 )
      return SCIP_OKAY;

   if( heurdata->maxvarsfac >= 0.0 && !SCIPisParamFixed(scip, "conflict/maxvarsfac") )
   {
      SCIP_CALL( SCIPgetRealParam(scip, "conflict/maxvarsfac", &maxvarsfac) );
      SCIP_CALL( SCIPsetRealParam(scip, "conflict/maxvarsfac", heurdata->maxvarsfac) );
   }
   if( heurdata->minmaxvars >= 0 && !SCIPisParamFixed(scip, "conflict/minmaxvars") )
   {
      SCIP_CALL( SCIPgetIntParam(scip, "conflict/minmaxvars", &minmaxvars) );
      SCIP_CALL( SCIPsetIntParam(scip, "conflict/minmaxvars", heurdata->minmaxvars) );
   }

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   if( heurdata->maxvarsfac >= 0.0 && !SCIPisParamFixed(scip, "conflict/maxvarsfac") )
   {
      assert(maxvarsfac != SCIP_INVALID); /*lint !e777*/
      SCIP_CALL( SCIPsetRealParam(scip, "conflict/maxvarsfac", maxvarsfac) );
   }
   if( heurdata->minmaxvars >= 0.0 && !SCIPisParamFixed(scip, "conflict/minmaxvars") )
   {
      assert(minmaxvars != INT_MAX);
      SCIP_CALL( SCIPsetIntParam(scip, "conflict/minmaxvars", minmaxvars) );
   }

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
   SCIP_Real conflictlocksum;
   SCIP_Real locksum;
   SCIP_Bool mayrounddown;
   SCIP_Bool mayroundup;
   int nlocksup;
   int nlocksdown;
   int nconflictlocksup;
   int nconflictlocksdown;

   rng = SCIPdivesetGetRandnumgen(diveset);
   assert(rng != NULL);

   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   nlocksup = SCIPvarGetNLocksUp(cand);
   nlocksdown = SCIPvarGetNLocksDown(cand);

   nconflictlocksup = SCIPvarGetNLocksUpType(cand, SCIP_LOCKTYPE_CONFLICT);
   nconflictlocksdown = SCIPvarGetNLocksDownType(cand, SCIP_LOCKTYPE_CONFLICT);

   conflictlocksum = nconflictlocksup + nconflictlocksdown;
   locksum = nlocksdown + nlocksup;

   mayrounddown = (nconflictlocksdown == 0);
   mayroundup = (nconflictlocksup == 0);

   /* variable can be rounded in exactly one direction and we try to go into the feasible direction */
   if( mayrounddown != mayroundup )
   {
      *roundup = mayroundup;
   }
   /* variable is locked in both directions */
   else if( !mayroundup )
   {
      assert(!mayrounddown);

      if( nconflictlocksup != nconflictlocksdown || nlocksup != nlocksdown )
      {
         if( nconflictlocksup != nconflictlocksdown )
         {
            *roundup = (nconflictlocksup > nconflictlocksdown);
         }
         else
         {
            assert(nlocksup != nlocksdown);
            *roundup = (nlocksup > nlocksdown);
         }

         if( !heurdata->maxviol )
            *roundup = !(*roundup);
      }
      else if( !SCIPisEQ(scip, candsfrac, 0.5) )
      {
         *roundup = (heurdata->maxviol ? candsfrac < 0.5 : candsfrac > 0.5);
      }
      else
      {
         *roundup = (SCIPrandomGetInt(rng, 0, 1) == 1);
      }
   }
   /* the variable is not locked by conflict constraints */
   else
   {
      assert(nconflictlocksdown == 0 && nconflictlocksup == 0);

      if( nlocksup != nlocksdown )
      {
         *roundup = (nconflictlocksup > nconflictlocksdown);

         if( !heurdata->maxviol )
            *roundup = !(*roundup);

      }
      else if( !SCIPisEQ(scip, candsfrac, 0.5) )
      {
         *roundup = (heurdata->maxviol ? candsfrac < 0.5 : candsfrac > 0.5);
      }
      else
      {
         *roundup = (SCIPrandomGetInt(rng, 0, 1) == 1);
      }
   }

   if( *roundup )
   {
      SCIP_Real scalefactor;

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

      scalefactor = (LOCKFRAC + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND));

      *score = candsfrac;

      if( nconflictlocksup > 0 )
         *score += 10 * nconflictlocksup / conflictlocksum;

      *score += (scalefactor * nlocksup / MAX(1.0, locksum));
   }
   else
   {
      SCIP_Real scalefactor;

      if ( divetype == SCIP_DIVETYPE_SOS1VARIABLE && SCIPisFeasNegative(scip, candsol) )
         candsfrac = 1.0 - candsfrac;

      scalefactor = (LOCKFRAC + SCIPrandomGetReal(rng, MIN_RAND, MAX_RAND));

      *score = candsfrac;

      if( nconflictlocksdown > 0 )
         *score += 10 * nconflictlocksdown / conflictlocksum;

      *score += (scalefactor * nlocksdown / MAX(1.0, locksum));
   }

   /* penalize too less conflict locks */
   if( conflictlocksum < heurdata->minconflictlocks )
      (*score) *= 0.1;

   if( conflictlocksum == 0 )
      (*score) *= 0.01;

   /* penalize too small fractions */
   if( candsfrac < 0.01 )
      (*score) *= 0.1;

   /* prefer decisions on binary variables */
   if( !SCIPvarIsBinary(cand) )
      (*score) *= 0.1;

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

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/maxviol", "try to maximize the violation",
         &heurdata->maxviol, TRUE, DEFAULT_MAXVIOL, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minconflictlocks",
         "minimal number of conflict locks per variable",
         &heurdata->minconflictlocks, TRUE, DEFAULT_MINCONFLICTLOCKS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/minmaxvars",
         "minimal absolute maximum of variables involved in a conflict constraint (-1: auto)",
         &heurdata->minmaxvars, TRUE, DEFAULT_MINMAXVARS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxvarsfac",
         "maximal fraction of variables involved in a conflict constraint (< 0: auto)",
         &heurdata->maxvarsfac, TRUE, DEFAULT_MAXVARSFAC, -1.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/maxnnzobjfac",
         "maximal portion of nonzero objective coeffcients.",
         &heurdata->maxnnzobjfac, TRUE, DEFAULT_MAXNNZOBJFAC
      , 0.0, 1.0, NULL, NULL) );


   return SCIP_OKAY;
}
