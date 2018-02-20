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

/**@file   heur_farkasdiving.c
 * @brief  LP diving heuristic that tries to construct a Farkas-proof
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_farkasdiving.h"

#define HEUR_NAME             "farkasdiving"
#define HEUR_DESC             "LP diving heuristic that tries to construct a Farkas-proof"
#define HEUR_DISPCHAR         'u'
#define HEUR_PRIORITY         -1001250
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */


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
#define DEFAULT_LPRESOLVEDOMCHGQUOT 0.15 /**< percentage of immediate domain changes during probing to trigger LP resolve */
#define DEFAULT_LPSOLVEFREQ           0 /**< LP solve frequency for diving heuristics */
#define DEFAULT_ONLYLPBRANCHCANDS FALSE /**< should only LP branching candidates be considered instead of the slower but
                                         *   more general constraint handler diving variable selection? */
#define DEFAULT_RANDSEED            151 /**< initial seed for random number generation */

#define DEFAULT_DIFFOBJFAC         0.15
#define DEFAULT_CHECKOBJ          FALSE
#define DEFAULT_CHECKOBJGLB        TRUE
#define DEFAULT_SCALESCORE         TRUE

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             diffobjfac;         /**< fraction of absolute different objective coefficients */
   SCIP_Bool             disabled;           /**< remember if the heuristic should not run at all */
   SCIP_Bool             checkobj;           /**< should objective function be checked before running? */
   SCIP_Bool             checkobjglb;        /**< check objective function only once w.r.t to the global problem */
   SCIP_Bool             objchecked;         /**< remember whether objection function was already checked */
   SCIP_Bool             scalescore;
};


/*
 * local methods
 */

/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyFarkasdiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurFarkasdiving(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeFarkasdiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitFarkasdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   heurdata->disabled = FALSE;
   heurdata->objchecked = FALSE;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitFarkasdiving) /*lint --e{715}*/
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
static
SCIP_DECL_HEURINITSOL(heurInitsolFarkasdiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   heurdata->disabled = FALSE;
   heurdata->objchecked = FALSE;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecFarkasdiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;
   int nnzobjvars;

   heurdata = SCIPheurGetData(heur);
   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);

   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   *result = SCIP_DIDNOTRUN;

   /* terminate if the heuristic has been disabled */
   if( heurdata->disabled )
      return SCIP_OKAY;

   /* terminate if there are no integer variables (note that, e.g., SOS1 variables may be present) */
   if( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0 )
      return SCIP_OKAY;

   nnzobjvars = SCIPgetNObjVars(scip);

   /* terminate if at most one variable has a non-zero objective coefficient */
   if( nnzobjvars == 0 )
   {
      heurdata->disabled = TRUE;
      return SCIP_OKAY;
   }

   /* initialize the hashing of all objective coefficients */
   if( heurdata->checkobj )
   {
      SCIP_VAR** divecandvars;
      SCIP_Real* objcoefs;
      SCIP_Real lastobjcoef;
      int nnzobjcoefs;
      int ndiffnnzobjs;
      int ndivecands;
      int i;

      divecandvars = NULL;
      ndivecands = -1;

      /* check w.r.t. all variables and decide whether we want to run or not */
      if( heurdata->checkobjglb )
      {
         if( !heurdata->objchecked )
         {
            divecandvars = SCIPgetVars(scip);
            ndivecands = SCIPgetNVars(scip);
         }
         else
         {
            goto PERFORMDIVING;
         }
      }
      else if( !heurdata->checkobjglb )
      {
         /* we can only access the branching candidates if the LP is solved to optimality */
         if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         {
            SCIP_CALL( SCIPgetLPBranchCands(scip, &divecandvars, NULL, NULL, &ndivecands, NULL, NULL) );
         }
         else
         {
            goto PERFORMDIVING;
         }
      }

      assert(divecandvars != NULL);
      assert(ndivecands >= 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &objcoefs, ndivecands) );

      /* collect all absolute values of objective coefficients */
      nnzobjcoefs = 0;
      for( i = 0; i < ndivecands; i++ )
      {
         SCIP_Real obj = SCIPvarGetObj(divecandvars[i]);

         if( SCIPisZero(scip, obj) )
            continue;

         objcoefs[nnzobjcoefs] = REALABS(obj);
         ++nnzobjcoefs;
      }

      /* sort in increasing order */
      SCIPsortReal(objcoefs, nnzobjcoefs);
      assert(!SCIPisZero(scip, objcoefs[0]));

      lastobjcoef = objcoefs[0];
      ndiffnnzobjs = 1;

      for( i = 1; i < nnzobjcoefs; i++ )
      {
         if( SCIPisGT(scip, objcoefs[i], lastobjcoef) )
         {
            lastobjcoef = objcoefs[i];
            ++ndiffnnzobjs;
         }
      }

      SCIPfreeBufferArray(scip, &objcoefs);

      SCIPdebugMsg(scip, "%d divecands; %d nnzobjs; %d diffnnzobjs", ndivecands, nnzobjcoefs, ndiffnnzobjs);

//      if( nnzobjcoefs == 0 || ndiffnnzobjs / (SCIP_Real)ndivecands < heurdata->diffobjfac )
      if( nnzobjcoefs == 0 || nnzobjcoefs / (SCIP_Real)ndivecands < heurdata->diffobjfac )
      {
         SCIPdebugMsg(scip, " ---> disable farkasdiving%s\n", heurdata->checkobjglb ? "" : " locally");

         /* disable the heuristic if we want to check the objective only once */
         if( heurdata->checkobjglb )
            heurdata->disabled = TRUE;

         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMsg(scip, "\n");
         heurdata->objchecked = TRUE;
      }
   }

  PERFORMDIVING:
   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible) );

   return SCIP_OKAY;
}

#define MIN_RAND 1e-06
#define MAX_RAND 1e-05

/** calculate score and preferred rounding direction for the candidate variable */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreFarkasdiving)
{
   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real obj;

   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   randnumgen = SCIPdivesetGetRandnumgen(diveset);
   assert(randnumgen != NULL);

   obj = SCIPvarGetObj(cand);

   /* for the reduced costs of a basic variable it holds c_i - y^TA_i = 0, so we approximate the reduced costs
    * within an infeasibility proof by r_i := -y^TA_i = -c_i.
    */
   if( SCIPisNegative(scip, obj) )
   {
      *roundup = TRUE;
   }
   else if( SCIPisPositive(scip, obj) )
   {
      *roundup = FALSE;
   }
   else
   {
      if( SCIPisEQ(scip, candsfrac, 0.5) )
         *roundup = !SCIPrandomGetInt(randnumgen, 0, 1);
      else
         *roundup = (candsfrac > 0.5);
   }

   /* larger score is better */
   *score = REALABS(obj) + SCIPrandomGetReal(randnumgen, MIN_RAND, MAX_RAND);

   if( heurdata->scalescore )
   {
      if( *roundup )
         *score *= (1.0 - candsfrac);
      else
         *score *= candsfrac;
   }

   /* prefer decisions on binary variables */
   if( SCIPvarGetType(cand) != SCIP_VARTYPE_BINARY )
      *score = -1.0 / *score;

   return SCIP_OKAY;

}

/*
 * heuristic specific interface methods
 */

/** creates the farkasdiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurFarkasdiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Farkasdiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecFarkasdiving, heurdata) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyFarkasdiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeFarkasdiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitFarkasdiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitFarkasdiving) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolFarkasdiving) );

   /* farkasdiving heuristic parameters */
   /* create a diveset (this will automatically install some additional parameters for the heuristic) */
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS, DIVESET_DIVETYPES, divesetGetScoreFarkasdiving) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/checkobj",
         "should objective function be check before running?",
         &heurdata->checkobj, TRUE, DEFAULT_CHECKOBJ, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/checkobjglb",
         "should objective function be check globally or locally?",
         &heurdata->checkobjglb, TRUE, DEFAULT_CHECKOBJGLB, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/scalescore",
         "should score be scaled by fractionality?",
         &heurdata->scalescore, TRUE, DEFAULT_SCALESCORE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/diffobjfac",
         "fraction of absolute different objective coefficients",
         &heurdata->diffobjfac, TRUE, DEFAULT_DIFFOBJFAC, 0.0, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}

