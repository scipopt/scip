/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_indicatordiving.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  indicator diving heuristic
 * @author Katrin Halbig
 * @author Alexander Hoen
 * #TODO: deactivate heuristic if model contains no indicator variable
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <assert.h>

#include "scip/cons_indicator.h"
#include "scip/heur_indicatordiving.h"
#include "scip/heuristics.h"
#include "scip/pub_cons.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_heur.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_sol.h"
#include "scip/scip_tree.h"
#include "struct_heur.h"
#include <string.h>

#define HEUR_NAME             "indicatordiving"
#define HEUR_DESC             "indicator diving heuristic"
#define HEUR_DISPCHAR         '?' /**< todo: change to SCIP_HEURDISPCHAR_DIVING */
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPPLUNGE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */
#define DIVESET_DIVETYPES     SCIP_DIVETYPE_INTEGRALITY /**< bit mask that represents all supported dive types */
#define DIVESET_ISPUBLIC      FALSE   /**< is this dive set publicly available (ie., can be used by other primal heuristics?) */


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
#define DEFAULT_RANDSEED             11  /**< initial seed for random number generation */
#define DEFAULT_ROUNDINGFRAC         50 /**< default parameter setting for parameter roundingfrac */
#define DEFAULT_MODE                  2 /**< default parameter setting for parameter mode */

/** locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Real             roundingfrac;       /**< in fractional case all fractional below this value are rounded up*/
   int                   mode;               /**< decides which mode is selected (0: rounding down; 1: rounding up; 2: fractional rounding (default))*/
};

/*
 * Local methods
 */

/** checks if variable is indicator variable and returns corresponding indicator constraint */
static
SCIP_RETCODE checkAndGetIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< candidate variable */
   SCIP_CONS**           cons,               /**< pointer to store indicator constraint */
   SCIP_Bool*            isindicator         /**< pointer to store whether candidate variable is indicator variable */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** indicatorconss;
   int c;

   assert(scip != NULL);
   assert(cand != NULL);
   assert(cons != NULL);
   assert(isindicator != NULL);

   if( SCIPvarGetType(cand) != 0 )
   {
      *cons = NULL;
      *isindicator = FALSE;
      return SCIP_OKAY;
   }

   conshdlr = SCIPfindConshdlr(scip, "indicator");
   indicatorconss = SCIPconshdlrGetConss(conshdlr);

   *isindicator = FALSE;
   for( c = 0; c < SCIPconshdlrGetNActiveConss(conshdlr); c++ )
   {
      SCIP_VAR* indicatorvar;
      indicatorvar = SCIPgetBinaryVarIndicator(indicatorconss[c]);

      if( cand == indicatorvar )
      {
         *cons = indicatorconss[c];
         *isindicator = TRUE;
         return SCIP_OKAY;
      }
   }

   *cons = NULL;
   *isindicator = FALSE;
   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyIndicatordiving)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurIndicatordiving(scip) );

   return SCIP_OKAY;
}


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeIndicatordiving) /*lint --e{715}*/
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
SCIP_DECL_HEURINIT(heurInitIndicatordiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitIndicatordiving) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecIndicatordiving)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_DIVESET* diveset;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   assert(SCIPheurGetNDivesets(heur) > 0);
   assert(SCIPheurGetDivesets(heur) != NULL);
   diveset = SCIPheurGetDivesets(heur)[0];
   assert(diveset != NULL);

   SCIPdebugMessage("call heurExecIndicatordiving at depth %d \n", SCIPgetDepth(scip));

   SCIP_CALL( SCIPperformGenericDivingAlgorithm(scip, diveset, heurdata->sol, heur, result, nodeinfeasible, -1L, SCIP_DIVECONTEXT_SINGLE) );

   SCIPdebugMessage("leave heurExecIndicatordiving\n");

   return SCIP_OKAY;
}

/** calculate score and preferred rounding direction for the candidate variable */
static
SCIP_DECL_DIVESETGETSCORE(divesetGetScoreIndicatordiving)
{
   /*input:
    * - scip : SCIP main data structure
    * - diveset : diving settings for scoring
    * - divetype : represents different methods for a dive set to explore the next children
    * - cand : candidate variable for which the score should be determined
    * - candsol : solution value of variable in LP relaxation solution
    * - candsfrac : fractional part of solution value of variable
    * - score : pointer for diving score value - the best candidate maximizes this score
    * - roundup : pointer to store whether the preferred rounding direction is upwards
    *
    */

   SCIP_HEUR* heur;
   SCIP_HEURDATA* heurdata;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_VAR** consvars;
   SCIP_CONS* indicatorcons;
   SCIP_CONS* lincons;
   SCIP_VAR* slackvar;
   SCIP_Real* consvals;
   SCIP_Real activity;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nconsvars;
   SCIP_Bool isindicatorvar;
   SCIP_Bool success;
   int v;

   /* check if cand variable is indicator variable */
   SCIP_CALL( checkAndGetIndicator(scip, cand, &indicatorcons, &isindicatorvar) );

   if( !isindicatorvar )
      return SCIP_OKAY;

   SCIPdebugMessage("cand: %s, candsol: %.2f, candobjcoeff: %f\n", SCIPvarGetName(cand), candsol, SCIPvarGetObj(cand));

   lincons = SCIPgetLinearConsIndicator(indicatorcons);
   slackvar = SCIPgetSlackVarIndicator(indicatorcons);
   lhs = SCIPconsGetLhs(scip, lincons, &success);
   rhs = SCIPconsGetRhs(scip, lincons, &success);
   assert(SCIPisLE(scip, lhs, rhs));

   SCIPdebugPrintCons(scip, indicatorcons, NULL);
   SCIPdebugPrintCons(scip, lincons, NULL);
   SCIPdebugMessage("slackvar has UB = %f\n", SCIPvarGetUbGlobal(slackvar));
   SCIPdebugMessage("cons lhs %f\n", SCIPconsGetLhs(scip, lincons, &success));

   SCIPgetConsNVars(scip, lincons, &nconsvars, &success);
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nconsvars) );
   SCIPgetConsVars(scip, lincons, consvars, nconsvars, &success);
   SCIPgetConsVals(scip, lincons, consvals, nconsvars, &success);

   /* get heuristic data */
   heur = SCIPdivesetGetHeur(diveset);
   assert(heur != NULL);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   randnumgen = SCIPdivesetGetRandnumgen(diveset);
   assert(randnumgen != NULL);

   activity = 0;
   for( v = 0; v < nconsvars - 1; v++ )
   {
      SCIPdebugMessage("%s lp sol %f %f\n", SCIPvarGetName(consvars[v]), SCIPvarGetLPSol(consvars[v]),
                       consvals[v]);
      activity += consvals[v] * SCIPvarGetLPSol(consvars[v]);
   }
   assert(consvars[v]==slackvar);
   SCIPdebugMessage("activity: %f\n", activity);

   if( heurdata->mode == 0 )
   {
      *roundup = FALSE;
      if( SCIPisGE(scip, activity, lhs) && SCIPisLE(scip, activity, rhs) )
      {
         /* indicator constraint is feasible */
         *score = SCIPrandomGetReal(randnumgen, -1, 0);
      }
      else if( SCIPisGT(scip, activity, rhs))
      {
         *score = 100 * (activity - rhs) / MAX(ABS(rhs),1);
         assert(score>=0);
      }
      else if( SCIPisLT(scip, activity, lhs))
      {
         *score = 100 * (lhs - activity) / MAX(ABS(lhs),1);
         assert(score>=0);
      }
      else
      {
         assert(FALSE);
      }
   }
   else if( heurdata->mode == 1)
   {
      *roundup = TRUE;
      if( SCIPisGE(scip, activity, lhs) && SCIPisLE(scip, activity, rhs) )
      {
         /* indicator constraint is feasible */
         *score = SCIPrandomGetReal(randnumgen, -1, 0);
      }
      else if( SCIPisGT(scip, activity, rhs))
      {
         *score = 100 * (activity - rhs) / MAX(ABS(rhs),1);
         assert(score>=0);
      }
      else if( SCIPisLT(scip, activity, lhs))
      {
         *score = 100 * (lhs - activity) / MAX(ABS(lhs),1);
         assert(score>=0);
      }
      else
      {
         assert(FALSE);
      }
   }
   else
   {
      if( SCIPisGE(scip, activity, lhs) && SCIPisLE(scip, activity, rhs) )
      {
         /* indicator constraint is feasible */
         *roundup = TRUE;
         *score = SCIPrandomGetReal(randnumgen, -1, 0);
      }
      else if( SCIPisGT(scip, activity, rhs))
      {
         *score = 100 * (activity - rhs) / MAX(ABS(rhs),1);
         assert(score>=0);
         *roundup = (*score < heurdata->roundingfrac);
      }
      else if( SCIPisLT(scip, activity, lhs))
      {
         *score = 100 * (lhs - activity) / MAX(ABS(lhs),1);
         assert(score>=0);
         *roundup = (*score < heurdata->roundingfrac);
      }
      else
      {
          assert(FALSE);
      }
   }


   /* free memory */
   SCIPfreeBufferArray(scip, &consvals);
   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}

#define divesetAvailableIndicatordiving NULL

/*
 * heuristic specific interface methods
 */

/** creates the indicatordiving heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurIndicatordiving(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create indicatordiving primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   heur = NULL;


   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecIndicatordiving, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyIndicatordiving) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeIndicatordiving) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitIndicatordiving) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitIndicatordiving) );

   /* create a diveset (this will automatically install some additional parameters for the heuristic)*/
   SCIP_CALL( SCIPcreateDiveset(scip, NULL, heur, HEUR_NAME, DEFAULT_MINRELDEPTH, DEFAULT_MAXRELDEPTH, DEFAULT_MAXLPITERQUOT,
         DEFAULT_MAXDIVEUBQUOT, DEFAULT_MAXDIVEAVGQUOT, DEFAULT_MAXDIVEUBQUOTNOSOL, DEFAULT_MAXDIVEAVGQUOTNOSOL, DEFAULT_LPRESOLVEDOMCHGQUOT,
         DEFAULT_LPSOLVEFREQ, DEFAULT_MAXLPITEROFS, DEFAULT_RANDSEED, DEFAULT_BACKTRACK, DEFAULT_ONLYLPBRANCHCANDS,
         DIVESET_ISPUBLIC, DIVESET_DIVETYPES, divesetGetScoreIndicatordiving, divesetAvailableIndicatordiving) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/" HEUR_NAME "/roundingfrac",
         "in fractional case all fractional below this value are rounded up",
         &heurdata->roundingfrac, FALSE, DEFAULT_ROUNDINGFRAC, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/mode",
         "decides which mode is selected (0: rounding down; 1: rounding up; 2: fractional rounding (default))",
         &heurdata->mode, FALSE, DEFAULT_MODE, 0, 2, NULL, NULL) );

   return SCIP_OKAY;
}
