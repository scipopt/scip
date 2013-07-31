/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define SCIP_STATISTIC

/**@file   heur_randrounding.c
 * @brief  rand and fast LP rounding heuristic
 * @author Gregor Hendel
 *
 * The heuristic also tries to round relaxation solutions if available.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_randrounding.h"


#define HEUR_NAME             "randrounding"
#define HEUR_DESC             "fast LP rounding heuristic"
#define HEUR_DISPCHAR         'G'
#define HEUR_PRIORITY         -200
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_DURINGPRICINGLOOP
#define HEUR_USESSUBSCIP      FALSE          /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_ONCEPERNODE   FALSE          /**< should the heuristic only be called once per node? */
#define DEFAULT_RANDSEED      14081986
#define DEFAULT_USESIMPLEROUNDING    FALSE
#define DEFAULT_MAXPROPROUNDS 1
#define DEFAULT_PROPAGATEONLYROOT FALSE

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Longint          lastlp;             /**< last LP number where the heuristic was applied */
   unsigned int          randseed;
   int                   maxproprounds;
   SCIP_Bool             oncepernode;        /**< should the heuristic only be called once per node? */
   SCIP_Bool             usesimplerounding;
   SCIP_Longint          ngeneratedconflicts;
   SCIP_Bool             propagateonlyroot;
};


/*
 * Local methods
 */

/** perform rounding */
static
SCIP_RETCODE performRandRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,
   SCIP_SOL*             sol,                /**< solution to round */
   SCIP_VAR**            cands,              /**< candidate variables */
   SCIP_Real*            candssol,           /**< solutions of candidate variables */
   int                   ncands,             /**< number of candidates */
   SCIP_Bool             propagate,          /**< should the rounding be propagated? */
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   int c;
   SCIP_Bool stored;
   SCIP_VAR** permutedcands;
   SCIP_Bool cutoff;

   assert(heurdata != NULL);

   /* round all roundable fractional columns in the corresponding direction as long as no unroundable column was found */
   if( propagate )
   {
      SCIP_CALL( SCIPstartProbing(scip) );
      SCIPenableVarHistory(scip);
   }

   SCIP_CALL( SCIPduplicateBufferArray(scip, &permutedcands, cands, ncands) );

   assert(permutedcands != NULL);

   SCIPpermuteArray((void **)permutedcands, 0, ncands, &heurdata->randseed);
   cutoff = FALSE;

   for (c = 0; c < ncands && !cutoff; ++c)
   {
      SCIP_VAR* var;
      SCIP_Real oldsolval;
      SCIP_Real newsolval;
      SCIP_Bool mayrounddown;
      SCIP_Bool mayroundup;
      SCIP_Longint ndomreds;
      SCIP_Bool lbadjust;
      SCIP_Bool ubadjust;
      SCIP_Real lb;
      SCIP_Real ub;


      var = permutedcands[c];
      oldsolval = SCIPgetSolVal(scip, sol, var);
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      assert( ! SCIPisFeasIntegral(scip, oldsolval) );
      assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );

      mayrounddown = SCIPvarMayRoundDown(var);
      mayroundup = SCIPvarMayRoundUp(var);

      SCIPdebugMessage("rand rounding heuristic: var <%s>, val=%g, rounddown=%u, roundup=%u\n",
         SCIPvarGetName(var), oldsolval, mayrounddown, mayroundup);

      if( lb > SCIPfeasCeil(scip, oldsolval) + 0.5 || ub < SCIPfeasFloor(scip, oldsolval) - 0.5 )
      {
         cutoff = TRUE;
         break;
      } else if( SCIPisFeasEQ(scip, lb, SCIPfeasCeil(scip, oldsolval)) )
      {
         assert(SCIPisFeasGE(scip, ub, SCIPfeasCeil(scip, oldsolval)));
         newsolval = SCIPfeasCeil(scip, oldsolval);
      }
      else if( SCIPisFeasEQ(scip, ub, SCIPfeasFloor(scip, oldsolval)) )
      {
         assert(SCIPisFeasLE(scip,lb, SCIPfeasFloor(scip, oldsolval)));
         newsolval = SCIPfeasFloor(scip, oldsolval);
      }
      else if( !heurdata->usesimplerounding || !(mayroundup && mayrounddown) )
      {
         SCIP_Real randnumber;

         randnumber = SCIPgetRandomReal(0, 1, &heurdata->randseed);
         if( randnumber <= oldsolval - SCIPfeasFloor(scip, oldsolval) )
            newsolval = SCIPfeasCeil(scip, oldsolval);
         else
            newsolval = SCIPfeasFloor(scip, oldsolval);
      }
      /* choose rounding direction */
      else if( mayrounddown && mayroundup )
      {
         /* we can round in both directions: round in objective function direction */
         if ( SCIPvarGetObj(var) >= 0.0 )
            newsolval = SCIPfeasFloor(scip, oldsolval);
         else
            newsolval = SCIPfeasCeil(scip, oldsolval);
      }
      else if( mayrounddown )
         newsolval = SCIPfeasFloor(scip, oldsolval);
      else
      {
         assert(mayroundup);
         newsolval = SCIPfeasCeil(scip, oldsolval);
      }

      assert(SCIPisFeasLE(scip, lb, newsolval));
      assert(SCIPisFeasGE(scip, ub, newsolval));

      if( propagate )
      {
         lbadjust = SCIPisGT(scip, newsolval, lb);
         ubadjust = SCIPisLT(scip, newsolval, ub);

         if( lbadjust || ubadjust )
         {
            SCIP_CALL( SCIPnewProbingNode(scip) );

            if( lbadjust )
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, newsolval) );
            }
            if( ubadjust )
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, newsolval) );
            }

            heurdata->ngeneratedconflicts -= SCIPgetNConflictConssFoundNode(scip);
            SCIP_CALL( SCIPpropagateProbing(scip, heurdata->maxproprounds, &cutoff, &ndomreds) );
            heurdata->ngeneratedconflicts += SCIPgetNConflictConssFoundNode(scip);
         }
      }
      /* store new solution value */
      SCIP_CALL( SCIPsetSolVal(scip, sol, var, newsolval) );
   }

   if( !cutoff )
   {
      /* check, if rounding was successful */
      if( SCIPallColsInLP(scip) )
      {
         /* check solution for feasibility, and add it to solution store if possible
          * neither integrality nor feasibility of LP rows has to be checked, because all fractional
          * variables were already moved in feasible direction to the next integer
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, FALSE, FALSE, TRUE, &stored) );
      }
      else
      {
         /* if there are variables which are not present in the LP, e.g., for
          * column generation, we need to check their bounds
          */
         SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, TRUE, &stored) );
      }

      if( stored )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("found feasible rounded solution:\n");
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif
         *result = SCIP_FOUNDSOL;
      }
   }

   assert( !propagate || SCIPinProbing(scip) );
   if( propagate )
   {
      SCIP_CALL( SCIPendProbing(scip) );
   }

   SCIPfreeBufferArray(scip, &permutedcands);

   return SCIP_OKAY;
}

/** perform LP-rounding */
static
SCIP_RETCODE performLPRandRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_HEURTIMING       heurtiming,         /**< heuristic timing mask */
   SCIP_Bool             propagate,
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandssol;
   SCIP_Longint nlps;
   int nlpcands;

   /* only call heuristic, if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get fractional variables, that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, &lpcandssol, NULL, &nlpcands, NULL) );

   /* only call heuristic, if LP solution is fractional; except we are called during pricing, in this case we
    * want to detect a (mixed) integer (LP) solution which is primal feasible */
   if ( nlpcands == 0  && heurtiming != SCIP_HEURTIMING_DURINGPRICINGLOOP )
      return SCIP_OKAY;

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert( sol != NULL );

   /* copy the current LP solution to the working solution */
   SCIP_CALL( SCIPlinkLPSol(scip, sol) );

   /* don't call heuristic, if we have already processed the current LP solution */
   nlps = SCIPgetNLPs(scip);
   if( nlps == heurdata->lastlp )
      return SCIP_OKAY;
   heurdata->lastlp = nlps;


   /* perform rand rounding */
   SCIPdebugMessage("executing rand LP-rounding heuristic: %d fractionals\n", nlpcands);
   SCIP_CALL( performRandRounding(scip, heurdata, sol, lpcands, lpcandssol, nlpcands, propagate, result) );

   return SCIP_OKAY;
}

/** perform relaxation solution rounding */
static
SCIP_RETCODE performRelaxRandRounding(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data */
   SCIP_Bool             propagate,
   SCIP_RESULT*          result              /**< pointer to store the result of the heuristic call */
   )
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_VAR** relaxcands;
   SCIP_Real* relaxcandssol;
   int nrelaxcands;
   int nbinvars;
   int nintvars;
   int nvars;
   int v;

   /* do not call heuristic if no relaxation solution is available */
   if ( ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* get variables */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );
   nvars = nbinvars + nintvars; /* consider integral variables (don't have to care for implicit ints) */

   /* get storage */
   SCIP_CALL( SCIPallocBufferArray(scip, &relaxcands, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &relaxcandssol, nvars) );

   /* get fractional variables, that should be integral */
   nrelaxcands = 0;
   for (v = 0; v < nvars; ++v)
   {
      SCIP_Real val;

      val = SCIPgetRelaxSolVal(scip, vars[v]);
      if ( ! SCIPisFeasZero(scip, val) )
      {
         relaxcands[nrelaxcands] = vars[v];
         relaxcandssol[nrelaxcands++] = val;
      }
   }

   /* get the working solution from heuristic's local data */
   sol = heurdata->sol;
   assert( sol != NULL );

   /* copy the current relaxation solution to the working solution */
   SCIP_CALL( SCIPlinkRelaxSol(scip, sol) );

   /* perform rand rounding */
   SCIPdebugMessage("executing rand rounding heuristic on relaxation solution: %d fractionals\n", nrelaxcands);
   SCIP_CALL( performRandRounding(scip, heurdata, sol, relaxcands, relaxcandssol, nrelaxcands, propagate, result) );

   /* free storage */
   SCIPfreeBufferArray(scip, &relaxcands);
   SCIPfreeBufferArray(scip, &relaxcandssol);

   return SCIP_OKAY;
}


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRandrounding)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRandrounding(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIPstatisticMessage("Random Rounding found a total of %"SCIP_LONGINT_FORMAT" conflicts \n", heurdata->ngeneratedconflicts);
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* create heuristic data */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );
   heurdata->lastlp = -1;
   heurdata->randseed = DEFAULT_RANDSEED;
   heurdata->ngeneratedconflicts = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolRandrounding)
{
   SCIP_HEURDATA* heurdata;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   heurdata->lastlp = -1;

   /* change the heuristic's timingmask, if it should be called only once per node */
   if( heurdata->oncepernode )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_AFTERLPNODE);

   return SCIP_OKAY;
}


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolRandrounding)
{
   /* reset the timing mask to its default value */
   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRandrounding) /*lint --e{715}*/
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool propagate;

   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(result != NULL);
   assert(SCIPhasCurrentNodeLP(scip));

   *result = SCIP_DIDNOTRUN;

   /* only call heuristic, if an optimal LP solution is at hand or if relaxation solution is available */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL && ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;

   /* only call heuristic, if the LP objective value is smaller than the cutoff bound */
   if( SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)) )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* don't call heuristic, if we have already processed the current LP solution but no relaxation solution is available */
   if ( SCIPgetNLPs(scip) == heurdata->lastlp && ! SCIPisRelaxSolValid(scip) )
      return SCIP_OKAY;


      *result = SCIP_DIDNOTFIND;

   propagate = !heurdata->propagateonlyroot || SCIPgetDepth(scip) == 0;

   /* try to round LP solution */
   SCIP_CALL( performLPRandRounding(scip, heurdata, heurtiming, propagate, result) );

   /* try to round relaxation solution */
   SCIP_CALL( performRelaxRandRounding(scip, heurdata, propagate, result) );

   return SCIP_OKAY;
}

/*
 * heuristic specific interface methods
 */

/** creates the rand rounding heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRandrounding(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecRandrounding, heurdata) );
   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyRandrounding) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitRandrounding) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitRandrounding) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolRandrounding) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolRandrounding) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeRandrounding) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/oncepernode",
         "should the heuristic only be called once per node?",
         &heurdata->oncepernode, TRUE, DEFAULT_ONCEPERNODE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/usesimplerounding", "bla?",
         &heurdata->usesimplerounding, TRUE, DEFAULT_USESIMPLEROUNDING, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/propagateonlyroot", "bla?",
            &heurdata->propagateonlyroot, TRUE, DEFAULT_PROPAGATEONLYROOT, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxproprounds", "bla", &heurdata->maxproprounds, TRUE, DEFAULT_MAXPROPROUNDS,
         -1, INT_MAX, NULL, NULL) );
   return SCIP_OKAY;
}
