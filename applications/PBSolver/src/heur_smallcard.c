/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/**@file   heur_smallcard.c
 * @ingroup HEURISTICS
 * @brief  smallcard heuristic
 * @author Marc Pfetsch
 *
 * Setup sub-SCIP and check whether there exists a solution with few indicator constraints not fulfilled.
 *
 * Based on sepa_rapidlearning.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#ifndef NDEBUG
#include <string.h>
#endif

#include "heur_smallcard.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_var.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"

#define HEUR_NAME             "smallcard"
#define HEUR_DESC             "search for small cardinality indicator solutions"
#define HEUR_DISPCHAR         '8'
#define HEUR_PRIORITY         100000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */


#define DEFAULT_MAXCARD               1 /**< maximum size of non-fulfilled indicator constraints */
#define DEFAULT_MAXVALUE            1.0 /**< maximum objective value allowed in subproblem */
#define DEFAULT_USECARD            TRUE /**< Use maximum size constraint (otherwise use maximum value)? */
#define DEFAULT_MAXNVARS          50000 /**< maximum problem size (variables) for which heuristic will be called */
#define DEFAULT_MAXNCONSS         50000 /**< maximum problem size (constraints) for which heuristic will be called */

#define DEFAULT_MINNODES           5000 /**< minimum number of nodes considered in heuristic run */
#define DEFAULT_MAXNODES          10000 /**< maximum number of nodes considered in heuristic run */

#define DEFAULT_CONTVARS          FALSE /**< Should the heuristic be be applied when there are continuous variables? */
#define DEFAULT_COPYCUTS          FALSE /**< Should all active cuts from the cutpool of the original scip be copied to constraints of the subscip? */
#define DEFAULT_SOLVELP           FALSE /**< Should intermediate lp's be solved? */

/*
 * Data structures
 */

/** heuristic data */
struct SCIP_HeurData
{
   int                   maxcard;            /**< maximum size of non-fulfilled indicator constraints */
   SCIP_Real             maxvalue;           /**< maximum objective value allowed in subproblem */
   SCIP_Bool             usecard;            /**< Use maximum size constraint (otherwise use maximum value)? */
   int                   maxnvars;           /**< maximum problem size (variables) for which heuristic will be called   */
   int                   maxnconss;          /**< maximum problem size (constraints) for which heuristic will be called */
   int                   minnodes;           /**< minimum number of nodes considered in heuristic run */
   int                   maxnodes;           /**< maximum number of nodes considered in heuristic run */
   SCIP_Bool             contvars;           /**< Should heuristic be applied when there are continuous variables? */
   SCIP_Bool             copycuts;           /**< Should all active cuts from cutpool be copied to constraints in subproblem? */
   SCIP_Bool             solvelp;            /**< Should intermediate lp's be solved? */
   SCIP_Bool             applied;            /**< whether the heuristic has been applied during this run */
};


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem */
   SCIP_HEUR*            heur,               /**< heuristic structure */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem */
   SCIP_SOL*             subsol,             /**< solution of the subproblem */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
)
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
   int nvars;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( heur != NULL );
   assert( subsol != NULL );
   assert( success != NULL );

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* subSCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert( nvars <= SCIPgetNOrigVars(subscip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* check feasibility of new solution and pass it to trysol heuristic */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/*
 * Callback methods of heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySmallcard)
{
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSmallcard(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSmallcard)
{   /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of heursitic (called after problem was transformed) */
#define heurInitSmallcard NULL

/** deinitialization method of heursitic (called before transformed problem is freed) */
#define heurExitSmallcard NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolSmallcard NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolSmallcard)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );
   heurdata->applied = FALSE;

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSmallcard)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* indicatorconshdlr;
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_RETCODE retcode;
   SCIP_Bool success;
   SCIP_Bool soladded;
   SCIP_CONS** indconss;
   SCIP_VAR** subvars;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
#ifdef SCIP_DISABLED_CODE
   SCIP_SOL* bestsol;
#endif
   SCIP* subscip;
   int ndiscvars;
   int nindconss;
   int nvars;
#ifdef SCIP_DISABLED_CODE
   int sum;
   int restarts;
   int restartnum;
#endif
   int i;

   SCIP_Longint nodelimit;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;

   soladded = FALSE;

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   ndiscvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip);

   /* only run when still not fixed discrete variables exists */
   if ( ndiscvars == 0 )
      return SCIP_OKAY;

   /* get heursitic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* only run for integer programs */
   if ( !heurdata->contvars && ndiscvars != SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* do not run if pricers are present */
   if ( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* do not go on if the probelm is too big */
   if ( SCIPgetNVars(scip) > heurdata->maxnvars || SCIPgetNConss(scip) > heurdata->maxnconss )
      return SCIP_OKAY;

   /* exit if heuristic has been run */
   if ( heurdata->applied )
      return SCIP_OKAY;
   heurdata->applied = TRUE;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if ( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if ( timelimit < 10.0 || memorylimit <= 0.0 )
      return SCIP_OKAY;

   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* get indicator constraint handler */
   indicatorconshdlr = SCIPfindConshdlr(scip, "indicator");

   /* exit if indicator constraint handler is not present */
   if ( indicatorconshdlr == NULL )
      return SCIP_OKAY;

   /* get data */
   nindconss = SCIPconshdlrGetNConss(indicatorconshdlr);
   indconss = SCIPconshdlrGetConss(indicatorconshdlr);

   /* exit if there are not indicator constraints */
   if ( nindconss == 0 )
   {
      SCIPdebugMessage("No indicator constraints present - skipping small cardinality heuristic.\n");
      return SCIP_OKAY;
   }

   SCIPdebugMessage("Small cardinality heuristic ...\n");
   *result = SCIP_DIDNOTFIND;

   /* initializing the subproblem */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcMultihashSize(5 * nvars)) );
   success = FALSE;

   /* copy the subproblem */
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "smallcard", FALSE, FALSE, FALSE, FALSE, &success) );

   if ( heurdata->copycuts )
   {
      /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, FALSE, NULL) );
   }

   /* set up variables of subproblem */
   for (i = 0; i < nvars; ++i)
   {
      subvars[i] = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmapfw, vars[i]);
      assert( subvars[i] != NULL );
   }

   SCIPdebugMessage("Copying SCIP was%s successful.\n", success ? "" : " not");

   /* add bounding constraint */
   if ( heurdata->usecard )
   {
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "boundcard", 0, NULL, NULL, -SCIPinfinity(subscip), (SCIP_Real)heurdata->maxcard,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "boundval", 0, NULL, NULL, -SCIPinfinity(subscip), heurdata->maxvalue,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   }

#ifdef SCIP_DISABLED_CODE /* get best solution in original problem and test whether it is feasible for the subproblem */
   bestsol = SCIPgetBestSol(scip);
   sum = 0;
#endif

   /* set up constraint */
   for (i = 0; i < nindconss; ++i)
   {
      SCIP_Bool negated;
      SCIP_VAR* binvar;
      SCIP_VAR* var;

      negated = FALSE;
      binvar = SCIPgetBinaryVarIndicator(indconss[i]);
      assert( binvar != NULL );

#ifdef SCIP_DISABLED_CODE
      assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, bestsol, binvar)) );
      if ( SCIPgetSolVal(scip, bestsol, binvar) > 0.5 )
         ++sum;
#endif

      if ( SCIPvarIsNegated(binvar) )
      {
         binvar = SCIPvarGetNegatedVar(binvar);
         negated = TRUE;
      }

      var = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmapfw, binvar);
      assert( var != NULL );
      if ( ! negated )
      {
         SCIP_CALL( SCIPgetNegatedVar(subscip, var, &var) );
      }

      if ( heurdata->usecard )
      {
         SCIP_CALL( SCIPaddCoefLinear(subscip, cons, var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPaddCoefLinear(subscip, cons, var, SCIPvarGetObj(var)) );
      }
   }
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

   SCIPhashmapFree(&varmapfw);

   if ( ! heurdata->solvelp )
   {
      /* mimic an FD solver: DFS, no LP solving, 1-FUIP instead of all-FUIP */
      SCIP_CALL( SCIPsetIntParam(subscip, "lp/solvefreq", -1) );
      SCIP_CALL( SCIPsetIntParam(subscip, "conflict/fuiplevels", 1) );
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/dfs/stdpriority", INT_MAX/4) );
      SCIP_CALL( SCIPsetIntParam(subscip, "propagating/pseudoobj/freq", -1) );
      SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/disableenfops", TRUE) );

      /* use inference branching */
      SCIP_CALL( SCIPsetBoolParam(subscip, "branching/inference/useweightedsum", FALSE) );
   }

   /* only create short conflicts */
   SCIP_CALL( SCIPsetRealParam(subscip, "conflict/maxvarsfac", 0.05) );

   /* set limits for the subproblem */
   nodelimit = SCIPgetNLPIterations(scip);
   nodelimit = MAX(heurdata->minnodes, nodelimit);
   nodelimit = MIN(heurdata->maxnodes, nodelimit);

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
#ifdef SCIP_DISABLED_CODE
   restarts = 0;
   restartnum = 1000;

   SCIP_CALL( SCIPsetIntParam(subscip, "limits/restarts", restarts) );
   SCIP_CALL( SCIPsetIntParam(subscip, "conflict/restartnum", restartnum) );
#endif

   /* forbid recursive call of heuristics and heuristics solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use only fast heuristics (it seems that heuristics cannot use the bouding constraint) */
   SCIP_CALL( SCIPsetHeuristics(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifndef SCIP_OUTPUT
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "smallcard.lp", "lp", FALSE) );
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "smallcard.cip", "cip", FALSE) );
#endif

   SCIPdebugMessage("Subproblem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
      SCIPgetNOrigVars(subscip), SCIPgetNOrigBinVars(subscip), SCIPgetNOrigIntVars(subscip), SCIPgetNOrigImplVars(subscip),
      SCIPgetNOrigContVars(subscip), SCIPgetNOrigConss(subscip));

   SCIPdebugMessage("Start solving subproblem ...\n");

   /* redundant constraints might be eliminated in presolving */
   SCIP_CALL( SCIPpresolve(subscip));

#ifdef SCIP_DISABLED_CODE
   /* if the best solution in the original problem is feasible for the subproblem */
   if ( sum <= 0 )
   {
      SCIP_SOL* sol;
      /* create solution for the subproblem */
      SCIP_CALL( SCIPcreateSol(subscip, &sol, NULL) );

      for (i = 0; i < nvars; ++i)
         SCIP_CALL( SCIPsetSolVal(subscip, sol, subvars[i], SCIPgetSolVal(scip, bestsol, vars[i])) );

      SCIP_CALL( SCIPtrySolFree(subscip, &sol, TRUE, TRUE, TRUE, TRUE, &success) );
      SCIPdebugMessage("Copied feasible solution of value %f from original problem to subproblem.\n", SCIPsolGetOrigObj(bestsol));
   }
#endif

   /* solve the subproblem */
   retcode = SCIPsolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process Hence, the return
    * code is catched and a warning is printed, only in debug mode, SCIP will stop.
    */
   if ( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in smallcard heuristic; subSCIP terminated with code <%d>\n",retcode);
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

   /* check, whether a solution was found */
   if ( SCIPgetNSols(subscip) > 0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found; due to numerics, it might happen that not all
       * solutions are feasible -> try all solutions until was declared to be feasible
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      soladded = FALSE;

      /* sequentially add solutions to trysol heuristic */
      for (i = 0; i < nsubsols && !soladded; ++i)
      {
	 SCIPdebugMessage("Try to create new solution by copying subscip solution.\n");
#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPprintSol(subscip, subsols[i], NULL, FALSE) );
#endif
         SCIP_CALL( createNewSol(scip, subscip, heur, subvars, subsols[i], &soladded) );
      }
   }

   SCIPdebugMessage("Smallcard %s primal solution.\n", soladded ? "found" : "found no");

   /* change result pointer */
   if ( soladded )
      *result = SCIP_FOUNDSOL;

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the twoopt primal heuristic for indicators and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSmallcard(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create smallcard heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->applied = FALSE;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopySmallcard, heurFreeSmallcard, heurInitSmallcard, heurExitSmallcard,
         heurInitsolSmallcard, heurExitsolSmallcard, heurExecSmallcard,
         heurdata) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxcard",
         "maximum size of non-fulfilled indicator constraints",
         &heurdata->maxcard, TRUE, DEFAULT_MAXCARD, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/maxvalue",
         "maximum objective value allowed in subproblem",
         &heurdata->maxvalue, TRUE, DEFAULT_MAXVALUE, -SCIP_REAL_MAX, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/usecard",
         "Use maximum size constraint (otherwise use maximum value)?",
         &heurdata->usecard, TRUE, DEFAULT_USECARD, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnvars",
         "maximum problem size (variables) for which heuristic will be called",
         &heurdata->maxnvars, TRUE, DEFAULT_MAXNVARS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnconss",
         "maximum problem size (constraints) for which heuristic will be called",
         &heurdata->maxnconss, TRUE, DEFAULT_MAXNCONSS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes considered in heuristic run",
         &heurdata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes considered in heuristic run",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/contvars",
         "Should the heuristic be applied when there are continuous variables?",
         &heurdata->contvars, TRUE, DEFAULT_CONTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "Should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/solvelp",
         "Should intermediate lp's be solved?",
         &heurdata->solvelp, TRUE, DEFAULT_SOLVELP, NULL, NULL) );

   return SCIP_OKAY;
}
