/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_activecons.c
 * @brief  activecons primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/heur_activecons.h"
#include "scip/scip.h"


#define HEUR_NAME             "activecons"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'q'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPLOOP
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/** primal heuristic data */
struct SCIP_HeurData
{
};


/*
 * Local methods
 */

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}


/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyActivecons)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of activecons primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyActivecons NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeActivecons)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of activecons primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeActivecons NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitActivecons)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of activecons primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitActivecons NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitActivecons)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of activecons primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitActivecons NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolActivecons)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of activecons primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolActivecons NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolActivecons)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of activecons primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolActivecons NULL
#endif

/** main procedure of the Activecons heuristic, creates and solves a sub-SCIP */
static
SCIP_RETCODE applyActivecons(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_RESULT*          result              /**< result data structure */
   )
{
   SCIP* subscip;                                 /* the subproblem created by OFINS */
   SCIP_HASHMAP* varmapfw;                        /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR** vars;                               /* source problem's variables */
   SCIP_VAR** subvars;                            /* subproblem's variables */
//   SCIP_EVENTHDLR* eventhdlr;                     /* event handler for LP events */

   SCIP_Real timelimit;                           /* time limit for OFINS subproblem */
   SCIP_Real memorylimit;                         /* memory limit for OFINS subproblem */

   int nvars;                                     /* number of source problem's variables */
   int i;

   SCIP_SOL** subsols;
   int nsubsols;

   SCIP_Bool valid;
   SCIP_Bool success;
   SCIP_RETCODE retcode;

   SCIP_ROW** rows; /* original scip rows */
   int nrows;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMessage("+---+ Start Activecons heuristic +---+\n");

   /* check whether there is enough time and memory left */
   timelimit = 0.0;
   memorylimit = 0.0;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );

   /* substract the memory already used by the main SCIP and the estimated memory usage of external software */
   if( !SCIPisInfinity(scip, memorylimit) )
   {
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
      memorylimit -= SCIPgetMemExternEstim(scip)/1048576.0;
   }

   /* abort if no time is left or not enough memory to create a copy of SCIP, including external memory usage */
   if( timelimit <= 0.0 || memorylimit <= 2.0*SCIPgetMemExternEstim(scip)/1048576.0 )
   {
      SCIPdebugMessage("-> not enough memory left\n");
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   /* get variable data */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

//   eventhdlr = NULL;

   valid = FALSE;

   /* copy all plugins */
   SCIP_CALL( SCIPcopyPlugins(scip, subscip, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, &valid) );

   SCIP_CALL( SCIPcopyProb(scip, subscip, varmapfw, NULL, TRUE, "activecons") );

   SCIP_CALL( SCIPcopyParamSettings(scip, subscip) );

   SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );

   for( i = 0; i < nvars; i++ )
   {
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
     assert(subvars[i] != NULL);
   }
   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* get the rows and their number */
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* copy all rows to linear constraints */
   for( i = 0; i < nrows; i++ )
   {
      SCIP_CONS* cons;
      SCIP_VAR** consvars;
      SCIP_COL** cols;
      SCIP_Real constant;
      SCIP_Real lhs;
      SCIP_Real rhs;
      SCIP_Real* vals;
      int nnonz;
      int j;

      /* ignore rows that are only locally valid */
      if( SCIProwGetAge(rows[i]) > 0 )
         continue;

      /* get the row's data */
      constant = SCIProwGetConstant(rows[i]);
      lhs = SCIProwGetLhs(rows[i]) - constant;
      rhs = SCIProwGetRhs(rows[i]) - constant;
      vals = SCIProwGetVals(rows[i]);
      nnonz = SCIProwGetNNonz(rows[i]);
      cols = SCIProwGetCols(rows[i]);

      assert(lhs <= rhs);

      /* allocate memory array to be filled with the corresponding subproblem variables */
      SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nnonz) );
      for( j = 0; j < nnonz; j++ )
         consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

      /* create a new linear constraint and add it to the subproblem */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

      /* free temporary memory */
      SCIPfreeBufferArray(subscip, &consvars);
   }

   SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");

   /* create event handler for LP events */
//   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecOfins, NULL) );
//   if( eventhdlr == NULL )
//   {
//      SCIPerrorMessage("event handler for " HEUR_NAME " heuristic not found.\n");
//      return SCIP_PLUGINNOTFOUND;
//   }

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
//   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(subscip, "estimate") != NULL && !SCIPisParamFixed(subscip, "nodeselection/estimate/stdpriority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use mostinference branching */
   if( SCIPfindBranchrule(subscip, "mostinference") != NULL && !SCIPisParamFixed(subscip, "branching/mostinference/priority") )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/mostinference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
//   if( !SCIPisParamFixed(subscip, "conflict/enable") )
//   {
//      SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/enable", FALSE) );
//   }

#ifdef SCIP_DEBUG
   /* for debugging OFINS, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

   /* presolve the subproblem */
   retcode = SCIPpresolve(subscip);

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while presolving subproblem in %s heuristic; sub-SCIP terminated with code <%d>\n", HEUR_NAME, retcode);

      /* free */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }

   SCIPdebugMessage("%s presolved subproblem: %d vars, %d cons\n", HEUR_NAME, SCIPgetNVars(subscip), SCIPgetNConss(subscip));

//   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPtransformProb(subscip) );
//   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   /* solve the subproblem */
//   SCIPdebugMessage("solving subproblem: nstallnodes=%" SCIP_LONGINT_FORMAT ", maxnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->maxnodes);
   retcode = SCIPsolve(subscip);

//   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_LPSOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in Activecons heuristic; sub-SCIP terminated with code <%d>\n", retcode);
   }

   /* print solving statistics of subproblem if we are in SCIP's debug mode */
   SCIPdebug( SCIP_CALL( SCIPprintStatistics(subscip, NULL) ) );

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for( i = 0; i < nsubsols && !success; i++ )
   {
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

   SCIPstatisticPrintf("%s statistic: fixed %6.3f integer variables, needed %6.1f seconds, %" SCIP_LONGINT_FORMAT " nodes, solution %10.4f found at node %" SCIP_LONGINT_FORMAT "\n",
      HEUR_NAME, 0.0, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
      nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 );

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecActivecons)
{  /*lint --e{715}*/
   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic of node was already detected to be infeasible */
   if( nodeinfeasible )
      return SCIP_OKAY;

   /* get heuristic data */
//   heurdata = SCIPheurGetData(heur);
//   assert( heurdata != NULL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
//   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));
//
//   /* reward OFINS if it succeeded often */
//   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
//   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
//   nstallnodes += heurdata->nodesofs;
//
//   /* determine the node limit for the current process */
//   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);
//
//   /* check whether we have enough nodes left to call subproblem solving */
//   if( nstallnodes < heurdata->minnodes )
//   {
//      SCIPdebugMessage("skipping OFINS: nstallnodes=%" SCIP_LONGINT_FORMAT ", minnodes=%" SCIP_LONGINT_FORMAT "\n", nstallnodes, heurdata->minnodes);
//      return SCIP_OKAY;
//   }

   SCIP_CALL( applyActivecons(scip, heur, result) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the activecons primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurActivecons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create activecons primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyActivecons, heurFreeActivecons, heurInitActivecons, heurExitActivecons, heurInitsolActivecons, heurExitsolActivecons, heurExecActivecons,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecActivecons, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyActivecons) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeActivecons) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitActivecons) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitActivecons) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolActivecons) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolActivecons) );
#endif

   /* add activecons primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
