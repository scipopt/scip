/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_hailmary.c
 * @brief  hailmary primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_hailmary.h"


#define HEUR_NAME             "hailmary"
#define HEUR_DESC             "heuristic trying to solve the problem without objective"
#define HEUR_DISPCHAR         'H'
#define HEUR_PRIORITY         -1000000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE | SCIP_HEURTIMING_BEFOREPRESOL
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */

/* event handler properties */
#define EVENTHDLR_NAME         "Hailmary"
#define EVENTHDLR_DESC         "LP event handler for "HEUR_NAME" heuristic"

/* default values for hailmary-specific plugins */
#define DEFAULT_MAXNODES      1000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which hailmary should at least improve the incumbent          */
#define DEFAULT_MINNODES      100LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_MAXLPITERS    5000LL   /* maximum number of LP iterations to be performed in the subproblem */
#define DEFAULT_NODESOFS      100LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_ADDALLSOLS   FALSE      /* should all subproblem solutions be added to the original SCIP?       */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          maxlpiters;         /**< maximum number of LP iterations to be performed in the subproblem   */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by hailmary in earlier calls                     */
   SCIP_Real             minimprove;         /**< factor by which hailmary should at least improve the incumbent      */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   SCIP_Bool             addallsols;         /**< should all subproblem solutions be added to the original SCIP?      */
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
   SCIP_HEUR*            heur,               /**< hailmary heuristic structure                            */
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

/* ---------------- Callback methods of event handler ---------------- */

/* exec the event handler
 *
 * we interrupt the solution process
 */
static
SCIP_DECL_EVENTEXEC(eventExecHailmary)
{
   SCIP_HEURDATA* heurdata;

   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_NODESOLVED);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata != NULL);

   /* interrupt solution process of sub-SCIP */
   if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_ITERLIMIT || SCIPgetNLPIterations(scip) >= heurdata->maxlpiters )
   {
      SCIP_CALL( SCIPinterruptSolve(scip) );
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyHailmary)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurHailmary(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeHailmary)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitHailmary)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitHailmary NULL


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolHailmary NULL


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolHailmary NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecHailmary)
{  /*lint --e{715}*/
  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Longint nnodes;                 /* number of stalling nodes for the subproblem */

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward hailmary if it succeeded often */
   nnodes = (SCIP_Longint)(nnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nnodes -= heurdata->usednodes;
   nnodes = MIN(nnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping hailmary: nnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   /* do not run hailmary, if the problem does not have an objective function anyway */
   if( SCIPgetNObjVars(scip) == 0 )
   {
      SCIPdebugMessage("skipping hailmary: pure feasibility problem anyway\n");
      return SCIP_OKAY;
   }

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPapplyHailmary(scip, heur, result, heurdata->minimprove, nnodes) );

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */


/** main procedure of the hailmary heuristic, creates and solves a sub-SCIP */
SCIP_RETCODE SCIPapplyHailmary(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minimprove,         /**< factor by which hailmary should at least improve the incumbent      */
   SCIP_Longint          nnodes              /**< node limit for the subproblem                                       */
   )
{
   SCIP*                 subscip;            /* the subproblem created by hailmary              */
   SCIP_HASHMAP*         varmapfw;           /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR**            vars;               /* original problem's variables                    */
   SCIP_VAR**            subvars;            /* subproblem's variables                          */
   SCIP_HEURDATA*        heurdata;           /* heuristic's private data structure              */
   SCIP_EVENTHDLR*       eventhdlr;          /* event handler for LP events                     */

   SCIP_Real cutoff;                         /* objective cutoff for the subproblem             */
   SCIP_Real timelimit;                      /* time limit for hailmary subproblem              */
   SCIP_Real memorylimit;                    /* memory limit for hailmary subproblem            */

   int nvars;                                /* number of original problem's variables          */
   int i;

   SCIP_Bool success;
   SCIP_Bool valid;
   SCIP_RETCODE retcode;
   SCIP_SOL** subsols;
   int nsubsols;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   assert(nnodes >= 0);
   assert(0.0 <= minimprove && minimprove <= 1.0);

   *result = SCIP_DIDNOTRUN;

   /* only call feaspump once at the root */
   if( SCIPgetDepth(scip) <= 0 && SCIPheurGetNCalls(heur) > 0 )
      return SCIP_OKAY;

   /* only call the heuristic if we do not have an incumbent  */
   if( SCIPgetNSolsFound(scip) >  0 )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

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
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   /* different methods to create sub-problem: either copy LP relaxation or the CIP with all constraints */
   valid = FALSE;

   /* copy complete SCIP instance */
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "hailmary", TRUE, FALSE, &valid) );
   SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");

   /* create event handler for LP events */
   eventhdlr = NULL;
   SCIP_CALL( SCIPincludeEventhdlrBasic(subscip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecHailmary, NULL) );
   if( eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for "HEUR_NAME" heuristic not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* get variable image and change to 0.0 in sub-SCIP */
   for( i = 0; i < nvars; i++ )
   {
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);
      SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );
   }

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", 1) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable expensive techniques that merely work on the dual bound */

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );
   SCIP_CALL( SCIPsetIntParam(subscip, "presolving/maxrounds", 50) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(scip, "dfs") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/dfs/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(scip, "inference") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/leastinf/priority", INT_MAX/4) );
   }

   /* disable feaspump and fracdiving */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/feaspump/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/fracdiving/freq", -1) );

   /* restrict LP iterations */
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/iterlim", 2*heurdata->maxlpiters / MAX(1,nnodes)) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "lp/rootiterlim", heurdata->maxlpiters) );

#ifdef SCIP_DEBUG
   /* for debugging hailmary, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      cutoff = SCIPinfinity(scip);
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound ( scip ) >= 0 )
            cutoff = ( 1 - minimprove ) * SCIPgetUpperbound ( scip );
         else
            cutoff = ( 1 + minimprove ) * SCIPgetUpperbound ( scip );
      }
      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   /* catch LP events of sub-SCIP */
   SCIP_CALL( SCIPtransformProb(subscip) );
   SCIP_CALL( SCIPcatchEvent(subscip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, NULL) );

   SCIPdebugMessage("solving subproblem: nnodes=%"SCIP_LONGINT_FORMAT"\n", nnodes);
   retcode = SCIPsolve(subscip);

   /* drop LP events of sub-SCIP */
   SCIP_CALL( SCIPdropEvent(subscip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, (SCIP_EVENTDATA*) heurdata, -1) );

   /* errors in solving the subproblem should not kill the overall solving process;
    * hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage(scip, "Error while solving subproblem in hailmary heuristic; sub-SCIP terminated with code <%d>\n",retcode);
   }

   /* check, whether a solution was found;
    * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
    */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for( i = 0; i < nsubsols && (!success || heurdata->addallsols); ++i )
   {
      SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      if( success )
         *result = SCIP_FOUNDSOL;
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

   SCIPstatistic(
      SCIPstatisticMessage("hailmary statistic: fixed %6.3f integer variables, %6.3f all variables, needed %6.1f seconds, %"SCIP_LONGINT_FORMAT" nodes, solution %10.4f found at node %"SCIP_LONGINT_FORMAT"\n",
         intfixingrate, allfixingrate, SCIPgetSolvingTime(subscip), SCIPgetNNodes(subscip), success ? SCIPgetPrimalbound(scip) : SCIPinfinity(scip),
         nsubsols > 0 ? SCIPsolGetNodenum(SCIPgetBestSol(subscip)) : -1 )
      );

   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/** creates the hailmary primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurHailmary(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyHailmary,
         heurFreeHailmary, heurInitHailmary, heurExitHailmary,
         heurInitsolHailmary, heurExitsolHailmary, heurExecHailmary,
         heurdata) );

   /* add hailmary primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes, TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxlpiters",
         "maximum number of LP iterations to be performed in the subproblem",
         &heurdata->maxlpiters, TRUE, DEFAULT_MAXLPITERS, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which hailmary should at least improve the incumbent",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/addallsols",
         "should all subproblem solutions be added to the original SCIP?",
         &heurdata->addallsols, TRUE, DEFAULT_ADDALLSOLS, NULL, NULL) );

   return SCIP_OKAY;
}
