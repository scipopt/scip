/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   concsolver_scip.c
 * @ingroup PARALLEL
 * @brief  implementation of concurrent solver interface for SCIP
 * @author Leona Gottwald
 * @author Marc Pfetsch
 */

/* activate the define below for a feasibility check of the solutions transferred to the main SCIP. */
/* #define SCIP_CHECK_MAINSCIP_SOLUTION */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/boundstore.h"
#include "scip/concsolver.h"
#include "scip/concsolver_scip.h"
#include "scip/concurrent.h"
#include "scip/heur_sync.h"
#include "scip/pub_disp.h"
#include "scip/pub_event.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_paramset.h"
#include "scip/pub_sol.h"
#include "scip/pub_var.h"
#include "scip/scip_concurrent.h"
#include "scip/scip_copy.h"
#include "scip/scip_event.h"
#include "scip/scip_general.h"
#include "scip/scip_heur.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solve.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_timing.h"
#include "scip/syncstore.h"

/* event handler for synchronization */
#define EVENTHDLR_NAME         "sync"
#define EVENTHDLR_DESC         "event handler for synchronization of concurrent scip solvers"

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   int             filterpos;
   int             filterpossol;
   SCIP_CONCSOLVER* concsolver;        /**< the concurrent solver this event handler belongs to */
};

/** data for a concurrent solver */
struct SCIP_ConcSolverData
{
   SCIP*                 solverscip;         /**< the concurrent solvers private SCIP data structure */
   SCIP_VAR**            vars;               /**< array of variables in the order of the main SCIP's variable array */
   int                   nvars;              /**< number of variables in the above arrays */
};

/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeSync)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, &eventhdlrdata);

   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitSync)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SYNCSTORE*  syncstore;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);

   if( eventhdlrdata->filterpos < 0 && SCIPsyncstoreIsInitialized(syncstore) )
   {
      /* notify SCIP that your event handler wants to react on synchronization events */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SYNC, eventhdlr, NULL, &eventhdlrdata->filterpos) );
   }

   /* in opportunistic mode new incumbent solutions are published in the solution pool
    * immediately instead of waiting for the next synchronization point; the event is also
    * caught when incumbents should be printed, which is supported in both parallel modes
    */
   if( eventhdlrdata->filterpossol < 0 && SCIPsyncstoreIsInitialized(syncstore) )
   {
      SCIP_Bool solpool;
      SCIP_Bool printincumbents;

      SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/solpool", &solpool) );
      SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/printincumbents", &printincumbents) );

      if( (SCIPsyncstoreGetMode(syncstore) == SCIP_PARA_OPPORTUNISTIC && solpool) || printincumbents )
      {
         SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, &eventhdlrdata->filterpossol) );
      }
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitSync)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* notify SCIP that your event handler wants to drop the event type synchronization found */
   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SYNC, eventhdlr, NULL, eventhdlrdata->filterpos) );
      eventhdlrdata->filterpos = -1;
   }

   if( eventhdlrdata->filterpossol >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, eventhdlrdata->filterpossol) );
      eventhdlrdata->filterpossol = -1;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSync)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(event != NULL);
   assert(scip != NULL);

   SCIP_STRINGEQ( SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME, SCIP_INVALIDCALL );

   /* communicate the objective value of new incumbent solutions to the other concurrent
    * solvers immediately; if the solution pool is enabled the solution values are published
    * there as well, otherwise they are exchanged at the synchronization points
    */
   if( SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND )
   {
      SCIP_EVENTHDLRDATA* eventhdlrdata;
      SCIP_SYNCSTORE* syncstore;
      SCIP_SOL* sol;
      SCIP_Real minobj;
      SCIP_Real* solvals = NULL;
      int nsolvals = 0;
      SCIP_Bool published = FALSE;

      eventhdlrdata = (SCIP_EVENTHDLRDATA*)SCIPeventhdlrGetData(eventhdlr);
      syncstore = SCIPgetSyncstore(scip);
      assert(syncstore != NULL);

      sol = SCIPeventGetSol(event);
      minobj = SCIPgetSolOrigObj(scip, sol);
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
         minobj = -minobj;

      /* extract the solution values only if this incumbent improves on the current global best;
       * the global best objective only decreases, so a solution filtered out here can never be
       * the one that wins the locked comparison inside SCIPsyncstoreUpdateBestMinObj
       */
      if( SCIPsyncstoreSolPoolEnabled(syncstore) )
      {
         SCIP_Real bestminobj;

         bestminobj = SCIPsyncstoreGetBestMinObj(syncstore);

         if( minobj < bestminobj )
         {
            SCIP_CONCSOLVERDATA* data;

            data = SCIPconcsolverGetData(eventhdlrdata->concsolver);
            assert(data != NULL);

            SCIP_CALL( SCIPallocBufferArray(scip, &solvals, data->nvars) );
            SCIP_CALL( SCIPgetSolVals(scip, sol, data->nvars, data->vars, solvals) );
            nsolvals = data->nvars;
         }
      }

      SCIPsyncstoreUpdateBestMinObj(syncstore, minobj, SCIPconcsolverGetName(eventhdlrdata->concsolver),
         SCIPconcsolverGetIdx(eventhdlrdata->concsolver), solvals, nsolvals, &published);

      /* the solution was published in the pool and will be drained by the other solvers, so count
       * it as shared here; in solution-pool mode no solutions are written at the synchronization
       * points, so this is the only place the shared-solution statistic is updated
       */
      if( published )
         SCIPconcsolverAddNSolsShared(eventhdlrdata->concsolver, 1LL);

      SCIPfreeBufferArrayNull(scip, &solvals);

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsynchronize(scip) );

   return SCIP_OKAY;
}


/** includes event handler for synchronization found */
static
SCIP_RETCODE includeEventHdlrSync(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONCSOLVER*      concsolver          /**< the concurrent solver this event handler belongs to */
   )
{
   SCIP_EVENTHDLR*     eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventhdlrdata) );
   eventhdlrdata->filterpos = -1;
   eventhdlrdata->filterpossol = -1;
   eventhdlrdata->concsolver = concsolver;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSync, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSync) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSync) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitSync) );

   return SCIP_OKAY;
}

/** data for a concurrent solver type */
struct SCIP_ConcSolverTypeData
{
   SCIP_Bool             loademphasis;       /**< should emphasis settings be loaded when creating an instance of this concurrent solver */
   SCIP_PARAMEMPHASIS    emphasis;           /**< parameter emphasis that will be loaded if loademphasis is true */
};

/** Disable dual reductions that might cut off optimal solutions. Although they keep at least
 *  one optimal solution intact, communicating these bounds may cut off all optimal solutions,
 *  if different optimal solutions were kept in different concurrent solvers. */
static
SCIP_RETCODE disableConflictingDualReductions(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool commvarbnds;

   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/commvarbnds", &commvarbnds) );

   if( !commvarbnds )
      return SCIP_OKAY;

   SCIP_CALL( SCIPsetBoolParam(scip, "misc/allowstrongdualreds", FALSE) );

   return SCIP_OKAY;
}

/** sets the child selection rule based on the index of the concurrent solver */
static
SCIP_RETCODE setChildSelRule(
   SCIP_CONCSOLVER*      concsolver          /**< the concurrent solver */
   )
{
   SCIP_CONCSOLVERDATA* data;
   static const char childsel[] = { 'h', 'i', 'p', 'r', 'l', 'd', 'u' };

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   SCIP_CALL( SCIPsetCharParam(data->solverscip, "nodeselection/childsel", childsel[SCIPconcsolverGetIdx(concsolver) % 7]) );

   return SCIP_OKAY;
}

/** number of configurations in the built-in racing parameter portfolio */
#define NRACINGCONFIGS 10

/** temporarily fixes the parameters of a concurrent solver's SCIP that must keep their current
 *  values while emphasis or portfolio settings are applied, i.e., the limit, numerics, memory,
 *  and synchronization parameters; the returned buffer array must be released with
 *  unfixProtectedParams()
 */
static
SCIP_RETCODE fixProtectedParams(
   SCIP*                 solverscip,         /**< the concurrent solver's SCIP datastructure */
   SCIP_PARAM***         fixedparams,        /**< buffer to store the array of fixed parameters */
   int*                  nfixedparams        /**< buffer to store the number of fixed parameters */
   )
{
   SCIP_PARAM** params;
   int nparams;
   int i;

   params = SCIPgetParams(solverscip);
   nparams = SCIPgetNParams(solverscip);
   SCIP_CALL( SCIPallocBufferArray(solverscip, fixedparams, nparams) );
   *nfixedparams = 0;

   for( i = 0; i < nparams; ++i )
   {
      const char* paramname;

      paramname = SCIPparamGetName(params[i]);

      if( strncmp(paramname, "limits/", 7) == 0 ||
          strncmp(paramname, "numerics/", 9) == 0 ||
          strncmp(paramname, "memory/", 7) == 0 ||
          strncmp(paramname, "concurrent/sync/", 16) == 0 ||
          strncmp(paramname, "heuristics/sync/", 16) == 0 ||
          strncmp(paramname, "propagating/sync/", 17) == 0 )
      {
         (*fixedparams)[(*nfixedparams)++] = params[i];
         SCIP_CALL( SCIPfixParam(solverscip, paramname) );
      }
   }

   return SCIP_OKAY;
}

/** unfixes the parameters that were fixed by fixProtectedParams() and releases the buffer array */
static
SCIP_RETCODE unfixProtectedParams(
   SCIP*                 solverscip,         /**< the concurrent solver's SCIP datastructure */
   SCIP_PARAM***         fixedparams,        /**< array of fixed parameters */
   int                   nfixedparams        /**< number of fixed parameters */
   )
{
   int i;

   for( i = 0; i < nfixedparams; ++i )
      SCIP_CALL( SCIPunfixParam(solverscip, SCIPparamGetName((*fixedparams)[i])) );

   SCIPfreeBufferArray(solverscip, fixedparams);

   return SCIP_OKAY;
}

/** applies one configuration of the built-in racing parameter portfolio to a concurrent solver's SCIP */
static
SCIP_RETCODE applyRacingSettings(
   SCIP*                 solverscip,         /**< the concurrent solver's SCIP datastructure */
   int                   config              /**< index of the racing configuration to apply */
   )
{
   assert(solverscip != NULL);
   assert(config >= 0 && config < NRACINGCONFIGS);

   switch( config )
   {
   /* default settings; differs from a sequential SCIP only by the racing seed */
   case 0:
      break;

   /* feasibility emphasis; used twice since the different racing seeds make the search distinct */
   case 1:
      SCIP_CALL( SCIPsetBoolParam(solverscip, "symmetries/enabled", FALSE) );
      SCIP_CALL( SCIPsetEmphasis(solverscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );
      break;
   case 2:
      SCIP_CALL( SCIPsetEmphasis(solverscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );
      break;

   /* no separation with fast presolving; trades bound quality for node throughput */
   case 3:
      SCIP_CALL( SCIPsetSeparating(solverscip, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPsetPresolving(solverscip, SCIP_PARAMSETTING_FAST, TRUE) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "symmetries/enabled", FALSE) );
      break;

   /* early feasibility jump; runs the feasibility jump heuristic before presolving to find feasible
    * solutions in the first seconds
    */
   case 4:
      SCIP_CALL( SCIPsetBoolParam(solverscip, "heuristics/feasjump/beforepresol", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "symmetries/enabled", FALSE) );
      break;

   /* cpsolver search with a geometric restart schedule and feasibility jump before presolving; never
    * solves an LP, so it cannot stall on instances with hard root LPs
    */
   case 5:
      SCIP_CALL( SCIPsetEmphasis(solverscip, SCIP_PARAMEMPHASIS_CPSOLVER, TRUE) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "conflict/restartnum", 100) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "conflict/restartfac", 1.5) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "presolving/maxrestarts", 100) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "heuristics/feasjump/beforepresol", TRUE) );
      break;

   /* large neighborhood search improver; fires the LNS heuristics aggressively and immediately on new
    * incumbents, so solutions from the other solvers are polished as soon as they arrive; separation is
    * throttled to leave the CPU to the LNS sub-MIPs
    */
   case 6:
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/rins/freq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/rins/nwaitingnodes", 0) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "heuristics/rins/minimprove", 0.001) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "heuristics/rins/nodesquot", 0.5) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/rins/maxnodes", 10000) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/crossover/freq", 10) );
      SCIP_CALL( SCIPsetLongintParam(solverscip, "heuristics/crossover/nwaitingnodes", 0LL) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "heuristics/crossover/dontwaitatroot", TRUE) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/crossover/nusedsols", 4) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "heuristics/crossover/minimprove", 0.001) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/dins/freq", 15) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/localbranching/freq", 25) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/mutation/freq", 25) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/trustregion/freq", 25) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/rens/freq", 15) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "heuristics/gins/freq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/maxrounds", 1) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/maxroundsroot", 5) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "symmetries/enabled", FALSE) );
      break;

   /* no presolving; starts the search immediately on the original formulation, providing early feasible
    * solutions on instances where presolving takes long
    */
   case 7:
      SCIP_CALL( SCIPsetPresolving(solverscip, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "symmetries/enabled", FALSE) );
      break;

   /* constraint programming style search; no LP relaxation, aggressive conflict analysis, depth first search */
   case 8:
      SCIP_CALL( SCIPsetEmphasis(solverscip, SCIP_PARAMEMPHASIS_CPSOLVER, TRUE) );
      break;

   /* intensified separation and reliability branching; strengthens the separators and presolvers that
    * are active by default (more rounds, more cuts, deeper probing, larger strong branching budgets) to
    * drive the dual bound
    */
   case 9:
      SCIP_CALL( SCIPsetRealParam(solverscip, "branching/relpscost/sbiterquot", 1.0) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "branching/relpscost/maxreliable", 10.0) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "presolving/restartfac", 0.0125) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "presolving/restartminred", 0.06) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/maxroundsrootsubrun", 5) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/maxaddrounds", 5) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/maxcutsroot", 5000) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/gomory/maxroundsroot", 15) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/gomory/maxsepacutsroot", 400) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/aggregation/maxfailsroot", 200) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/aggregation/maxsepacutsroot", 1000) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/zerohalf/maxroundsroot", 30) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/zerohalf/maxsepacutsroot", 200) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/rlt/freq", 20) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/rlt/maxroundsroot", 15) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/disjunctive/freq", 20) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/disjunctive/maxroundsroot", 150) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/clique/freq", 20) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/mcf/freq", 20) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/mcf/maxtestdelta", -1) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "separating/mcf/trynegscaling", TRUE) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "separating/mcf/maxsepacutsroot", 400) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/linear/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/linear/maxsepacutsroot", 500) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "constraints/linear/separateall", TRUE) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/knapsack/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/knapsack/maxsepacutsroot", 500) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/logicor/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/or/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/setppc/sepafreq", 10) );
      SCIP_CALL( SCIPsetBoolParam(solverscip, "constraints/setppc/cliquelifting", TRUE) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/SOS2/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/varbound/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "constraints/xor/sepafreq", 10) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "propagating/probing/maxuseless", 1500) );
      SCIP_CALL( SCIPsetIntParam(solverscip, "propagating/probing/maxtotaluseless", 75) );
      SCIP_CALL( SCIPsetRealParam(solverscip, "cutselection/hybrid/minorthoroot", 0.1) );
      break;

   default:
      SCIPerrorMessage("invalid racing configuration index <%d>\n", config);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** initialize the concurrent SCIP solver, i.e., setup the copy of the problem and the
 *  mapping of the variables */
static
SCIP_RETCODE initConcsolver(
   SCIP*                 scip,               /**< the main SCIP instance */
   SCIP_CONCSOLVER*      concsolver          /**< the concurrent solver to set up */
   )
{
   SCIP_HASHMAP* varmapfw;
   SCIP_CONCSOLVERDATA* data;
   SCIP_VAR** mainvars;
   SCIP_VAR** mainfixedvars;
   SCIP_VAR** mainallvars;
   SCIP_Bool symmetrybefore;
   SCIP_Bool valid;
   SCIP_Bool usesolpool;
   int paramode;
   int nmainvars;
   int nmainfixedvars;
   int* varperm;
   int cnt;
   int v;

   assert(scip != NULL);
   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   /* we force the copying of symmetry constraints that may have been detected during a central presolving step;
    * otherwise, the copy may become invalid */
   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/symmetrybefore", &symmetrybefore) );
   if( symmetrybefore
      && ( SCIPsetBoolParam(scip, "constraints/orbitope_full/forceconscopy", TRUE) != SCIP_OKAY
      || SCIPsetBoolParam(scip, "constraints/orbitope_pp/forceconscopy", TRUE) != SCIP_OKAY
      || SCIPsetBoolParam(scip, "constraints/orbisack/forceconscopy", TRUE) != SCIP_OKAY
      || SCIPsetBoolParam(scip, "constraints/symresack/forceconscopy", TRUE) != SCIP_OKAY ) )
   {
      SCIPdebugMessage("Could not force copying of symmetry constraints\n");
   }

   /* get number of active variables in main SCIP */
   nmainvars = SCIPgetNVars(scip);
   mainvars = SCIPgetVars(scip);

   /* create the concurrent solver's SCIP instance, set up the problem, and copy symmetry handlers
    * if symmetry wasn't computed before and user wants symmetry */
   SCIP_CALL( SCIPcreate(&data->solverscip) );
   SCIPsetMessagehdlrQuiet(data->solverscip, SCIPmessagehdlrIsQuiet(SCIPgetMessagehdlr(scip)));
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(data->solverscip), nmainvars + SCIPgetNFixedVars(scip)) );
   SCIP_CALL( SCIPcopyConsCompression(scip, data->solverscip, varmapfw, NULL, SCIPconcsolverGetName(concsolver),
         NULL, NULL, 0, TRUE, FALSE, !symmetrybefore, FALSE, FALSE, &valid) );
   assert(valid);

   /* Note that because some aggregations or fixed variables cannot be resolved by some constraint handlers (in
    * particular cons_sos1, cons_sos2, cons_and), the copied problem may contain more variables than the original
    * problem has active variables. */
   data->nvars = SCIPgetNOrigVars(data->solverscip);
   assert( nmainvars <= data->nvars );
   assert(data->nvars <= SCIPgetNVars(scip) + SCIPgetNFixedVars(scip));

   /* allocate memory for the arrays to store the variable mapping */
   SCIP_CALL( SCIPallocBlockMemoryArray(data->solverscip, &data->vars, data->nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(data->solverscip, &varperm, data->nvars) );

   /* In the following, we create a variable mapping between the solver and main SCIP variables. The mapping is created
    * by first retrieving the active variables, then the variables that are "fixed" in the main SCIP. This order is
    * taken because when performing SCIPcopyConsCompression, the active variables are copied first. This is followed by
    * the variables, which might involve (multi-)aggregated/fixed variables and coupling linear constraints. The
    * latter variables appear in the main SCIP as "fixed" variables. */

   /* set up the arrays for the variable mapping */
   SCIP_CALL( SCIPallocBufferArray(data->solverscip, &mainallvars, data->nvars) );
   for( v = 0; v < nmainvars; v++ )
   {
      SCIP_VAR* var;
      int idx;

      var = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, mainvars[v]);
      assert(var != NULL);
      idx = SCIPvarGetProbindex(var);
      assert(0 <= idx && idx < data->nvars);

      data->vars[v] = var;
      assert(varperm[idx] == 0);
      varperm[idx] = v;

      /* for copying solutions below */
      mainallvars[v] = mainvars[v];
   }

   nmainfixedvars = SCIPgetNFixedVars(scip);
   mainfixedvars = SCIPgetFixedVars(scip);
   cnt = nmainvars;
   for( v = 0; v < nmainfixedvars; v++ )
   {
      SCIP_VAR* var;
      int idx;

      var = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, mainfixedvars[v]);
      if( var != NULL )
      {
         idx = SCIPvarGetProbindex(var);
         if( idx >= 0 )
         {
            assert(idx < data->nvars);

            data->vars[cnt] = var;
            assert(varperm[idx] == 0);
            varperm[idx] = cnt;

            /* for copying solutions below */
            mainallvars[cnt] = mainfixedvars[v];
            ++cnt;
         }
      }
   }
   assert( cnt == data->nvars );

   /* transfer solutions from original problem to concurrent instances */
   if( SCIPgetNSols(scip) != 0 )
   {
      SCIP_Bool stored;
      SCIP_SOL* mainsol;
      SCIP_SOL* solversol;

      mainsol = SCIPgetBestSol(scip);
      SCIP_CALL( SCIPcreateSol(data->solverscip, &solversol, NULL) );
      for( v = 0; v < data->nvars; ++v )
      {
         SCIP_Real val;

         val = SCIPgetSolVal(scip, mainsol, mainallvars[v]);
         assert(data->vars[v] != NULL);
         SCIP_CALL( SCIPsetSolVal(data->solverscip, solversol, data->vars[v], val) );
      }
      SCIP_CALL( SCIPaddSolFree(data->solverscip, &solversol, &stored) );
      assert(stored);
   }

   /* create the concurrent data structure for the concurrent solver's SCIP */
   SCIP_CALL( SCIPcreateConcurrent(data->solverscip, concsolver, varperm, data->nvars) );

   /* enable the solution-pool drain in the sync heuristic, so that solutions published by the
    * other concurrent solvers are tried at every node instead of only at synchronization points
    */
   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/solpool", &usesolpool) );
   SCIP_CALL( SCIPgetIntParam(scip, "parallel/mode", &paramode) );

   if( usesolpool && paramode == (int) SCIP_PARA_OPPORTUNISTIC )
   {
      SCIP_HEUR* heursync;

      heursync = SCIPfindHeur(data->solverscip, "sync");
      assert(heursync != NULL);

      SCIPheurSyncEnableSolPool(data->solverscip, heursync, concsolver, data->vars, data->nvars);
   }
   SCIPfreeBufferArray(data->solverscip, &mainallvars);
   SCIPfreeBufferArray(data->solverscip, &varperm);

   /* free the hashmap */
   SCIPhashmapFree(&varmapfw);

   return SCIP_OKAY;
}

/** creates an instance of a concurrent SCIP solver */
static
SCIP_DECL_CONCSOLVERCREATEINST(concsolverScipCreateInstance)
{
   char filename[SCIP_MAXSTRLEN];
   SCIP_CONCSOLVERDATA* data;
   SCIP_CONCSOLVERTYPEDATA* typedata;
   SCIP_Bool changechildsel;
   char* prefix;

   assert(scip != NULL);
   assert(concsolvertype != NULL);
   assert(concsolver != NULL);

   typedata = SCIPconcsolverTypeGetData(concsolvertype);

   SCIP_ALLOC( BMSallocMemory(&data) );
   SCIPconcsolverSetData(concsolver, data);

   SCIP_CALL( initConcsolver(scip, concsolver) );

   /* check if emphasis setting should be loaded; fix certain parameters before loading the
    * emphasis to avoid setting them to default values
    */
   if( typedata->loademphasis )
   {
      SCIP_PARAM** fixedparams;
      int nfixedparams;

      SCIP_CALL( fixProtectedParams(data->solverscip, &fixedparams, &nfixedparams) );
      SCIP_CALL( SCIPsetEmphasis(data->solverscip, typedata->emphasis, TRUE) );
      SCIP_CALL( unfixProtectedParams(data->solverscip, &fixedparams, nfixedparams) );
   }
   else
   {
      SCIP_Bool racingportfolio;

      /* diversify the concurrent solvers with the built-in racing portfolio; settings files
       * loaded below through concurrent/paramsetprefix may override the portfolio settings
       */
      SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/racingportfolio", &racingportfolio) );

      if( racingportfolio )
      {
         SCIP_PARAM** fixedparams;
         int config;
         int nfixedparams;

         config = SCIPconcsolverGetIdx(concsolver) % NRACINGCONFIGS;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "applying racing configuration <%d> to concurrent solver <%s>\n",
            config, SCIPconcsolverGetName(concsolver));

         SCIP_CALL( fixProtectedParams(data->solverscip, &fixedparams, &nfixedparams) );
         SCIP_CALL( applyRacingSettings(data->solverscip, config) );
         SCIP_CALL( unfixProtectedParams(data->solverscip, &fixedparams, nfixedparams) );
      }
   }

   /* load settings file if it exists */
   SCIP_CALL( SCIPgetStringParam(scip, "concurrent/paramsetprefix", &prefix) );
   (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s%s.set", prefix, SCIPconcsolverGetName(concsolver));

   if( SCIPfileExists(filename) )
   {
      /* load settings file and print info message */
      SCIPinfoMessage(scip, NULL, "reading parameter file <%s> for concurrent solver <%s>\n", filename, SCIPconcsolverGetName(concsolver));
      SCIP_CALL( SCIPreadParams(data->solverscip, filename) );
   }
   else
   {
      /* print message about missing setting files only in verblevel full */
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "skipping non existent parameter file <%s> for concurrent solver <%s>\n",
         filename, SCIPconcsolverGetName(concsolver));
   }

   /* include eventhandler for synchronization */
   SCIP_CALL( includeEventHdlrSync(data->solverscip, concsolver) );

   /* disable output for subscip */
   SCIP_CALL( SCIPsetIntParam(data->solverscip, "display/verblevel", 0) );

   /* use wall clock time in subscips */
   SCIP_CALL( SCIPsetIntParam(data->solverscip, "timing/clocktype", (int)SCIP_CLOCKTYPE_WALL) );

   /* don't catch ctrlc since already caught in main SCIP */
   SCIP_CALL( SCIPsetBoolParam(data->solverscip, "misc/catchctrlc", FALSE) );

   /* one solver can do all dual reductions and share them with the other solvers */
   if( SCIPconcsolverGetIdx(concsolver) != 0 )
   {
      SCIP_CALL( disableConflictingDualReductions(data->solverscip) );
   }

   /* set different child selection rules if corresponding parameter is TRUE */
   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/changechildsel", &changechildsel) );
   if( changechildsel )
   {
      SCIP_CALL( setChildSelRule(concsolver) );
   }

   return SCIP_OKAY;
}

/** destroys an instance of a concurrent SCIP solver */
static
SCIP_DECL_CONCSOLVERDESTROYINST(concsolverScipDestroyInstance)
{
   SCIP_CONCSOLVERDATA* data;

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);
   assert(data->solverscip != NULL);

   /* free the array with the variable mapping */
   SCIPfreeBlockMemoryArray(data->solverscip, &data->vars, data->nvars);

   /* free subscip */
   SCIP_CALL( SCIPfree(&data->solverscip) );
   BMSfreeMemory(&data);
   SCIPconcsolverSetData(concsolver, NULL);

   return SCIP_OKAY;
}

/** frees the data of a concurrent solver type */
static
SCIP_DECL_CONCSOLVERTYPEFREEDATA(concsolverTypeScipFreeData)
{
   assert(data != NULL);
   BMSfreeMemory(data);
}

/** initializes the random and permutation seeds of the concurrent solver
 *
 *  Uses the mild FiberSCIP-style racing scheme (seed shift idx+1, permutation seed idx)
 *  instead of large random seeds: solver 0 keeps an unpermuted problem, the others get
 *  small distinct permutation seeds, and the permutevars/permuteconss defaults are left
 *  untouched. Large random seeds combined with forced variable permutation measurably
 *  slowed every solver relative to a sequential SCIP with the same settings.
 */
static
SCIP_DECL_CONCSOLVERINITSEEDS(concsolverScipInitSeeds)
{  /*lint --e{715}*/
   SCIP_CONCSOLVERDATA* data;
   int idx;

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   idx = SCIPconcsolverGetIdx(concsolver);

   SCIPinfoMessage(data->solverscip, NULL, "initializing seeds to shift %d, permutation %d in concurrent solver '%s'\n",
      idx + 1, idx, SCIPconcsolverGetName(concsolver));

   SCIP_CALL( SCIPsetIntParam(data->solverscip, "randomization/randomseedshift", idx + 1) );
   SCIP_CALL( SCIPsetIntParam(data->solverscip, "randomization/permutationseed", idx) );

   return SCIP_OKAY;
}

/** extracts solving status of this concurrent solver and the solving statistics
 *  into the given SCIP instance
 */
static
SCIP_DECL_CONCSOLVERCOPYSOLVINGDATA(concsolverGetSolvingData)
{
   SCIP_CONCSOLVERDATA* data;
   int nsols;

   assert(scip != NULL);
   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);
   assert(data->solverscip != NULL);

   nsols = SCIPgetNSols(data->solverscip);
   if( nsols > 0 )
   {
      SCIP_VAR** mainvars;
      SCIP_SOL** solversols;
      SCIP_Real* solvals;
      int nmainvars;
      int i;

      mainvars = SCIPgetVars(scip);
      nmainvars = SCIPgetNVars(scip);
      assert(nmainvars <= data->nvars);

      solversols = SCIPgetSols(data->solverscip);

      /* allocate buffer array used for translating the solution to the given SCIP */
      SCIP_CALL( SCIPallocBufferArray(scip, &solvals, nmainvars) );

      /* add the solutions to the given SCIP */
      for( i = 0; i < nsols; ++i )
      {
         SCIP_SOL* mainsol;
         SCIP_HEUR* heur;
         SCIP_Bool stored;

         /* only get the first nmainvars, which correspond to the active variables */
         SCIP_CALL( SCIPgetSolVals(data->solverscip, solversols[i], nmainvars, data->vars, solvals) );

         heur = SCIPsolGetHeur(solversols[i]);
         if( heur != NULL )
            heur = SCIPfindHeur(scip, SCIPheurGetName(heur));

         SCIP_CALL( SCIPcreateSol(scip, &mainsol, heur) );
         SCIP_CALL( SCIPsetSolVals(scip, mainsol, nmainvars, mainvars, solvals) );
         SCIP_CALL( SCIPcopySolStats(solversols[i], mainsol) );

#ifdef SCIP_CHECK_MAINSCIP_SOLUTION
         /* The following sometimes fails because we do not copy aggregations and cons_fixedvar can reject solutions in
          * mainscip, because of these. */
         {
            SCIP_Bool feasible;
            SCIP_CALL( SCIPcheckSol(scip, mainsol, TRUE, TRUE, TRUE, TRUE, FALSE, &feasible) );
            assert( feasible );
         }
#endif

         SCIP_CALL( SCIPaddSolFree(scip, &mainsol, &stored) );
      }

      /* free the buffer array */
      SCIPfreeBufferArray(scip, &solvals);
   }

   /* copy solving statistics and status from the solver SCIP to the given SCIP */
   SCIP_CALL( SCIPcopyConcurrentSolvingStats(data->solverscip, scip) );

   return SCIP_OKAY;
}

/** execution method of SCIP concsolver solver
 *
 *  Start solving the problem until the solving reaches a limit, gets interrupted, or just finished successfully.
 */
static
SCIP_DECL_CONCSOLVEREXEC(concsolverScipExec)
{
   SCIP_CONCSOLVERDATA* data;

   assert(concsolver != NULL);
   assert(solvingtime != NULL);
   assert(nlpiterations != NULL);
   assert(nnodes != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   /* print info message that solving has started */
   SCIPinfoMessage(data->solverscip, NULL, "starting solve in concurrent solver '%s'\n", SCIPconcsolverGetName(concsolver));

   /* solve */
   SCIP_CALL( SCIPsolve(data->solverscip) );

   /* first output time */
   SCIPinfoMessage(data->solverscip, NULL, " ");
   SCIPdispTime(SCIPgetMessagehdlr(data->solverscip), NULL, SCIPgetSolvingTime(data->solverscip), 5);
   /* print info message with status */
   SCIPinfoMessage(data->solverscip, NULL, ": concurrent solver '%s' stopped with status ", SCIPconcsolverGetName(concsolver));
   SCIP_CALL( SCIPprintStatus(data->solverscip, NULL) );
   SCIPinfoMessage(data->solverscip, NULL, "\n");

   /* set solving statistics */
   *solvingtime = SCIPgetSolvingTime(data->solverscip);
   *nlpiterations = SCIPgetNLPIterations(data->solverscip);
   *nnodes = SCIPgetNNodes(data->solverscip);

   return SCIP_OKAY;
}

/** stops the concurrent solver as soon as possible */
static
SCIP_DECL_CONCSOLVERSTOP(concsolverScipStop)
{
   SCIP_CONCSOLVERDATA* data;

   assert(concsolver != NULL);

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   /* Interrupts both the LP solve and sets the user interrupt flag */
   SCIP_CALL( SCIPinterruptLP(data->solverscip, TRUE) );

   return SCIP_OKAY;
}

/** writes new solutions and global boundchanges to the given synchronization data */
static
SCIP_DECL_CONCSOLVERSYNCWRITE(concsolverScipSyncWrite)
{
   SCIP_SOL** sols;
   SCIP_CONCSOLVERDATA* data;
   SCIP_BOUNDSTORE* boundstore;
   SCIP_STATUS solverstatus;
   int concsolverid;
   int nsols;
   int i;

   assert(concsolver != NULL);
   assert(syncstore != NULL);
   assert(syncdata != NULL);
   assert(nsolsshared != NULL);

   *nsolsshared = 0;

   if ( maxcandsols <= 0 )
      return SCIP_OKAY;

   if( SCIPsyncdataGetStatus(syncdata) != SCIP_STATUS_UNKNOWN )
      return SCIP_OKAY;

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);
   assert(data->solverscip != NULL);
   concsolverid = SCIPconcsolverGetIdx(concsolver);
   solverstatus = SCIPgetStatus(data->solverscip);

   SCIPsyncdataSetStatus(syncdata, solverstatus, concsolverid);
   SCIPsyncdataSetLowerbound(syncdata, SCIPgetDualbound(data->solverscip));
   SCIPsyncdataSetUpperbound(syncdata, SCIPgetPrimalbound(data->solverscip));

   /* push the finalized dual bound into bestdualbound; the sync heuristic samples it per node and
    * misses the terminal bound that this sync write (incl. the post-solve one) sees finalized */
   if( SCIPsyncstoreSolPoolEnabled(syncstore) )
   {
      SCIP_Real mindual;

      mindual = SCIPgetDualbound(data->solverscip);
      if( SCIPgetObjsense(data->solverscip) == SCIP_OBJSENSE_MAXIMIZE )
         mindual = -mindual;

      SCIPsyncstoreUpdateBestDualbound(syncstore, mindual);
   }

   SCIPdebugMessage("syncing in concurrent solver %s\n", SCIPconcsolverGetName(concsolver));

   /* in solution-pool mode new incumbents are shared immediately through the pool, so the
    * synchronization data carries only bound changes; in the regular protocol the best
    * solutions are collected into the synchronization data here
    */
   if( !SCIPsyncstoreSolPoolEnabled(syncstore) )
   {
      /* consider at most maxcandsols many solutions, and since the solution array is sorted, consider best solutions */
      nsols = SCIPgetNSols(data->solverscip);
      nsols = MIN(nsols, maxcandsols);
      sols = SCIPgetSols(data->solverscip);

      for( i = 0; i < nsols; ++i )
      {
         if( SCIPIsConcurrentSolNew(data->solverscip, sols[i]) )
         {
            SCIP_Real solobj;
            SCIP_Real* solvals;

            solobj = SCIPgetSolOrigObj(data->solverscip, sols[i]);

            SCIPdebugMessage("adding sol in concurrent solver %s\n", SCIPconcsolverGetName(concsolver));
            SCIPsyncdataGetSolutionBuffer(syncstore, syncdata, solobj, concsolverid, &solvals);

            /* if syncstore has no place for this solution, we can stop, since the next solution will have
             * a worse objective value and thus won't be accepted either
             */
            if( solvals == NULL )
               break;

            ++(*nsolsshared);
            SCIP_CALL( SCIPgetSolVals(data->solverscip, sols[i], data->nvars, data->vars, solvals) );

            /* if we have added the maximum number of solutions we can also stop */
            if( *nsolsshared == maxsharedsols )
               break;
         }
      }
   }

   boundstore = SCIPgetConcurrentGlobalBoundChanges(data->solverscip);

   if( boundstore != NULL )
   {
      SCIP_CALL( SCIPsyncdataAddBoundChanges(syncstore, syncdata, boundstore) );
   }

   SCIPsyncdataAddMemTotal(syncdata, SCIPgetMemTotal(data->solverscip));

   return SCIP_OKAY;
}

/** reads the solutions and bounds from the given synchronization data */
static
SCIP_DECL_CONCSOLVERSYNCREAD(concsolverScipSyncRead)
{  /*lint --e{715}*/
   SCIP_Real** solvals;
   SCIP_CONCSOLVERDATA* data;
   SCIP_BOUNDSTORE* boundstore;
   int* concsolverids;
   int concsolverid;
   int nbndchgs;
   int nsols;
   int i;

   assert(concsolver != NULL);
   assert(syncstore != NULL);
   assert(syncdata != NULL);
   assert(nsolsrecvd != NULL);
   assert(ntighterbnds != NULL);
   assert(ntighterintbnds != NULL);

   *nsolsrecvd = 0;

   data = SCIPconcsolverGetData(concsolver);
   assert(data != NULL);

   concsolverid = SCIPconcsolverGetIdx(concsolver);

   /* get solutions from synchronization data; when the solution pool is enabled the solutions
    * are exchanged through the pool already, so the synchronization data only carries bounds
    */
   if( SCIPsyncstoreSolPoolEnabled(syncstore) )
   {
      solvals = NULL;
      concsolverids = NULL;
      nsols = 0;
   }
   else
      SCIPsyncdataGetSolutions(syncdata, &solvals, &concsolverids, &nsols);

   for( i = 0; i < nsols; ++i )
   {
      SCIP_SOL* newsol;
      SCIP_Bool feasible;

      /* do not add own solutions */
      if( concsolverids[i] == concsolverid )
         continue;

      /* solution is from another solver, so translate to this solver's variable space and add it to SCIP */
      ++(*nsolsrecvd);
      SCIP_CALL( SCIPcreateOrigSol(data->solverscip, &newsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(data->solverscip, newsol, data->nvars, data->vars, solvals[i]) );

      /* check whether solution is feasible */
      SCIP_CALL( SCIPcheckSol(data->solverscip, newsol, FALSE, FALSE, TRUE, TRUE, FALSE, &feasible) );

      if( feasible )
      {
         SCIPdebugMessage("adding solution in concurrent solver %s\n", SCIPconcsolverGetName(concsolver));
         SCIP_CALL( SCIPaddConcurrentSol(data->solverscip, newsol) );
      }
   }

   /* get bound changes from the synchronization data and add it to this concurrent solvers SCIP */
   *ntighterbnds = 0;
   *ntighterintbnds = 0;
   boundstore = SCIPsyncdataGetBoundChgs(syncdata);
   nbndchgs = SCIPboundstoreGetNChgs(boundstore);

   for( i = 0; i < nbndchgs; ++i )
   {
      SCIP_VAR* var;
      SCIP_BOUNDTYPE boundtype;
      SCIP_Real newbound;
      int idx;

      idx = SCIPboundstoreGetChgVaridx(boundstore, i);
      assert(0 <= idx && idx < data->nvars);
      var = data->vars[idx];
      boundtype = SCIPboundstoreGetChgType(boundstore, i);
      newbound = SCIPboundstoreGetChgVal(boundstore, i);

      SCIP_CALL( SCIPvarGetProbvarBound(&var, &newbound, &boundtype) );

      /* cannot change bounds of multi-aggregated variables so do not pass this bound-change to the propagator */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      /* if bound is not better then do not pass this bound and do not waste memory for storing this boundchange */
      if( boundtype == SCIP_BOUNDTYPE_LOWER && SCIPisGE(data->solverscip, SCIPvarGetLbGlobal(var), newbound) )
         continue;

      if( boundtype == SCIP_BOUNDTYPE_UPPER && SCIPisLE(data->solverscip, SCIPvarGetUbGlobal(var), newbound) )
         continue;

      /* bound is better so incremented counters for statistics and pass it to the sync propagator */
      ++(*ntighterbnds);

      if( SCIPvarIsNonimpliedIntegral(var) )
         ++(*ntighterintbnds);

      SCIP_CALL( SCIPaddConcurrentBndchg(data->solverscip, var, newbound, boundtype) );
   }

   return SCIP_OKAY;
}


/** creates the concurrent SCIP solver plugins and includes them in SCIP */
SCIP_RETCODE SCIPincludeConcurrentScipSolvers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONCSOLVERTYPEDATA* data;

   assert(scip != NULL);

   /* Include concurrent solvers for SCIP for all emphasis settings and without an emphasis setting.
    * For the SCIP without an emphasis setting we set the default preferred priority to 1 and for the other types to 0
    * so that the default concurrent solve will use multiple SCIP's using settings as specified by the user in the main SCIP.
    */
   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = FALSE;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip", 1.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_DEFAULT;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-default", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_CPSOLVER;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-cpsolver", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_EASYCIP;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-easycip", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_FEASIBILITY;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-feas", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_HARDLP;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-hardlp", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_OPTIMALITY;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-opti", 0.0, concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   SCIP_CALL( SCIPallocMemory(scip, &data) );
   data->loademphasis = TRUE;
   data->emphasis = SCIP_PARAMEMPHASIS_COUNTER;
   SCIP_CALL( SCIPincludeConcsolverType(scip, "scip-counter", 0.0,  concsolverScipCreateInstance, concsolverScipDestroyInstance, concsolverScipInitSeeds,
         concsolverScipExec, concsolverGetSolvingData, concsolverScipStop, concsolverScipSyncWrite,
         concsolverScipSyncRead, concsolverTypeScipFreeData, data) );

   return SCIP_OKAY;
}
