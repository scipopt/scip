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

/**@file   event_estimation.c
 * @brief  eventhdlr for estimation event
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define FILEOUTPUT
#include "scip/event_estimation.h"
#include "scip/scip.h"
#include <string.h>

#define EVENTHDLR_NAME         "estimation"
#define EVENTHDLR_DESC         "event handler for estimation event"
#define EVENT_TYPE_ESTIMATION SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE
#define DEFAULT_ENABLED        FALSE
#define HEADER "NUMBER  PRIMAL LOWER ESTIM ROOTLPPSESTIM DEPTH NCANDS NODELOWER LPSTAT"
/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_VAR**           rootlpcands;
   int                  nrootlpcands;
   SCIP_Real*           rootlpcandsfrac;
   SCIP_Bool            enabled;            /**< is the estimation event handler enabled? */
   int                  eventfilterpos;     /**< the event filter position, or -1, if event has not (yet) been caught */
   char*                filename;           /**< file to keep node results */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** free root LP solution from event handler data */
static
void dataFreeRootLPSol(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   assert(eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands != NULL);

   if( eventhdlrdata->rootlpcands != NULL )
   {
      SCIPfreeMemoryArray(scip, &eventhdlrdata->rootlpcands);
      SCIPfreeMemoryArray(scip, &eventhdlrdata->rootlpcandsfrac);
   }
   eventhdlrdata->nrootlpcands = 0;
}

/** stores root LP solution (after separation) in eventhdlrdata */
static
SCIP_RETCODE storeRootLPSol(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*     eventhdlrdata
   )
{
   SCIP_VAR** lpcandsroot;
   SCIP_Real* lpcandsfrac;

   assert(scip != NULL);
   assert(eventhdlrdata != NULL);
   assert(SCIPgetDepth(scip) == 0);

   /* free previously stored root LP solutions */
   if( eventhdlrdata->rootlpcands != NULL )
   {
      dataFreeRootLPSol(scip, eventhdlrdata);
   }
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OBJLIMIT )
   {
      SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcandsroot, NULL, &lpcandsfrac, &eventhdlrdata->nrootlpcands, NULL, NULL) );
   }

   if( eventhdlrdata->nrootlpcands > 0 )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &eventhdlrdata->rootlpcands, lpcandsroot, eventhdlrdata->nrootlpcands) );
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &eventhdlrdata->rootlpcandsfrac, lpcandsfrac, eventhdlrdata->nrootlpcands) );
   }
   return SCIP_OKAY;
}

static
SCIP_Real recalcClassicEstimate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real estimate;
   int i;
   estimate = SCIPgetLowerboundRoot(scip);

   if( eventhdlrdata->nrootlpcands == 0 || eventhdlrdata->rootlpcands == NULL )
      return SCIP_INVALID;

   for( i = 0; i < eventhdlrdata->nrootlpcands; ++i )
   {
      SCIP_Real upscore;
      SCIP_Real downscore;

      upscore = SCIPgetVarPseudocost(scip, eventhdlrdata->rootlpcands[i], SCIP_BRANCHDIR_UPWARDS);
      downscore = SCIPgetVarPseudocost(scip, eventhdlrdata->rootlpcands[i], SCIP_BRANCHDIR_DOWNWARDS);

      upscore *= (1.0 - eventhdlrdata->rootlpcandsfrac[i]);
      downscore *= eventhdlrdata->rootlpcandsfrac[i];
      estimate += MIN(upscore, downscore);
   }

   return estimate;
}

SCIP_Real SCIPgetRootLPSolPscostEstimate(
   SCIP*                 scip
   )
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   return recalcClassicEstimate(scip, eventhdlrdata);
}

/*
 * Callback methods of event handler
 */

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolEstimation)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos == -1);

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENT_TYPE_ESTIMATION, eventhdlr, NULL, &eventhdlrdata->eventfilterpos) );
      assert(eventhdlrdata->eventfilterpos >= 0);
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTFREE(eventFreeEstimation)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->rootlpcands == NULL && eventhdlrdata->rootlpcandsfrac == NULL);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->filename != NULL )
      SCIPfreeMemoryArray(scip, &eventhdlrdata->filename);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTINIT(eventInitEstimation)
{
#ifdef FILEOUTPUT
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   FILE* file;
   char filename[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   if( eventhdlrdata->filename != NULL )
      SCIPfreeMemory(scip, &eventhdlrdata->filename);

   sprintf(filename, "/OPTI/bzfhende/nodefiles/%s.npf",SCIPstringGetBasename(SCIPgetProbName(scip)));

   SCIPduplicateMemoryArray(scip, &eventhdlrdata->filename, filename, strlen(filename) + 1);

   /* open the file and overwrite its content */
   if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      file = fopen(eventhdlrdata->filename, "w+");
      assert(file != NULL);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, file, HEADER"\n");
      fclose(file);
   }
#endif

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitEstimation)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   dataFreeRootLPSol(scip, eventhdlrdata);

   assert(eventhdlrdata->eventfilterpos == -1 || eventhdlrdata->enabled );
   if( eventhdlrdata->eventfilterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENT_TYPE_ESTIMATION, eventhdlr, NULL, eventhdlrdata->eventfilterpos) );
      eventhdlrdata->eventfilterpos = -1;
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolEstimation)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   dataFreeRootLPSol(scip, eventhdlrdata);

   assert(eventhdlrdata->eventfilterpos == -1 || eventhdlrdata->enabled );
   if( eventhdlrdata->eventfilterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENT_TYPE_ESTIMATION, eventhdlr, NULL, eventhdlrdata->eventfilterpos) );
      eventhdlrdata->eventfilterpos = -1;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecEstimation)
{  /*lint --e{715}*/
   SCIP_NODE* focusnode;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos >= 0);

   focusnode = SCIPgetCurrentNode(scip);

   if( focusnode == NULL )
      return SCIP_OKAY;

   /* store the last root LP solution, if it is not already integer feasible */
   if( SCIPnodeGetDepth(focusnode) == 0 && SCIPeventGetType(event) == SCIP_EVENTTYPE_NODEBRANCHED )
   {
      SCIP_CALL( storeRootLPSol(scip, eventhdlrdata) );
   }

#ifdef FILEOUTPUT
   if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_NORMAL )
   {
      FILE* file;
      SCIP_Real lowerbound;
      SCIP_Real estimate;
      SCIP_Real upperbound;
      SCIP_NODE* lowerboundnode;
      int depth;
      int nlpcands;

      nlpcands = -1;

      if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
         nlpcands = SCIPgetNLPBranchCands(scip);


      lowerboundnode = SCIPgetBestboundNode(scip);

      lowerbound = lowerboundnode != NULL ? SCIPnodeGetLowerbound(lowerboundnode) : SCIPinfinity(scip);
      estimate = SCIPnodeGetEstimate(focusnode);
      upperbound = SCIPgetUpperbound(scip);
      depth = SCIPnodeGetDepth(focusnode);
      file = fopen(eventhdlrdata->filename, "a");
      assert(file != NULL);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, file, "%lld %15g %15g %15g %15g %d %d %15g %d\n",
            SCIPnodeGetNumber(focusnode),
            SCIPgetExternalValue(scip, upperbound),
            SCIPgetExternalValue(scip, lowerbound),
            SCIPgetExternalValue(scip, estimate),
            SCIPgetExternalValue(scip, SCIPgetRootLPSolPscostEstimate(scip)),
            depth, nlpcands, SCIPgetExternalValue(scip, SCIPnodeGetLowerbound(focusnode)), SCIPgetLPSolstat(scip));
      fclose(file);
   }
#endif


   return SCIP_OKAY;
}

/** creates event handler for estimation event */
SCIP_RETCODE SCIPincludeEventHdlrEstimation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create estimation event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);
   eventhdlrdata->rootlpcands = NULL;
   eventhdlrdata->rootlpcandsfrac = NULL;
   eventhdlrdata->eventfilterpos = -1;
   eventhdlrdata->filename = NULL;

   eventhdlr = NULL;
   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecEstimation, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitEstimation) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolEstimation) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolEstimation) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeEstimation) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitEstimation) );

   /* add estimation event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/estimation/enabled","enable event handler to perform estimation",
               &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   return SCIP_OKAY;
}
