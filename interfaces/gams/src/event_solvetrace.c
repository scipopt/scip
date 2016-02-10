/* Copyright (C) GAMS Development and others 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

/**@file   event_solvetrace.c
 * @ingroup EVENTS 
 * @brief  event handler to write GAMS solve trace file
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "event_solvetrace.h"
#include "GamsSolveTrace.h"
#include "gmomcc.h"

#include <assert.h>

#define EVENTHDLR_NAME         "solvetrace"
#define EVENTHDLR_DESC         "event handler to write GAMS solve trace file"


/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   GAMS_SOLVETRACE*      solvetrace;         /**< GAMS solve trace object */
   gmoHandle_t           gmo;                /**< GAMS model object */
   char*                 filename;           /**< name of trace file */
   int                   nodefreq;           /**< node frequency of writing to trace file */
   SCIP_Real             timefreq;           /**< time frequency of writing to trace file */
   int                   filterpos;          /**< position of node event in SCIPs eventfilter */
   SCIP_Real             dualbound;          /**< dual bound at last exit solve stage */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */



/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_EVENTCOPY(eventCopySolveTrace)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of solvetrace dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventCopySolveTrace NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_EVENTFREE(eventFreeSolveTrace)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->solvetrace == NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}
#else
#define eventFreeSolveTrace NULL
#endif

/** initialization method of event handler (called after problem was transformed) */
#if 1
static
SCIP_DECL_EVENTINIT(eventInitSolveTrace)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   char buffer[GMS_SSSIZE];
   int rc;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->gmo != NULL);

   if( eventhdlrdata->filename[0] == '\0' )
      return SCIP_OKAY;

   rc = GAMSsolvetraceCreate(&eventhdlrdata->solvetrace, eventhdlrdata->filename, "SCIPDEV", gmoOptFile(eventhdlrdata->gmo), gmoNameInput(eventhdlrdata->gmo, buffer), SCIPinfinity(scip), eventhdlrdata->nodefreq, eventhdlrdata->timefreq);
   if( rc != 0 )
   {
      SCIPerrorMessage("GAMSsolvetraceCreate returned with error %d, trace file name = %s\n", rc, eventhdlrdata->filename);
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED | SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_PRESOLVEROUND | SCIP_EVENTTYPE_LPEVENT, eventhdlr, (SCIP_EVENTDATA*)eventhdlrdata->solvetrace, &eventhdlrdata->filterpos) );

   return SCIP_OKAY;
}
#else
#define eventInitSolveTrace NULL
#endif

/** deinitialization method of event handler (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_EVENTEXIT(eventExitSolveTrace)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->solvetrace == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_NODESOLVED | SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_PRESOLVEROUND | SCIP_EVENTTYPE_LPEVENT, eventhdlr, (SCIP_EVENTDATA*)eventhdlrdata->solvetrace, eventhdlrdata->filterpos) );

   GAMSsolvetraceAddEndLine(eventhdlrdata->solvetrace, SCIPgetNTotalNodes(scip),
      eventhdlrdata->dualbound == SCIP_INVALID ? -(int)SCIPgetObjsense(scip) * SCIPinfinity(scip) : eventhdlrdata->dualbound,   /*lint !e777*/
      SCIPgetNSols(scip) > 0 ? SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)) : (int)SCIPgetObjsense(scip) * SCIPinfinity(scip));

   GAMSsolvetraceFree(&eventhdlrdata->solvetrace);
   assert(eventhdlrdata->solvetrace == NULL);

   return SCIP_OKAY;
}
#else
#define eventExitSolveTrace NULL
#endif

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_EVENTINITSOL(eventInitsolSolveTrace)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of solvetrace event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitsolSolveTrace NULL
#endif

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_EVENTEXITSOL(eventExitsolSolveTrace)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of solvetrace event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitsolSolveTrace NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteSolveTrace)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);

   return SCIP_OKAY;
}
#else
#define eventDeleteSolveTrace NULL
#endif

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSolveTrace)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(event != NULL);
   assert(eventdata != NULL);

   GAMSsolvetraceAddLine((GAMS_SOLVETRACE*)eventdata, SCIPgetNTotalNodes(scip),
      SCIPgetDualbound(scip),
      SCIPgetNSols(scip) > 0 ? SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)) : (int) SCIPgetObjsense(scip) * SCIPinfinity(scip));

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->dualbound = SCIPgetDualbound(scip);

   return SCIP_OKAY;
}

/** creates event handler for solvetrace event */
SCIP_RETCODE SCIPincludeEventHdlrSolveTrace(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoHandle_t           gmo                 /**< GAMS model object */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create solvetrace event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->solvetrace = NULL;
   eventhdlrdata->gmo = gmo;
   eventhdlrdata->filename = NULL;
   eventhdlrdata->dualbound = SCIP_INVALID;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopySolveTrace,
         eventFreeSolveTrace, eventInitSolveTrace, eventExitSolveTrace,
         eventInitsolSolveTrace, eventExitsolSolveTrace, eventDeleteSolveTrace, eventExecSolveTrace,
         eventhdlrdata) );

   /* add solvetrace event handler parameters */
   SCIP_CALL( SCIPaddStringParam(scip, "gams/solvetrace/file",
      "name of file where to write branch-and-bound trace information too",
      &eventhdlrdata->filename, FALSE, "", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "gams/solvetrace/nodefreq",
      "frequency in number of nodes when to write branch-and-bound trace information, 0 to disable",
      &eventhdlrdata->nodefreq, FALSE, 100, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "gams/solvetrace/timefreq",
      "frequency in seconds when to write branch-and-bound trace information, 0.0 to disable",
      &eventhdlrdata->timefreq, FALSE, 5.0, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
