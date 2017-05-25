/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_treesizeprediction
 * @brief  eventhdlr for tree-size prediction related events
 * @author Pierre Le Bodic
 */


/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG 1

#include "scip/event_treesizeprediction.h"

#include "string.h"
#define EVENTHDLR_NAME         "treesizeprediction"
#define EVENTHDLR_DESC         "event handler for tree-size prediction related events"


/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** event handler data */
struct SCIP_EventhdlrData
{
   unsigned int nodesfound;
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
SCIP_DECL_EVENTCOPY(eventCopyTreeSizePrediction)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tree-size prediction dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventCopyTreeSizePrediction NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_EVENTFREE(eventFreeTreeSizePrediction)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   /* SCIPdebugMsg(scip, "eventfree method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);
   return SCIP_OKAY;
}
#else
#define eventFreeTreeSizePrediction NULL
#endif

/** initialization method of event handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_EVENTINIT(eventInitTreeSizePrediction)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tree-size prediction event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitTreeSizePrediction NULL
#endif

/** deinitialization method of event handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_EVENTEXIT(eventExitTreeSizePrediction)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tree-size prediction event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitTreeSizePrediction NULL
#endif

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#if 1
static
SCIP_DECL_EVENTINITSOL(eventInitsolTreeSizePrediction)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   /* SCIPdebugMsg(scip, "initsol method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->nodesfound = 0;

   /* We catch node solved */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODESOLVED, eventhdlr, NULL, NULL) );

   /* We catch updates to the primal bound */
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}
#else
#define eventInitsolTreeSizePrediction NULL
#endif

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 1
static
SCIP_DECL_EVENTEXITSOL(eventExitsolTreeSizePrediction)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   /* SCIPdebugMsg(scip, "exitsol method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPdebugMessage("Found %u nodes in the B&B tree\n", eventhdlrdata->nodesfound);
   return SCIP_OKAY;
}
#else
#define eventExitsolTreeSizePrediction NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteTreeSizePrediction)
{  /*lint --e{715}*/
}
#else
#define eventDeleteTreeSizePrediction NULL
#endif

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecTreeSizePrediction)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_NODE* foundnode; /* The node found by the event (if any) */
   SCIP_Real  newlowerbound; /* The new lower bound found by the event (if any) */
   /* SCIPdebugMsg(scip, "exec method of eventhdlr "EVENTHDLR_NAME"\n"); */
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   switch(SCIPeventGetType(event)) {
      case SCIP_EVENTTYPE_NODEINFEASIBLE:
      case SCIP_EVENTTYPE_NODEBRANCHED:
      case SCIP_EVENTTYPE_NODEFEASIBLE:
         eventhdlrdata->nodesfound += 1;
         foundnode = SCIPeventGetNode(event);
         assert(foundnode != NULL);
         break;
      case SCIP_EVENTTYPE_BESTSOLFOUND:
      {
         newlowerbound = SCIPgetLowerbound(scip);
         SCIPdebugMsg(scip, "New best solution found\n");
         foundnode = NULL;
         break;
      }
      default:
         SCIPerrorMessage("Missing case in this switch.\n");
         SCIPABORT();
   }

   return SCIP_OKAY;
}

/** creates event handler for tree-size prediction event */
SCIP_RETCODE SCIPincludeEventHdlrTreeSizePrediction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create tree-size prediction event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   /* TODO: (optional) create event handler specific data here */

   eventhdlr = NULL;

   /* include event handler into SCIP */
#if 0
   /* use SCIPincludeEventhdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopyTreeSizePrediction,
         eventFreeTreeSizePrediction, eventInitTreeSizePrediction, eventExitTreeSizePrediction, 
         eventInitsolTreeSizePrediction, eventExitsolTreeSizePrediction, eventDeleteTreeSizePrediction, eventExecTreeSizePrediction,
         eventhdlrdata) );
#else
   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecTreeSizePrediction, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolTreeSizePrediction) );
   SCIP_CALL( SCIPsetEventhdlrDelete(scip, eventhdlr, eventDeleteTreeSizePrediction) );
#endif

   /* add tree-size prediction event handler parameters */
   /* TODO: (optional) add event handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
