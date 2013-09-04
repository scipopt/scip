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

/**@file   event_nodeevent.c
 * @brief  eventhdlr for nodeevent event
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_nodeevent.h"
#include "string.h"

#define EVENTHDLR_NAME         "nodeevent"
#define EVENTHDLR_DESC         "event handler for nodeevent event"
#define EVENTTOCATCH SCIP_EVENTTYPE_BESTSOLFOUND | SCIP_EVENTTYPE_NODEINFEASIBLE

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Real avgcutoffdepth;
   SCIP_Longint ncutoffnodes;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyNodeevent)
{
   SCIP_CALL( SCIPincludeEventHdlrNodeevent(scip) );

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeNodeevent)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Average bounding depth: %5.2f (%"SCIP_LONGINT_FORMAT" cutoffs)\n", eventhdlrdata->ncutoffnodes == 0 ? -1 : eventhdlrdata->avgcutoffdepth, eventhdlrdata->ncutoffnodes);

   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitNodeevent)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   eventhdlrdata->avgcutoffdepth = .0;
   eventhdlrdata->ncutoffnodes = 0L;
/*
   SCIPcatchEvent(scip, EVENTTOCATCH, eventhdlr, NULL, NULL);
*/
   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecNodeevent)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_BESTSOLFOUND )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Average bounding depth: %5.2f (%"SCIP_LONGINT_FORMAT" cutoffs)\n", eventhdlrdata->ncutoffnodes == 0 ? -1 : eventhdlrdata->avgcutoffdepth, eventhdlrdata->ncutoffnodes);
   }
   else if( (SCIPeventGetType(event) & SCIP_EVENTTYPE_NODEINFEASIBLE) && !(SCIPinProbing(scip) || SCIPinDive(scip)) )
   {
      ++(eventhdlrdata->ncutoffnodes);

      eventhdlrdata->avgcutoffdepth = (eventhdlrdata->ncutoffnodes - 1)/(SCIP_Real)(eventhdlrdata->ncutoffnodes) * eventhdlrdata->avgcutoffdepth;
      eventhdlrdata->avgcutoffdepth += SCIPgetFocusDepth(scip)/(SCIP_Real)(eventhdlrdata->ncutoffnodes);
   }



   return SCIP_OKAY;
}

/** creates event handler for nodeevent event */
SCIP_RETCODE SCIPincludeEventHdlrNodeevent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create nodeevent event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   /* TODO: (optional) create event handler specific data here */

   eventhdlr = NULL;

   /* include event handler into SCIP */
   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecNodeevent, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyNodeevent) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeNodeevent) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitNodeevent) );

   /* add nodeevent event handler parameters */
   /* TODO: (optional) add event handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
