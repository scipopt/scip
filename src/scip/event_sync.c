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

/**@file   event_sync.c
 * @brief  eventhandler for synchronization
 * @author Robert Lion Gottwald
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_sync.h"
#include "scip/event.h"
#include "scip/scip.h"
#include "scip/concurrent.h"
#include "scip/syncstore.h"
#include <string.h>

#define EVENTHDLR_NAME         "sync"
#define EVENTHDLR_DESC         "event handler for synchronization point"

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   int             filterpos;
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
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

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
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);

   if( eventhdlrdata->filterpos < 0 && SCIPsyncstoreIsInitialized(syncstore) )
   {
      /* notify SCIP that your event handler wants to react on synchronization events */
      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SYNC, eventhdlr, NULL, &eventhdlrdata->filterpos) );
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
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* notify SCIP that your event handler wants to drop the event type synchronization found */
   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SYNC, eventhdlr, NULL, eventhdlrdata->filterpos) );
      eventhdlrdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecSync)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);

   SCIP_CALL( SCIPsynchronize(scip) );

   return SCIP_OKAY;
}


/** includes event handler for synchronization found */
SCIP_RETCODE SCIPincludeEventHdlrSync(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLR*     eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->filterpos = -1;

   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecSync, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeSync) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitSync) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitSync) );

   return SCIP_OKAY;
}
