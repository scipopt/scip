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

/**@file   event_focusroot.c
 * @brief  eventhdlr for focusroot event
 * @author Jakob Witzig
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <string.h>

#include "scip/event_focusroot.h"
#include "scip/branch_nodereopt.h"

#define EVENTHDLR_NAME         "focusroot"
#define EVENTHDLR_DESC         "event handler for focusroot event"

/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool init; /** is already initialized? */
   SCIP_Bool reopt; /** is reoptimization enabled? */
   SCIP_Bool savelpbasis; /** should the LP basis saved? */
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
SCIP_DECL_EVENTCOPY(eventCopyFocusroot)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of focusroot dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventCopyFocusroot NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
static SCIP_DECL_EVENTFREE(eventFreeFocusroot)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static SCIP_DECL_EVENTINIT(eventInitFocusroot)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   /** disable the event handler if there are integer variables
    * todo: extend reoptimization to integer variables
    * */
   if (!eventhdlrdata->init)
   {
      SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/enable", &eventhdlrdata->reopt));

      if (eventhdlrdata->reopt && SCIPgetNIntVars(scip) > 0)
         eventhdlrdata->reopt = FALSE;

      if (!eventhdlrdata->reopt)
         eventhdlrdata->savelpbasis = FALSE;
      else
      {
         SCIP_CALL(SCIPgetBoolParam(scip, "reoptimization/savelpbasis", &eventhdlrdata->savelpbasis));
      }

      eventhdlrdata->init = TRUE;
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_EVENTEXIT(eventExitFocusroot)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of focusroot event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitFocusroot NULL
#endif

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static SCIP_DECL_EVENTINITSOL(eventInitsolFocusroot)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );

   if (eventhdlrdata->init && eventhdlrdata->savelpbasis)
   {
      SCIP_CALL(SCIPcatchEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, NULL));
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_EVENTEXITSOL(eventExitsolFocusroot)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of focusroot event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitsolFocusroot NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteFocusroot)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of focusroot event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventDeleteFocusroot NULL
#endif

/** execution method of event handler */
static SCIP_DECL_EVENTEXEC(eventExecFocusroot)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL );
   assert(eventhdlr != NULL );
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL );
   assert(SCIPgetCurrentNode(scip) == SCIPgetRootNode(scip));

   /** set the lpi (if some exists) */
   if (eventhdlrdata->savelpbasis)
   {
      /* @todo: Set the root LPI */
      SCIPdebugMessage(">> need to implement setting root LPI.\n");
   }

   /** disable this eventhandler for the current round */
   SCIP_CALL(SCIPdropEvent(scip, SCIP_EVENTTYPE_NODEFOCUSED, eventhdlr, NULL, -1));

   return SCIP_OKAY;
}

/** creates event handler for focusroot event */
SCIP_RETCODE
SCIPincludeEventHdlrFocusroot(SCIP* scip /**< SCIP data structure */
)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create focusroot event handler data */
   eventhdlr = NULL;
   SCIP_CALL(SCIPallocMemory(scip, &eventhdlrdata));

   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL(
         SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecFocusroot, eventhdlrdata));
   assert(eventhdlr != NULL );

   eventhdlrdata->init = FALSE;
   eventhdlrdata->reopt = FALSE;
   eventhdlrdata->savelpbasis = FALSE;

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeFocusroot));
   SCIP_CALL(SCIPsetEventhdlrInit(scip, eventhdlr, eventInitFocusroot));
   SCIP_CALL(SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolFocusroot));

   return SCIP_OKAY;
}
