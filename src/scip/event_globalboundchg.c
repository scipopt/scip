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

/**@file   event_globalboundchg.c
 * @brief  eventhdlr for globalboundchg event
 * @author Jakob Witzig
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include "scip/branch_pseudo.h"
#include "scip/event_globalboundchg.h"

#define EVENTHDLR_NAME         "globalboundchg"
#define EVENTHDLR_DESC         "event handler for globalboundchg event"

/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool             init;                    /** is initialized */
   SCIP_Bool             reopt;                   /** is reoptimization enabled */
};

/** copy method for event handler plugins (called when SCIP copies plugins) */
#define eventCopyGlobalboundchg NULL;

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeGlobalboundchg)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolGlobalboundchg)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR** vars;
   int varnr;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   vars = SCIPgetVars(scip);

   if( eventhdlrdata->reopt )
   {
      for(varnr = 0; varnr < SCIPgetNVars(scip); ++varnr)
      {
         SCIP_CALL(SCIPcatchVarEvent(scip, vars[varnr], SCIP_EVENTTYPE_GBDCHANGED, eventhdlr, NULL, NULL));
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolGlobalboundchg)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_VAR** vars;
   int varnr;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   vars = SCIPgetVars(scip);

   if( eventhdlrdata->reopt )
   {
      for(varnr = 0; varnr < SCIPgetNVars(scip); ++varnr)
      {
         SCIP_CALL(SCIPdropVarEvent(scip, vars[varnr], SCIP_EVENTTYPE_GBDCHANGED , eventhdlr, NULL, -1));
      }
   }
   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitGlobalboundchg)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /** HACK */
   /** check if all variable are binary, if not, disable reoptimization */
   if( !eventhdlrdata->init )
   {
      int maxsavednodes;

      SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/enable", &eventhdlrdata->reopt) );
      if( eventhdlrdata->reopt && SCIPgetNImplVars(scip) + SCIPgetNIntVars(scip) > 0 )
         eventhdlrdata->reopt = FALSE;

      SCIP_CALL( SCIPgetIntParam(scip, "reoptimization/maxsavednodes", &maxsavednodes) );

      if( maxsavednodes == 0 )
         eventhdlrdata->reopt = FALSE;
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
#define eventExitGlobalboundchg NULL;

/** frees specific event data */
#define eventDeleteGlobalboundchg NULL;

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecGlobalboundchg)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_NODE*          eventnode;
   SCIP_BOUNDTYPE      boundtype;
   int                 oldbound;
   int                 newbound;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( !eventhdlrdata->reopt )
      return SCIP_OKAY;

   eventnode = SCIPgetCurrentNode(scip);
   oldbound = SCIPeventGetOldbound(event);
   newbound = SCIPeventGetNewbound(event);
   boundtype = (SCIPisFeasGE(scip, oldbound, newbound) ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);

   assert( eventnode != NULL );

   /* skip if the node is not the focus nodes */
   if( SCIPnodeGetType(eventnode) != SCIP_NODETYPE_FOCUSNODE )
      return SCIP_OKAY;

   SCIPdebugMessage("catch event %x for node %lld\n", SCIPeventGetType(event), SCIPnodeGetNumber(SCIPeventGetNode(event)));
   SCIPdebugMessage(" -> depth: %d\n");
   SCIPdebugMessage(" -> change bound for <%s>: %g -> %g\n", SCIPvarGetName(SCIPeventGetVar(event)), SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

   SCIP_CALL(SCIPbranchrulePseudoAddPseudoVar(scip, eventnode, SCIPeventGetVar(event), boundtype, newbound));

   return SCIP_OKAY;
}

/** creates event handler for globalboundchg event */
SCIP_RETCODE SCIPincludeEventHdlrGlobalboundchg(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create globalboundchg event handler data */
   eventhdlr = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecGlobalboundchg, eventhdlrdata) );

   assert(eventhdlr != NULL);
   eventhdlrdata->init = FALSE;

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolGlobalboundchg) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolGlobalboundchg) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeGlobalboundchg) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitGlobalboundchg) );

   return SCIP_OKAY;
}
