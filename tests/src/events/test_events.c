/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   test_events.c
 * @brief  unit test for checking events
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/scip_event.h"


/** GLOBAL VARIABLES **/
#define EVENTHDLR_NAME         "dualboundimproved"
#define EVENTHDLR_DESC         "event handler for catching dual bound improvements"
int ndualboundimprovements = 0;
static SCIP* scip_test;

/** TEST SUITES **/
static
void setup(void)
{
   scip_test = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip_test) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip_test) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip_test) );

   cr_assert_null(scip_test);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyDualBoundImproved)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitDualBoundImproved)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* notify SCIP that your event handler wants to react on the event type best solution found */
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_DUALBOUNDIMPROVED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitDualBoundImproved)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* notify SCIP that your event handler wants to drop the event type best solution found */
   SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_DUALBOUNDIMPROVED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecDualBoundImproved)
{  /*lint --e{715}*/
   SCIP_Real dualvalue;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_DUALBOUNDIMPROVED);

   SCIPdebugMsg(scip, "exec method of event handler for dual bound improvement\n");

   dualvalue = SCIPgetDualbound(scip);

   /* print best solution value */
   SCIPinfoMessage(scip, NULL, "found new best dual bound with value <%g> in SCIP <%s>\n",
      dualvalue, SCIPgetProbName(scip) );

   ndualboundimprovements++;

   return SCIP_OKAY;
}

/** includes event handler for best solution found */
static
SCIP_RETCODE SCIPincludeEventHdlrDualBoundImproved(void)
{
   assert(scip_test != NULL);
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
   eventhdlrdata = NULL;
   eventhdlr = NULL;

   /* create event handler for dual bound improved */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip_test, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecDualBoundImproved, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrCopy(scip_test, eventhdlr, eventCopyDualBoundImproved) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip_test, eventhdlr, eventInitDualBoundImproved) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip_test, eventhdlr, eventExitDualBoundImproved) );

   return SCIP_OKAY;
}


/* test that we correctly catch the following events:
 *  SCIP_EVENTTYPE_DUALBOUNDIMPROVED
 */
TestSuite(events, .init = setup, .fini = teardown);

Test(events, dualboundimproved)
{
   char testfile[SCIP_MAXSTRLEN];
   strcpy(testfile, __FILE__);
   testfile[strlen(testfile) - 13] = '\0';  /* cutoff "test_events.c" */
   strcat(testfile, "../../../check/instances/MIP/bell5.mps");
   SCIP_CALL( SCIPreadProb(scip_test, testfile, NULL) );

   SCIP_CALL( SCIPincludeEventHdlrDualBoundImproved() );
   SCIP_CALL( SCIPsetEventhdlrExit(scip_test, SCIPfindEventhdlr(scip_test, EVENTHDLR_NAME), eventExitDualBoundImproved) );

   SCIP_CALL( SCIPsolve(scip_test) );

   cr_expect_neq(ndualboundimprovements, 0, "No dual bound improvements were caught");
}
