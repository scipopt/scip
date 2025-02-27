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

/**@file   event_bestsol.c
 * @brief  eventhdlr to print the best solution value
 * @author Alexander Hoen
 * @author Gioni Mexi
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "message_pb.h"
#include "event_bestsol.h"

#include <string.h>

#define EVENTHDLR_NAME         "bestsol"
#define EVENTHDLR_DESC         "event handler for best solutions found"


/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitBestsol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBestsol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** frees specific event data */
static
SCIP_DECL_EVENTDELETE(eventDeleteBestsol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(eventdata != NULL);
   assert(*eventdata != NULL);

   SCIPfreeMemory(scip, eventdata);

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBestsol)
{  /*lint --e{715}*/
   SCIP_MESSAGEHDLR* messagehdlr;
   SCIP_MESSAGEHDLRDATA* messagehdlrdata;
   SCIP_SOL* bestsol;
   SCIP_Real solvalue;
   SCIP_Bool quiet;
   SCIP_Bool comment;
   int solval;

   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);

   SCIPdebugMessage("exec method of best solution event of pbsolver\n");

   bestsol = SCIPgetBestSol(scip);
   assert(bestsol != NULL);
   solvalue = SCIPgetSolOrigObj(scip, bestsol);
   solval = (solvalue > 0.0) ? (solvalue + 0.5) : (solvalue - 0.5); /*lint !e524*/

   messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr != NULL);
   messagehdlrdata = SCIPmessagehdlrGetData(messagehdlr);
   assert(messagehdlrdata != NULL);

   /* store old quiet and comment value of message handler */
   quiet = SCIPmessagehdlrIsQuiet(messagehdlr);
   comment = messagehdlrdata->comment;

   /* turn on message handler */
   SCIPmessagehdlrSetQuiet(messagehdlr, FALSE);
   messagehdlrdata->comment = FALSE;

   /* print best solution value */
   SCIPinfoMessage(scip, NULL, "o %d\n", solval);

   /* reset old values of message handler data */
   SCIPmessagehdlrSetQuiet(messagehdlr, quiet);
   messagehdlrdata->comment = comment;

   return SCIP_OKAY;
}

/** creates event handler for best solution found */
SCIP_RETCODE SCIPcreateEventHdlrBestsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* create event handler for events on best solutions */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL, NULL, eventInitBestsol, eventExitBestsol, NULL, NULL, eventDeleteBestsol, eventExecBestsol, NULL) );

   return SCIP_OKAY;
}
