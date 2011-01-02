/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: event_xxx.c,v 1.4 2011/01/02 11:10:47 bzfheinz Exp $"

/**@file   event_xxx.c
 * @ingroup EVENTS 
 * @brief  eventhdlr for xxx event
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_xxx.h"

#define EVENTHDLR_NAME         "xxx"
#define EVENTHDLR_DESC         "event handler for xxx event"


/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** event handler data */
struct SCIP_EventhdlrData
{
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
SCIP_DECL_EVENTCOPY(eventCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventCopyXxx NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_EVENTFREE(eventFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventFreeXxx NULL
#endif

/** initialization method of event handler (called after problem was transformed) */
#if 0
static
SCIP_DECL_EVENTINIT(eventInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitXxx NULL
#endif

/** deinitialization method of event handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_EVENTEXIT(eventExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitXxx NULL
#endif

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_EVENTINITSOL(eventInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitsolXxx NULL
#endif

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_EVENTEXITSOL(eventExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitsolXxx NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventDeleteXxx NULL
#endif

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** creates event handler for xxx event */
SCIP_RETCODE SCIPincludeEventHdlrXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create xxx event handler data */
   eventhdlrdata = NULL;
   /* TODO: (optional) create event handler specific data here */

   /* imclude event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventCopyXxx,
         eventFreeXxx, eventInitXxx, eventExitXxx, 
         eventInitsolXxx, eventExitsolXxx, eventDeleteXxx, eventExecXxx,
         eventhdlrdata) );

   /* add xxx event handler parameters */
   /* TODO: (optional) add event handler specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
