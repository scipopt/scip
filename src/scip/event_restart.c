/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_restart.c
 * @brief  eventhdlr for restart event
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_restart.h"
#include "scip/event_treesizeprediction.h"
#include "scip/event_treeprofile.h"
#include "pub_event.h"
#include "pub_message.h"
#include "scip_event.h"
#include "scip_mem.h"
#include "scip_message.h"
#include "scip_param.h"
#include "scip_solve.h"
#include "scip_solvingstats.h"
#include "type_event.h"
#include "type_message.h"
#include "type_retcode.h"

#define EVENTHDLR_NAME         "restart"
#define EVENTHDLR_DESC         "event handler for restart event"
#define EVENTTYPE_RESTART      SCIP_EVENTTYPE_NODESOLVED

/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */
enum RestartPolicy
{
   RESTARTPOLICY_NEVER      = 0,             /**< never restart (disable this event handler) */
   RESTARTPOLICY_ALWAYS     = 1,             /**< always restart (can be fine tuned by using minimum number of nodes and restart limit) */
   RESTARTPOLICY_ESTIMATION = 2,             /**< base restart on the estimation method */
   RESTARTPOLICY_PROGRESS   = 3              /**< use progress measure to trigger restart */
};

typedef enum RestartPolicy RESTARTPOLICY;

#define RESTARTPOLICY_CHAR_NEVER 'n'
#define RESTARTPOLICY_CHAR_ALWAYS 'a'
#define RESTARTPOLICY_CHAR_ESTIMATION 'e'
#define RESTARTPOLICY_CHAR_PROGRESS 'p'

#define ESTIMATION_CHAR_TREESIZE         't' /**< should estimation use probability based tree size prediction? */
#define ESTIMATION_CHAR_PROFILE          'p'  /**< should estimation use profile based prediction a la Cornuejols? */


/** event handler data */
struct SCIP_EventhdlrData
{
   char                  restartpolicyparam; /**< restart policy parameter */
   char                  estimationparam;    /**< parameter to select the estimation method */
   char                  progressparam;      /**< progress method to use */
   int                   restartlimit;       /**< how often should a restart be triggered? (-1 for no limit) */
   int                   nrestartsperformed; /**< number of restarts performed so far */
   SCIP_Longint          minnodes;           /**< minimum number of nodes in a run before restart is triggered */
   SCIP_Bool             countonlyleaves;    /**< should only leaves count for the minnodes parameter? */
   SCIP_Real             estim_factor;       /**< factor by which the estimated number of nodes should exceed the current number of nodes */
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
SCIP_DECL_EVENTCOPY(eventCopyRestart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of restart dialog not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventCopyRestart NULL
#endif

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeRestart)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitRestart)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPcatchEvent(scip, EVENTTYPE_RESTART, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_EVENTEXIT(eventExitRestart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of restart event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitRestart NULL
#endif

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_EVENTINITSOL(eventInitsolRestart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of restart event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventInitsolRestart NULL
#endif

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_EVENTEXITSOL(eventExitsolRestart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of restart event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventExitsolRestart NULL
#endif

/** frees specific event data */
#if 0
static
SCIP_DECL_EVENTDELETE(eventDeleteRestart)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of restart event handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define eventDeleteRestart NULL
#endif

/** todo get restartpolicy based on the value of the restart parameter */
static
RESTARTPOLICY getRestartPolicy(
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   switch (eventhdlrdata->restartpolicyparam) {
      case RESTARTPOLICY_CHAR_ALWAYS:
         return RESTARTPOLICY_ALWAYS;
      case RESTARTPOLICY_CHAR_NEVER:
         return RESTARTPOLICY_NEVER;
      case RESTARTPOLICY_CHAR_ESTIMATION:
         return RESTARTPOLICY_ESTIMATION;
      case RESTARTPOLICY_CHAR_PROGRESS:
         return RESTARTPOLICY_PROGRESS;
      default:
         SCIPerrorMessage("Unknown restart policy %c\n", eventhdlrdata->restartpolicyparam);
         SCIPABORT();
         break;
   }
}

/** check conditions before applying restart policy */
static
SCIP_Bool checkConditions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Longint nnodes;
   /* check if max number of restarts has been reached */
   if( eventhdlrdata->restartlimit != -1 &&
         eventhdlrdata->nrestartsperformed >= eventhdlrdata->restartlimit )

      return FALSE;

   /* check if number of nodes exceeds the minimum number of nodes */
   if( eventhdlrdata->countonlyleaves )
      nnodes = SCIPgetNFeasibleLeaves(scip) + SCIPgetNInfeasibleLeaves(scip) + SCIPgetNObjlimLeaves(scip);
   else
      nnodes = SCIPgetNNodes(scip);

   if( nnodes < eventhdlrdata->minnodes )
      return FALSE;

   return TRUE;
}

/** should a restart be applied based on the current tree size estimation? */
static
SCIP_Bool shouldApplyRestartEstimation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_Real estimation = -1.0;
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   /* query a tree size estimation based on the user parameter */
   switch (eventhdlrdata->estimationparam) {
      case ESTIMATION_CHAR_TREESIZE:
         /* use the probability based tree size prediction */
         estimation = SCIPtreeSizeGetEstimateTotal(scip);
         break;
      case ESTIMATION_CHAR_PROFILE:
         /* use the estimation based on the tree profile */
         estimation = SCIPpredictTotalSizeTreeprofile(scip);
         break;
      default:
         break;
   }

   /* no estimation is available yet */
   if( estimation < 0.0 )
      return FALSE;

   /* if the estimation exceeds the current number of nodes by a dramatic factor, restart */
   if( estimation > SCIPgetNNodes(scip) * eventhdlrdata->estim_factor )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "Estimation %g exceeds current number of nodes %lld by a factor of %.1f\n",
               estimation, SCIPgetNNodes(scip), estimation / SCIPgetNNodes(scip));
      return TRUE;
   }

   return FALSE;
}

/** check if a restart should be performed based on the given restart policy */
static
SCIP_Bool shouldApplyRestart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   switch (getRestartPolicy(eventhdlrdata)) {
      case RESTARTPOLICY_ALWAYS:
         return TRUE;
      case RESTARTPOLICY_NEVER:
         return FALSE;
      case RESTARTPOLICY_ESTIMATION:
         return shouldApplyRestartEstimation(scip, eventhdlrdata);
      case RESTARTPOLICY_PROGRESS:
         /* author bzfhende
          *
          * TODO both need to be implemented
          */

         break;
      default:
         break;
   }

   return FALSE;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecRestart)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA*   eventhdlrdata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* check if all conditions are met such that the event handler should run */
   if( ! checkConditions(scip, eventhdlrdata) )
      return SCIP_OKAY;

   /* author bzfhende
    *
    * TODO comment
    */

   if( shouldApplyRestart(scip, eventhdlrdata) )
   {
      eventhdlrdata->nrestartsperformed++;


      SCIPrestartSolve(scip);
   }

   return SCIP_OKAY;
}

/** creates event handler for restart event */
SCIP_RETCODE SCIPincludeEventHdlrRestart(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create restart event handler data */
   eventhdlrdata = NULL;
   /* TODO: (optional) create event handler specific data here */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   /* author bzfhende
    *
    * TODO add parameters
    */

   BMSclearMemory(eventhdlrdata);

   eventhdlr = NULL;

   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecRestart, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyRestart) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeRestart) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitRestart) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitRestart) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolRestart) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolRestart) );
   SCIP_CALL( SCIPsetEventhdlrDelete(scip, eventhdlr, eventDeleteRestart) );

   /* add restart event handler parameters */
   /* TODO: (optional) add event handler specific parameters with SCIPaddTypeParam() here */
   SCIP_CALL( SCIPaddCharParam(scip, "restarts/restartpolicy", "restart policy: aenp",
            &eventhdlrdata->restartpolicyparam, FALSE, 'n', "aenp", NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "restarts/estimationmethod", "select estimation method",
               &eventhdlrdata->estimationparam, FALSE, 't', "t", NULL, NULL) );
   SCIP_CALL( SCIPaddCharParam(scip, "restarts/progressmeasure", "select progress measure",
               &eventhdlrdata->progressparam, FALSE, 'p', "p", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "restarts/restartlimit", "restart limit",
      &eventhdlrdata->restartlimit, FALSE, 1, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "restarts/minnodes", "minimum number of nodes before restart",
         &eventhdlrdata->minnodes, FALSE, 1000, -1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "restarts/countonlyleaves", "should only leaves count for the minnodes parameter?",
         &eventhdlrdata->countonlyleaves, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "restarts/estimation/factor",
         "factor by which the estimated number of nodes should exceed the current number of nodes",
         &eventhdlrdata->estim_factor, FALSE, 2.0, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
