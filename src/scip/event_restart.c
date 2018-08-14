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
 * @brief  event handler for restart event
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

/** enumerator for available restart policies */
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

#define PROGRESS_CHAR_RATIO               'r' /**< should the search progress be measured using ratio-based probabilities? */
#define PROGRESS_CHAR_UNIFORM             'u' /**< should the search progress be measured using even probabilities? */
#define PROGRESS_CHAR_GAP                 'g' /**< should the search progress be measured in terms of the gap? */
#define DEFAULT_WINDOWSIZE               100  /**< window size for search progress */
#define DEFAULT_DES_ALPHA               0.15  /**< default level smoothing constant for double exponential smoothing */
#define DEFAULT_DES_BETA                0.15   /**< default trend smoothing constant for double exponential smoothing */
#define DEFAULT_DES_USETRENDINLEVEL      TRUE /**< should the trend be used in the level update? */

/** double exponential smoothing data structure */
struct DoubleExpSmooth
{
   SCIP_Real             alpha;              /**< level smoothing constant */
   SCIP_Real             beta;               /**< trend smoothing constant */
   SCIP_Real             level;              /**< estimation of the current level used for smoothing */
   SCIP_Real             trend;              /**< estimation of the current trend (slope) */
   SCIP_Bool             usetrendinlevel;    /**< should the trend be used in the level update? */
   int                   n;                  /**< number of observations */
};
typedef struct DoubleExpSmooth DOUBLEEXPSMOOTH;

/** data structure to hold the search progress */
struct SearchProgress
{
   SCIP_Real*            progressarray;       /**< captures the current search progress in an array */
   SCIP_Real*            resourcearray;       /**< captures the resource measurements, e.g., nodes */
   int                   curr;                /**< index of current element */
   int                   nobservations;       /**< total number of training observations */
   DOUBLEEXPSMOOTH       desprogress;         /**< double exponential smoothing data structure for progress */
   DOUBLEEXPSMOOTH       desresources;        /**< double exponential smoothing data structure for resources */
};

typedef struct SearchProgress SEARCHPROGRESS;

/** event handler data */
struct SCIP_EventhdlrData
{
   SEARCHPROGRESS*       ratioprogress;      /**< ratio progress data structure */
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

/** initialize a double exponential smoothing data structure */
static
void initDoubleexpsmooth(
   DOUBLEEXPSMOOTH*      des,                /**< double exponential smoothing data structure */
   SCIP_Real             x1                  /**< the first sample value */
   )
{
   assert(des != NULL);

   des->n = 1;
   des->level = x1;
   des->trend = x1;

   /* author bzfhende
    *
    * TODO make these user parameters
    */
   des->alpha = DEFAULT_DES_ALPHA;
   des->beta = DEFAULT_DES_BETA;
   des->usetrendinlevel = DEFAULT_DES_USETRENDINLEVEL;

   return;
}

/** update a double exponential smoothing data structure */
static
void updateDoubleexpsmooth(
   DOUBLEEXPSMOOTH*      des,                /**< double exponential smoothing data structure */
   SCIP_Real             xnew                /**< new sample value */
   )
{
   if( des->n == 0 )
      initDoubleexpsmooth(des, xnew);
   else
   {
      SCIP_Real newlevel;
      SCIP_Real newtrend;

      newlevel = des->alpha * xnew + (1.0 - des->alpha) * (des->level + des->usetrendinlevel ? des->trend : 0.0);
      newtrend = des->beta * (newlevel - des->level) + (1.0 - des->beta) * des->trend;

      des->level = newlevel;
      des->trend = newtrend;
   }
}

/** get the current trend (slope) computed by this double exponential smoothing */
static
SCIP_Real getTrendDoubleexpsmooth(
   DOUBLEEXPSMOOTH*      des                 /**< double exponential smoothing data structure */
   )
{
   assert(des != NULL);

   if( des->n == 0 )
      return SCIP_INVALID;

   return des->trend;
}

/** reset search progress */
static
void resetSearchprogress(
   SEARCHPROGRESS*       progress            /**< search progress data structure */
   )
{
   progress->curr = -1;
   progress->nobservations = 0;
   progress->desprogress.n = 0;
   progress->desresources.n = 0;
}

/** create a search progress */
static
SCIP_RETCODE createSearchprogress(
   SEARCHPROGRESS**      progress            /**< pointer to store search progress data structure */
   )
{
   assert(progress != NULL);

   SCIP_ALLOC( BMSallocMemory(progress) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*progress)->progressarray, DEFAULT_WINDOWSIZE) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*progress)->resourcearray, DEFAULT_WINDOWSIZE) );

   resetSearchprogress(*progress);

   return SCIP_OKAY;
}

/** free search progress */
static
void freeSearchprogress(
   SEARCHPROGRESS**      progress            /**< pointer to search progress data structure */
   )
{
   assert(progress != NULL);

   if( *progress == NULL )
      return;

   BMSfreeMemoryArray(&(*progress)->progressarray);
   BMSfreeMemoryArray(&(*progress)->resourcearray);
   BMSfreeMemory(progress);
}

/** add a new sample to the search progress */
static
void addSampleSearchprogress(
   SEARCHPROGRESS*       progress,           /**< search progress data structure */
   SCIP_Real             obs,                /**< new observation */
   SCIP_Real             res                 /**< total resources, e.g., nodes, to reach observation */
   )
{
   assert(progress != NULL);
   progress->nobservations++;
   progress->curr = (progress->curr + 1) % DEFAULT_WINDOWSIZE;
   progress->progressarray[progress->curr] = obs;
   progress->resourcearray[progress->curr] = res;

   updateDoubleexpsmooth(&progress->desprogress, obs);
   updateDoubleexpsmooth(&progress->desresources, res);
}

/** get the current search progress */
static
SCIP_Real getCurrentProgress(
   SEARCHPROGRESS*       progress            /**< search progress data structure */
   )
{
   assert(progress != NULL);
   if( progress->curr == -1 )
      return 0.0;

   assert(0 <= progress->curr && progress->curr <= DEFAULT_WINDOWSIZE - 1);

   return progress->progressarray[progress->curr];
}

/** get the current resource measurement */
static
SCIP_Real getCurrentResources(
   SEARCHPROGRESS*       progress            /**< search progress data structure */
   )
{
   if( progress->curr == -1 )
      return 0.0;

   return progress->resourcearray[progress->curr];
}

/** forecast how many additional resources are necessary to reach a certain level of progress */
static
SCIP_Real forecastRemainingResources(
   SEARCHPROGRESS*       progress,           /**< search progress data structure */
   SCIP_Real             targetlevel         /**< targeted progress level, e.g., 1.0 to finish the search */
   )
{
   SCIP_Real progresstrend;
   SCIP_Real resourcetrend;
   SCIP_Real remprogress;

   assert(progress != NULL);

   remprogress = targetlevel - getCurrentProgress(progress);

   /* we have already reached the target level */
   if( remprogress <= 0.0 )
      return 0.0;

   /* no observation available yet */
   if( progress->nobservations == 0 )
      return SCIP_REAL_MAX;

   progresstrend = getTrendDoubleexpsmooth(&progress->desprogress);
   resourcetrend = getTrendDoubleexpsmooth(&progress->desresources);

   if( progresstrend == 0.0 )
      return SCIP_REAL_MAX;

   /* the remaining progress to the target level will be reached in approximately remprogress /progresstrend
    * many samples. The corresponding resource trend per time step yields the remaining ressources
    */
   return remprogress * resourcetrend / progresstrend;
}

/** measure the velocity between the indices at t1 and t2 */
static
SCIP_Real measureVelocity(
   SEARCHPROGRESS*       progress,           /**< search progress data structure */
   int                   t1,                 /**< the earlier time index */
   int                   t2                  /**< the later time index */
   )
{
   return (progress->progressarray[t2] - progress->progressarray[t1]) / (progress->resourcearray[t2] - progress->resourcearray[t1]);
}

/** forecast how many additional resources are needed to reach a target level by using a moving window */
static
SCIP_Real forecastRollingAverageWindow(
   SEARCHPROGRESS*       progress,           /**< search progress data structure */
   SCIP_Real             targetlevel,        /**< targeted progress level, e.g., 1.0 to finish the search */
   int                   windowsize,         /**< the size of the moving window */
   SCIP_Bool             useacceleration     /**< should the acceleration within the window in speed be taken into account? */
   )
{
   assert(progress != NULL);
   SCIP_Real remprogress;
   int windowstart;
   int windowend;

   windowsize = MIN(windowsize, progress->nobservations);

   remprogress = targetlevel - getCurrentProgress(progress);
   windowend = progress->curr;
   assert(progress->curr == (progress->nobservations - 1) % DEFAULT_WINDOWSIZE);
   if( progress->nobservations > DEFAULT_WINDOWSIZE )
      windowstart = (progress->curr - windowsize + 1 + DEFAULT_WINDOWSIZE) % DEFAULT_WINDOWSIZE;
   else
      windowstart = progress->curr - windowsize + 1;

   assert(windowstart >= 0);

   if( useacceleration )
   {
      SCIP_Real rootdiscriminant;
      SCIP_Real remres1;
      SCIP_Real remres2;
      int windowmid = (windowstart + windowsize / 2) % DEFAULT_WINDOWSIZE;
      SCIP_Real vel1 = measureVelocity(progress, windowstart, windowmid);
      SCIP_Real vel2 = measureVelocity(progress, windowmid, windowend);
      SCIP_Real acceleration = (vel2 - vel1) / progress->resourcearray[windowend] - progress->resourcearray[windowmid];
      SCIP_Real discriminant = vel2 * vel2 + 2 * acceleration * remprogress;
      discriminant = MAX(0, discriminant);
      rootdiscriminant = sqrt(discriminant);
      remres1 = (-vel2 + rootdiscriminant) / acceleration;
      remres2 = (-vel2 - rootdiscriminant) / acceleration;

      return MAX(remres1, remres2);
   }
   else
   {
      SCIP_Real velocitywindow = measureVelocity(progress, windowstart, windowend);

      return remprogress / velocitywindow;
   }



   return 0.0;
}

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

   freeSearchprogress(&eventhdlrdata->ratioprogress);

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
static
SCIP_DECL_EVENTINITSOL(eventInitsolRestart)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   resetSearchprogress(eventhdlrdata->ratioprogress);

   return SCIP_OKAY;
}

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

/** should a restart be applied based on the current progress? */
static
SCIP_Bool shouldApplyRestartProgress(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{



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
          * TODO still needs to be implemented
          */
         return shouldApplyRestartProgress(scip, eventhdlrdata);
         break;
      default:
         break;
   }

   return FALSE;
}

/** update the search progress after a new leaf has been reached */
static
SCIP_RETCODE updateSearchProgress(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   SCIP_NODE*            leafnode            /**< current leaf node of SCIP */
   )
{
   SEARCHPROGRESS* searchprogress;
   SCIP_Real currentprogress;

   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   searchprogress = eventhdlrdata->ratioprogress;

   switch (eventhdlrdata->progressparam) {
      case PROGRESS_CHAR_GAP:
         currentprogress = 1.0 - MIN(SCIPgetGap(scip), 1.0);
         break;
      case PROGRESS_CHAR_UNIFORM:
         currentprogress = getCurrentProgress(searchprogress) + pow(2.0, -SCIPnodeGetDepth(leafnode));

         break;
      case PROGRESS_CHAR_RATIO:
         SCIP_CALL( SCIPgetNodeProbability(scip, leafnode, &currentprogress) );
         currentprogress += getCurrentProgress(searchprogress);
         break;
      default:
         break;
   }

   /* author bzfhende
    *
    * TODO add different resource types than nodes
    */
   addSampleSearchprogress(searchprogress, currentprogress, SCIPgetNNodes(scip));

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecRestart)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA*   eventhdlrdata;
   SCIP_Bool isleaf;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* update the search progress at leaf nodes */
   isleaf = (SCIPeventGetType(event) & (SCIP_EVENTTYPE_NODEFEASIBLE | SCIP_EVENTTYPE_NODEINFEASIBLE));
   if( isleaf )
   {
      SCIP_Real remainnodes;
      SCIP_CALL( updateSearchProgress(scip, eventhdlrdata, SCIPeventGetNode(event)) );

      remainnodes = forecastRemainingResources(eventhdlrdata->ratioprogress, 1.0);
      SCIPdebugMsg(scip, "Updated search progress to %.8f tree size estimation %g (%lld + %g)\n",
         getCurrentProgress(eventhdlrdata->ratioprogress),
         SCIPgetNNodes(scip) + remainnodes,
         SCIPgetNNodes(scip), remainnodes);
      SCIPdebugMsg(scip, "Remaining nodes forecast at node %lld: %g %g %g (window sizes 2, 10, 100)\n",
         SCIPgetNNodes(scip),
         forecastRollingAverageWindow(eventhdlrdata->ratioprogress, 1.0, 2, FALSE),
         forecastRollingAverageWindow(eventhdlrdata->ratioprogress, 1.0, 10, FALSE),
         forecastRollingAverageWindow(eventhdlrdata->ratioprogress, 1.0, 100, FALSE)
         );
   }
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

   SCIP_CALL( createSearchprogress(&eventhdlrdata->ratioprogress) );

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
               &eventhdlrdata->progressparam, FALSE, 'r', "gru", NULL, NULL) );

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
