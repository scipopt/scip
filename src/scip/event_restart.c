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

#include <string.h>
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
#define EVENTTYPE_RESTART      (SCIP_EVENTTYPE_PQNODEINFEASIBLE | SCIP_EVENTTYPE_NODEBRANCHED)

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
#define NREPORTS                    100      /**< maximum number of reports that should be generated */

#define ESTIMATION_CHAR_TREESIZE         't' /**< should estimation use probability based tree size prediction? */
#define ESTIMATION_CHAR_PROFILE          'p'  /**< should estimation use profile based prediction a la Cornuejols? */

#define PROGRESS_CHAR_RATIO               'r' /**< should the search progress be measured using ratio-based probabilities? */
#define PROGRESS_CHAR_UNIFORM             'u' /**< should the search progress be measured using even probabilities? */
#define PROGRESS_CHAR_GAP                 'g' /**< should the search progress be measured in terms of the gap? */
#define PROGRESS_CHAR_FIXED               'f' /**< should the search progress be measured using fixed, ratio based probabilities? */

#define DEFAULT_WINDOWSIZE               100  /**< window size for search progress */
#define MAX_WINDOWSIZE                   500  /**< window size for search progress */
#define DEFAULT_DES_ALPHA               0.95  /**< default level smoothing constant for double exponential smoothing */
#define DEFAULT_DES_BETA                0.10   /**< default trend smoothing constant for double exponential smoothing */
#define DEFAULT_DES_USETRENDINLEVEL      TRUE /**< should the trend be used in the level update? */

#define FORECAST_BACKTRACKESTIM           'b' /**< use backtrack estimation for forecasting */
#define FORECAST_LINEAR                   'l' /**< use linear trends based on double exponential smoothing for forecasting */
#define FORECAST_WINDOW                   'w' /**< use either linear or quadratic trends within window for forecasting */

#define TABLE_NAME              "restart"
#define TABLE_DESC              "restart statistics table"
#define TABLE_POSITION          22000           /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_INIT /**< output of the statistics table is only printed from this stage onwards */


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

/** estimation of tree size that is updated at every leaf node */
struct BacktrackEstim
{
   SCIP_Real             numerator;          /**< weighted sample sizes based on the path probability */
   SCIP_Real             denominator;        /**< sum of weights (aka progress) */
   char                  progressmethod;     /**< 'f'ixed or 'u'niform? */
};
typedef struct BacktrackEstim BACKTRACKESTIM;

/** time series data structure for leaf time series
 *
 *  these time series are the basic ingredient for tree size estimation via forecasting.
 *
 *  This general class represents concrete time series such as the closed gap, progress, and leaf frequency.
 *  Through callbacks for data (de-)initialization and value queries, it provides a common interface
 *  to which double exponential smoothing or window forecasts can be applied.
 *  */
typedef struct TimeSeries TIMESERIES;

/** data structure for convenient access of tree information */
typedef struct TreeData TREEDATA;
#define NTIMESERIES 3

/** event handler data */
struct SCIP_EventhdlrData
{
   SEARCHPROGRESS*       ratioprogress;      /**< ratio progress data structure */
   BACKTRACKESTIM*       backtrackestim;     /**< backtrack estimator for tree size */
   TIMESERIES*           timeseries[NTIMESERIES]; /**< array of time series slots */
   TREEDATA*             treedata;           /**< tree data */
   char                  restartpolicyparam; /**< restart policy parameter */
   char                  estimationparam;    /**< parameter to select the estimation method */
   char                  progressparam;      /**< progress method to use */
   char                  forecastparam;      /**< method used for forecasting */
   int                   windowsize;         /**< the window size used */
   SCIP_Bool             useacceleration;    /**< consider also acceleration within window? */
   int                   restartlimit;       /**< how often should a restart be triggered? (-1 for no limit) */
   int                   nrestartsperformed; /**< number of restarts performed so far */
   int                   restarthitcounter;  /**< the number of successive samples that would trigger a restart */
   int                   hitcounterlim;      /**< limit on the number of successive samples to really trigger a restart */
   SCIP_Longint          minnodes;           /**< minimum number of nodes in a run before restart is triggered */
   SCIP_Bool             countonlyleaves;    /**< should only leaves count for the minnodes parameter? */
   SCIP_Real             estim_factor;       /**< factor by which the estimated number of nodes should exceed the current number of nodes */
   SCIP_Real             proglastreport;     /**< progress at which last report was printed */
   int                   nreports;           /**< the number of reports already printed */
};

struct TreeData
{
   SCIP_Longint          nnodes;             /**< the total number of nodes */
   SCIP_Longint          nopen;              /**< the current number of open nodes */
   SCIP_Longint          ninner;             /**< the number of inner nodes */
   SCIP_Longint          nleaves;            /**< the number of final leaf nodes */
   SCIP_Longint          nvisited;           /**< the number of visited nodes */
   SCIP_Real             progress;           /**< the current progress (sum of leaf weights) */
};

/** data initialization callback of time series */
#define DECL_TIMESERIESINIT(x) SCIP_RETCODE x ( \
   SCIP* scip,                                  \
   TIMESERIES* ts                               \
   )

/** data deinitialization callback of time series */
#define DECL_TIMESERIESEXIT(x) SCIP_RETCODE x ( \
   SCIP* scip,                                  \
   TIMESERIES* ts                               \
   )

/** update callback of time series */
#define DECL_TIMESERIESUPDATE(x) SCIP_RETCODE x (\
   SCIP*                 scip,                   \
   TIMESERIES*           ts,                     \
   TREEDATA*             treedata,               \
   SCIP_Real*            value                   \
   )

/** time series data structure for leaf time series */
struct TimeSeries
{
   DOUBLEEXPSMOOTH       des;                /**< double exponential smoothing data structure */
   char*                 name;               /**< name of this time series */
   SCIP_Real*            vals;               /**< value array of this time series */
   SCIP_Real             targetvalue;        /**< target value of this time series */
   SCIP_Real             currentvalue;       /**< current value of time series */
   SCIP_Longint          nobs;               /**< total number of observations */
   int                   valssize;           /**< size of value array */
   int                   nvals;              /**< number of values */
   int                   resolution;         /**< current (inverse of) resolution */
   DECL_TIMESERIESINIT((*timeseriesinit));   /**< data initialization callback (at the beginning of search) */
   DECL_TIMESERIESEXIT((*timeseriesexit));   /**< data deinitialization callback */
   DECL_TIMESERIESUPDATE((*timeseriesupdate));/**< update callback at nodes */
   union {
   }                     data;               /**< time series data pointer */
};

/*
 * Local methods
 */

/** reset tree data */
static
void treedataReset(
   TREEDATA*             treedata            /**< tree data */
   )
{
   /* simply set everything to 0 */
   BMSclearMemory(treedata);

   /* set up root node */
   treedata->nnodes = 1;
   treedata->nopen = 1;
}

/** create tree data structure */
static
SCIP_RETCODE treedataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA**            treedata            /**< pointer to store tree data */
   )
{
   assert(treedata != NULL);
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, treedata) );

   treedataReset(*treedata);

   return SCIP_OKAY;
}

/** free tree data structure */
static
void treedataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA**            treedata            /**< pointer to tree data */
   )
{
   assert(scip != NULL);
   assert(treedata != NULL);

   SCIPfreeMemory(scip, treedata);
   *treedata = NULL;
}

/** todo update tree data structure after a node has been solved/is about to be deleted */
static
void treedataUpdate(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEDATA*             treedata,           /**< tree data */
   int                   depth,              /**< the depth of the node */
   int                   nchildren           /**< the number of children */
   )
{
   ++treedata->nvisited;
   treedata->nopen--;

   if( nchildren == 0 )
   {
      treedata->nleaves++;
      treedata->progress += pow(0.5, depth);
   }
   else
   {
      treedata->nnodes += nchildren;
      treedata->nopen += nchildren;
      ++treedata->ninner;
   }
}

/* print method for tree data */
static
char* treedataPrint(
   TREEDATA*             treedata,           /**< tree data */
   char*                 strbuf              /**< string buffer */
   )
{
   sprintf(strbuf,
      "Tree Data: %lld nodes (%lld visited, %lld inner, %lld leaves, %lld open), progress: %.4f",
      treedata->nnodes,
      treedata->nvisited,
      treedata->ninner,
      treedata->nleaves,
      treedata->nopen,
      treedata->progress
      );
   return strbuf;
}

/** reset double exponential smoothing */
static
void doubleexpsmoothReset(
   DOUBLEEXPSMOOTH*      des                 /**< double exponential smoothing data structure */
   )
{
  des->n = 0;
  des->level = SCIP_INVALID;
  des->trend = SCIP_INVALID;
}

/** initialize a double exponential smoothing data structure */
static
void doubleexpsmoothInit(
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
void doubleexpsmoothUpdate(
   DOUBLEEXPSMOOTH*      des,                /**< double exponential smoothing data structure */
   SCIP_Real             xnew                /**< new sample value */
   )
{
   if( des->n == 0 )
      doubleexpsmoothInit(des, xnew);
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
SCIP_Real doubleexpsmoothGetTrend(
   DOUBLEEXPSMOOTH*      des                 /**< double exponential smoothing data structure */
   )
{
   assert(des != NULL);

   if( des->n == 0 )
      return SCIP_INVALID;

   return des->trend;
}

/** reset time series */
static
void timeseriesReset(
   TIMESERIES*           timeseries          /**< pointer to store time series */
   )
{
   timeseries->resolution = 1;
   timeseries->nvals = 0;
   timeseries->nobs = 0L;
   timeseries->currentvalue = 0;

   doubleexpsmoothReset(&timeseries->des);
}



/** create a time series object */
static
SCIP_RETCODE timeseriesCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMESERIES**          timeseries,         /**< pointer to store time series */
   const char*           name,               /**< name of this time series */
   SCIP_Real             targetvalue,        /**< target value of this time series */
   DECL_TIMESERIESINIT   ((*timeseriesinit)),/**< data initialization callback (at the beginning of search), or NULL */
   DECL_TIMESERIESEXIT   ((*timeseriesexit)),/**< data deinitialization callback, or NULL */
   DECL_TIMESERIESUPDATE ((*timeseriesupdate)) /**< update callback at nodes, or NULL */
   )
{
   TIMESERIES* timeseriesptr;
   assert(scip != NULL);
   assert(timeseries != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPallocMemory(scip, timeseries) );

   timeseriesptr = *timeseries;
   assert(timeseriesptr != NULL);

   /* copy name */
   SCIP_ALLOC( BMSduplicateMemoryArray(&timeseriesptr->name, name, strlen(name)+1) );

   /* copy callbacks */
   assert(timeseriesupdate != NULL);
   timeseriesptr->timeseriesinit = timeseriesinit;
   timeseriesptr->timeseriesexit = timeseriesexit;
   timeseriesptr->timeseriesupdate = timeseriesupdate;

   timeseriesptr->targetvalue = targetvalue;
   timeseriesptr->valssize = 1024;

   SCIP_CALL( SCIPallocMemoryArray(scip, &timeseriesptr->vals, timeseriesptr->valssize) );

   timeseriesReset(timeseriesptr);

   SCIPdebugMsg(scip, "Finished creation of time series '%s'\n", timeseriesptr->name);

   return SCIP_OKAY;
}

/** free a time series */
static
void timeseriesFree(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMESERIES**          timeseries          /**< pointer to time series */
   )
{
   assert(scip != NULL);
   assert(timeseries != NULL);

   BMSfreeMemoryArray(&(*timeseries)->name);

   SCIPfreeMemoryArray(scip, &(*timeseries)->vals);

   SCIPfreeMemory(scip, timeseries);

   *timeseries = NULL;
}


/** resample to lower resolution */
static
void timeseriesResample(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   int i;

   assert(timeseries->nvals % 2 == 0);

   doubleexpsmoothReset(&timeseries->des);

   /* compress vals array to store only every second entry */
   for( i = 0; i < timeseries->nvals / 2; ++i )
   {
      timeseries->vals[i] = timeseries->vals[2 * i];
      doubleexpsmoothUpdate(&timeseries->des, timeseries->vals[i]);
   }

   timeseries->resolution *= 2;
   timeseries->nvals = timeseries->nvals / 2;
}

/** update time series */
static
SCIP_RETCODE timeseriesUpdate(
   SCIP*                 scip,               /**< SCIP data structure */
   TIMESERIES*           timeseries,         /**< time series */
   TREEDATA*             treedata,           /**< tree data */
   SCIP_Bool             isleaf              /**< are we at a leaf node? */
   )
{

   SCIP_Real value;

   assert(scip != NULL);
   assert(timeseries != NULL);
   assert(treedata != NULL);

   /* call update callback */
   assert(timeseries->timeseriesupdate != NULL);
   SCIP_CALL( timeseries->timeseriesupdate(scip, timeseries, treedata, &value) );

   if( !isleaf )
      return SCIP_OKAY;

   /* store the value as current value */
   timeseries->currentvalue = value;
   timeseries->nobs++;


   /* if this is a leaf that matches the time series resolution, store the value */
   if( timeseries->nobs % timeseries->resolution == 0 )
   {
      assert(timeseries->nvals < timeseries->valssize);
      timeseries->vals[timeseries->nvals++] = value;
      doubleexpsmoothUpdate(&timeseries->des, value);
   }

   /* if the time series has reached its capacity, resample and increase the resolution */
   if( timeseries->nvals == timeseries->valssize )
      timeseriesResample(timeseries);

   return SCIP_OKAY;
}

/** get current value of time series */
static
SCIP_Real timeseriesGet(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   assert(timeseries != NULL);

   return timeseries->currentvalue;
}

/** get name of time series */
static
char* timeseriesGetName(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->name;
}

/** get target value (which this time series reaches at the end of the solution process) */
static
SCIP_Real timeseriesGetTargetValue(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->targetvalue;
}

/** get resolution of time series */
static
int timeseriesGetResolution(
   TIMESERIES*           timeseries          /**< time series */
   )
{
   return timeseries->resolution;
}

/** estimate tree size at which time series reaches target value */
static
SCIP_Real timeseriesEstimate(
   TIMESERIES*           timeseries,         /**< time series */
   TREEDATA*             treedata            /**< tree data for fallback estimation */
   )
{
   SCIP_Real val;
   SCIP_Real targetval;
   SCIP_Real trend;

   /* if no observations have been made yet, return infinity */
   if( timeseries->nobs == 0L )
      return -1.0;

   val = timeseriesGet(timeseries);
   targetval = timeseriesGetTargetValue(timeseries);

   /* if the value has reached the target value already, return the number of observations */
   if( EPSZ(val - targetval, 1e-6) )
      return 2.0 * timeseries->nobs - 1;


   trend = doubleexpsmoothGetTrend(&timeseries->des);
   /* get current value and trend. The linear trend estimation may point into the wrong direction
    * In this case, we use the fallback mechanism that we will need twice as many nodes.
    */
   if( (targetval > val && trend < 1e-6) || (targetval < val && trend > -1e-6) )
   {
      return 2.0 * treedata->nvisited;
   }

   /* compute after how many additional steps the current trend reaches the target value; multiply by resolution */
   return 2.0 * timeseriesGetResolution(timeseries) * (timeseries->nvals + (targetval - val) / (SCIP_Real)trend) - 1.0;
}


/* put your local methods here, and declare them static */



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
   SCIP_ALLOC( BMSallocMemoryArray(&(*progress)->progressarray, MAX_WINDOWSIZE) );
   SCIP_ALLOC( BMSallocMemoryArray(&(*progress)->resourcearray, MAX_WINDOWSIZE) );

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
   progress->curr = (progress->curr + 1) % MAX_WINDOWSIZE;
   progress->progressarray[progress->curr] = obs;
   progress->resourcearray[progress->curr] = res;

   doubleexpsmoothUpdate(&progress->desprogress, obs);
   doubleexpsmoothUpdate(&progress->desresources, res);
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

   assert(0 <= progress->curr && progress->curr <= MAX_WINDOWSIZE - 1);

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
   /*SCIP_Real resourcetrend;*/
   SCIP_Real remprogress;
   SCIP_Real remleaves;
   SCIP_Real totalleaves;

   assert(progress != NULL);

   remprogress = targetlevel - getCurrentProgress(progress);

   /* we have already reached the target level */
   if( remprogress <= 0.0 )
      return 0.0;

   /* no observation available yet */
   if( progress->nobservations == 0 )
      return SCIP_REAL_MAX;

   progresstrend = doubleexpsmoothGetTrend(&progress->desprogress);
   /*resourcetrend = getTrendDoubleexpsmooth(&progress->desresources);*/

   if( progresstrend == 0.0 )
      return SCIP_REAL_MAX;

   /* the remaining progress to the target level will be reached in approximately remprogress /progresstrend
    * many samples. The corresponding resource trend per time step yields the remaining ressources
    */
   remleaves = remprogress / progresstrend;
   totalleaves = remleaves + progress->nobservations;

   /* the total number of nodes is the 2 * N (leave number) - 1 */
   return 2 * totalleaves - 1 - getCurrentResources(progress);
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
   SCIP_Real remprogress;
   int windowstart;
   int windowend;

   assert(progress != NULL);
   windowsize = MIN(windowsize, progress->nobservations);

   /* we need at least 3 observations in our window to compute the acceleration */
   useacceleration = useacceleration && windowsize >= 3;

   remprogress = targetlevel - getCurrentProgress(progress);

   /* nothing to forecast anymore */
   if( remprogress <= 0.0 )
      return 0.0;

   windowend = progress->curr;
   assert(progress->curr == (progress->nobservations - 1) % MAX_WINDOWSIZE);

   /* compute the start index of the window */
   if( progress->nobservations > MAX_WINDOWSIZE )
      windowstart = (progress->curr - windowsize + 1 + MAX_WINDOWSIZE) % MAX_WINDOWSIZE;
   else
      windowstart = progress->curr - windowsize + 1;

   assert(windowstart >= 0);

   /* try to compute remaining ressources as the root of a quadratic function
    *
    * s(r) = s_0 + v * r + .5 a * r^2
    *
    * where s_0, v, and a are computed by using the start, end, and midpoint of the current
    * window.
    * */
   if( useacceleration )
   {
      SCIP_Real rootdiscriminant;
      SCIP_Real remres1;
      SCIP_Real remres2;
      SCIP_Real v;
      SCIP_Real s0;
      int windowmid = ((windowstart + windowsize) / 2) % MAX_WINDOWSIZE;
      SCIP_Real w1 = progress->resourcearray[windowstart];
      SCIP_Real w3 = progress->resourcearray[windowend];
      SCIP_Real w2 = progress->resourcearray[windowmid];
      SCIP_Real vel1 = measureVelocity(progress, windowstart, windowmid);
      SCIP_Real velwindow = measureVelocity(progress, windowstart, windowend);
      SCIP_Real discriminant;

      /* coefficient a, the acceleration, in the above formula */
      SCIP_Real acceleration = (velwindow - vel1) / (w3 - w2) * 2.0;

      /* coefficient v, the velocity, and s_0, the y intercept in the quadratic function */
      v = vel1 - .5 * acceleration * (w1 + w2);
      s0 = progress->progressarray[windowstart] - v * w1 - .5 * acceleration * w1 * w1;

      if( ! EPSZ(acceleration, 1e-9) )
      {
         /* solve the quadratic equation s(r) = targetlevel = s_0 + v * r + 0.5 * a * r^2
          *
          * r1/2 = 2 / a * (-v +/- sqrt(v^2 - 2 * a * (s_0 - targetlevel)))
          * */
         discriminant = v * v - 2 * acceleration * (s0 - targetlevel);
         discriminant = MAX(0, discriminant);
         rootdiscriminant = sqrt(discriminant);
         remres1 = (-v + rootdiscriminant) / acceleration;
         remres2 = (-v - rootdiscriminant) / acceleration;

         return MAX(remres1, remres2);
      }
      else
      {
         /* solve the linear displacement formula because the acceleration is 0 */
         return remprogress / v;
      }
   }
   else
   {
      SCIP_Real velocitywindow = measureVelocity(progress, windowstart, windowend);

      return remprogress / velocitywindow;
   }
}

/** reset a backtrack estimator */
static
void resetBacktrackestim(
   BACKTRACKESTIM*       backtrackestim      /**< backtrack estimator */
   )
{
   BMSclearMemory(backtrackestim);

   return;
}

/** create a backtrack estimator */
static
SCIP_RETCODE createBacktrackestim(
   BACKTRACKESTIM**      backtrackestim,     /**< pointer to store backtrack estimator */
   char                  progressmethod      /**< 'f'ixed or 'u'niform? */
   )
{
   assert(backtrackestim != NULL);

   SCIP_ALLOC( BMSallocMemory(backtrackestim) );

   resetBacktrackestim(*backtrackestim);

   return SCIP_OKAY;
}

/** free a backtrack estimator */
static
void freeBacktrackestim(
   BACKTRACKESTIM**      backtrackestim      /**< pointer to store backtrack estimator */
   )
{
   assert(backtrackestim != NULL);

   if( *backtrackestim == NULL )
      return;

   BMSfreeMemory(backtrackestim);
}

/** update backtrack estimator by a new leaf node */
static
void updateBacktrackestim(
   BACKTRACKESTIM*       backtrackestim,     /**< backtrack estimator */
   SCIP_NODE*            leafnode            /**< new, previously unseen leaf node */
   )
{
   SCIP_Real probability;
   SCIP_Real num;
   SCIP_Real arcprobability;
   SCIP_Real pathprobability;
   SCIP_NODE* parent;
   SCIP_NODE* current;


   assert(backtrackestim != NULL);
   assert(leafnode != NULL);

   switch (backtrackestim->progressmethod) {
      case PROGRESS_CHAR_FIXED:
         probability = SCIPnodeGetFixedProbability(leafnode);
         pathprobability = probability;

         current = leafnode;
         num = 1.0;
         /* loop back along all arcs along the path */
         while( (parent = SCIPnodeGetParent(current)) != NULL )
         {
            arcprobability = SCIPnodeGetFixedProbability(current) / SCIPnodeGetFixedProbability(parent);
            num += probability / pathprobability;
            pathprobability /= arcprobability;

            current = parent;
         }
         break;
      case PROGRESS_CHAR_UNIFORM:
         probability = pow(0.5, SCIPnodeGetDepth(leafnode));
         num = 2 - probability;
         break;
      default:
         SCIPerrorMessage("Unsupported progress type <%c> for backtrack estimation\n");
         SCIPABORT();

         break;
   }

   backtrackestim->numerator += num;
   backtrackestim->denominator += probability;
}

/** estimate the total tree size using the backtrack estimation */
static
SCIP_Real estimateTreesizeBacktrackestim(
   BACKTRACKESTIM*       backtrackestim      /**< backtrack estimator */
   )
{
   assert(backtrackestim != NULL);

   if( backtrackestim->denominator == 0.0 )
      return -1;

   return backtrackestim->numerator / backtrackestim->denominator;
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


/** free all time series */
static
void freeTimeseries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   TIMESERIES** tss = eventhdlrdata->timeseries;
   int t;

   /* loop over time series and reset them */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      assert(tss[t] != NULL);
      timeseriesFree(scip, &tss[t]);
   }
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeRestart)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   freeSearchprogress(&eventhdlrdata->ratioprogress);

   freeBacktrackestim(&eventhdlrdata->backtrackestim);

   treedataFree(scip, &eventhdlrdata->treedata);

   freeTimeseries(scip, eventhdlrdata);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitRestart)
{  /*lint --e{715}*/


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


/** todo reset all time series */
static
void resetTimeseries(
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   TIMESERIES** tss = eventhdlrdata->timeseries;
   int t;

   /* loop over time series and reset them */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      assert(tss[t] != NULL);
      timeseriesReset(tss[t]);
   }
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolRestart)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   resetSearchprogress(eventhdlrdata->ratioprogress);

   resetBacktrackestim(eventhdlrdata->backtrackestim);

   /* backtrack estimator only allows fixed or uniform progress */
   if( eventhdlrdata->progressparam == PROGRESS_CHAR_FIXED )
      eventhdlrdata->backtrackestim->progressmethod = PROGRESS_CHAR_FIXED;
   else
      eventhdlrdata->backtrackestim->progressmethod = PROGRESS_CHAR_UNIFORM;

   eventhdlrdata->restarthitcounter = 0;
   eventhdlrdata->proglastreport = 0.0;
   eventhdlrdata->nreports = 0;

   /* reset tree data */
   treedataReset(eventhdlrdata->treedata);

   resetTimeseries(eventhdlrdata);

   SCIP_CALL( SCIPcatchEvent(scip, EVENTTYPE_RESTART, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolRestart)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPdropEvent(scip, EVENTTYPE_RESTART, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

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

   return RESTARTPOLICY_NEVER;
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

/** forecast the number of remaining nodes depending on the selected user parameters */
static
SCIP_Real forecastRemainingNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{

   switch (eventhdlrdata->forecastparam) {
      case FORECAST_BACKTRACKESTIM:
         return MAX(0.0, estimateTreesizeBacktrackestim(eventhdlrdata->backtrackestim) - SCIPgetNNodes(scip));
         break;
      case FORECAST_LINEAR:
         return forecastRemainingResources(eventhdlrdata->ratioprogress, 1.0);
         break;
      case FORECAST_WINDOW:
         return forecastRollingAverageWindow(eventhdlrdata->ratioprogress, 1.0,
                  eventhdlrdata->windowsize,
                  eventhdlrdata->useacceleration);
         break;
      default:
         break;
   }

   return -1.0;
}

/** should a restart be applied based on the current progress? */
static
SCIP_Bool shouldApplyRestartProgress(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{

   SCIP_Real estimation;
   SCIP_Real remainnodes;

   remainnodes = forecastRemainingNodes(scip, eventhdlrdata);

   if( remainnodes < 0.0 )
      return FALSE;

   estimation = SCIPgetNNodes(scip) + remainnodes;

   /* if the estimation exceeds the current number of nodes by a dramatic factor, restart */
   if( estimation > SCIPgetNNodes(scip) * eventhdlrdata->estim_factor )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "Estimation %g exceeds current number of nodes %lld by a factor of %.1f\n",
               estimation, SCIPgetNNodes(scip), estimation / SCIPgetNNodes(scip));
      return TRUE;
   }

   return FALSE;



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
         currentprogress = getCurrentProgress(searchprogress) + pow(0.5, SCIPnodeGetDepth(leafnode));

         break;
      case PROGRESS_CHAR_RATIO:
         SCIP_CALL( SCIPgetNodeProbability(scip, leafnode, &currentprogress) );
         currentprogress += getCurrentProgress(searchprogress);
         break;
      case PROGRESS_CHAR_FIXED:
         currentprogress = getCurrentProgress(searchprogress) + SCIPnodeGetFixedProbability(leafnode);
         break;
      default:
         break;
   }

   /* author bzfhende
    *
    * TODO add different resource types than nodes
    */
   addSampleSearchprogress(searchprogress, currentprogress, SCIPgetNNodes(scip));

   updateBacktrackestim(eventhdlrdata->backtrackestim, leafnode);

   SCIPdebugMsg(scip, "Update search progress by leaf %lld at depth %d: %g\n",
      SCIPnodeGetNumber(leafnode), SCIPnodeGetDepth(leafnode), pow(0.5, SCIPnodeGetDepth(leafnode)));

   return SCIP_OKAY;
}

/** update all time series */
static
SCIP_RETCODE updateTimeseries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   TREEDATA*             treedata,           /**< tree data */
   SCIP_Bool             isleaf              /**< are we at a leaf node? */
   )
{
   TIMESERIES** tss = eventhdlrdata->timeseries;
   int t;

   /* loop over time series */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      assert(tss[t] != NULL);
      timeseriesUpdate(scip, tss[t], treedata, isleaf);

#if 0
      SCIPdebugMsg(scip,
         "Update of time series '%s', current value %.4f (%lld observations)\n",
         timeseriesGetName(tss[t]), timeseriesGet(tss[t]), tss[t]->nobs);
#endif
   }

   return SCIP_OKAY;
}

/** print a treesize estimation report into the string buffer */
static
char* printReport(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata,      /**< event handler data */
   char*                 strbuf,             /**< string buffer */
   int                   reportnum           /**< report number, or 0 to omit number */
   )
{
   TREEDATA* treedata = eventhdlrdata->treedata;
   char* ptr = strbuf;
   int t;

   /* print report number */
   if( reportnum > 0 )
      ptr += sprintf(ptr, "Report %d\n", reportnum);

   /* print tree data */
   ptr += sprintf(ptr,
         "  %-17s: %lld nodes (%lld visited, %lld inner, %lld leaves, %lld open), progress: %.4f\n",
         "Tree Data",
         treedata->nnodes,
         treedata->nvisited,
         treedata->ninner,
         treedata->nleaves,
         treedata->nopen,
         treedata->progress
         );

   /* print estimations */
   ptr += sprintf(ptr, "Tree Estimation    : %11s %11s %11s",
            "estim",
            "value",
            "trend");

   ptr += sprintf(ptr, "\n");

   ptr += sprintf(ptr, "  wbe              : %11.0f %11s %11s\n",
            estimateTreesizeBacktrackestim(eventhdlrdata->backtrackestim), "-", "-");
   ptr += sprintf(ptr, "  tree profile     : %11.0f %11s %11s\n",
            SCIPpredictTotalSizeTreeprofile(scip),
            "-", "-");

   /* print time series forecasts */
   for( t = 0; t < NTIMESERIES; ++t )
   {
      TIMESERIES* ts = eventhdlrdata->timeseries[t];
      ptr += sprintf(ptr, "  %-17s: %11.0f %11.5f %11.5f\n",
            timeseriesGetName(ts), timeseriesEstimate(ts, eventhdlrdata->treedata), timeseriesGet(ts), doubleexpsmoothGetTrend(&ts->des));
   }

   return strbuf;
}


/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecRestart)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA*   eventhdlrdata;
   SCIP_Bool isleaf;
   SCIP_Bool isleafbit;
   SCIP_EVENTTYPE eventtype;
   TREEDATA* treedata;

   assert(scip != NULL);
   assert(eventhdlr != NULL);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);
   eventtype = SCIPeventGetType(event);
   treedata = eventhdlrdata->treedata;

   if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED || eventtype == SCIP_EVENTTYPE_PQNODEINFEASIBLE )
   {
      int nchildren = 0;

      if( eventtype == SCIP_EVENTTYPE_NODEBRANCHED )
         nchildren = SCIPgetNChildren(scip);

      treedataUpdate(scip, treedata, SCIPnodeGetDepth(SCIPeventGetNode(event)), nchildren);

#ifdef SCIP_DEBUG
      {
         char strbuf[SCIP_MAXSTRLEN];
         SCIPdebugMsg(scip, "%s\n", treedataPrint(treedata, strbuf));
      }
#endif

      SCIP_CALL( updateTimeseries(scip, eventhdlrdata, treedata, nchildren == 0) );

      /* should a new report be printed? */
      if( treedata->progress >= eventhdlrdata->proglastreport + 1.0 / (SCIP_Real)NREPORTS )
      {
         char strbuf[SCIP_MAXSTRLEN];

         SCIPinfoMessage(scip, NULL, "%s\n", printReport(scip, eventhdlrdata, strbuf, ++eventhdlrdata->nreports));

         eventhdlrdata->proglastreport = 1 / (SCIP_Real)NREPORTS * (int)(treedata->progress * NREPORTS);
      }
   }

   if( eventtype != SCIP_EVENTTYPE_PQNODEINFEASIBLE )
      return SCIP_OKAY;

   /* update the search progress at leaf nodes */
   isleaf = (eventtype == SCIP_EVENTTYPE_PQNODEINFEASIBLE);

   isleafbit = 0 != (eventtype & (SCIP_EVENTTYPE_PQNODEINFEASIBLE));

   if( isleaf != isleafbit )
   {
      SCIPABORT();
   }



   if( isleaf )
   {
      SCIP_Real remainnodes;
      SCIP_CALL( updateSearchProgress(scip, eventhdlrdata, SCIPeventGetNode(event)) );

      remainnodes = forecastRemainingNodes(scip, eventhdlrdata);
      SCIPdebugMsg(scip, "Updated search progress to %.8f tree size estimation %g (%lld + %g)\n",
         getCurrentProgress(eventhdlrdata->ratioprogress),
         SCIPgetNNodes(scip) + remainnodes,
         SCIPgetNNodes(scip), remainnodes);
   }

   /* if nodes have been pruned, this is usually an indication that things are progressing, don't restart */
   if( eventtype & SCIP_EVENTTYPE_PQNODEINFEASIBLE )
   {
      SCIPdebugMsg(scip, "PQ node %lld (depth %d) infeasible, isleaf: %u\n",
               SCIPnodeGetNumber(SCIPeventGetNode(event)),
               SCIPnodeGetDepth(SCIPeventGetNode(event)),
               isleaf);
   }

   /* check if all conditions are met such that the event handler should run */
   if( ! checkConditions(scip, eventhdlrdata) )
      return SCIP_OKAY;


   /* author bzfhende
    *
    * test if a restart should be applied
    */

   if( shouldApplyRestart(scip, eventhdlrdata) )
   {
      eventhdlrdata->restarthitcounter++;

      if( eventhdlrdata->restarthitcounter >= eventhdlrdata->hitcounterlim )
      {
         eventhdlrdata->nrestartsperformed++;

         SCIP_CALL( SCIPrestartSolve(scip) );
      }
   }
   else
   {
      eventhdlrdata->restarthitcounter = 0;
   }

   return SCIP_OKAY;
}

/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputRestart)
{  /*lint --e{715}*/
   SCIP_EVENTHDLR* eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   char strbuf[SCIP_MAXSTRLEN];

   assert(eventhdlr != NULL);

   SCIPinfoMessage(scip, file, "%s\n", printReport(scip, eventhdlrdata, strbuf, 0));

   return SCIP_OKAY;
}

/** query callback for current value */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateGap)
{
   SCIP_Real primalbound;
   SCIP_Real dualbound;

   assert(scip != NULL);
   assert(ts != NULL);
   assert(value != NULL);

   /* avoid to call SCIPgetDualbound during a restart where the queue is simply emptied */
   if( SCIPisInRestart(scip) )
   {
      *value = timeseriesGet(ts);

      return SCIP_OKAY;
   }

   primalbound = SCIPgetPrimalbound(scip);
   dualbound = SCIPgetDualbound(scip);
   if( SCIPisInfinity(scip, REALABS(primalbound)) || SCIPisInfinity(scip, REALABS(dualbound)) )
      *value = 0;
   else if( SCIPisEQ(scip, primalbound, dualbound) )
      *value = 1.0;
   else
      *value = 1.0 - REALABS(primalbound - dualbound)/MAX(REALABS(primalbound), REALABS(dualbound));

   /* using this max, we set the closed gap to 0 in the case where the primal and dual bound differ in their sign */
   *value = MAX(*value, 0.0);

   return SCIP_OKAY;
}

/** query callback for current value */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateProgress)
{
   *value = treedata->progress;

   return SCIP_OKAY;
}

/** query callback for current value */
static
DECL_TIMESERIESUPDATE(timeseriesUpdateLeaffreq)
{
   if( treedata->nvisited == 0 )
      *value = -0.5;
   else
      *value = (treedata->nleaves - 0.5)/(SCIP_Real)treedata->nvisited;

   return SCIP_OKAY;
}




/** include time series to forecast into event handler */
static
SCIP_RETCODE includeTimeseries(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);

   /* include gap time series */
   SCIP_CALL( timeseriesCreate(scip, &eventhdlrdata->timeseries[0], "gap", 1.0,
         NULL, NULL, timeseriesUpdateGap) );


   /* include progress time series */
   SCIP_CALL( timeseriesCreate(scip, &eventhdlrdata->timeseries[1], "progress", 1.0,
      NULL, NULL, timeseriesUpdateProgress) );

   /* include leaf time series */
   SCIP_CALL( timeseriesCreate(scip, &eventhdlrdata->timeseries[2], "leaf-frequency", 0.5,
         NULL, NULL, timeseriesUpdateLeaffreq) );

   /* author bzfhende
    *
    * TODO include SSG time series
    */




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

   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   /* author bzfhende
    *
    * TODO add parameters
    */

   BMSclearMemory(eventhdlrdata);

   SCIP_CALL( createSearchprogress(&eventhdlrdata->ratioprogress) );

   SCIP_CALL( createBacktrackestim(&eventhdlrdata->backtrackestim, PROGRESS_CHAR_UNIFORM) );

   SCIP_CALL( treedataCreate(scip, &eventhdlrdata->treedata) );

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
               &eventhdlrdata->progressparam, FALSE, 'u', "fgru", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "restarts/restartlimit", "restart limit",
      &eventhdlrdata->restartlimit, FALSE, 1, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip, "restarts/minnodes", "minimum number of nodes before restart",
         &eventhdlrdata->minnodes, FALSE, 1000, -1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "restarts/countonlyleaves", "should only leaves count for the minnodes parameter?",
         &eventhdlrdata->countonlyleaves, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "restarts/estimation/factor",
         "factor by which the estimated number of nodes should exceed the current number of nodes",
         &eventhdlrdata->estim_factor, FALSE, 2.0, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "restarts/forecast", "method used for forecasting",
         &eventhdlrdata->forecastparam, FALSE, FORECAST_LINEAR, "blw", NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "restarts/windowsize", "the window size for window forecasting",
         &eventhdlrdata->windowsize, FALSE, DEFAULT_WINDOWSIZE, 2, MAX_WINDOWSIZE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "restarts/useacceleration", "consider also acceleration within window?",
         &eventhdlrdata->useacceleration, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "restarts/hitcounterlim", "limit on the number of successive samples to really trigger a restart",
         &eventhdlrdata->hitcounterlim, FALSE, 50, 1, INT_MAX, NULL, NULL) );


   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         NULL, NULL, NULL, NULL,
         NULL, NULL, tableOutputRestart,
         NULL, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   /* author bzfhende
    *
    * TODO include time series into event handler
    */
   SCIP_CALL( includeTimeseries(scip, eventhdlrdata) );

   return SCIP_OKAY;
}
