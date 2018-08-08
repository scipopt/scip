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

/**@file   event_treeprofile.c
 * @brief  event handler to manage the tree profile and use it for predictions
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_treeprofile.h"

#include <assert.h>
#include <stddef.h>

#include "scip/def.h"
#include "scip/pub_event.h"
#include "scip/pub_message.h"
#include "scip/pub_tree.h"
#include "scip/scip_event.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_param.h"
#include "scip/scip_solvingstats.h"
#include "scip/type_event.h"
#include "scip/type_tree.h"

#define EVENTHDLR_NAME         "treeprofile"
#define EVENTHDLR_DESC         "event handler for treeprofile event"

#define EVENTTYPE_TREEPROFILE SCIP_EVENTTYPE_NODESOLVED
#define TREEPRROFILE_MINSIZE    512 /**< minimum size (depth) that tree profile can hold */
#define DEFAULT_ENABLED       FALSE /**< should the event handler collect data? */
#define DEFAULT_MAXDEPTHFACTOR 20.0 /**< factor by which the number of nodes exceeds the maximum width before producing estimations */
#define DEFAULT_FREQ           -1   /**< frequency for periodic output of estimation, or -1 for no output */
/*
 * Data structures
 */

/** statistics collected from profile used for prediction */
struct TreeProfileStats
{
   int                   maxdepth;           /**< maximum node depth encountered */
   int                   lastfulldepth;      /**< deepest layer for which all nodes have been explored */
   int                   minwaistdepth;      /**< minimum depth of the waist, ie the widest part of the tree */
   int                   maxwaistdepth;      /**< maximum depth of the waist, ie the widest part of the tree */
};

typedef struct TreeProfileStats TREEPROFILESTATS;


/** profile data structure for tree */
struct TreeProfile
{
   SCIP_Longint*         profile;            /**< array to store the tree profile */
   int                   profilesize;        /**< size of the profile array */
   TREEPROFILESTATS      stats;              /**< statistics collected from profile used for prediction */
};

typedef struct TreeProfile TREEPROFILE;

/** event handler data */
struct SCIP_EventhdlrData
{
   TREEPROFILE*          treeprofile;       /**< tree profile data structure */
   SCIP_Bool             enabled;           /**< should the event handler collect data? */
   SCIP_Real             maxdepthfactor;    /**< factor by which the number of nodes exceeds the maximum width before producing estimations */
   SCIP_Real             lastestimate;      /**< the last estimate predicted by SCIPpredictTotalSizeTreeprofile() */
   TREEPROFILESTATS      lastestimatestats; /**< tree profile statistics at last estimation */
   int                   freq;              /**< frequency for periodic output of estimation, or -1 for no output */
   SCIP_Longint          nextoutputnode;    /**< next output node for geometric output frequency */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** compare two tree profile statistics for equality */
static
SCIP_Bool isEqualTreeprofilestats(
   TREEPROFILESTATS*     stats,              /**< first tree profile statistics */
   TREEPROFILESTATS*     other               /**< other tree profile statistics */
   )
{
   assert(stats != NULL);
   assert(other != NULL);

   return  stats->maxdepth == other->maxdepth &&
      stats->lastfulldepth == other->lastfulldepth &&
      stats->minwaistdepth == other->minwaistdepth &&
      stats->maxwaistdepth == other->maxwaistdepth;
}

/** copy source tree profile into destination */
static
void copyTreeprofilestats(
   TREEPROFILESTATS*     dest,               /**< destination tree profile statistics */
   TREEPROFILESTATS*     src                 /**< source tree profile statistics */
   )
{
   assert(dest != NULL);
   assert(src != NULL);

   dest->maxdepth = src->maxdepth;
   dest->lastfulldepth = src->lastfulldepth;
   dest->minwaistdepth = src->minwaistdepth;
   dest->maxwaistdepth = src->maxwaistdepth;
}

/** reset tree profile statistics */
static
void resetTreeprofilestats(
   TREEPROFILESTATS*     treeprofilestats   /**< tree profile statistics */
   )
{
   assert(treeprofilestats != NULL);

   BMSclearMemory(treeprofilestats);
}


/** extend tree profile to deeper tree */
static
SCIP_RETCODE extendMemoryTreeprofile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE*          treeprofile,        /**< tree profile data structure */
   int                   mindepth            /**< minimum depth that the tree profile should hold */
   )
{
   if( mindepth < treeprofile->profilesize )
      return SCIP_OKAY;

   if( treeprofile->profile == NULL )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &treeprofile->profile, mindepth) );
      treeprofile->profilesize = mindepth;
   }
   else
   {
      int newsize = SCIPcalcMemGrowSize(scip, mindepth);

      assert(newsize > treeprofile->profilesize);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &treeprofile->profile, newsize) );
      BMSclearMemoryArray(&treeprofile->profile[treeprofile->profilesize], newsize - treeprofile->profilesize);
      treeprofile->profilesize = newsize;
   }

   return SCIP_OKAY;
}

/** create a tree profile */
static
SCIP_RETCODE createTreeprofile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE**         treeprofile         /**< pointer to store tree profile data structure */
   )
{
   assert(scip != NULL);
   assert(treeprofile != NULL);

   SCIP_CALL( SCIPallocMemory(scip, treeprofile) );

   (*treeprofile)->profile = NULL;
   (*treeprofile)->profilesize = 0;
   SCIP_CALL( extendMemoryTreeprofile(scip, *treeprofile, TREEPRROFILE_MINSIZE) );

   resetTreeprofilestats(&(*treeprofile)->stats);

   return SCIP_OKAY;
}

/** free a tree profile */
static
void freeTreeprofile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE**         treeprofile         /**< pointer to tree profile data structure */
   )
{
   assert(scip != NULL);
   assert(treeprofile != NULL);

   if( *treeprofile == NULL )
      return;

   SCIPfreeMemoryArray(scip, &(*treeprofile)->profile);

   SCIPfreeMemory(scip, treeprofile);

   *treeprofile = NULL;
}

/** update tree profile */
static
SCIP_RETCODE updateTreeprofile(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE*          treeprofile,        /**< tree profile data structure */
   SCIP_NODE*            node                /**< node that should be added to the profile */
   )
{
   int nodedepth;
   SCIP_Longint nodedepthcnt;
   SCIP_Longint maxnodes;

   assert(scip != NULL);
   assert(treeprofile != NULL);
   assert(node != NULL);

   nodedepth = SCIPnodeGetDepth(node);
   assert(nodedepth >= 0);
   maxnodes = treeprofile->profile[treeprofile->stats.minwaistdepth];
   assert(treeprofile->stats.minwaistdepth == treeprofile->stats.maxwaistdepth ||
      maxnodes == treeprofile->profile[treeprofile->stats.maxwaistdepth]);

   /* ensure that the memory can hold at least this depth */
   SCIP_CALL( extendMemoryTreeprofile(scip, treeprofile, nodedepth) );

   nodedepthcnt = ++treeprofile->profile[nodedepth];

   /* is this level full explored? We assume binary branching */
   if( (unsigned int)nodedepth < 8*sizeof(int) && nodedepthcnt == (1U << nodedepth) )
   {
      SCIPdebugMsg(scip, "Level %d fully explored: %lld nodes\n", nodedepth, nodedepthcnt);

      treeprofile->stats.lastfulldepth = nodedepth;
   }

   /* update maximum depth */
   if( treeprofile->stats.maxdepth < nodedepth )
   {
      assert(treeprofile->stats.maxdepth == nodedepth - 1);
      treeprofile->stats.maxdepth = nodedepth;
      SCIPdebugMsg(scip, "Maximum depth increased to %d\n", treeprofile->stats.maxdepth);
   }

   /* minimum and maximum waist now coincide */
   if( nodedepthcnt > maxnodes )
   {
      treeprofile->stats.minwaistdepth = treeprofile->stats.maxwaistdepth = nodedepth;
      SCIPdebugMsg(scip, "Updating depth of tree waist: %d (%lld nodes)\n", treeprofile->stats.minwaistdepth, nodedepthcnt);
   }
   else if( nodedepthcnt == maxnodes )
   {
      /* enlarge the interval in which the waist lies */
      if( treeprofile->stats.minwaistdepth > nodedepth )
         treeprofile->stats.minwaistdepth = nodedepth;
      else if( treeprofile->stats.maxwaistdepth < nodedepth )
         treeprofile->stats.maxwaistdepth = nodedepth;
   }
   assert(treeprofile->stats.minwaistdepth <= treeprofile->stats.maxwaistdepth);



   return SCIP_OKAY;
}

/** make a prediction of the total tree size based on the current tree profile */
static
SCIP_Real predictTotalSize(
   SCIP*                 scip,               /**< SCIP data structure */
   TREEPROFILE*          treeprofile         /**< tree profile data structure */
   )
{
   SCIP_Real estimate;
   SCIP_Real gamma_prod;
   int d;
   int waist;

   assert(scip != NULL);
   assert(treeprofile != NULL);

   waist = (2 * treeprofile->stats.maxwaistdepth + treeprofile->stats.minwaistdepth) / 3;

   gamma_prod = 2;
   estimate = 1;

   /* loop over all full levels */
   for( d = 1; d < treeprofile->stats.lastfulldepth; ++d )
   {
      SCIP_Real gamma_d = 2.0;

      estimate += gamma_prod;
      gamma_prod *= gamma_d;
   }

   /* loop until the waist is reached */
   for( ; d < waist; ++d )
   {
      SCIP_Real gamma_d = 2.0 - (d - treeprofile->stats.lastfulldepth + 1.0)/(waist - treeprofile->stats.lastfulldepth + 1.0);

      assert(1.0 <= gamma_d && gamma_d <= 2.0);
      estimate += gamma_prod;
      gamma_prod *= gamma_d;
   }

   /* loop over the remaining levels */
   for( ; d <= treeprofile->stats.maxdepth; ++d )
   {
      SCIP_Real gamma_d = (1.0 - (d - waist + 1.0)/(treeprofile->stats.maxdepth - waist + 1.0));
      assert(0.0 <= gamma_d && gamma_d <= 1.0);

      estimate += gamma_prod;
      gamma_prod *= gamma_d;
   }

   return estimate;
}



/*
 * Callback methods of event handler
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeTreeprofile)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitTreeprofile)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->enabled )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Activating tree profile data collection\n");

      SCIP_CALL( SCIPcatchEvent(scip, EVENTTYPE_TREEPROFILE, eventhdlr, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitTreeprofile)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(eventhdlrdata != NULL);


   if( eventhdlrdata->enabled )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Disabling tree profile data collection\n");

      SCIP_CALL( SCIPdropEvent(scip, EVENTTYPE_TREEPROFILE, eventhdlr, NULL, -1) );
   }

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolTreeprofile)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->treeprofile == NULL);

   eventhdlrdata->lastestimate = -1.0;
   eventhdlrdata->nextoutputnode = 1L;
   resetTreeprofilestats(&eventhdlrdata->lastestimatestats);

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( createTreeprofile(scip, &eventhdlrdata->treeprofile) );
   }


   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolTreeprofile)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   if( eventhdlrdata->treeprofile != NULL )
      freeTreeprofile(scip, &eventhdlrdata->treeprofile);

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecTreeprofile)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   SCIP_NODE* node;
   assert(eventhdlrdata != NULL);

   node = SCIPeventGetNode(event);
   assert(node != NULL);

   SCIP_CALL( updateTreeprofile(scip, eventhdlrdata->treeprofile, node) );

   if( (eventhdlrdata->freq > 0 && SCIPgetNNodes(scip) % eventhdlrdata->freq == 0)
            || (eventhdlrdata->freq == 0 && SCIPgetNNodes(scip) == eventhdlrdata->nextoutputnode) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "SCIP nodes: %lld Estimation: %g\n", SCIPgetNNodes(scip), SCIPpredictTotalSizeTreeprofile(scip));

      if( eventhdlrdata->freq == 0 )
         eventhdlrdata->nextoutputnode *= 2;
   }

   return SCIP_OKAY;
}

/** returns a prediction of the total size of the final B&B tree, or -1 if no suitable prediction can be computed (yet) */
SCIP_Real SCIPpredictTotalSizeTreeprofile(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);

   /* if the event handler has not been included with scip, return here */
   if( eventhdlr == NULL )
      return -1.0;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   /* if the event handler is deactivated or it is too early during the solving process, return */
   if( eventhdlrdata->treeprofile == NULL )
      return -1;

   /* two few nodes to make a prediction */
   if( eventhdlrdata->maxdepthfactor * eventhdlrdata->treeprofile->stats.maxdepth > SCIPgetNNodes(scip) )
      return -1;

   /* reuse previous estimation if tree profile hasn't changed */
   if( isEqualTreeprofilestats(&eventhdlrdata->lastestimatestats, &eventhdlrdata->treeprofile->stats) )
   {
      SCIPdebugMsg(scip, "Reusing previous estimation result %g\n", eventhdlrdata->lastestimate);

      return eventhdlrdata->lastestimate;
   }

   /* copy tree profile statistics */
   copyTreeprofilestats(&eventhdlrdata->lastestimatestats, &eventhdlrdata->treeprofile->stats);

   eventhdlrdata->lastestimate = predictTotalSize(scip, eventhdlrdata->treeprofile);

   return eventhdlrdata->lastestimate;
}

/** creates event handler for treeprofile event */
SCIP_RETCODE SCIPincludeEventHdlrTreeprofile(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create treeprofile event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   eventhdlrdata->treeprofile = NULL;

   eventhdlr = NULL;

   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecTreeprofile, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeTreeprofile) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitTreeprofile) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitTreeprofile) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolTreeprofile) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolTreeprofile) );

   /* add treeprofile event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "estimates/profile/enabled",
         "should the event handler collect data?", &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "estimates/profile/maxdepthfactor",
         "factor by which the number of nodes exceeds the maximum width before producing estimations",
         &eventhdlrdata->maxdepthfactor, FALSE, DEFAULT_MAXDEPTHFACTOR, 1.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "estimates/profile/freq",
         "positive frequency for periodic output of estimation, or -1 for no output",
         &eventhdlrdata->freq, FALSE, DEFAULT_FREQ, -1, INT_MAX, NULL, NULL) );


   return SCIP_OKAY;
}
