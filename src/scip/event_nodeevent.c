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

/**@file   event_nodeevent.c
 * @brief  eventhdlr for nodeevent event
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/event_nodeevent.h"
#include "string.h"

#define EVENTHDLR_NAME         "nodeevent"
#define EVENTHDLR_DESC         "event handler for nodeevent event"
#define EVENTTOCATCH SCIP_EVENTTYPE_NODEFOCUSED
#define DEFAULT_ENABLED FALSE
#define DEFAULT_ARRAYSIZE 150
#define DEFAULT_DISPFREQ 10000
#define DEFAULT_RANDSEED 314156
/*
 * Data structures
 */

/** event handler data */
struct SCIP_EventhdlrData
{
   int*                  nnodesdepth;        /**< the number of nodes at a particular depth */
   int                   nnodesdepthsize;    /**< size of the array */
   int                   filterpos;          /**< filter position of this event for quicker dropping */
   SCIP_Bool             enabled;            /**< user parameter to enable the event handler */
   int                   dispfreq;
   int                   maxnodes;           /**< maximum number of nodes at single depth */
   int                   waist;              /**< depth of maxnodes (not necessarily unique)*/
   unsigned int          randseed;           /**< random seed */

};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */
static
SCIP_RETCODE ensureArraySize(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata,
   int                   newsize
   )
{
   assert(scip != NULL);
   assert(eventhdlrdata != NULL);
   assert(eventhdlrdata->nnodesdepthsize == 0 || eventhdlrdata->nnodesdepth != NULL);
   if( newsize > eventhdlrdata->nnodesdepthsize )
   {
      int neededsize;
      int i;

      neededsize = MAX(newsize, DEFAULT_ARRAYSIZE);
      if( eventhdlrdata->nnodesdepthsize == 0 )
      {
         SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->nnodesdepth, neededsize) );
      }
      else
      {
         SCIP_CALL( SCIPreallocMemoryArray(scip, &eventhdlrdata->nnodesdepth, neededsize) );
      }

      for( i = eventhdlrdata->nnodesdepthsize; i < neededsize; ++i )
         eventhdlrdata->nnodesdepth[i] = 0;

      eventhdlrdata->nnodesdepthsize = neededsize;
   }

   return SCIP_OKAY;
}

static
void freeArray(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata
   )
{
   assert(eventhdlrdata->nnodesdepthsize == 0 || eventhdlrdata->nnodesdepth != NULL);
   if( eventhdlrdata->nnodesdepthsize > 0 )
      SCIPfreeMemoryArray(scip, &eventhdlrdata->nnodesdepth);

   eventhdlrdata->nnodesdepth = NULL;
   eventhdlrdata->nnodesdepthsize = 0;
}

static
void copyDepthData(
   int*                  nodesdepth,         /**< array to copy nodes data to */
   int                   nodesdepthsize,     /**< size of the array */
   SCIP_NODE**           nodes,              /**< nodes to count depths */
   int                   nodessize,          /**< size of nodes array */
   int*                  maxnodes            /**< maximum number of nodes in single depth */
   )
{
   int n;

   assert(maxnodes != NULL);
   assert(*maxnodes >= 0);

   for( n = 0; n < nodessize; ++n )
   {
      int depth;

      depth = SCIPnodeGetDepth(nodes[n]);
      assert(depth < nodesdepthsize);
      assert(0 <= depth);
      ++nodesdepth[depth];

      if( nodesdepth[depth] > *maxnodes )
         *maxnodes = nodesdepth[depth];
   }
}

static
SCIP_Longint calcEstimation(
   int                   maxdepth,
   int                   lastfulldepth,
   int                   waist
)
{
   int i;
   SCIP_Real lastlevelnodes;
   SCIP_Real estimatednodes;

   lastlevelnodes = 1.0;
   estimatednodes = 1.0;

   for( i = 1; i <= maxdepth; ++i )
    {
       SCIP_Real estimatedgamma;

       if( i  < lastfulldepth )
          estimatedgamma = 2;
       else if( i < waist )
          estimatedgamma = 2 - (i - lastfulldepth + 1)/(SCIP_Real)(waist - lastfulldepth + 1);
       else
          estimatedgamma = 1 - (i - waist + 1)/(SCIP_Real)(maxdepth - waist + 1);

       lastlevelnodes = estimatedgamma * lastlevelnodes;
       estimatednodes += lastlevelnodes;
    }

   return (SCIP_Longint)estimatednodes;
}

static
void determineTreeCharacteristics(
   SCIP*              scip,            /**< SCIP data structure */
   int*               nnodesdepth,
   int                nnodesdepthsize,
   int*               lastfulldepth,
   int*               waist,
   int*               maxdepth,
   int                maxnodes
   )
{
   int i;
   int minwaistlevel;
   int maxwaistlevel;

   *lastfulldepth = -1;
   minwaistlevel = -1;
   maxwaistlevel = 0;

   /* loop over all depths with more than zero nodes. Determine waist, last full level, and maximum depth */
   for( i = 0; i < nnodesdepthsize && nnodesdepth[i] > 0; ++i )
   {
      if( minwaistlevel == -1 && nnodesdepth[i] > 0.5 * maxnodes )
      {
         minwaistlevel = i;

      }
      if( minwaistlevel != -1 && nnodesdepth[i] > 0.5 * maxnodes )
         maxwaistlevel = i;

      if( i > 0 && *lastfulldepth == -1 && (nnodesdepth[i] / nnodesdepth[i - 1]) < 2 )
      {
         *lastfulldepth = i - 1;
      }
   }

   assert(i == 0 || nnodesdepth[i - 1] > 0);
   *maxdepth = i - 1;
   assert(minwaistlevel <=maxwaistlevel);
   *waist = SCIPfeasCeil(scip, (maxwaistlevel + minwaistlevel) / 2.0);
}

static
SCIP_Longint calcTreeSizeActiveNodes(
   SCIP_NODE**           opennodes,
   SCIP_EVENTHDLRDATA*   eventhdlrdata,
   int                   nopennodes,
   int                   maxdepth,
   int                   lastfulldepth,
   int                   waist
   )
{
   int n;
   SCIP_Longint estimatednodes;

   estimatednodes = 0L;
   for( n = 0;  n < nopennodes; ++n )
   {
      SCIP_NODE* node;
      int depth;
      SCIP_Real factor;

      node = opennodes[n];
      depth = SCIPnodeGetDepth(node);
      factor = (maxdepth - depth) / (SCIP_Real)(maxdepth);
      factor = SCIPgetRandomReal(0, factor, &eventhdlrdata->randseed);
      estimatednodes += calcEstimation((int)((maxdepth) * factor), (int)(lastfulldepth * factor), (int)(waist * factor));
   }

   return estimatednodes;
}

static
SCIP_Longint estimateTreeSizeActiveNodes(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata,
   int*                  lastfulldepth,
   int*                  waist,
   int*                  maxdepth
   )
{

   int* nodesdepth;
   SCIP_Longint estimatednodes;

   if( eventhdlrdata->nnodesdepth == NULL )
      return 0;

   assert(lastfulldepth != NULL);
   assert(waist != NULL);
   assert(maxdepth != NULL);

   nodesdepth = eventhdlrdata->nnodesdepth;
   estimatednodes = SCIPgetNNodes(scip);

   /* loop over all depths with more than zero nodes. Determine waist, last full level, and maximum depth */
   determineTreeCharacteristics(scip, nodesdepth, eventhdlrdata->nnodesdepthsize, lastfulldepth, waist, maxdepth, eventhdlrdata->maxnodes);

   /* during solving stage, the tree estimation gets refined by open nodes data */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_NODE** leaves;
      SCIP_NODE** children;
      SCIP_NODE** siblings;
      int nsiblings;
      int nchildren;
      int nleaves;

      leaves = NULL;
      children = NULL;
      siblings = NULL;
      nchildren = nleaves = nsiblings = 0;

      assert(eventhdlrdata->nnodesdepthsize > 0);

      SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );

      if( nsiblings > 0 )
      {
         estimatednodes += calcTreeSizeActiveNodes(siblings, eventhdlrdata, nsiblings, *maxdepth, *lastfulldepth, *waist);
      }
      if( nchildren > 0 )
         estimatednodes += calcTreeSizeActiveNodes(children, eventhdlrdata, nchildren, *maxdepth, *lastfulldepth, *waist);
      if( nleaves > 0 )
         estimatednodes += calcTreeSizeActiveNodes(leaves, eventhdlrdata, nleaves, *maxdepth, *lastfulldepth, *waist);
   }
   return estimatednodes;

}

static
SCIP_Longint estimateTreeSize(
   SCIP*                 scip,
   SCIP_EVENTHDLRDATA*   eventhdlrdata,
   int*                  lastfulldepth,
   int*                  waist,
   int*                  maxdepth
   )
{

   int* nodesdepth;
   SCIP_Longint estimatednodes;
   int maxnodes;

   if( eventhdlrdata->nnodesdepth == NULL )
      return 0;

   assert(lastfulldepth != NULL);
   assert(waist != NULL);
   assert(maxdepth != NULL);

   maxnodes = eventhdlrdata->maxnodes;
   /* during solving stage, the tree estimation gets refined by open nodes data */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_NODE** leaves;
      SCIP_NODE** children;
      SCIP_NODE** siblings;
      int nsiblings;
      int nchildren;
      int nleaves;

      leaves = NULL;
      children = NULL;
      siblings = NULL;
      nchildren = nleaves = nsiblings = 0;

      assert(eventhdlrdata->nnodesdepthsize > 0);

      SCIP_CALL( SCIPallocBufferArray(scip, &nodesdepth, eventhdlrdata->nnodesdepthsize) );

      BMScopyMemoryArray(nodesdepth, eventhdlrdata->nnodesdepth, eventhdlrdata->nnodesdepthsize);

      SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );

      if( nsiblings > 0 )
         copyDepthData(nodesdepth, eventhdlrdata->nnodesdepthsize, siblings, nsiblings, &maxnodes);
      if( nchildren > 0 )
         copyDepthData(nodesdepth, eventhdlrdata->nnodesdepthsize, children, nchildren, &maxnodes);
      if( nleaves > 0 )
         copyDepthData(nodesdepth, eventhdlrdata->nnodesdepthsize, leaves, nleaves, &maxnodes);
   }
   else
      /* if scip is not in solving stage, the entire tree is in the eventhdlrdata */
      nodesdepth = eventhdlrdata->nnodesdepth;

   /* loop over all depths with more than zero nodes. Determine waist, last full level, and maximum depth */
    determineTreeCharacteristics(scip, nodesdepth, eventhdlrdata->nnodesdepthsize, lastfulldepth, waist, maxdepth, eventhdlrdata->maxnodes);

   estimatednodes = calcEstimation(*maxdepth, *lastfulldepth, *waist);

   /* nodesdepth is a buffer array only in solving stage */
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      SCIPfreeBufferArray(scip, &nodesdepth);

   return estimatednodes;

}
/*
 * Callback methods of event handler
 */

/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyNodeevent)
{
   SCIP_CALL( SCIPincludeEventHdlrNodeevent(scip) );

   return SCIP_OKAY;
}

/** destructor of event handler to free user data (called when SCIP is exiting) */
static
SCIP_DECL_EVENTFREE(eventFreeNodeevent)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   assert(eventhdlrdata->nnodesdepth == NULL);
   SCIPfreeMemory(scip, &eventhdlrdata);
   SCIPeventhdlrSetData(eventhdlr, NULL);

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before branch and bound terminates) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolNodeevent)
{

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   assert(!eventhdlrdata->enabled || eventhdlrdata->filterpos >= 0);

   if( eventhdlrdata->filterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENTTOCATCH, eventhdlr, NULL, eventhdlrdata->filterpos) );
      eventhdlrdata->filterpos = -1;
   }

   if( eventhdlrdata->enabled )
   {
      int lastfulldepth;
      int waist;
      int maxdepth;
      SCIP_Longint estimatednodes;

      lastfulldepth = 0;
      waist = 0;
      maxdepth = 0;

      estimatednodes = estimateTreeSizeActiveNodes(scip, eventhdlrdata, &lastfulldepth, &waist, &maxdepth);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Explored nodes: %"SCIP_LONGINT_FORMAT", estimated number of nodes: %"SCIP_LONGINT_FORMAT" (%d,%d,%d)\n",
            SCIPgetNNodes(scip), estimatednodes, lastfulldepth, waist, maxdepth);

   }
   freeArray(scip, eventhdlrdata);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolNodeevent)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   if( eventhdlrdata->enabled )
   {
      assert(eventhdlrdata->filterpos == -1);
      SCIP_CALL( SCIPcatchEvent(scip, EVENTTOCATCH, eventhdlr, NULL, &eventhdlrdata->filterpos) );
   }
   eventhdlrdata->waist = 0;
   eventhdlrdata->maxnodes = 1;
   eventhdlrdata->randseed = DEFAULT_RANDSEED;

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecNodeevent)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   int focusnodedepth;

   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   focusnodedepth = SCIPgetFocusDepth(scip);
   SCIP_CALL( ensureArraySize(scip, eventhdlrdata, SCIPgetMaxDepth(scip) + 1) );
   assert(eventhdlrdata->nnodesdepth[focusnodedepth] >= 0);
   ++eventhdlrdata->nnodesdepth[focusnodedepth];
   if( eventhdlrdata->nnodesdepth[focusnodedepth] > eventhdlrdata->maxnodes )
   {
      eventhdlrdata->waist = focusnodedepth;
      eventhdlrdata->maxnodes = eventhdlrdata->nnodesdepth[focusnodedepth];
   }

   if( SCIPgetNNodes(scip) >= eventhdlrdata->dispfreq && SCIPgetNNodes(scip) % eventhdlrdata->dispfreq == 0 )
   {
      int lastfulldepth;
      int waist;
      int maxdepth;
      SCIP_Longint estimatednodes;

      lastfulldepth = 0;
      waist = 0;
      maxdepth = 0;

      estimatednodes = estimateTreeSizeActiveNodes(scip, eventhdlrdata, &lastfulldepth, &waist, &maxdepth);

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Explored nodes: %"SCIP_LONGINT_FORMAT", estimated number of nodes: %"SCIP_LONGINT_FORMAT" (%d,%d,%d)\n",
            SCIPgetNNodes(scip), estimatednodes, lastfulldepth, waist, maxdepth);
   }

   return SCIP_OKAY;
}

/** creates event handler for nodeevent event */
SCIP_RETCODE SCIPincludeEventHdlrNodeevent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create nodeevent event handler data */
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );

   eventhdlrdata->nnodesdepth = NULL;
   eventhdlrdata->nnodesdepthsize = 0;
   eventhdlrdata->filterpos = -1;
   eventhdlr = NULL;

   /* include event handler into SCIP */
   /* use SCIPincludeEventhdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecNodeevent, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyNodeevent) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeNodeevent) );
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolNodeevent) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolNodeevent) );

   /* add nodeevent event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/"EVENTHDLR_NAME"/enabled", "bla", &eventhdlrdata->enabled, FALSE,
               DEFAULT_ENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "eventhdlr/"EVENTHDLR_NAME"/dispfreq", "bla", &eventhdlrdata->dispfreq, FALSE, DEFAULT_DISPFREQ,
               1, INT_MAX / 4, NULL, NULL) );
   return SCIP_OKAY;
}
