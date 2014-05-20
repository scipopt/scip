/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_treeinfos.c
 * @brief  eventhdlr for treeinfos event
 * @author Gregor Hendel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/event_treeinfos.h"
#include "scip/scip.h"
#include <string.h>

#define EVENTHDLR_NAME         "treeinfos"
#define EVENTHDLR_DESC         "event handler for treeinfos event"
#define EVENT_TYPE_TREEINFOS  SCIP_EVENTTYPE_NODEFOCUSED | SCIP_EVENTTYPE_NODEBRANCHED | SCIP_EVENTTYPE_BESTSOLFOUND
#define DEFAULT_ENABLED        FALSE
#define HEADER "NUMBER  PRIMAL LOWER ESTIM ROOTLPPSESTIM DEPTH NCANDS NODELOWER LPSTAT"
#define DEFAULT_FILEOUTPUT FALSE
/*
 * Data structures
 */

/* TODO: fill in the necessary event handler data */

/** depth information structure */
struct DepthInfo
{
   int          nsolvednodes;  /**< nodes that were solved so far at this depth */
   SCIP_Real    minestimate;   /**< the minimum estimate of a solved node */
   SCIP_NODE**  minnodes;      /**< points to the rank1nodes at this depth (open nodes whose estimate is lower than current
                                    minimum estimate over solved nodes) */
   int          nminnodes;     /**< the number of minimum nodes */
   int          minnodescapacity; /**< the capacity of the min nodes array */
};

typedef struct DepthInfo DEPTHINFO;

/** event handler data */
struct SCIP_EventhdlrData
{
   SCIP_Bool            enabled;            /**< is the treeinfos event handler enabled? */
   int                  eventfilterpos;     /**< the event filter position, or -1, if event has not (yet) been caught */
   DEPTHINFO**          depthinfos;
   int                  maxdepth;
   int                  nrank1nodes;        /**< number of rank1 nodes */
};

/*
 * Local methods
 */

/** nodes are sorted first by their estimates, and if estimates are equal, by their number */
static
SCIP_DECL_SORTPTRCOMP(sortCompTreeinfo)
{
   SCIP_NODE* node1;
   SCIP_NODE* node2;
   SCIP_Real estim1, estim2;
   SCIP_Longint number1, number2;
   node1 = (SCIP_NODE*)elem1;
   node2 = (SCIP_NODE*)elem2;

   estim1 = SCIPnodeGetEstimate(node1);
   estim2 = SCIPnodeGetEstimate(node2);
   number1 = SCIPnodeGetNumber(node1);
   number2 = SCIPnodeGetNumber(node2);
   if( estim1 < estim2 )
      return -1;
   else if( estim1 > estim2 )
      return 1;
   else if( number1 < number2 )
      return -1;
   else if( number1 > number2 )
      return 1;

   return 0;
}

static
SCIP_RETCODE nodesUpdateRank1Nodes(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata,
   SCIP_NODE** nodes,
   int nnodes
   )
{
   int n;

   assert(nnodes == 0 || nodes != NULL);

   for( n = 0; n < nnodes; ++n )
   {
      SCIP_NODE* node = nodes[n];
      DEPTHINFO* depthinfo = eventhdlrdata->depthinfos[SCIPnodeGetDepth(node)];

      assert(SCIPnodeGetType(node) == SCIP_NODETYPE_CHILD || SCIPnodeGetType(node) == SCIP_NODETYPE_LEAF
            || SCIPnodeGetType(node) == SCIP_NODETYPE_SIBLING);
      /* an open node has rank 1 if it has an estimate at least as small as the best solved node
       * at this depth
       */
      if( depthinfo->nsolvednodes == 0 || SCIPisGE(scip, depthinfo->minestimate, SCIPnodeGetEstimate(node)) )
      {
         int pos;

         /* allocate additional memory to hold new node */
         if( depthinfo->nminnodes == depthinfo->minnodescapacity )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &depthinfo->minnodes, depthinfo->minnodescapacity, 2 * depthinfo->minnodescapacity) );
            depthinfo->minnodescapacity *= 2;
         }

         /* find correct insert position */
         SCIPsortedvecInsertPtr((void **)depthinfo->minnodes, sortCompTreeinfo, (void*)node, &depthinfo->nminnodes, &pos);
         assert(pos >= 0 && pos < depthinfo->nminnodes);
         assert(depthinfo->minnodes[pos] == node);
         ++eventhdlrdata->nrank1nodes;
      }
   }

   return SCIP_OKAY;
}

static
void removeNode(
   SCIP_NODE* node,
   SCIP_EVENTHDLRDATA* eventhdlrdata
   )
{
   DEPTHINFO* depthinfo;
   int pos;
   SCIP_Bool contained;

   depthinfo = eventhdlrdata->depthinfos[SCIPnodeGetDepth(node)];

   if( depthinfo->nminnodes == 0 )
      return;

   contained = SCIPsortedvecFindPtr((void **)depthinfo->minnodes, sortCompTreeinfo, (void *)node, depthinfo->nminnodes, &pos);

   if( contained )
   {
      SCIPsortedvecDelPosPtr((void **)depthinfo->minnodes, sortCompTreeinfo, pos, &(depthinfo->nminnodes));
      --eventhdlrdata->nrank1nodes;
   }
}

/** returns the current number of rank 1 nodes in the tree */
int SCIPgetNRank1Nodes(
   SCIP* scip
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   eventhdlrdata = SCIPeventhdlrGetData(SCIPfindEventhdlr(scip, EVENTHDLR_NAME));
   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
      return eventhdlrdata->nrank1nodes;
   else
      return -1;
}

/** discards all previous depth information and renews it afterwards */
static
SCIP_RETCODE storeRank1Nodes(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_NODE** leaves;
   SCIP_NODE** children;
   SCIP_NODE** siblings;

   int nleaves;
   int nchildren;
   int nsiblings;
   int d;
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return 0;

   for( d = 0; d < eventhdlrdata->maxdepth; ++d )
      eventhdlrdata->depthinfos[d]->nminnodes = 0;

   eventhdlrdata->nrank1nodes = 0;

   assert(eventhdlrdata != NULL);

   nleaves = nchildren = nsiblings = 0;

   SCIP_CALL( SCIPgetOpenNodesData(scip, &leaves, &children, &siblings, &nleaves, &nchildren, &nsiblings) );

   SCIP_CALL ( nodesUpdateRank1Nodes(scip, eventhdlrdata, children, nchildren) );

   SCIP_CALL ( nodesUpdateRank1Nodes(scip, eventhdlrdata, siblings, nsiblings) );

   SCIP_CALL ( nodesUpdateRank1Nodes(scip, eventhdlrdata, leaves, nleaves) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE createDepthinfo(
   SCIP* scip,
   DEPTHINFO** depthinfo
   )
{
   SCIP_CALL( SCIPallocMemory(scip, depthinfo) );

   (*depthinfo)->minestimate = SCIPinfinity(scip);
   (*depthinfo)->nsolvednodes = 0;
   (*depthinfo)->nminnodes = 0;
   (*depthinfo)->minnodescapacity = 2;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*depthinfo)->minnodes, (*depthinfo)->minnodescapacity) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE freeDepthinfo(
   SCIP* scip,
   DEPTHINFO** depthinfo
   )
{
   SCIPfreeBlockMemoryArray(scip, &(*depthinfo)->minnodes, (*depthinfo)->minnodescapacity);
   SCIPfreeMemory(scip, depthinfo);

   return SCIP_OKAY;
}

static
void updateDepthinfo(
   SCIP* scip, /**< SCIP data structure */
   SCIP_EVENTHDLRDATA* eventhdlrdata,
   SCIP_NODE* node
   )
{
   DEPTHINFO* depthinfo;

   depthinfo = eventhdlrdata->depthinfos[SCIPnodeGetDepth(node)];

   removeNode(node, eventhdlrdata);

   if( SCIPisLT(scip, SCIPnodeGetEstimate(node), depthinfo->minestimate) )
      depthinfo->minestimate = SCIPnodeGetEstimate(node);

   /* loop over remaining, unsolved nodes and decide whether they are still rank1 nodes */
   while(  depthinfo->nminnodes > 0 && SCIPisGT(scip, SCIPnodeGetEstimate(depthinfo->minnodes[depthinfo->nminnodes - 1]), depthinfo->minestimate) )
   {
      /* forget about node */
      --(depthinfo->nminnodes);
      --(eventhdlrdata->nrank1nodes);
   }

   ++(depthinfo->nsolvednodes);
}
static
SCIP_RETCODE storeDepthInfo(
   SCIP* scip,
   SCIP_EVENTHDLRDATA* eventhdlrdata,
   SCIP_NODE* node
   )
{
   int nodedepth;
   int newsize;
   int oldsize;

   nodedepth = SCIPnodeGetDepth(node);
   oldsize = eventhdlrdata->maxdepth;
   newsize = oldsize;
   if( oldsize == 0 )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &eventhdlrdata->depthinfos, 10) );
      newsize = 10;
   }
   else if( nodedepth  + 1 >= eventhdlrdata->maxdepth )
   {
      assert(nodedepth > 0);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &eventhdlrdata->depthinfos, 2 * nodedepth) );
      newsize = 2 * nodedepth;
   }

   if( newsize > oldsize )
   {
      int c;
      for( c = oldsize; c < newsize; ++c )
      {
         SCIP_CALL( createDepthinfo(scip, &(eventhdlrdata->depthinfos[c])) );

      }

      eventhdlrdata->maxdepth = newsize;
   }

   assert(newsize > nodedepth);

   updateDepthinfo(scip, eventhdlrdata, node);

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

static
SCIP_DECL_EVENTCOPY(eventCopyTreeinfos)
{
   SCIP_CALL( SCIPincludeEventHdlrTreeinfos(scip) );

   return SCIP_OKAY;
}
/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINITSOL(eventInitsolTreeinfos)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos == -1);
   eventhdlrdata->depthinfos = NULL;
   eventhdlrdata->maxdepth = 0;

   if( eventhdlrdata->enabled )
   {
      SCIP_CALL( SCIPcatchEvent(scip, EVENT_TYPE_TREEINFOS, eventhdlr, NULL, &eventhdlrdata->eventfilterpos) );
      assert(eventhdlrdata->eventfilterpos >= 0);
   }
   return SCIP_OKAY;
}

static
SCIP_DECL_EVENTFREE(eventFreeTreeinfos)
{  /*lint --e{715}*/

   SCIP_EVENTHDLRDATA* eventhdlrdata;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata != NULL);

   SCIPfreeMemory(scip, &eventhdlrdata);

   return SCIP_OKAY;
}

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
static
SCIP_DECL_EVENTEXITSOL(eventExitsolTreeinfos)
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);

   assert(eventhdlrdata->eventfilterpos == -1 || eventhdlrdata->enabled );
   if( eventhdlrdata->eventfilterpos >= 0 )
   {
      SCIP_CALL( SCIPdropEvent(scip, EVENT_TYPE_TREEINFOS, eventhdlr, NULL, eventhdlrdata->eventfilterpos) );
      eventhdlrdata->eventfilterpos = -1;
   }

   if( eventhdlrdata->maxdepth > 0 )
   {
      int c;
      for( c = 0; c < eventhdlrdata->maxdepth; ++c )
      {
         SCIP_CALL( freeDepthinfo(scip, &(eventhdlrdata->depthinfos[c])) );
      }
      SCIPfreeMemoryArray(scip, &eventhdlrdata->depthinfos);
      eventhdlrdata->maxdepth = 0;
   }

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecTreeinfos)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTTYPE eventtype;
   assert(scip != NULL);
   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   assert(eventhdlrdata->eventfilterpos >= 0);

   eventtype = SCIPeventGetType(event);

   assert(eventtype & (EVENT_TYPE_TREEINFOS));

   if( eventtype & SCIP_EVENTTYPE_BESTSOLFOUND )
   {
      SCIP_CALL( storeRank1Nodes(scip, eventhdlrdata) );
   }
   else if( eventtype & SCIP_EVENTTYPE_NODEFOCUSED )
   {
      SCIP_NODE* focusnode;
      focusnode = SCIPgetCurrentNode(scip);
      storeDepthInfo(scip, eventhdlrdata, focusnode);
   }

   else if( eventtype & SCIP_EVENTTYPE_NODEBRANCHED )
   {
      SCIP_NODE** children;
      int nchildren;
      SCIP_CALL( SCIPgetChildren(scip, &children, &nchildren) );
      SCIP_CALL ( nodesUpdateRank1Nodes(scip, eventhdlrdata, children, nchildren) );
   }

   return SCIP_OKAY;
}

/** creates event handler for treeinfos event */
SCIP_RETCODE SCIPincludeEventHdlrTreeinfos(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;

   /* create treeinfos event handler data */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &eventhdlrdata) );
   assert(eventhdlrdata != NULL);
   eventhdlrdata->eventfilterpos = -1;
   eventhdlrdata->depthinfos = NULL;
   eventhdlrdata->maxdepth = 0;

   eventhdlr = NULL;
   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecTreeinfos, eventhdlrdata) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInitsol(scip, eventhdlr, eventInitsolTreeinfos) );
   SCIP_CALL( SCIPsetEventhdlrExitsol(scip, eventhdlr, eventExitsolTreeinfos) );
   SCIP_CALL( SCIPsetEventhdlrFree(scip, eventhdlr, eventFreeTreeinfos) );
   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyTreeinfos) );

   /* add treeinfos event handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "eventhdlr/treeinfos/enabled","enable event handler to perform treeinfos",
               &eventhdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   return SCIP_OKAY;
}
