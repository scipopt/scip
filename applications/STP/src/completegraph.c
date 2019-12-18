/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   completegraph.c
 * @brief  includes complete graph methods, in particular for MST calculation
 * @author Daniel Rehfeldt
 *
 * Complete graph methods, in particular for MST calculation. Assumes a maximum size of the graph.
 * Only sensible for small maximum sizes.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

//#define SCIP_DEBUG


#include "completegraph.h"
#include "portab.h"

#define NODE_ID_UNDEFINED -2
#define EDGECOST_UNDEFINED_VALUE -1.0


/** gets current number of nodes */
static inline
int getNnodesCurr(
   const CGRAPH*         cgraph              /**< the graph */
)
{
   assert(cgraph);
   return cgraph->nnodes_curr;
}


/** gets maximum number of nodes */
static inline
int getNnodesMax(
   const CGRAPH*         cgraph              /**< the graph */
)
{
   assert(cgraph);
   return cgraph->nnodes_max;
}


/** gets start position in edge cost array */
static inline
int getEdgeStart(
   int                   nodepos,            /**< node position */
   int                   nnodes_max          /**< maximum number of nodes */
)
{
   assert(nodepos >= 0);
   assert(nodepos < nnodes_max);

   return (nodepos * nnodes_max);
}


/** gets end position in edge cost array; NOT included! */
static inline
int getEdgeEnd(
   int                   nodepos,            /**< node position */
   int                   nnodes_curr,        /**< current number of nodes */
   int                   nnodes_max          /**< maximum number of nodes */
)
{
   assert(nodepos >= 0);
   assert(nnodes_curr >= 0);
   assert(nnodes_curr <= nnodes_max);
   assert(nodepos < nnodes_curr);

   return (nodepos * nnodes_max + nnodes_curr);
}


/** is the graph valid? */
SCIP_Bool cgraph_valid(
   const CGRAPH*         cgraph              /**< the graph */
)
{
   int start;
   const int nnodes_curr = getNnodesCurr(cgraph);
   const int nnodes_max = getNnodesMax(cgraph);
   const SCIP_Real* const edgecosts = cgraph->edgecosts;
   const int* const nodeids = cgraph->nodeids;

   assert(edgecosts && nodeids);
   assert(nnodes_max > 1);
   assert(nnodes_curr <= nnodes_max);
   assert(nnodes_curr >= 0);

   for( int i = 0; i < nnodes_curr - 1; i++ )
   {
      if( getEdgeEnd(i, nnodes_curr, nnodes_max) > getEdgeStart(i + 1, nnodes_max) )
      {
         SCIPdebugMessage("positions are broken \n");
         return FALSE;
      }
   }

   for( int i = 0; i < nnodes_curr; i++ )
   {
      if( getEdgeStart(i, nnodes_max) > getEdgeStart(i + 1, nnodes_max) )
      {
         SCIPdebugMessage("positions are broken \n");
         return FALSE;
      }
   }

   start = nnodes_curr > 0 ? getEdgeEnd(nnodes_curr - 1, nnodes_curr, nnodes_max) : 0;

   for( int i = start; i < nnodes_max * nnodes_max; i++ )
   {
      if( !EQ(EDGECOST_UNDEFINED_VALUE, edgecosts[i]) )
      {
         SCIPdebugMessage("unused edge value is wrong: %f \n", edgecosts[i]);
         return FALSE;
      }
   }

   for( int i = 0; i < nnodes_curr; i++ )
   {
      if( NODE_ID_UNDEFINED == nodeids[i] )
      {
         SCIPdebugMessage("node %d is wrongly set to NODE_ID_UNDEFINED \n", i);
         return FALSE;
      }
   }

   for( int i = nnodes_curr; i < nnodes_max; i++ )
   {
      if( NODE_ID_UNDEFINED != nodeids[i] )
      {
         SCIPdebugMessage("node %d is not set to NODE_ID_UNDEFINED \n", i);
         return FALSE;
      }
   }

   return TRUE;
}


/** initialize complete, undirected graph */
SCIP_RETCODE cgraph_init(
   SCIP*                 scip,               /**< SCIP data structure */
   CGRAPH**              cgraph,             /**< new graph */
   int                   maxnnodes           /**< maximum number of nodes */
   )
{
   CGRAPH* g;

   assert(scip && cgraph);
   assert(maxnnodes > 1);

   SCIP_CALL( SCIPallocMemory(scip, cgraph) );

   g = *cgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->edgecosts), maxnnodes * maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->nodeids), maxnnodes) );

   g->nnodes_curr = 0;
   g->nnodes_max = maxnnodes;

#ifndef NDEBUG
   for( int i = 0; i < maxnnodes * maxnnodes; i++ )
      g->edgecosts[i] = EDGECOST_UNDEFINED_VALUE;

   for( int i = 0; i < maxnnodes; i++ )
      g->nodeids[i] = NODE_ID_UNDEFINED;
#endif

   assert(cgraph_valid(g));

   SCIPdebugMessage("cgraph has been successfully built \n");

   return SCIP_OKAY;
}


/** free complete graph */
void cgraph_free(
   SCIP*                 scip,               /**< SCIP data structure */
   CGRAPH**              cgraph              /**< new graph */
   )
{
   CGRAPH* g;

   assert(scip && cgraph);

   g = *cgraph;

   SCIPfreeMemoryArray(scip, &(g->nodeids));
   SCIPfreeMemoryArray(scip, &(g->edgecosts));

   SCIPfreeMemory(scip, cgraph);
}



/** adds node */
void cgraph_node_append(
   CGRAPH*               cgraph,             /**< new graph */
   const SCIP_Real*      adjcosts,           /**< array with edge costs to neigbors of size nnodes_curr (use -1.0 as debug sentinel) */
   int                   nodeid              /**< the node id */
   )
{
   int lastedge;
   SCIP_Real* const edgecosts = cgraph->edgecosts;
   const int nnodes_curr_org = cgraph->nnodes_curr;
   const int nnodes_curr_new = nnodes_curr_org + 1;
   const int nnodes_max = cgraph->nnodes_max;
   const int start_new = getEdgeStart(nnodes_curr_org, nnodes_max);

#ifndef NDEBUG
   assert(cgraph_valid(cgraph));
   assert(adjcosts && edgecosts && cgraph->nodeids);
   assert(nodeid >= 0);
   assert(nnodes_curr_org < nnodes_max);
   assert(NODE_ID_UNDEFINED == cgraph->nodeids[nnodes_curr_org]);

   for( int j = 0; j < nnodes_curr_org; j++ )
      assert(GE(adjcosts[j], 0.0));
#endif

   BMScopyMemoryArray(edgecosts + start_new, adjcosts, nnodes_curr_org);

   cgraph->nnodes_curr++;

   /* adapt all other edges (going to the new node) */
   for( int i = 0; i < nnodes_curr_org; i++ )
   {
      const int end = getEdgeEnd(i, nnodes_curr_org, nnodes_max);

      assert(EQ(edgecosts[end], EDGECOST_UNDEFINED_VALUE));

      edgecosts[end] = adjcosts[i];
   }

   lastedge = getEdgeEnd(nnodes_curr_org, nnodes_curr_new, nnodes_max) - 1;
   assert(lastedge >= 0);
   assert(EQ(EDGECOST_UNDEFINED_VALUE, edgecosts[lastedge]));

   edgecosts[lastedge] = FARAWAY;
   cgraph->nodeids[nnodes_curr_org] = nodeid;

   assert(cgraph_valid(cgraph));
}


/** deletes node at the top */
void cgraph_node_deleteTop(
   CGRAPH*               cgraph              /**< new graph */
   )
{
   assert(cgraph_valid(cgraph));
   assert(cgraph->nnodes_curr <= cgraph->nnodes_max);
   assert(cgraph->nnodes_curr > 0);

   cgraph->nnodes_curr--;

#ifndef NDEBUG
   {
      int last_start;
      int last_end;

      assert(NODE_ID_UNDEFINED != cgraph->nodeids[cgraph->nnodes_curr]);

      cgraph->nodeids[cgraph->nnodes_curr] = NODE_ID_UNDEFINED;

      /* remove edge entries going to deleted vertex */
      for( int i = 0; i < cgraph->nnodes_curr; i++ )
      {
         const int end = getEdgeEnd(i, cgraph->nnodes_curr, cgraph->nnodes_max);

         assert(!EQ(cgraph->edgecosts[end], EDGECOST_UNDEFINED_VALUE));

         cgraph->edgecosts[end] = EDGECOST_UNDEFINED_VALUE;
      }

      last_start = getEdgeStart(cgraph->nnodes_curr, cgraph->nnodes_max);
      last_end = getEdgeEnd(cgraph->nnodes_curr, cgraph->nnodes_curr + 1, cgraph->nnodes_max);

      /* remove entries of deleted vertex */
      for( int i = last_start; i < last_end; i++ )
      {
         assert(!EQ(cgraph->edgecosts[i], EDGECOST_UNDEFINED_VALUE));
         cgraph->edgecosts[i] = EDGECOST_UNDEFINED_VALUE;
      }
   }

   assert(cgraph_valid(cgraph));

#endif
}


/** deletes node */
void cgraph_node_exchange(
   CGRAPH*               cgraph,             /**< new graph */
   const SCIP_Real*      adjcosts,           /**< array with edge costs to neigbors of size nnodes_curr (FARAWAY at own position) */
   int                   nodepos,            /**< the node position */
   int                   nodeid_old,         /**< node id (for debugging) */
   int                   nodeid_new          /**< node id of new entry */
   )
{
   int todo; // what do we really need? nodeid_old, nodeid_new?

   assert(cgraph);
   assert(cgraph_valid(cgraph));
   assert(nodeid_new >= 0);
   assert(cgraph->nnodes_curr <= cgraph->nnodes_max);

}

/** Get edge cost.
 * Note: quite slow, only use for debugging! */
SCIP_Real cgraph_edge_getCost(
   const CGRAPH*         cgraph,             /**< new graph */
   int                   nodeid_tail,        /**< the node id */
   int                   nodeid_head         /**< the node id */
   )
{
   assert(cgraph);
   assert(cgraph_valid(cgraph));

   return 0.0;
}

/** is the MST struct valid? */
SCIP_Bool cmst_valid(
   const CMST*           cmst                /**< the MST */
)
{

   return TRUE;
}


/** initializes MST */
SCIP_RETCODE cmst_init(
   SCIP*                 scip,               /**< SCIP data structure */
   CMST**                cmst,               /**< the MST */
   int                   maxnnodes           /**< maximum number of nodes */
   )
{
   CMST* mst;

   assert(scip && cmst);
   assert(maxnnodes > 1);

   SCIP_CALL( SCIPallocMemory(scip, cmst) );

   mst = *cmst;

   graph_heap_create(scip, maxnnodes, NULL, NULL, &(mst->heap));

   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->dist), maxnnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mst->predecessors), maxnnodes) );

   mst->nnodes_max = maxnnodes;

   assert(cmst_valid(mst));

   SCIPdebugMessage("cmst has been successfully built \n");

   return SCIP_OKAY;
}


/** frees MST */
void cmst_free(
   SCIP*                 scip,               /**< SCIP data structure */
   CMST**                cmst                /**< MST */
   )
{
   CMST* mst;

   assert(scip && cmst);

   mst = *cmst;

   SCIPfreeMemoryArray(scip, &(mst->predecessors));
   SCIPfreeMemoryArray(scip, &(mst->dist));

   graph_heap_free(scip, TRUE, TRUE, &(mst->heap));

   SCIPfreeMemory(scip, cmst);
}


/** compute MST on given graph */
void cmst_computeMst(
   const CGRAPH*         cgraph,             /**< the graph to run on */
   int                   mstroot,            /**< root for the MST */
   CMST*                 cmst,               /**< the MST */
   SCIP_Real*            mstobj              /**< objective of computed MST */
)
{
   const int nnodes_curr = getNnodesCurr(cgraph);
   const int nnodes_max = getNnodesMax(cgraph);
   SCIP_Real mstcost = 0.0;
   const SCIP_Real* const edgecosts = cgraph->edgecosts;
   int* const preds = cmst->predecessors;
   SCIP_Real* const dist = cmst->dist;
   DHEAP* const dheap = cmst->heap;
   int* const state = dheap->position;

   assert(mstobj && dheap);
   assert(cgraph_valid(cgraph));
   assert(nnodes_curr <= cmst->nnodes_max);
   assert(mstroot >= 0 && mstroot < nnodes_curr);

   for( int i = 0; i < nnodes_curr; i++ )
   {
      preds[i] = -1;
      state[i] = UNKNOWN;
      dist[i] = FARAWAY;
   }

   assert(graph_heap_isClean(dheap));

   preds[mstroot] = -1;
   dist[mstroot] = 0.0;
   graph_heap_correct(mstroot, 0.0, dheap);

   assert(dheap->size > 0);

   /* build MST */
   while( dheap->size > 0 )
   {
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int start = getEdgeStart(k, nnodes_max);

      mstcost += dist[k];

      assert(k >= 0 && k < nnodes_curr);
      assert(state[k] == CONNECT);

      for( int m = 0; m < nnodes_curr; m++ )
      {
         if( state[m] != CONNECT  )
         {
            const SCIP_Real ecost = edgecosts[start + m];

            if( ecost < dist[m] )
            {
               assert(m != k);

               dist[m] = ecost;
               preds[m] = k;

               graph_heap_correct(m, ecost, dheap);
            }
         }
      }
   }

   *mstobj = mstcost;

#ifndef NDEBUG
   for( int i = 0; i < nnodes_curr; i++ )
   {
      assert(LT(dist[i], FARAWAY));
      assert(preds[i] != -1 || i == mstroot);
   }
#endif
}

