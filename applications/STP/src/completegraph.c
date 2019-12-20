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
      if( i != (nnodes_max - 1) && getEdgeStart(i, nnodes_max) > getEdgeStart(i + 1, nnodes_max) )
      {
         SCIPdebugMessage("positions are broken \n");
         return FALSE;
      }
   }

   start = nnodes_curr > 0 ? getEdgeEnd(nnodes_curr - 1, nnodes_curr, nnodes_max) : 0;

   for( int i = start; i < nnodes_max * nnodes_max; i++ )
   {
      if( !EQ(CGRAPH_EDGECOST_UNDEFINED_VALUE, edgecosts[i]) )
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

   /* check if the edge costs are symmetric */
   for( int i = 0; i < nnodes_curr; i++ )
   {
      const int start_i = getEdgeStart(i, nnodes_max);

      for( int j = 0; j < nnodes_curr; j++ )
      {
         const int start_j = getEdgeStart(j, nnodes_max);
         const SCIP_Real cost_ij = edgecosts[start_i + j];
         const SCIP_Real cost_ji = edgecosts[start_j + i];

         assert(i != j || EQ(cost_ij, FARAWAY));

         if( !EQ(cost_ji, cost_ij) )
         {
            SCIPdebugMessage("wrong edge costs between %d and %d (%f != %f) \n", i, j, cost_ij, cost_ji);
            return FALSE;
         }
      }
   }

   return TRUE;
}


/** is the graph in sync with given node id list? */
SCIP_Bool cgraph_idsInSync(
   const CGRAPH*         cgraph,             /**< the graph */
   const int*            ids,                /**< node ids */
   int                   nids                /**< number of ids */
)
{
   assert(cgraph && ids);
   assert(cgraph->nodeids);


   if( nids != cgraph->nnodes_curr )
   {
      SCIPdebugMessage("wrong number of ids \n");
      return FALSE;
   }

   for( int i = 0; i < nids; i++ )
   {
      if( cgraph->nodeids[i] != ids[i] )
         return FALSE;
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
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->adjedgecosts), maxnnodes + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g->nodeids), maxnnodes) );

   g->nnodes_curr = 0;
   g->nnodes_max = maxnnodes;

#ifndef NDEBUG
   for( int i = 0; i < maxnnodes * maxnnodes; i++ )
      g->edgecosts[i] = CGRAPH_EDGECOST_UNDEFINED_VALUE;

   for( int i = 0; i < maxnnodes; i++ )
      g->nodeids[i] = NODE_ID_UNDEFINED;

   for( int i = 0; i < maxnnodes + 1; i++ )
      g->adjedgecosts[i] = CGRAPH_EDGECOST_UNDEFINED_VALUE;
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
   SCIPfreeMemoryArray(scip, &(g->adjedgecosts));
   SCIPfreeMemoryArray(scip, &(g->edgecosts));

   SCIPfreeMemory(scip, cgraph);
}


/** applies adjacency costs to node, but only use if smaller than existing ones. */
void cgraph_node_applyMinAdjCosts(
   CGRAPH*               cgraph,             /**< new graph */
   int                   nodepos,            /**< the node position */
   int                   nodeid              /**< the node id (for debugging only) */
   )
{
   const int nnodes_curr = getNnodesCurr(cgraph);
   const int nnodes_max = getNnodesMax(cgraph);
   const int start = getEdgeStart(nodepos, nnodes_max);
   SCIP_Real* const edgecosts = cgraph->edgecosts;
   const SCIP_Real* const adjcosts = cgraph->adjedgecosts;

   assert(nodepos >= 0 && nodepos < nnodes_curr);
   assert(cgraph->nodeids[nodepos] == nodeid);

   for( int i = start, j = 0; i != start + nnodes_curr; i++, j++ )
   {
      const SCIP_Real newcost = adjcosts[j];
      assert(GE(newcost, 0.0));
      assert(!EQ(edgecosts[i], CGRAPH_EDGECOST_UNDEFINED_VALUE));

      if( newcost < edgecosts[i] )
         edgecosts[i] = newcost;
   }

   for( int i = 0; i < nnodes_curr; i++ )
   {
      const SCIP_Real newcost = adjcosts[i];
      const int pos = getEdgeStart(i, nnodes_max) + nodepos;

      assert(!EQ(edgecosts[pos], CGRAPH_EDGECOST_UNDEFINED_VALUE));
      assert(GE(newcost, 0.0));

      if( newcost < edgecosts[pos] )
         edgecosts[pos] = newcost;
   }

   assert(EQ(adjcosts[nnodes_curr], CGRAPH_EDGECOST_UNDEFINED_VALUE));
   assert(cgraph_valid(cgraph));
}

/** add node (at the end, so at position cgraph->nnodes_curr) */
void cgraph_node_append(
   CGRAPH*               cgraph,             /**< new graph */
   int                   nodeid              /**< the node id */
   )
{
   int lastedge;
   SCIP_Real* const edgecosts = cgraph->edgecosts;
   const int nnodes_curr_org = cgraph->nnodes_curr;
   const int nnodes_curr_new = nnodes_curr_org + 1;
   const int nnodes_max = cgraph->nnodes_max;
   const int start_new = getEdgeStart(nnodes_curr_org, nnodes_max);

   assert(cgraph_valid(cgraph));
   assert(edgecosts && cgraph->nodeids);
   assert(nodeid >= 0);
   assert(nnodes_curr_org < nnodes_max);
   assert(NODE_ID_UNDEFINED == cgraph->nodeids[nnodes_curr_org]);

   cgraph->nnodes_curr++;

   for( int i = start_new; i < start_new + nnodes_curr_org; i++ )
   {
      edgecosts[i] = FARAWAY;
   }

   for( int i = 0; i < nnodes_curr_org; i++ )
   {
      const int end = getEdgeEnd(i, nnodes_curr_org, nnodes_max);
      assert(EQ(edgecosts[end], CGRAPH_EDGECOST_UNDEFINED_VALUE));

      edgecosts[end] = FARAWAY;
   }

   lastedge = getEdgeEnd(nnodes_curr_org, nnodes_curr_new, nnodes_max) - 1;
   assert(lastedge >= 0 && lastedge == start_new + nnodes_curr_org);
   assert(EQ(CGRAPH_EDGECOST_UNDEFINED_VALUE, edgecosts[lastedge]));

   edgecosts[lastedge] = FARAWAY;
   cgraph->nodeids[nnodes_curr_org] = nodeid;

   assert(cgraph_valid(cgraph));
}


/** replaces node at nodepos_new with top node */
void cgraph_node_repositionTop(
   CGRAPH*               cgraph,             /**< new graph */
   int                   nodepos_new         /**< the new node position */
   )
{
   const int nnodes_curr = getNnodesCurr(cgraph);
   const int nnodes_max = getNnodesMax(cgraph);
   const int nodepos_top = nnodes_curr - 1;
   const int start_new = getEdgeStart(nodepos_new, nnodes_max);
   const int start_top = getEdgeStart(nodepos_top, nnodes_max);

   SCIP_Real* const edgecosts = cgraph->edgecosts;

   assert(cgraph_valid(cgraph));
   assert(nodepos_new >= 0 && nodepos_new < nodepos_top);

   for( int i = 0; i < nodepos_top; i++ )
   {
      if( i != nodepos_new )
      {
         const int edgepos = getEdgeStart(i, nnodes_max) + nodepos_new;

         assert(EQ(edgecosts[edgepos], edgecosts[start_new + i]));

         edgecosts[edgepos] = edgecosts[start_top + i];
      }
   }

   assert(start_new + nodepos_top < start_top);

   BMScopyMemoryArray(edgecosts + start_new, edgecosts + start_top, nodepos_top);

   /* adapt diagonal entry */
   edgecosts[start_new + nodepos_new] = FARAWAY;

   cgraph->nodeids[nodepos_new] = cgraph->nodeids[nodepos_top];

   cgraph_node_deleteTop(cgraph);

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

         assert(!EQ(cgraph->edgecosts[end], CGRAPH_EDGECOST_UNDEFINED_VALUE));

         cgraph->edgecosts[end] = CGRAPH_EDGECOST_UNDEFINED_VALUE;
      }

      last_start = getEdgeStart(cgraph->nnodes_curr, cgraph->nnodes_max);
      last_end = getEdgeEnd(cgraph->nnodes_curr, cgraph->nnodes_curr + 1, cgraph->nnodes_max);

      /* remove entries of deleted vertex */
      for( int i = last_start; i < last_end; i++ )
      {
         assert(!EQ(cgraph->edgecosts[i], CGRAPH_EDGECOST_UNDEFINED_VALUE));
         cgraph->edgecosts[i] = CGRAPH_EDGECOST_UNDEFINED_VALUE;
      }
   }

   assert(cgraph_valid(cgraph));

#endif
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
   assert(0 && "not yet implemented");

   return 0.0;
}


/** is the MST struct valid? */
SCIP_Bool cmst_isSync(
   const CGRAPH*         cgraph,             /**< new graph */
   const CMST*           cmst                /**< the MST */
)
{
   SCIP_Real obj;
   const int nnodes = getNnodesCurr(cgraph);
   const int nnodes_max = getNnodesMax(cgraph);
   const SCIP_Real* const edgecosts = cgraph->edgecosts;
   const int* const preds = cmst->predecessors;

   assert(edgecosts && preds);

   if( cmst->nnodes_max < nnodes )
   {
      SCIPdebugMessage("mst has not enough nodes \n");
      return FALSE;
   }

   /* now check the objective */

   obj = 0.0;

   for( int i = 0; i < nnodes; i++ )
   {
      const int start = getEdgeStart(i, nnodes_max);
      const int pred = preds[i];

      if( pred == -1 )
         continue;

      assert(pred >= 0 && pred < nnodes);
      assert(pred != i);

      obj += edgecosts[start + pred];
   }

   if( !EQ(obj, cmst->mstobj) )
   {
      SCIPdebugMessage("wrong objective: %f != %f \n", obj, cmst->mstobj);
      return FALSE;
   }

   SCIPdebugMessage("graph and mst are in sync! \n");

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
   mst->mstobj = 0.0;

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
   CMST*                 cmst                /**< the MST */
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

   assert(edgecosts && preds && dist && dheap && state);
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

   cmst->mstobj = mstcost;

#ifndef NDEBUG
   for( int i = 0; i < nnodes_curr; i++ )
   {
      assert(LT(dist[i], FARAWAY)); // todo that might actually happen if the costs are too high...
      assert(preds[i] != -1 || i == mstroot);
   }

   assert(cmst_isSync(cgraph, cmst));
#endif
}

