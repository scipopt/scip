/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   extreduce_dist.c
 * @brief  (special) distance computation methods for Steiner tree extended reductions
 * @author Daniel Rehfeldt
 *
 * This file implements methods for Steiner tree problem extended reduction techniques
 * to compute and recompute (special) distances between vertices.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


// #define SCIP_DEBUG
#include "extreduce.h"
#include "misc_stp.h"
#include "portab.h"



#ifdef RED_UTIL_TIME
#include <time.h>
#endif

#define EDGE_FORBIDDEN_NONE -2

/** auxiliary data for (single) closenodes run*/
typedef struct close_nodes_run
{
   int*                  prededge;           /**< predecessors */
   SCIP_Bool*            edgemark;           /**< debug only */
   int                   startvertex;        /**< start vertex */
   SCIP_Bool             is_buildphase;      /**< in build up phase? */
} CNODESRUN;


/**@name Local methods
 *
 * @{
 */

/** returns entry of element within sorted array of size arraysize, or -1 if element could not be found
 *  NOTE: optimized binary search */
static inline
int findEntryFromSorted(
   const int*            array,              /**< array */
   unsigned int          arraysize,          /**< size */
   int                   element             /**< element to look for */
   )
{
   unsigned int lower = 0;
   unsigned int shift = arraysize;

#ifndef NDEBUG
   assert(arraysize > 0);
   for( unsigned int i = 1; i < arraysize; i++ )
      assert(array[i - 1] < array[i] );
#endif

   while( shift > 1 )
   {
      const unsigned int middle = shift / 2;

      if( element >= array[lower + middle] )
         lower += middle;

      shift -= middle;
   }

   if( element == array[lower] )
      return lower;

   return -1;
}


/** it the path between node and the close node forbidden?
 * todo better have an additional closenodes_ancestor that saves the previous node! */
static inline
SCIP_Bool closeNodesPathIsForbidden(
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< to be initialized */
   int                   edge_forbidden,     /**< forbidden edge */
   const SCIP_Bool*      edges_isForbidden,  /**< blocked edges marked (half) */
   int                   node,               /**< the node */
   int                   closenode,          /**< the close node */
   int                   closenode_pos       /**< the close node position */

)
{
   const int* const indices = distdata->closenodes_indices;
   const int* const prededges = distdata->closenodes_prededges;
   const RANGE* const range = distdata->closenodes_range;
   const int node_startpos = range[node].start;
   const int size =  range[node].end - node_startpos;
   int pred = closenode;
   int pred_pos = closenode_pos;
   SCIP_Bool isForbidden = FALSE;

   assert(distdata->closenodes_indices[closenode_pos] == closenode);
   assert(graph_edge_isInRange(g, prededges[pred_pos]));
   assert(g->head[prededges[pred_pos]] == closenode);
   assert(pred != node);
   assert(size > 0);


  // printf("check path %d->%d \n", node, closenode);

   while( pred != node )
   {
      int pred_edge;
      int pred_edgeHalf;

      /* not in first iteration? */
      if( pred != closenode )
         pred_pos = node_startpos + findEntryFromSorted(&indices[node_startpos], size, pred);

      assert(pred_pos >= node_startpos);

      pred_edge = prededges[pred_pos];
      assert(graph_edge_isInRange(g, pred_edge));

    //  graph_edge_printInfo(g, pred_edge);

      pred_edgeHalf = pred_edge / 2;

      if( edges_isForbidden && edges_isForbidden[pred_edgeHalf] )
      {
         isForbidden = TRUE;
         break;
      }

      if( pred_edgeHalf == (edge_forbidden / 2) )
      {
         isForbidden = TRUE;
         break;
      }

      pred = g->tail[pred_edge];
   }

   return isForbidden;
}


/** gets distance and path nodes */
static inline
SCIP_Real getCloseNodePath(
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< distance data */
   int                   node,               /**< the node */
   int                   closenode,          /**< the close node whose position is to be found */
   int*                  pathnodes,          /**< inner nodes */
   int*                  npathnodes         /**< number of inner nodes */
)
{
   const int* const indices = distdata->closenodes_indices;
   const RANGE* const range = distdata->closenodes_range;
   const int start = range[node].start;
   const int end = range[node].end;
   const int size = end - start;
   int position;
   SCIP_Real dist = -1.0;
   assert(size > 0);

   *npathnodes = 0;

   position = findEntryFromSorted(&indices[start], size, closenode);

   if( position >= 0 )
   {
      const int* const prededges = distdata->closenodes_prededges;
      int pred = closenode;
      int pred_pos = start + position;

      assert(indices[pred_pos] == closenode);
      assert(graph_edge_isInRange(g, prededges[pred_pos]));
      assert(g->head[prededges[pred_pos]] == closenode);
      assert(pred != node && closenode != node);

      dist = distdata->closenodes_distances[start + position];

      while( pred != node )
      {
         int pred_edge;

         /* not in first iteration? */
         if( pred != closenode )
         {
            pred_pos = start + findEntryFromSorted(&indices[start], size, pred);
            pathnodes[(*npathnodes)++] = pred;
         }

         assert(pred_pos >= start);

         pred_edge = prededges[pred_pos];
         assert(graph_edge_isInRange(g, pred_edge));

         pred = g->tail[pred_edge];
      }
   }

   return dist;
}


/** returns distance of closenode from node, or -1.0 if this distance is not stored in close nodes list of node */
static inline
SCIP_Real getCloseNodeDistance(
   const DISTDATA*       distdata,           /**< to be initialized */
   int                   node,               /**< the node */
   int                   closenode           /**< the close node whose position is to be found */
)
{
   const int* const indices = distdata->closenodes_indices;
   const RANGE* const range = distdata->closenodes_range;
   const int start = range[node].start;
   const int end = range[node].end;
   const int size = end - start;
   int position;
   SCIP_Real dist = -1.0;

   assert(size >= 0);

   if( size == 0 )
   {
      assert(distdata->hasPathReplacement);
      return dist;
   }

   position = findEntryFromSorted(&indices[start], size, closenode);

   if( position >= 0 )
   {
      assert(indices[start + position] == closenode);
      dist = distdata->closenodes_distances[start + position];
   }

   return dist;
}


/** as above, but with forbidden edge */
static inline
SCIP_Real getCloseNodeDistanceForbidden(
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< to be initialized */
   int                   edge_forbidden,     /**< forbidden edge */
   const SCIP_Bool*      edges_isForbidden,  /**< blocked edges marked or NULL */
   int                   node,               /**< the node */
   int                   closenode           /**< the close node whose position is to be found */
)
{
   const int* const indices = distdata->closenodes_indices;
   const RANGE* const range = distdata->closenodes_range;
   const int start = range[node].start;
   const int end = range[node].end;
   const int size = end - start;
   int position;
   SCIP_Real dist = -1.0;

   assert(size >= 0);

   if( size == 0 )
   {
      assert(distdata->hasPathReplacement);
      return dist;
   }

   position = findEntryFromSorted(&indices[start], size, closenode);

   if( position >= 0 )
   {
      assert(indices[start + position] == closenode);

      if( !closeNodesPathIsForbidden(g, distdata, edge_forbidden, edges_isForbidden, node, closenode, start + position) )
      {
         dist = distdata->closenodes_distances[start + position];
      }
   }

   return dist;
}


/** as above, but with forbidden last edge */
static inline
SCIP_Real getCloseNodeDistanceForbiddenLast(
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< to be initialized */
   int                   node,               /**< the node */
   int                   closenode,          /**< the close node whose position is to be found */
   int                   lastedge_node2close /**< last edge */
)
{
   const int* const prededges = distdata->closenodes_prededges;
   const int* const indices = distdata->closenodes_indices;
   const RANGE* const range = distdata->closenodes_range;
   const int start = range[node].start;
   const int end = range[node].end;
   const int size = end - start;
   int position_rel;
   SCIP_Real dist = -1.0;

   if( size == 0 )
   {
      assert(distdata->hasPathReplacement);
      return dist;
   }

   assert(size > 0);
   assert(graph_edge_isInRange(g, lastedge_node2close));
   assert(g->head[lastedge_node2close] == closenode);

   position_rel = findEntryFromSorted(&indices[start], size, closenode);

   if( position_rel >= 0 )
   {
      const int position_abs = start + position_rel;
      const int prededge = prededges[position_abs];

      assert(indices[position_abs] == closenode);
      assert(graph_edge_isInRange(g, prededge));

      if( prededge != lastedge_node2close )
      {
         assert((prededge / 2) != (lastedge_node2close / 2));

         dist = distdata->closenodes_distances[position_abs];
      }
   }

   return dist;
}


/** as above, but with forbidden edge/edges and known equality value */
static inline
SCIP_Real getCloseNodeDistanceForbiddenEq(
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< to be initialized */
   SCIP_Real             dist_eq,            /**< critical distance or -1.0 if not known */
   int                   edge_forbidden,     /**< forbidden edge */
   const SCIP_Bool*      edges_isForbidden,  /**< blocked edges marked */
   int                   node,               /**< the node */
   int                   closenode           /**< the close node whose position is to be found */
)
{
   const int* const indices = distdata->closenodes_indices;
   const RANGE* const range = distdata->closenodes_range;
   const int start = range[node].start;
   const int end = range[node].end;
   const int size = end - start;
   int position;
   SCIP_Real dist = -1.0;

   if( size == 0 )
   {
      assert(distdata->hasPathReplacement);
      return dist;
   }

   assert(size > 0);

   position = findEntryFromSorted(&indices[start], size, closenode);

   if( position >= 0 )
   {
      assert(indices[start + position] == closenode);
      dist = distdata->closenodes_distances[start + position];

      if( GT(dist, dist_eq) )
      {
         dist = -1.0;
      }
      else if( EQ(dist, dist_eq)  )
      {
         if( closeNodesPathIsForbidden(g, distdata, edge_forbidden, edges_isForbidden, node,
             closenode, start + position) )
         {
            dist = -1.0;
         }
      }
   }

   return dist;
}


/** compute paths root list */
static inline
SCIP_RETCODE distDataPathRootsInsertRoot(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   halfedge,           /**< halfedge to insert at */
   int                   root,               /**< root to insert */
   DISTDATA*             distdata            /**< distance data */
)
{
   int* const pathroot_blocksizes = distdata->pathroot_blocksizes;
   int* const pathroot_blocksizesmax = distdata->pathroot_blocksizesmax;

   assert(scip && g && distdata && pathroot_blocksizes && pathroot_blocksizesmax);
   assert(halfedge >= 0 && halfedge < g->edges / 2);
   assert(root >= 0 && root < g->knots);
   assert(pathroot_blocksizes[halfedge] <= pathroot_blocksizesmax[halfedge]);

   /* need to reallocate? */
   if( pathroot_blocksizes[halfedge] == pathroot_blocksizesmax[halfedge] )
   {
      const int oldsize = pathroot_blocksizesmax[halfedge];
      int newsize;

      assert(oldsize >= 0);

      if( oldsize == 0 )
      {
         assert(distdata->pathroot_blocks[halfedge] == NULL);
         newsize = 2;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(distdata->pathroot_blocks[halfedge]), newsize) );
      }
      else
      {
         assert(distdata->pathroot_blocks[halfedge] != NULL);
         newsize = 2 * pathroot_blocksizesmax[halfedge];
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(distdata->pathroot_blocks[halfedge]), oldsize, newsize) );
      }

      pathroot_blocksizesmax[halfedge] = newsize;
   }

   assert(pathroot_blocksizes[halfedge] < pathroot_blocksizesmax[halfedge] );

   /* now add the root */
   distdata->pathroot_blocks[halfedge][pathroot_blocksizes[halfedge]].pathroot_id = root;
   distdata->pathroot_blocks[halfedge][pathroot_blocksizes[halfedge]++].pathroot_nrecomps = distdata->pathroot_nrecomps[root];

   return SCIP_OKAY;
}


/** compute paths root list */
static
SCIP_RETCODE distDataPathRootsInitialize(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* pathroot_blocksizes;
   int* pathroot_blocksizesmax;
   int* pathroot_blockcount;
   PRSTATE** pathroot_blocks;
   const int* const closenodes_prededges = distdata->closenodes_prededges;
   const int nnodes = g->knots;
   const int halfnedges = g->edges / 2;
   const RANGE* const range_closenodes = distdata->closenodes_range;

   assert(scip && g && closenodes_prededges && distdata);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->pathroot_nrecomps), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->pathroot_isdirty), nnodes) );

   for( int i = 0; i < nnodes; i++ )
      distdata->pathroot_isdirty[i] = FALSE;

   for( int i = 0; i < nnodes; i++ )
      distdata->pathroot_nrecomps[i] = 0;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(pathroot_blocks), halfnedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pathroot_blocksizes), halfnedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pathroot_blocksizesmax), halfnedges) );

   SCIP_CALL( SCIPallocBufferArray(scip, &(pathroot_blockcount), halfnedges) );

   for( int k = 0; k < halfnedges; k++ )
      pathroot_blocksizes[k] = 0;

   /* compute the edge range sizes */
   for( int j = 0; j < range_closenodes[nnodes - 1].end; j++ )
   {
      const int edge = closenodes_prededges[j] / 2;
      assert(edge >= 0 && edge < halfnedges);
      assert(g->oeat[2 * edge] != EAT_FREE);

      pathroot_blocksizes[edge]++;
   }

   for( int e = 0; e < halfnedges; e++ )
   {
      const int size = pathroot_blocksizes[e];

      /* is edge used? */
      if( size > 0 )
      {
         assert(g->oeat[2 * e] != EAT_FREE);

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pathroot_blocks[e]), size) );
      }
      else
      {
         pathroot_blocks[e] = NULL;
      }
   }

   /* fill the path roots in */

   for( int k = 0; k < halfnedges; k++ )
      pathroot_blockcount[k] = 0;

   for( int k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 )
         continue;

      for( int j = range_closenodes[k].start; j < range_closenodes[k].end; j++ )
      {
         const int edge = closenodes_prededges[j] / 2;
         const int blockcount = pathroot_blockcount[edge];

         assert(edge >= 0 && edge < halfnedges);
         assert(g->oeat[2 * edge] != EAT_FREE);
         assert(blockcount < pathroot_blocksizes[edge]);

         pathroot_blocks[edge][blockcount].pathroot_id = k;
         pathroot_blocks[edge][blockcount].pathroot_nrecomps = 0;

         pathroot_blockcount[edge]++;
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < halfnedges; e++ )
      assert(pathroot_blockcount[e] == pathroot_blocksizes[e]);
#endif

   for( int k = 0; k < halfnedges; k++ )
      pathroot_blocksizesmax[k] = pathroot_blocksizes[k];

   distdata->pathroot_blocks = pathroot_blocks;
   distdata->pathroot_blocksizes = pathroot_blocksizes;
   distdata->pathroot_blocksizesmax = pathroot_blocksizesmax;

   SCIPfreeBufferArray(scip, &pathroot_blockcount);

   return SCIP_OKAY;
}


/** initializes run from 'startvertex' */
static inline
SCIP_RETCODE closeNodesRunInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   DISTDATA*             distdata,           /**< distance data */
   CNODESRUN*            cnodesrun           /**< auxiliary data */
)
{
   DIJK* const dijkdata = distdata->dijkdata;
   SCIP_Real* const dist = dijkdata->node_distance;
   DHEAP* const dheap = dijkdata->dheap;
   const int nnodes = g->knots;
   const int startvertex = cnodesrun->startvertex;

   assert(dheap->size == 0 && distdata->closenodes_maxnumber > 0);
   assert(graph_knot_isInRange(g, startvertex));
   assert(distdata->closenodes_range[startvertex].start == distdata->closenodes_range[startvertex].end);

   SCIP_CALL( SCIPallocBufferArray(scip, &(cnodesrun->prededge), nnodes) );

#ifndef NDEBUG
   SCIP_CALL( SCIPallocMemoryArray(scip, &(cnodesrun->edgemark), g->edges / 2) );

   for( int e = 0; e < g->edges / 2; e++ )
      (cnodesrun->edgemark)[e] = FALSE;

   for( int k = 0; k < nnodes; k++ )
   {
      (cnodesrun->prededge)[k] = -1;
      assert(dist[k] == FARAWAY && dheap->position[k] == UNKNOWN);
   }
#endif

   dist[startvertex] = 0.0;
   dijkdata->visitlist[0] = startvertex;
   graph_heap_correct(startvertex, 0.0, dheap);

   dijkdata->nvisits = 1;

   assert(dheap->size == 1);

   return SCIP_OKAY;
}


/** computes sorted shortest path list to constant number of neighbors */
static inline
SCIP_RETCODE closeNodesRunCompute(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   CNODESRUN*            cnodesrun,          /**< auxiliary data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   const SDPROFIT* sdprofit = NULL;
   DIJK* const dijkdata = distdata->dijkdata;
   int* const visitlist = dijkdata->visitlist;
   SCIP_Real* const dist = dijkdata->node_distance;
   DHEAP* const dheap = dijkdata->dheap;
   STP_Bool* const visited = dijkdata->node_visited;
   int* const prededge = cnodesrun->prededge;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const int* const edgeid = dcsr->edgeid;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const SCIP_Real* const pc_costshifts = dijkdata->node_bias;
   RANGE* const range_closenodes = distdata->closenodes_range;
   int* const closenodes_indices = distdata->closenodes_indices;
   int* const closenodes_prededges = distdata->closenodes_prededges;
   SCIP_Real* const closenodes_distances = distdata->closenodes_distances;
   int nvisits = dijkdata->nvisits;
   int nclosenodes = 0;
   const int startvertex = cnodesrun->startvertex;
   const int closenodes_limit = distdata->closenodes_maxnumber;
   const SCIP_Bool isPc = graph_pc_isPc(g);
   const SCIP_Bool withProfit = (distdata->sdistdata != NULL && distdata->sdistdata->sdprofit != NULL);
   const SCIP_Bool is_buildphase = cnodesrun->is_buildphase;

   assert(dcsr && dist && visitlist && visited && dheap && prededge && cnodesrun->edgemark);
   assert(range_closenodes && closenodes_indices && closenodes_prededges);
   assert(dheap->size == 1);
   assert(!isPc || pc_costshifts);

   if( withProfit )
      sdprofit = distdata->sdistdata->sdprofit;

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;
      const SCIP_Real k_dist = dist[k];
      int k_pred = -1;

      if( k != startvertex )
      {
         const int closenodes_pos = range_closenodes[startvertex].end;

#ifndef NDEBUG
         assert(graph_edge_isInRange(g, prededge[k]));
         assert(cnodesrun->edgemark[prededge[k] / 2] == FALSE);  /* make sure that the edge is not marked twice */
         assert(closenodes_pos < distdata->closenodes_totalsize && state[k] == CONNECT);

         cnodesrun->edgemark[prededge[k] / 2] = TRUE;
#endif

         closenodes_indices[closenodes_pos] = k;
         closenodes_distances[closenodes_pos] = isPc ? (k_dist + pc_costshifts[k]) : (k_dist);

         if( is_buildphase )
         {
            closenodes_prededges[closenodes_pos] = prededge[k];
         }
         else
         {
            closenodes_prededges[closenodes_pos] = prededge[k];
            SCIP_CALL( distDataPathRootsInsertRoot(scip, g, prededge[k] / 2, startvertex, distdata) );
         }

         range_closenodes[startvertex].end++;

         if( ++nclosenodes >= closenodes_limit )
            break;

         if( withProfit )
            k_pred = g->tail[prededge[k]];
      }

      /* correct adjacent nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];
         assert(g->mark[m]);

         if( state[m] != CONNECT )
         {
            SCIP_Real distnew = isPc ?
               k_dist + cost_csr[e] - pc_costshifts[m]
             : k_dist + cost_csr[e];

            assert(GE(distnew, 0.0));

            if( withProfit && m != k_pred )
            {
               SCIP_Real profitBias = reduce_sdprofitGetProfit(sdprofit, k, m, k_pred);
               profitBias = MIN(profitBias, cost_csr[e]);
               profitBias = MIN(profitBias, k_dist);
               distnew -= profitBias;
            }

            if( LT(distnew, dist[m]) )
            {
               if( !visited[m] )
               {
                  visitlist[nvisits++] = m;
                  visited[m] = TRUE;
               }

               dist[m] = distnew;
               prededge[m] = edgeid[e];
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   dijkdata->nvisits = nvisits;

   return SCIP_OKAY;
}


/** sort close nodes list of node */
static inline
void closeNodesRunSort(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to sort for */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const RANGE* const range_closenodes = distdata->closenodes_range;
   int* const closenodes_indices = distdata->closenodes_indices;
   SCIP_Real* const closenodes_distances = distdata->closenodes_distances;
   const int start = range_closenodes[node].start;
   const int length = range_closenodes[node].end - start;

   assert(g && distdata);
   assert(node >= 0 && node < g->knots);
   assert(g->grad[node] > 0);
   assert(length >= 0);
   assert(distdata->closenodes_prededges);

   SCIPsortIntIntReal(&closenodes_indices[start], &(distdata->closenodes_prededges[start]), &closenodes_distances[start], length);

#ifndef NDEBUG
   for( int i = 1; i < length; i++ )
      assert(closenodes_indices[start] < closenodes_indices[start + i]);
#endif
}


/** exits */
static inline
void closeNodesRunExit(
   SCIP*                 scip,               /**< SCIP */
   CNODESRUN*            cnodesrun           /**< auxiliary data */
)
{
   SCIPfreeMemoryArrayNull(scip, &(cnodesrun ->edgemark));
   SCIPfreeBufferArray(scip, &(cnodesrun ->prededge));
}


#ifdef RED_UTIL_TIME

typedef struct pathroot_info
{
   int pathroot_id;
} PRINFO;

//#define USE_STRUCT
/** compute paths root list */
static
SCIP_RETCODE distDataPathRootsInitializeBench(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* pathroot_blocksizes;
#ifdef USE_STRUCT
   PRINFO** pathroot_blocks;
#else
   int** pathroot_blocks;

#endif
   const int halfnedges = 1000000;


   clock_t start, end;
     double cpu_time_used;

     start = clock();


   SCIP_CALL( SCIPallocMemoryArray(scip, &(pathroot_blocks), halfnedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(pathroot_blocksizes), halfnedges) );


   for( int e = 0; e < halfnedges; e++ )
   {
      const int size = 1 + halfnedges % 32;

      /* is edge used? */
      if( size > 0 )
      {
         assert(g->oeat[2 * e] != EAT_FREE);

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pathroot_blocks[e]), size) );
#ifdef USE_STRUCT
         pathroot_blocks[e][0].pathroot_id = size;
#else
         pathroot_blocks[e][0] = size;
#endif
      }
      else
      {
         pathroot_blocks[e] = NULL;
      }
   }


   for( int e = halfnedges - 1; e >= 0 ; e-- )
   {
#ifdef USE_STRUCT
      const int size = pathroot_blocks[e][0].pathroot_id;
#else
      const int size = pathroot_blocks[e][0];
#endif

      /* is edge used? */
      if( size > 0 )
      {
         assert(pathroot_blocks[e] != NULL);

         SCIPfreeBlockMemoryArray(scip, &(pathroot_blocks[e]), size);
      }
      else
      {
         assert(pathroot_blocks[e] == NULL);
      }
   }

   SCIPfreeMemoryArray(scip, &pathroot_blocksizes);
   SCIPfreeMemoryArray(scip, &pathroot_blocks);

   end = clock();
   cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

   printf("time %f \n", cpu_time_used);

   exit(1);

   return SCIP_OKAY;
}
#endif

/** frees paths root list */
static
void distDataPathRootsFree(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* pathroot_blocksizesmax;
   PRSTATE** pathroot_blocks;
   const int halfnedges = g->edges / 2;

   assert(scip && g && distdata);

   pathroot_blocksizesmax = distdata->pathroot_blocksizesmax;
   pathroot_blocks = distdata->pathroot_blocks;

   for( int e = halfnedges - 1; e >= 0 ; e-- )
   {
      const int maxsize = pathroot_blocksizesmax[e];

      /* is edge used? */
      if( maxsize > 0 )
      {
         assert(pathroot_blocks[e] != NULL);

         SCIPfreeBlockMemoryArray(scip, &(pathroot_blocks[e]), maxsize);
      }
      else
      {
         assert(pathroot_blocks[e] == NULL);
      }
   }

   SCIPfreeMemoryArray(scip, &distdata->pathroot_blocksizesmax);
   SCIPfreeMemoryArray(scip, &distdata->pathroot_blocksizes);
   SCIPfreeMemoryArray(scip, &distdata->pathroot_blocks);
   SCIPfreeMemoryArray(scip, &distdata->pathroot_isdirty);
   SCIPfreeMemoryArray(scip, &distdata->pathroot_nrecomps);
}


/** computes sorted shortest path list to constant number of neighbors */
static
SCIP_RETCODE distDataComputeCloseNodes(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   CNODESRUN*            cnodesrun,          /**< auxiliary data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   assert(scip && g && distdata && cnodesrun);

   SCIP_CALL( closeNodesRunInit(scip, g, distdata, cnodesrun) );

   SCIP_CALL( closeNodesRunCompute(scip, g, cnodesrun, distdata) );

   closeNodesRunExit(scip, cnodesrun);

   /* sort close nodes according to their index */
   closeNodesRunSort(g, cnodesrun->startvertex, distdata);

   return SCIP_OKAY;
}


/** sets sizes */
static
void distDataInitSizes(
   const GRAPH*          g,                  /**< graph data structure */
   int                   maxnclosenodes,     /**< maximum number of close nodes to each node */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   int nnodes_undeleted;
   int closenodes_totalsize;

   graph_get_nVET(g, &nnodes_undeleted, NULL, NULL);
   assert(nnodes_undeleted >= 1 && maxnclosenodes >= 1);

   closenodes_totalsize = nnodes_undeleted * maxnclosenodes;

   assert(closenodes_totalsize >= 1);

   distdata->closenodes_totalsize = closenodes_totalsize;
   distdata->closenodes_maxnumber = maxnclosenodes;
}


/** allocates memory for some distance data members */
static
SCIP_RETCODE distDataAllocateNodesArrays(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = g->knots;
   const int closenodes_totalsize = distdata->closenodes_totalsize;

   assert(scip && g && distdata);
   assert(distdata->closenodes_totalsize > 0);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_range), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_indices), closenodes_totalsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_distances), closenodes_totalsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_prededges), closenodes_totalsize) );

   return SCIP_OKAY;
}


/** re-compute close nodes for default distance from vertex1 */
static inline
void distDataRecomputeNormalDist(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   CNODESRUN closenodes = { .edgemark = NULL, .prededge = NULL,
                            .startvertex = vertex1, .is_buildphase = FALSE };

   assert(distdata->pathroot_isdirty[vertex1]);
   assert(!distdata->sdistdata || !distdata->sdistdata->sdprofit);

   distdata->pathroot_nrecomps[vertex1]++;
   SCIP_CALL_ABORT( distDataComputeCloseNodes(scip, g, &closenodes, distdata) );

   graph_dijkLimited_reset(g, distdata->dijkdata);

   SCIPdebugMessage("vertex %d is dirty, recompute \n", vertex1);

   distdata->pathroot_isdirty[vertex1] = FALSE;
}


/** returns (normal) shortest path distance between vertex1 and vertex2
 *  and provides inner shortest path vertices. Returns -1.0 if no shortest path was found */
static inline
SCIP_Real distDataGetSp(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   int*                  pathnodes,          /**< inner nodes */
   int*                  npathnodes,         /**< number of inner nodes */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   /* neighbors list not valid anymore? */
   if( distdata->pathroot_isdirty[vertex1] )
   {
      distDataRecomputeNormalDist(scip, g, vertex1, distdata);
   }

   dist = getCloseNodePath(g, distdata, vertex1, vertex2, pathnodes, npathnodes);

   return dist;
}


/** Gets shortest v1->v2 (standard) distance.
 *  Returns -1.0 if the distance is not known. */
static inline
SCIP_Real distDataGetNormalDist(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);

   /* neighbors list not valid anymore? */
   if( distdata->pathroot_isdirty[vertex1] )
   {
      /* recompute */
      distDataRecomputeNormalDist(scip, g, vertex1, distdata);
   }

   /* look in neighbors list of vertex1 */
   dist = getCloseNodeDistance(distdata, vertex1, vertex2);

   return dist;
}


/** as above, but with forbidden edge */
static inline
SCIP_Real distDataGetNormalDistForbidden(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   edge_forbidden,     /**< forbidden edge */
   const SCIP_Bool*      edges_isForbidden,  /**< blocked edges marked or NULL */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);
   assert(edge_forbidden == EDGE_FORBIDDEN_NONE || graph_edge_isInRange(g, edge_forbidden));

   /* neighbors list not valid anymore? */
   if( distdata->pathroot_isdirty[vertex1] )
   {
      distDataRecomputeNormalDist(scip, g, vertex1, distdata);
   }

   dist = getCloseNodeDistanceForbidden(g, distdata, edge_forbidden, edges_isForbidden, vertex1, vertex2);

   return dist;
}


/** as above, but with forbidden last edge */
static inline
SCIP_Real distDataGetNormalDistForbiddenLast(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   int                   lastedge12,         /**< forbidden last edge for v1->v2 path */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);
   assert(graph_edge_isInRange(g, lastedge12));

   /* neighbors list not valid anymore? */
   if( distdata->pathroot_isdirty[vertex1] )
   {
      distDataRecomputeNormalDist(scip, g, vertex1, distdata);
   }

   dist = getCloseNodeDistanceForbiddenLast(g, distdata, vertex1, vertex2, lastedge12);

   return dist;
}

/** as above, but with forbidden edges and known equality value */
static inline
SCIP_Real distDataGetNormalDistForbiddenEq(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real             dist_eq,            /**< critical distance */
   int                   edge_forbidden,     /**< forbidden edge */
   const SCIP_Bool*      edges_isForbidden,  /**< blocked edges marked, or NULL */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);
   assert(graph_edge_isInRange(g, edge_forbidden));

   /* neighbors list not valid anymore? */
   if( distdata->pathroot_isdirty[vertex1] )
   {
      /* recompute */
      distDataRecomputeNormalDist(scip, g, vertex1, distdata);
   }

   /* look in neighbors list of vertex1 */
   dist = getCloseNodeDistanceForbiddenEq(g, distdata, dist_eq, edge_forbidden, edges_isForbidden, vertex1, vertex2);

   return dist;
}


/** returns v1->v2 special distance */
static inline
SCIP_Real distDataGetSpecialDist(
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);
   assert(distdata->sdistdata);

   dist = reduce_sdGetSd(g, vertex1, vertex2, FARAWAY, 0.0, distdata->sdistdata);

   return dist;
}


/** returns v1->v2 special distance, but only for SDs with intermediate terms */
static inline
SCIP_Real distDataGetSpecialDistIntermedTerms(
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);
   assert(distdata->sdistdata);

   dist = reduce_sdGetSdIntermedTerms(g, vertex1, vertex2, FARAWAY, 0.0, distdata->sdistdata);

   return dist;
}



/**@} */

/**@name Interface methods
 *
 * @{
 */



/** initializes distance data */
SCIP_RETCODE extreduce_distDataInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   int                   maxnclosenodes,     /**< maximum number of close nodes to each node */
   SCIP_Bool             computeSD,          /**< also compute special distances? */
   SCIP_Bool             useBias,            /**< use bias? */
   DISTDATA**            distdata            /**< to be initialized */
)
{
   CNODESRUN closenodes = { NULL, NULL, -1, FALSE };
   const int nnodes = g->knots;
   RANGE* range_closenodes;
   DHEAP* dheap;
   DISTDATA* dist;

   assert(distdata && g && scip && g->dcsr_storage);
   assert(graph_valid_dcsr(g, FALSE));
   assert(!graph_pc_isPcMw(g) || !g->extended);

   SCIP_CALL( SCIPallocMemory(scip, distdata) );
   dist = *distdata;
   dist->hasPathReplacement = FALSE;

   distDataInitSizes(g, maxnclosenodes, dist);
   SCIP_CALL( distDataAllocateNodesArrays(scip, g, dist) );

   SCIP_CALL( graph_dijkLimited_init(scip, g, &(dist->dijkdata)) );

   if( graph_pc_isPc(g) )
   {
      SCIP_CALL( graph_dijkLimited_initPcShifts(scip, g, dist->dijkdata) );
   }

   if( computeSD )
   {
      if( useBias )
      {
         SCIP_CALL( reduce_sdInitBiasedBottleneck(scip, g, &(dist->sdistdata)) );
      }
      else
      {
         SCIP_CALL( reduce_sdInit(scip, g, &(dist->sdistdata)) );
      }
   }
   else
   {
      assert(!useBias);
      dist->sdistdata = NULL;
   }

   range_closenodes = dist->closenodes_range;

   /* compute close nodes to each not yet deleted node */
   for( int k = 0; k < nnodes; k++ )
   {
      range_closenodes[k].start = (k == 0) ? 0 : range_closenodes[k - 1].end;
      range_closenodes[k].end = range_closenodes[k].start;

      if( !g->mark[k] )
      {
         assert(g->grad[k] == 0 || graph_pc_isPcMw(g));
      //   assert(g->dcsr_storage->range[k].end == g->dcsr_storage->range[k].start);

         continue;
      }

      if( g->grad[k] == 0 )
      {
         assert(graph_pc_isPcMw(g));
         assert(graph_pc_knotIsNonLeafTerm(g, k));

         continue;
      }

      assert(g->grad[k] > 0);

      closenodes.startvertex = k;
      closenodes.is_buildphase = TRUE;
      SCIP_CALL( distDataComputeCloseNodes(scip, g, &closenodes, dist) );

      graph_dijkLimited_reset(g, dist->dijkdata);
   }

   /* store for each edge the roots of all paths it is used for */
   SCIP_CALL( distDataPathRootsInitialize(scip, g, dist) );

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &dheap) );
   dist->dheap = dheap;

   assert(dist->closenodes_prededges);
   assert(useBias || extreduce_distCloseNodesAreValid(scip, g, dist));

   return SCIP_OKAY;
}


/** updates data for edge deletion */
void extreduce_distDataDeleteEdge(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   edge,               /**< edge to be deleted */
   DISTDATA*             distdata            /**< distance data */
)
{
   const int halfedge = edge / 2;
   PRSTATE* const pathrootblock = distdata->pathroot_blocks[halfedge];
   SCIP_Bool* const node_isdirty = distdata->pathroot_isdirty;
   RANGE* const range_closenodes = distdata->closenodes_range;
   const int* const node_nrecomps = distdata->pathroot_nrecomps;
   const int npathroots = distdata->pathroot_blocksizes[halfedge];

   assert(scip && distdata && distdata->pathroot_blocks);
   assert(edge >= 0);
 //  assert(g->oeat[edge] == EAT_FREE && g->ieat[edge] == EAT_FREE);

   for( int k = 0; k < npathroots; k++ )
   {
      const int node = pathrootblock[k].pathroot_id;
      const int nrecomps = pathrootblock[k].pathroot_nrecomps;

      if( nrecomps == node_nrecomps[node] && !node_isdirty[node] )
      {
         node_isdirty[node] = TRUE;

#ifndef NDEBUG
         for( int j = range_closenodes[node].start; j < range_closenodes[node].end; j++ )
         {
            distdata->closenodes_indices[j] = -1;
            distdata->closenodes_distances[j] = -FARAWAY;
         }
#endif
         range_closenodes[node].end = range_closenodes[node].start;
      }
   }

   if( distdata->pathroot_blocksizesmax[halfedge] > 0  )
   {
      assert(distdata->pathroot_blocksizes[halfedge] <= distdata->pathroot_blocksizesmax[halfedge]);

      SCIPfreeBlockMemoryArray(scip, &(distdata->pathroot_blocks[halfedge]), distdata->pathroot_blocksizesmax[halfedge]);
   }
   else
   {
      assert(distdata->pathroot_blocks[halfedge] == NULL);
   }

   distdata->pathroot_blocksizes[halfedge] = 0;
   distdata->pathroot_blocksizesmax[halfedge] = 0;
}


/** recomputes shortest paths for dirty nodes */
void extreduce_distDataRecomputeDirtyPaths(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Bool* pathroot_isdirty;
   const int* isMarked;
   const int nnodes = graph_get_nNodes(g);

   assert(scip && distdata);
   assert(graph_isMarked(g));

   pathroot_isdirty = distdata->pathroot_isdirty;
   isMarked = g->mark;

   for( int i = 0; i < nnodes; i++ )
   {
      if( pathroot_isdirty[i] && isMarked[i] )
      {
         assert(g->grad[i] > 0);
         distDataRecomputeNormalDist(scip, g, i, distdata);
         pathroot_isdirty[i] = FALSE;
      }
   }
}


/** returns (normal) shortest path distance between vertex1 and vertex2
 *  and provides inner shortest path vertices. Returns -1.0 if no shortest path was found */
SCIP_Real extreduce_distDataGetSp(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   int*                  pathnodes,          /**< inner nodes */
   int*                  npathnodes,         /**< number of inner nodes */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(scip && g && distdata && pathnodes && npathnodes);
   assert(graph_knot_isInRange(g, vertex1) && graph_knot_isInRange(g, vertex2));

   dist = distDataGetSp(scip, g, vertex1, vertex2, pathnodes, npathnodes, distdata);

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}


/** Gets bottleneck (or special) distance between v1 and v2.
 *  Will check shortest known v1->v2 path, but NOT shortest known v2->v1 path.
 *  Returns -1.0 if no distance is known. (might only happen for RPC or PC) */
SCIP_Real extreduce_distDataGetSd(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(distdata);

   dist = distDataGetNormalDist(scip, g, vertex1, vertex2, distdata);

   assert(EQ(dist, -1.0) || dist >= 0.0);

   if( distdata->sdistdata )
   {
      const SCIP_Real dist_sd = distDataGetSpecialDist(g, vertex1, vertex2, distdata);

      if( EQ(dist, -1.0) || dist_sd < dist )
      {
     //    printf("%f->%f \n", dist, dist_sd);

         dist = dist_sd;
      }

      assert(GE(dist, 0.0));
   }

   return dist;
}


/** Gets bottleneck (or special) distance between v1 and v2.
 *  Will check shortest known v1->v2 path, and also shortest known v2->v1 path if no v1-v2 path is known.
 *  Returns -1.0 if no distance is known. (might only happen for RPC or PC) */
SCIP_Real extreduce_distDataGetSdDouble(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist = distDataGetNormalDist(scip, g, vertex1, vertex2, distdata);

   /* no distance found? */
   if( dist < -0.5 )
   {
      assert(EQ(dist, -1.0));
      dist = distDataGetNormalDist(scip, g, vertex2, vertex1, distdata);
   }
   else
   {
      const SCIP_Real distrev = distDataGetNormalDist(scip, g, vertex2, vertex1, distdata);

      if( distrev > -0.5 && distrev < dist  )
         dist = distrev;

      assert(GE(dist, 0.0));
   }

   if( distdata->sdistdata )
   {
      const SCIP_Real dist_sd = distDataGetSpecialDist(g, vertex1, vertex2, distdata);

      if( EQ(dist, -1.0) || dist_sd < dist )
         dist = dist_sd;

      assert(GE(dist, 0.0));
   }

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}


/** As 'extreduce_distDataGetSdDouble', but with critial value for early abort  */
SCIP_Real extreduce_distDataGetSdDoubleEq(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real             dist_eq,            /**< critical distance */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist = distDataGetNormalDist(scip, g, vertex1, vertex2, distdata);

   /* no distance found? */
   if( dist < -0.5 )
   {
      assert(EQ(dist, -1.0));
      dist = distDataGetNormalDist(scip, g, vertex2, vertex1, distdata);

      if( dist > -0.5 && LT(dist, dist_eq) )
         return dist;
   }
   else
   {
      const SCIP_Real distrev = distDataGetNormalDist(scip, g, vertex2, vertex1, distdata);

      if( distrev > -0.5 && distrev < dist  )
         dist = distrev;

      if( dist > -0.5 && LT(dist, dist_eq) )
         return dist;

      assert(GE(dist, 0.0));
   }

   if( distdata->sdistdata )
   {
      const SCIP_Real dist_sd = distDataGetSpecialDist(g, vertex1, vertex2, distdata);

      if( EQ(dist, -1.0) || dist_sd < dist )
         dist = dist_sd;

      assert(GE(dist, 0.0));
   }

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}



/** Same as extreduce_distDataGetSdDouble, but only takes paths that do not include any edges marked as blocked. */
SCIP_Real extreduce_distDataGetSdDoubleForbidden(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   DISTDATA* distdata = extdata->distdata;
   SCIP_Real dist = -1.0;

   assert(scip && g && distdata);

   /* NOTE: does not work for pseudo-elimination, because paths might get destroyed */
   if( extInitialCompIsEdge(extdata) || extInitialCompIsGenStar(extdata) )
   {
      const SCIP_Bool* const edges_isForbidden = extdata->sdeq_edgesIsForbidden;
      dist = distDataGetNormalDistForbidden(scip, g, EDGE_FORBIDDEN_NONE, edges_isForbidden, vertex1, vertex2, distdata);

      /* no distance found? */
      if( dist < -0.5 )
      {
         assert(EQ(dist, -1.0));
         dist = distDataGetNormalDistForbidden(scip, g, EDGE_FORBIDDEN_NONE, edges_isForbidden, vertex1, vertex2, distdata);
      }
      else
      {
         const SCIP_Real distrev = distDataGetNormalDistForbidden(scip, g, EDGE_FORBIDDEN_NONE, edges_isForbidden, vertex1, vertex2, distdata);

         if( distrev > -0.5 && distrev < dist  )
            dist = distrev;

         assert(GE(dist, 0.0));
      }
   }

   if( distdata->sdistdata )
   {
      if( !Is_term(g->term[vertex1]) || !Is_term(g->term[vertex2]) )
      {
         dist = distDataGetSpecialDistIntermedTerms(g, vertex1, vertex2, distdata);
      }
   }

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}

/** Same as extreduce_distDataGetSdDouble, but only takes paths that do not include
 *  given edge or any edges marked as blocked.
 *  User needs to provide (known) equality value */
SCIP_Real extreduce_distDataGetSdDoubleForbiddenEq(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real             dist_eq,            /**< critical distance */
   int                   edge_forbidden,     /**< forbidden edge */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   DISTDATA* distdata = extdata->distdata;
   SCIP_Real dist = -1.0;

   assert(scip && g && distdata);
   assert(graph_edge_isInRange(g, edge_forbidden));

   /* NOTE: does not work for pseudo-elimination, because paths might get destroyed */
   if( extInitialCompIsEdge(extdata) || extInitialCompIsGenStar(extdata) )
   {
      const SCIP_Bool* const edges_isForbidden = extdata->sdeq_edgesIsForbidden;
      dist = distDataGetNormalDistForbiddenEq(scip, g, dist_eq, edge_forbidden, edges_isForbidden, vertex1, vertex2, distdata);

      if( EQ(dist_eq, dist) )
         return dist;

      /* no distance found? */
      if( dist < -0.5 )
      {
         assert(EQ(dist, -1.0));
         dist = distDataGetNormalDistForbiddenEq(scip, g, dist_eq, edge_forbidden, edges_isForbidden, vertex2, vertex1, distdata);
      }
      else
      {
         const SCIP_Real distrev = distDataGetNormalDistForbiddenEq(scip, g, dist_eq, edge_forbidden, edges_isForbidden, vertex2, vertex1, distdata);

         if( distrev > -0.5 && distrev < dist  )
            dist = distrev;

         assert(GE(dist, 0.0));
      }

      if( EQ(dist_eq, dist) )
         return dist;
   }

   if( distdata->sdistdata )
   {
      if( !Is_term(g->term[vertex1]) || !Is_term(g->term[vertex2]) )
      {
         dist = distDataGetSpecialDistIntermedTerms(g, vertex1, vertex2, distdata);
      }
   }

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}


/** Same as extreduce_distDataGetSdDouble, but only takes paths that do not contain given edge */
SCIP_Real extreduce_distDataGetSdDoubleForbiddenSingle(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   edge_forbidden,     /**< edge */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist = -1.0;

   assert(scip && g && distdata);
   assert(graph_edge_isInRange(g, edge_forbidden));

   /* NOTE: does not work for pseudo-elimination, because paths might get destroyed */
   {
      dist = distDataGetNormalDistForbidden(scip, g, edge_forbidden, NULL, vertex1, vertex2, distdata);

      /* no distance found? */
      if( dist < -0.5 )
      {
         assert(EQ(dist, -1.0));

         dist = distDataGetNormalDistForbidden(scip, g, edge_forbidden, NULL, vertex2, vertex1, distdata);
      }
      else
      {
         const SCIP_Real distrev = distDataGetNormalDistForbidden(scip, g, edge_forbidden, NULL, vertex2, vertex1, distdata);

         if( distrev > -0.5 && distrev < dist  )
            dist = distrev;

         assert(GE(dist, 0.0));
      }
   }

   if( distdata->sdistdata )
   {
      if( !Is_term(g->term[vertex1]) || !Is_term(g->term[vertex2]) )
      {
         const SCIP_Real dist_sd = distDataGetSpecialDistIntermedTerms(g, vertex1, vertex2, distdata);

         assert(GE(dist_sd, 0.0));

         if( LT(dist_sd, dist) || dist < -0.5 )
            dist = dist_sd;
      }
   }

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}



/** Same as extreduce_distDataGetSdDouble, but only takes paths that do not contain given last edges.
 *  NOTE: Behaves peculiarly for terminal-paths! */
SCIP_Real extreduce_distDataGetSdDoubleForbiddenLast(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   int                   lastedge12,         /**< forbidden last edge on v1->v2 path; needs to be in-edge of v2 */
   int                   lastedge21,         /**< forbidden last edge on v2->v1 path; needs to be in-edge of v1 */
   DISTDATA*             distdata            /**< distance data */
)
{
   SCIP_Real dist;

   assert(scip && g && distdata);
   assert(graph_edge_isInRange(g, lastedge12));
   assert(graph_edge_isInRange(g, lastedge21));

   dist = distDataGetNormalDistForbiddenLast(scip, g, vertex1, vertex2, lastedge12, distdata);

   /* no distance found? */
   if( dist < -0.5 )
   {
      assert(EQ(dist, -1.0));

      dist = distDataGetNormalDistForbiddenLast(scip, g, vertex2, vertex1, lastedge21, distdata);
   }
   else
   {
      const SCIP_Real distrev = distDataGetNormalDistForbiddenLast(scip, g, vertex2, vertex1, lastedge21, distdata);

      if( distrev > -0.5 && distrev < dist  )
         dist = distrev;

      assert(GE(dist, 0.0));
   }

   if( distdata->sdistdata )
   {
#if 1
      if( !Is_term(g->term[vertex1]) && !Is_term(g->term[vertex2]) )
      {
         const SCIP_Real dist_sd = distDataGetSpecialDistIntermedTerms(g, vertex1, vertex2, distdata);

         assert(GE(dist_sd, 0.0));

         if( LT(dist_sd, dist) || dist < -0.5 )
            dist = dist_sd;
      }
#endif
   }

   assert(EQ(dist, -1.0) || dist >= 0.0);

   return dist;
}


/** frees distance data */
void extreduce_distDataFree(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   DISTDATA**            distdata            /**< to be freed */
)
{
   DISTDATA* dist;

   assert(scip && graph && distdata);

   dist = *distdata;

   graph_heap_free(scip, TRUE, TRUE, &(dist->dheap));
   SCIPfreeMemoryArray(scip, &(dist->closenodes_range));
   SCIPfreeMemoryArray(scip, &(dist->closenodes_indices));
   SCIPfreeMemoryArray(scip, &(dist->closenodes_distances));
   SCIPfreeMemoryArrayNull(scip, &(dist->closenodes_prededges));

   if( dist->sdistdata )
      reduce_sdFree(scip, &(dist->sdistdata));

   distDataPathRootsFree(scip, graph, dist);
   graph_dijkLimited_free(scip, &(dist->dijkdata));

   SCIPfreeMemory(scip, distdata);
}

/**@} */
