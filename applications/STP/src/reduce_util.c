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

/**@file   reduce_util.c
 * @brief  utility methods for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem reduction techniques.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "reduce.h"
#include "misc_stp.h"
#include "portab.h"

#ifdef RED_UTIL_TIME
#include <time.h>
#endif

/** returns entry of element within sorted array of size arraysize, or -1 if element could not be found */
static
int findEntryFromSorted(
   const int* array,
   int arraysize,
   int element)
{
   int l = 0;
   int u = arraysize - 1;

#ifndef NDEBUG
   for( int i = 1; i < arraysize; i++ )
      assert(array[i - 1] < array[i] );
#endif

   while( l <= u )
   {
      const int m = (u - l) / 2 + l;

      if( array[m] < element )
         l = m + 1;
      else if( array[m] == element )
         return m;
      else
         u = m - 1;
   }
   // try http://eigenjoy.com/2011/09/09/binary-search-revisited/ ?

   return -1;
}


/** sort close nodes list of node */
static
void distDataSortCloseNodes(
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

   SCIPsortIntReal(&closenodes_indices[start], &closenodes_distances[start], length);

#ifndef NDEBUG
   for( int i = 1; i < length; i++ )
      assert(closenodes_indices[start] < closenodes_indices[start + i]);
#endif

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

   assert(size > 0);

   position = findEntryFromSorted(&indices[start], size, closenode);

   if( position >= 0 )
   {
      assert(indices[start + position] == closenode);
      dist = distdata->closenodes_distances[start + position];
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
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* pathroot_blocksizes;
   int* pathroot_blocksizesmax;
   int* pathroot_blockcount;
   PRSTATE** pathroot_blocks;

   const int nnodes = g->knots;
   const int halfnedges = g->edges / 2;
   const RANGE* const range_closenodes = distdata->closenodes_range;

   assert(scip && g && closenodes_edges && distdata);

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
      const int edge = closenodes_edges[j];
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
         const int edge = closenodes_edges[j];
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


/** limited Dijkstra to constant number of neighbors, taking SD distances into account */
static
SCIP_RETCODE distDataComputeCloseNodesSD(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   startvertex,        /**< start vertex */
   SCIP_Bool             is_buildphase,      /**< in build up phase? */
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DIJK*                 dijkdata,           /**< limited Dijkstra data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* edgemark = NULL;

#ifndef NDEBUG
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgemark, g->edges / 2) );
    for( int e = 0; e < g->edges / 2; e++ )
       edgemark[e] = FALSE;

#endif

   assert(0);


#if 0
   const int closenodes_limit = distdata->closenodes_maxnumber;

      for( int v = prededge[k]; v != startvertex; v = g->tail[prededge[v]] )
      {
         /* edge already used? */
         if( prededge[v] ) // todo
         assert(prededge[v] >= 0);

      }
#endif
   SCIPfreeMemoryArrayNull(scip, &edgemark);

   return SCIP_OKAY;
}


/** computes sorted shortest path list to constant number of neighbors */
static
SCIP_RETCODE distDataComputeCloseNodes(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   startvertex,        /**< start vertex */
   SCIP_Bool             is_buildphase,      /**< in build up phase? */
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DIJK*                 dijkdata,           /**< limited Dijkstra data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* const visitlist = dijkdata->visitlist;
   SCIP_Real* const dist = dijkdata->distance;
   DHEAP* const dheap = dijkdata->dheap;
   STP_Bool* const visited = dijkdata->visited;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const int* const edgeid = dcsr->edgeid;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   RANGE* const range_closenodes = distdata->closenodes_range;
   int* const closenodes_indices = distdata->closenodes_indices;
   SCIP_Real* const closenodes_distances = distdata->closenodes_distances;
   const int nnodes = g->knots;
   const int closenodes_limit = distdata->closenodes_maxnumber;
   int* prededge;
   SCIP_Bool* edgemark = NULL;
   int nvisits;
   int nclosenodes;

   assert(dcsr && g && dist && visitlist && visited && dheap && dijkdata && distdata);
   assert(dheap->size == 0 && closenodes_limit > 0);
   assert(startvertex >= 0 && startvertex < g->knots);
   assert(range_closenodes[startvertex].start == range_closenodes[startvertex].end);

   SCIP_CALL( SCIPallocBufferArray(scip, &prededge, nnodes) );

#ifndef NDEBUG
   SCIP_CALL( SCIPallocMemoryArray(scip, &edgemark, g->edges / 2) );
   for( int e = 0; e < g->edges / 2; e++ )
      edgemark[e] = FALSE;

   for( int k = 0; k < g->knots; k++ )
   {
      prededge[k] = -1;
      assert(dist[k] == FARAWAY && state[k] == UNKNOWN);
   }
#endif

   nvisits = 0;
   nclosenodes = 0;
   dist[startvertex] = 0.0;
   visitlist[nvisits++] = startvertex;
   graph_heap_correct(startvertex, 0.0, dheap);

   assert(dheap->size == 1);

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      if( k != startvertex )
      {
         const int closenodes_pos = range_closenodes[startvertex].end;

#ifndef NDEBUG
         assert((prededge[k] >= 0 && prededge[k] < g->edges));
         assert(edgemark[prededge[k] / 2] == FALSE);  /* make sure that the edge is not marked twice */
         assert(closenodes_pos < distdata->closenodes_totalsize && state[k] == CONNECT);

         edgemark[prededge[k] / 2] = TRUE;
#endif

         closenodes_indices[closenodes_pos] = k;
         closenodes_distances[closenodes_pos] = dist[k];

         if( is_buildphase )
         {
            assert(closenodes_edges != NULL);
            closenodes_edges[closenodes_pos] = prededge[k] / 2;
         }
         else
         {
            assert(closenodes_edges == NULL);
            SCIP_CALL( distDataPathRootsInsertRoot(scip, g, prededge[k] / 2, k, distdata) );
         }

         range_closenodes[startvertex].end++;

         if( ++nclosenodes >= closenodes_limit )
            break;
      }

      /* correct adjacent nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];
         assert(g->mark[m]);

         if( state[m] != CONNECT )
         {
            const SCIP_Real distnew = dist[k] + cost_csr[e];

            if( distnew < dist[m] )
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

   SCIPfreeMemoryArrayNull(scip, &edgemark);
   SCIPfreeBufferArray(scip, &prededge);

   /* sort close nodes according to their index */
   distDataSortCloseNodes(g, startvertex, distdata);

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
   SCIP_Bool             computeSD,          /**< also compute special distances? */
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

   return SCIP_OKAY;
}


/** initializes distance data */
SCIP_RETCODE reduce_distDataInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   maxnclosenodes,     /**< maximum number of close nodes to each node */
   SCIP_Bool             computeSD,          /**< also compute special distances? */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = g->knots;
   int* closenodes_edges;
   RANGE* range_closenodes;
   DIJK* dijkdata;
   DHEAP* dheap;

   assert(distdata && g && scip && g->dcsr_storage);
   assert(graph_valid_dcsr(g, FALSE));
   assert(!graph_pc_isPcMw(g) || !g->extended);

   distDataInitSizes(g, maxnclosenodes, distdata);
   SCIP_CALL( distDataAllocateNodesArrays(scip, g, computeSD, distdata) );

   SCIP_CALL( SCIPallocMemory(scip, &distdata->dijkdata) );
   dijkdata = distdata->dijkdata;
   SCIP_CALL( graph_dijkLimited_init(scip, g, dijkdata) );

   range_closenodes = distdata->closenodes_range;

   SCIP_CALL( SCIPallocMemoryArray(scip, &closenodes_edges, distdata->closenodes_totalsize) );

   /* compute close nodes to each not yet deleted node */
   for( int k = 0; k < nnodes; k++ )
   {
      range_closenodes[k].start = (k == 0) ? 0 : range_closenodes[k - 1].end;
      range_closenodes[k].end = range_closenodes[k].start;

      if( !g->mark[k] )
      {
         assert(g->grad[k] == 0 || graph_pc_isPcMw(g));
         assert(g->dcsr_storage->range[k].end == g->dcsr_storage->range[k].start);

         continue;
      }


      if( g->grad[k] == 0 )
      {
         assert(graph_pc_isPcMw(g));
         assert(graph_pc_knotIsNonLeafTerm(g, k));

         continue;
      }

      assert(g->grad[k] > 0);

      if( computeSD )
         SCIP_CALL( distDataComputeCloseNodesSD(scip, g, k, TRUE, closenodes_edges, dijkdata, distdata) );
      else
         SCIP_CALL( distDataComputeCloseNodes(scip, g, k, TRUE, closenodes_edges, dijkdata, distdata) );

      graph_dijkLimited_reset(g, dijkdata);
   }

   /* store for each edge the roots of all paths it is used for */
   SCIP_CALL( distDataPathRootsInitialize(scip, g, closenodes_edges, distdata) );

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &dheap) );
   distdata->dheap = dheap;

   SCIPfreeMemoryArray(scip, &closenodes_edges);

   return SCIP_OKAY;
}


/** updates data for edge deletion */
void reduce_distDataDeleteEdge(
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
   assert(g->oeat[edge] == EAT_FREE && g->ieat[edge] == EAT_FREE);

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

/** gets bottleneck (or special) distance between v1 and v2; -1.0 if no distance is known */
SCIP_Real reduce_distDataGetSD(
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

   /* try to find SD via Duin's approximation todo */
   // if( distdata->nodeSDpaths_dirty[vertex1] && !pcmw
   // if( distdata->nodeSDpaths_dirty[vertex2] )


   /* neighbors list not valid anymore? */
   if( distdata->pathroot_isdirty[vertex1] )
   {
      /* recompute */
      distdata->pathroot_nrecomps[vertex1]++;
      SCIP_CALL_ABORT( distDataComputeCloseNodes(scip, g, vertex1, FALSE, NULL, distdata->dijkdata, distdata) );

      graph_dijkLimited_reset(g, distdata->dijkdata);

      SCIPdebugMessage("vertex %d is dirty, recompute \n", vertex1);

      distdata->pathroot_isdirty[vertex1] = FALSE;
   }

   /* look in neighbors list of vertex1 */
   dist = getCloseNodeDistance(distdata, vertex1, vertex2);

   /* if no success, binary search on neighbors of vertex2? todo too expensive? */

   //   if( distdata->pathroot_isdirty[vertex2] )


   return dist;
}

/** frees members of distance data */
void reduce_distDataFreeMembers(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   DISTDATA*             distdata            /**< to be freed */
)
{
   graph_heap_free(scip, TRUE, TRUE, &(distdata->dheap));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_range));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_indices));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_distances));

   distDataPathRootsFree(scip, graph, distdata);
   graph_dijkLimited_freeMembers(scip, distdata->dijkdata);

   SCIPfreeMemory(scip, &(distdata->dijkdata));
}


/** initialize permanent extension data struct */
SCIP_RETCODE reduce_extPermaInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   STP_Bool*             edgedeleted,        /**< edge array to mark which directed edge can be removed */
   EXTPERMA*             extperm             /**< (uninitialized) extension data */
)
{
   SCIP_Bool* isterm = NULL;
   SCIP_Real* bottleneckDistNode = NULL;
   SCIP_Real* pcSdToNode = NULL;
   int* tree_deg = NULL;
   const int nnodes = graph_get_nNodes(graph);
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(scip && extperm);

   SCIP_CALL( SCIPallocMemoryArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tree_deg, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bottleneckDistNode, nnodes) );

   if( pcmw )
      SCIP_CALL( SCIPallocMemoryArray(scip, &pcSdToNode, nnodes) );

   extperm->edgedeleted = edgedeleted;
   extperm->isterm = isterm;
   extperm->bottleneckDistNode = bottleneckDistNode;
   extperm->pcSdToNode = pcSdToNode;
   extperm->tree_deg = tree_deg;
   extperm->nnodes = nnodes;

   if( pcmw )
   {
      assert(pcSdToNode);

      for( int k = 0; k < nnodes; k++ )
         pcSdToNode[k] = -1.0;
   }

   graph_get_isTerm(graph, isterm);

   for( int k = 0; k < nnodes; k++ )
   {
      bottleneckDistNode[k] = -1.0;

      if( graph->mark[k] )
         tree_deg[k] = 0;
      else
         tree_deg[k] = -1;
   }

   assert(reduce_extPermaIsClean(graph, extperm));

   return SCIP_OKAY;
}


/** initialize permanent extension data struct */
SCIP_Bool reduce_extPermaIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperm             /**< extension data */
)
{
   const SCIP_Real* bottleneckDistNode = NULL;
   const SCIP_Real* pcSdToNode = NULL;
   const int* tree_deg = NULL;
   const int nnodes = graph_get_nNodes(graph);

   assert(extperm);

   bottleneckDistNode = extperm->bottleneckDistNode;
   pcSdToNode = extperm->pcSdToNode;
   tree_deg = extperm->tree_deg;

   if( nnodes != extperm->nnodes )
      return FALSE;

   for( int i = 0; i < nnodes; i++ )
   {
      if( !(tree_deg[i] == 0 || tree_deg[i] == -1) )
         return FALSE;

      if( bottleneckDistNode[i] > -0.5 )
         return FALSE;

      if( pcSdToNode && !EQ(pcSdToNode[i], -1.0) )
         return FALSE;
   }

   return TRUE;
}


/** frees members of extension data */
void reduce_extPermaFreeMembers(
   SCIP*                 scip,               /**< SCIP */
   EXTPERMA*             extperm             /**< extension data */
)
{
   assert(scip && extperm);

   SCIPfreeMemoryArrayNull(scip, &(extperm->pcSdToNode));
   SCIPfreeMemoryArray(scip, &(extperm->bottleneckDistNode));
   SCIPfreeMemoryArray(scip, &(extperm->tree_deg));
   SCIPfreeMemoryArray(scip, &(extperm->isterm));

   extperm->nnodes = -1;
}
