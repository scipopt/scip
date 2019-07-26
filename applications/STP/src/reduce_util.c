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

#define EXT_PATHROOTS_BUFFER 20

#include "reduce.h"
#include "misc_stp.h"

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
static
SCIP_RETCODE distDataComputePathRootsList(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int halfnedges_undeleted;
   int* edgehits;
   const int buffer_per_edge = EXT_PATHROOTS_BUFFER;
   const int nnodes = g->knots;
   const int halfnedges = g->edges / 2;
   const RANGE* const range_closenodes = distdata->closenodes_range;
   int* pathroots;
   RANGE* range_pathroots;

   assert(scip && g && closenodes_edges && distdata);

   graph_get_NVET(g, NULL, &halfnedges_undeleted, NULL);

   halfnedges_undeleted /= 2;
   distdata->pathroots_totalsize = range_closenodes[g->knots - 1].end + halfnedges_undeleted * buffer_per_edge;

   assert(halfnedges_undeleted >= 1);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->pathroots_range), halfnedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->pathroots), distdata->pathroots_totalsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(edgehits), halfnedges) );

   for( int k = 0; k < halfnedges; k++ )
      edgehits[k] = 0;

   /* compute the edge range sizes */
   for( int j = 0; j < range_closenodes[nnodes - 1].end; j++ )
   {
      const int edge = closenodes_edges[j];
      assert(edge >= 0 && edge < halfnedges);
      edgehits[edge]++;
   }

   /* compute the edge ranges */
   range_pathroots = distdata->pathroots_range;
   range_pathroots[0].start = 0;
   range_pathroots[0].end = 0;

   for( int e = 1; e < halfnedges; e++ )
   {
      const int prevedge = e - 1;

      /* has previous edge been deleted? */
      if( g->oeat[2 * prevedge] == EAT_FREE )
      {
         assert(edgehits[prevedge] == 0);
         range_pathroots[e].start = range_pathroots[prevedge].start;
      }
      else
      {
         range_pathroots[e].start = range_pathroots[prevedge].start + edgehits[prevedge] + buffer_per_edge;
      }

      range_pathroots[e].end = range_pathroots[e].start;
   }

   assert(range_pathroots[halfnedges - 1].end <= distdata->pathroots_totalsize);

   pathroots = distdata->pathroots;

   /* fill the path roots in */
   for( int k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 )
         continue;

      for( int j = range_closenodes[k].start; j < range_closenodes[k].end; j++ )
      {
         const int edge = closenodes_edges[j];

         assert(edge >= 0 && edge < halfnedges);
         assert(g->oeat[2 * edge] != EAT_FREE);

         pathroots[range_pathroots[edge].end++] = k;

         assert(edge == halfnedges - 1 ||
               range_pathroots[edge].end + buffer_per_edge <= range_pathroots[edge + 1].start);
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < halfnedges; e++ )
   {
      assert(range_pathroots[e].end - range_pathroots[e].start == edgehits[e]);
      assert(e == halfnedges - 1 || g->oeat[2 * e] == EAT_FREE
            || range_pathroots[e].end + buffer_per_edge == range_pathroots[e + 1].start );
   }
#endif

   SCIPfreeBufferArray(scip, &edgehits);

   return SCIP_OKAY;
}

/** limited Dijkstra to constant number of neighbors, taking SD distances into account */
static
SCIP_RETCODE distDataComputeCloseNodesSD(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   startvertex,        /**< start vertex */
   int                   closenodes_limit,   /**< close nodes limit */
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DIJK*                 dijkdata,           /**< limited Dijkstra data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* edgemark = NULL;

#ifndef NDEBUG
   SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, g->edges / 2) );
    for( int e = 0; e < g->edges / 2; e++ )
       edgemark[e] = FALSE;

#endif

   assert(0);


#if 0
      for( int v = prededge[k]; v != startvertex; v = g->tail[prededge[v]] )
      {
         /* edge already used? */
         if( prededge[v] ) // todo
         assert(prededge[v] >= 0);

      }
#endif
   SCIPfreeBufferArrayNull(scip, &edgemark);

   return SCIP_OKAY;
}


/** limited Dijkstra to constant number of neighbors */
static
SCIP_RETCODE distDataComputeCloseNodes(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   startvertex,        /**< start vertex */
   int                   closenodes_limit,   /**< close nodes limit */
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
   int* prededge;
   SCIP_Bool* edgemark = NULL;
   int nvisits;
   int nclosenodes;

   assert(dcsr && g && dist && visitlist && visited && dheap && dijkdata && distdata);
   assert(dheap->size == 0);
   assert(startvertex >= 0 && startvertex < g->knots);
   assert(range_closenodes[startvertex].start == range_closenodes[startvertex].end);

   SCIP_CALL( SCIPallocBufferArray(scip, &prededge, nnodes) );

#ifndef NDEBUG
   SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, g->edges / 2) );
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
         closenodes_edges[closenodes_pos] = prededge[k] / 2;
         closenodes_distances[closenodes_pos] = dist[k];

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

   SCIPfreeBufferArrayNull(scip, &edgemark);
   SCIPfreeBufferArray(scip, &prededge);

   return SCIP_OKAY;
}

/** returns maximum total number of close nodes that should be computed */
static
int distDataGetCloseNodesTotalSize(
   const GRAPH*          g,                  /**< graph data structure */
   int                   maxnclosenodes      /**< maximum number of close nodes to each node */
)
{
   int nnodes_undeleted;
   int closenodes_totalsize;

   graph_get_NVET(g, &nnodes_undeleted, NULL, NULL);
   assert(nnodes_undeleted >= 1 && maxnclosenodes >= 1);

   closenodes_totalsize = nnodes_undeleted * maxnclosenodes;

   assert(closenodes_totalsize >= 1);

   return closenodes_totalsize;
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

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->nodepaths_dirty), nnodes) );
   BMSclearMemoryArray(distdata->nodepaths_dirty, nnodes);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_range), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_indices), closenodes_totalsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_distances), closenodes_totalsize) );

   return SCIP_OKAY;
}

static
void distDataSortCloseNodes(
   const GRAPH*          g,                  /**< graph data structure */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = g->knots;
   const RANGE* const range_closenodes = distdata->closenodes_range;
   int* const closenodes_indices = distdata->closenodes_indices;
   SCIP_Real* const closenodes_distances = distdata->closenodes_distances;

   assert(g && distdata);

   for( int k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] != 0 )
      {
         const int start = range_closenodes[k].start;
         const int length = range_closenodes[k].end - start;
         assert(length > 0);

         SCIPsortIntReal(&closenodes_indices[start], &closenodes_distances[start], length);

#ifndef NDEBUG
         for( int i = 1; i < length; i++ )
            assert(closenodes_indices[start] < closenodes_indices[start + i]);
#endif
      }
   }
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
   DIJK dijkdata;
   DHEAP* dheap;

   assert(distdata && g && scip && g->dcsr_storage);
   assert(maxnclosenodes >= 1);
   assert(graph_valid_dcsr(g, FALSE));

   distdata->closenodes_totalsize = distDataGetCloseNodesTotalSize(g, maxnclosenodes);

   SCIP_CALL( distDataAllocateNodesArrays(scip, g, computeSD, distdata) );

   SCIP_CALL( SCIPallocBufferArray(scip, &closenodes_edges, distdata->closenodes_totalsize) );

   /* build auxiliary data */
   SCIP_CALL( graph_dijkLimited_init(scip, g, &dijkdata) );


   range_closenodes = distdata->closenodes_range;

   /* compute close nodes to each not yet deleted node */
   for( int k = 0; k < nnodes; k++ )
   {
      range_closenodes[k].start = (k == 0) ? 0 : range_closenodes[k - 1].end;
      range_closenodes[k].end = range_closenodes[k].start;

      if( g->grad[k] == 0 )
         continue;

      if( computeSD )
         SCIP_CALL( distDataComputeCloseNodesSD(scip, g, k, maxnclosenodes, closenodes_edges, &dijkdata, distdata) );
      else
         SCIP_CALL( distDataComputeCloseNodes(scip, g, k, maxnclosenodes, closenodes_edges, &dijkdata, distdata) );

      graph_dijkLimited_reset(g, &dijkdata);
   }

   /* sort close nodes according to their index */
   distDataSortCloseNodes(g, distdata);

   /* store for each edge the roots of all paths it is used for */
   SCIP_CALL( distDataComputePathRootsList(scip, g, closenodes_edges, distdata) );

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &dheap) );
   distdata->dheap = dheap;

   graph_dijkLimited_free(scip, &dijkdata);

   SCIPfreeBufferArray(scip, &closenodes_edges);

   return SCIP_OKAY;
}

/** gets bottleneck (or special) distance between v1 and v2; -1.0 if no distance is known */
SCIP_Real reduce_distDataGetSD(
   DISTDATA*             distdata,           /**< distance data */
   int                   vertex1,            /**< first vertex */
   int                   vertex2             /**< second vertex */
)
{
   SCIP_Real dist;

   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);

   /* try to find SD via Duin's approximation todo */
   // if( distdata->nodeSDpaths_dirty[vertex1]
   // if( distdata->nodeSDpaths_dirty[vertex2] )


   /* neighbors list not valid anymore? */
   if( distdata->nodepaths_dirty[vertex1] )
   {
      /* recompute */
      assert(0);


      distdata->nodepaths_dirty[vertex1] = FALSE;
   }

   /* look in neighbors list of vertex1 */
   dist = getCloseNodeDistance(distdata, vertex1, vertex2);


   /* if no success, binary search on neighbors of vertex2? todo too expensive? */

   //   if( distdata->nodepaths_dirty[vertex2] )


   return dist;
}

/** frees members of distance data */
void reduce_distDataFree(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   DISTDATA*             distdata            /**< to be freed */
)
{
   graph_heap_free(scip, TRUE, TRUE, &(distdata->dheap));
   SCIPfreeMemoryArray(scip, &(distdata->nodepaths_dirty));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_range));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_indices));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_distances));
   SCIPfreeMemoryArray(scip, &(distdata->pathroots_range));
   SCIPfreeMemoryArray(scip, &(distdata->pathroots));
}
