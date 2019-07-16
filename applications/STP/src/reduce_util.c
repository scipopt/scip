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

#define EXT_CLOSENODES_MAXN 64
#define EXT_PATHROOTS_BUFFER 20

#include "reduce.h"

#if 0
static
int findEntry(const int* array, int arraysize, int value)
{
   int l = 0;
   int u = arraysize - 1;

   while( l <= u )
   {
      const int m = (u - l) / 2 + l;

      if( array[m] < value )
         l = m + 1;
      else if( array[m] == value )
         return m;
      else
         u = m - 1;
   }
   // try http://eigenjoy.com/2011/09/09/binary-search-revisited/ ?

   return -1;
}
#endif

/** compute paths root list */
static
SCIP_RETCODE distDataComputePathRootsList(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int*                  closenodes_edges,   /**< edges used to reach close nodes */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int nedges_undelete;
   RANGE* const range_closenodes = distdata->closenodes_range;

   assert(0);
   assert(scip && g && closenodes_edges && distdata);

   graph_get_NVET(g, NULL, &nedges_undelete, NULL);

   nedges_undelete /= 2;
   distdata->pathroots_totalsize = range_closenodes[g->knots - 1].end + nedges_undelete * EXT_PATHROOTS_BUFFER;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->pathroots_range), distdata->pathroots_totalsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->pathroots), distdata->pathroots_totalsize) );

   // todo: fill the stuff in


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
   int*                  edgemark,           /**< edge mark array */
   DLIMITED*             dijkdata,           /**< limited Dijkstra data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(edgemark[k] < startvertex);
   }
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
   int*                  edgemark,           /**< edge mark array */
   DLIMITED*             dijkdata,           /**< limited Dijkstra data */
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

   int nvisits;
   int nclosenodes;

   assert(dcsr && g && dist && visitlist && visited && dheap && dijkdata && distdata);
   assert(dheap->size == 0);
   assert(startvertex >= 0 && startvertex < g->knots);
   assert(range_closenodes[startvertex].start == range_closenodes[startvertex].end);

   SCIP_CALL( SCIPallocBufferArray(scip, &prededge, nnodes) );

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      prededge[k] = -1;
      assert(dist[k] == FARAWAY && state[k] == UNKNOWN);
   }

   for( int e = 0; e < g->edges / 2; e++ )
      assert(edgemark[e] < startvertex);
#endif

   nvisits = 0;
   nclosenodes = 0;
   dist[startvertex] = 0.0;
   state[startvertex] = CONNECT;
   visitlist[nvisits++] = startvertex;

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;
      const int closenodes_pos = range_closenodes[startvertex].end;

#ifndef NDEBUG
      assert(prededge[k] >= 0 && prededge[k] < g->edges);
      assert(edgemark[prededge[k] / 2] < startvertex);  /* make sure that the edge is not marked twice */
      assert(closenodes_pos < distdata->closenodes_totalsize);
      edgemark[prededge[k] / 2] = startvertex;
#endif

      closenodes_indices[closenodes_pos] = k;
      closenodes_edges[closenodes_pos] = prededge[k] / 2;
      closenodes_distances[closenodes_pos] = dist[k];

      range_closenodes[startvertex].end++;

      if( ++nclosenodes >= closenodes_limit )
         break;

      assert(state[k] == CONNECT);

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

   SCIPfreeBufferArray(scip, &prededge);

   return SCIP_OKAY;
}



/** allocates memory for some distance data members */
static
SCIP_RETCODE distDataAllocateClosenodesArrays(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             computeSD,          /**< also compute special distances? */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = g->knots;
   int nnodes_undeleted;
   int closenodes_totalsize;

   graph_get_NVET(g, &nnodes_undeleted, NULL, NULL);
   assert(nnodes_undeleted >= 1);

   closenodes_totalsize = nnodes_undeleted * EXT_CLOSENODES_MAXN;
   distdata->closenodes_totalsize = closenodes_totalsize;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->nodepaths_dirty), nnodes) );
   BMSclearMemoryArray(distdata->nodepaths_dirty, nnodes);

   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_range), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_indices), closenodes_totalsize) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_distances), closenodes_totalsize) );

   return SCIP_OKAY;
}

/** initializes distance data */
SCIP_RETCODE reduce_distDataInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             computeSD,          /**< also compute special distances? */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = g->knots;
   const int halfnedges = g->edges / 2;
   int* closenodes_edges;
   int* edgemark = NULL;
   RANGE* range_closenodes;
   DLIMITED dijkdata;

   assert(distdata && g && scip);
   assert(halfnedges >= 1);
   assert(distdata->nodepaths_dirty == NULL && distdata->closenodes_distances == NULL && distdata->closenodes_range == NULL);
   assert(distdata->closenodes_indices == NULL && distdata->pathroots == NULL && distdata->pathroots_range == NULL);

   SCIP_CALL( distDataAllocateClosenodesArrays(scip, g, computeSD, distdata) );

   SCIP_CALL( SCIPallocBufferArray(scip, &closenodes_edges, distdata->closenodes_totalsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgemark, halfnedges) );

   for( int e = 0; e < halfnedges; e++ )
      edgemark[e] = -1;

   range_closenodes = distdata->closenodes_range;

   SCIP_CALL( graph_dijkLimited_init(scip, g, &dijkdata) ); // todo keep the heap for recomputation in distdata?

   /* compute close nodes to each not yet deleted node */
   for( int k = 0; k < nnodes; k++ )
   {
      range_closenodes[k].start = (k == 0) ? 0 : range_closenodes[k - 1].end;
      range_closenodes[k].end = range_closenodes[k].start;

      if( g->grad[k] == 0 )
         continue;

      if( computeSD )
         SCIP_CALL( distDataComputeCloseNodesSD(scip, g, k, EXT_CLOSENODES_MAXN, closenodes_edges, edgemark, &dijkdata, distdata) );
      else
         SCIP_CALL( distDataComputeCloseNodes(scip, g, k, EXT_CLOSENODES_MAXN, closenodes_edges, edgemark, &dijkdata, distdata) );

      graph_dijkLimited_reset(g, &dijkdata);
   }


   /* sort */

   SCIP_CALL( distDataComputePathRootsList(scip, g, closenodes_edges, distdata) );

   graph_dijkLimited_free(scip, &dijkdata);

   SCIPfreeBufferArray(scip, &edgemark);
   SCIPfreeBufferArray(scip, &closenodes_edges);

   return SCIP_OKAY;
}

/** gets bottleneck (or special) distance between v1 and v2; -1.0 if no distance is known */
SCIP_Real reduce_distDataGetSD(
   const DISTDATA*       distdata,           /**< to be initialized */
   int                   vertex1,            /**< first vertex */
   int                   vertex2             /**< second vertex */
)
{
   assert(distdata);
   assert(vertex1 >= 0 && vertex2 >= 0);

   /* try to find SD via Duin's approximation todo */
   // if( distdata->nodeSDpaths_dirty[vertex1]
   // if( distdata->nodeSDpaths_dirty[vertex2] )



   /* neighbors list not valid anymore? */
   if( distdata->nodepaths_dirty[vertex1] )
   {
      /* recompute */

      distdata->nodepaths_dirty[vertex1] = FALSE;
   }

   /* binary search on neighbors of vertex1 */


   /* binary search on neighbors of vertex2? todo too expensive? */

   return -1.0;
}

/** frees members of distance data */
void reduce_distDataFree(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   DISTDATA*             distdata            /**< to be freed */
)
{
   SCIPfreeMemoryArray(scip, &(distdata->nodepaths_dirty));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_range));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_indices));
   SCIPfreeMemoryArray(scip, &(distdata->closenodes_distances));
   SCIPfreeMemoryArray(scip, &(distdata->pathroots_range));
   SCIPfreeMemoryArray(scip, &(distdata->pathroots));
}

