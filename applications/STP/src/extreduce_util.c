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

/**@file   extreduce_util.c
 * @brief  utility methods for Steiner tree extended reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem extended reduction techniques.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include "extreduce.h"
#include "misc_stp.h"
#include "portab.h"

#ifdef RED_UTIL_TIME
#include <time.h>
#endif

#define MLDISTS_EMPTYSLOT_NONE -1
#define MLDISTS_ID_UNSET      -1
#define MLDISTS_DIST_UNSET     -1.0

/** Structure for storing distances in the extension tree.
 *  Organized in slots that can be filled by the user.
 *  On each level there are number of slots available.
 *  Each slots consists of a base (id) and a number of targets (distances, ids).
 *  Each slot on a level has the same number of targets, namely level_ntargets[level].  */
struct multi_level_distances_storage
{
#ifndef NDEBUG
   int*                  target_ids;        /**< target ids only in DEBUG mode! */
#endif
   SCIP_Real*            target_dists;      /**< target ids */
   int*                  base_ids;          /**< bases ids */
   int*                  level_basestart;   /**< start of bases for given level */
   int*                  level_targetstart; /**< start of targets for given level */
   int*                  level_ntargets;    /**< number of targets per base on given level */
   int                   level_maxntargets; /**< maximum number of targets per level */
   int                   level_maxnslots;   /**< maximum number of bases per level */
   int                   nlevels;           /**< number of levels */
   int                   maxnlevels;        /**< maximum number of levels */
   int                   maxntargets;       /**< total maximum number of targets */
   int                   maxnslots;         /**< total maximum number of bases */
   int                   emptyslot_number;  /**< number (0,...) of current empty slot, or EMPTYSLOT_NONE if none exists */
};


/** returns entry of element within sorted array of size arraysize, or -1 if element could not be found */
static inline
int findEntryFromSorted(
   const int*            array,              /**< array */
   int                   arraysize,          /**< size */
   int                   element             /**< element to look for */
   )
{
   int l = 0;
   int u = arraysize - 1;

#ifndef NDEBUG
   assert(u >= 0);

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


/** gets current level */
static inline
int mldistsGetTopLevel(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(mldists->nlevels > 0 && mldists->nlevels <= mldists->maxnlevels);

   return (mldists->nlevels - 1);
}


/** gets start of bases */
static inline
int mldistsGetPosBasesStart(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< the level */
)
{
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_basestart[level] >= 0 && mldists->level_basestart[level] < mldists->maxnslots);

   return mldists->level_basestart[level];
}


/** gets start of bases */
static inline
int mldistsGetPosBasesEnd(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< the level */
)
{
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_basestart[level + 1] >= 0 && mldists->level_basestart[level] <= mldists->maxnslots);

   return mldists->level_basestart[level + 1];
}


/** Gets (internal) position of given base id,
 *  or -1 if it could not be found. */
static
int mldistsGetPosBase(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base id */
)
{
   const int start = mldistsGetPosBasesStart(mldists, level);
   const int end = mldistsGetPosBasesEnd(mldists, level);
   const int* const base_ids = mldists->base_ids;

   assert(baseid >= 0);

   /* top level and with empty slots? */
   if( level == mldistsGetTopLevel(mldists) && extreduce_mldistsEmptySlotExists(mldists) )
   {
      const int emptyslotnumber = mldists->emptyslot_number;

      for( int i = start, j = 0; i < end && j < emptyslotnumber; ++i, ++j )
      {
         const int id = base_ids[i];

         assert(id >= 0);

         if( id == baseid )
            return i;
      }
   }
   else
   {
      for( int i = start; i < end; ++i )
      {
         const int id = base_ids[i];

         assert(id >= 0);

         if( id == baseid )
            return i;
      }
   }

   return -1;
}


/** gets (internal) position of targets for given base id */
static inline
int mldistsGetPosTargetsStart(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base id */
)
{
   int targetpos;
   int offset = 0;
   const int ntargets = mldists->level_ntargets[level];

   assert(extreduce_mldistsLevelContainsBase(mldists, level, baseid));
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(baseid >= 0);

   offset = mldistsGetPosBase(mldists, level, baseid) - mldists->level_basestart[level];

   assert(offset >= 0);
   assert(ntargets > 0);

   targetpos = mldists->level_targetstart[level] + offset * ntargets;

   assert(targetpos >= 0 && targetpos < mldists->maxntargets);

   return targetpos;
}


/** gets targets start of current empty slot */
static inline
int mldistsGetPosEmptyTargetsStart(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int ntargets = mldists->level_ntargets[level];
   const int start = mldists->level_targetstart[level] + mldists->emptyslot_number * ntargets;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(start >= 0 && start < mldists->maxntargets);
   assert(mldists->emptyslot_number >= 0);
   assert(ntargets > 0 && ntargets <= mldists->level_maxntargets);

   return start;
}


/** only for debugging */
static inline
void mldistsTopLevelUnset(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
#ifndef NDEBUG
   const int nlevels = mldists->nlevels;
   const int start_base = mldists->level_basestart[nlevels - 1];
   const int end_base = mldists->level_basestart[nlevels];
   const int start_target = mldists->level_targetstart[nlevels - 1];
   const int end_target = mldists->level_targetstart[nlevels];

   assert(start_base >= 0);
   assert(start_base <= end_base);
   assert(end_base <= mldists->maxnslots);

   for( int i = start_base; i < end_base; ++i )
   {
      mldists->base_ids[i] = MLDISTS_ID_UNSET;
   }

   assert(start_target >= 0);
   assert(start_target <= end_target);
   assert(end_target <= mldists->maxntargets);

   for( int i = start_target; i < end_target; ++i )
   {
      mldists->target_ids[i] = MLDISTS_ID_UNSET;
      mldists->target_dists[i] = MLDISTS_DIST_UNSET;
   }
#endif
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
            SCIP_CALL( distDataPathRootsInsertRoot(scip, g, prededge[k] / 2, startvertex, distdata) );
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

/*
 * Interface methods
 */

/** is the edge valid? */
SCIP_Bool extreduce_edgeIsValid(
   const GRAPH*          graph,              /**< graph data structure */
   int                   e                   /**< edge to be checked */
)
{
   if( EAT_FREE == graph->oeat[e] )
   {
      return FALSE;
   }
   else if( graph_pc_isPcMw(graph) )
   {
      const int tail = graph->tail[e];
      const int head = graph->head[e];

      if( (!graph->mark[tail] || !graph->mark[head]) )
      {
         assert(graph_pc_knotIsDummyTerm(graph, tail) || graph_pc_knotIsDummyTerm(graph, head));

         return FALSE;
      }

      assert(!graph_pc_knotIsDummyTerm(graph, tail));
      assert(!graph_pc_knotIsDummyTerm(graph, head));
   }

   return TRUE;
}


/** initializes distance data */
SCIP_RETCODE extreduce_distDataInit(
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

   assert(extreduce_distCloseNodesAreValid(scip, g, distdata));

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
SCIP_Real extreduce_distDataGetSD(
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
void extreduce_distDataFreeMembers(
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
SCIP_RETCODE extreduce_extPermaInit(
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
   const int msts_datasize = STP_EXT_MAXDFSDEPTH * STP_EXTTREE_MAXNLEAVES_GUARD * 2;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);

   assert(scip && extperm);

   SCIP_CALL( SCIPallocMemoryArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &tree_deg, nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &bottleneckDistNode, nnodes) );

   if( pcmw )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &pcSdToNode, nnodes) );
   }

   SCIP_CALL( reduce_dcmstInit(scip, STP_EXTTREE_MAXNLEAVES_GUARD, &(extperm->dcmst)) );
   SCIP_CALL( graph_csrdepo_init(scip, STP_EXT_MAXDFSDEPTH, msts_datasize, &(extperm->msts)) );
   SCIP_CALL( graph_csrdepo_init(scip, STP_EXT_MAXDFSDEPTH, msts_datasize, &(extperm->msts_reduced)) );

   SCIP_CALL( extreduce_mldistsInit(scip, STP_EXT_MAXDFSDEPTH, STP_EXT_MAXGRAD,
         STP_EXTTREE_MAXNLEAVES_GUARD, &(extperm->sds_vertical)) );
   SCIP_CALL( extreduce_mldistsInit(scip, STP_EXT_MAXDFSDEPTH, STP_EXT_MAXGRAD,
         STP_EXT_MAXGRAD - 1, &(extperm->sds_horizontal)) );

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

   assert(extreduce_extPermaIsClean(graph, extperm));

   return SCIP_OKAY;
}


/** initialize permanent extension data struct */
SCIP_Bool extreduce_extPermaIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperm             /**< extension data */
)
{
   const SCIP_Real* bottleneckDistNode = NULL;
   const SCIP_Real* pcSdToNode = NULL;
   const int* tree_deg = NULL;
   const int nnodes = graph_get_nNodes(graph);

   assert(extperm);

   // todo maybe need to add to extreduce_reddataClean

   bottleneckDistNode = extperm->bottleneckDistNode;
   pcSdToNode = extperm->pcSdToNode;
   tree_deg = extperm->tree_deg;

   if( nnodes != extperm->nnodes )
   {
      SCIPdebugMessage("inconsistent number of nodes \n");
      return FALSE;
   }

   if( !graph_csrdepo_isEmpty(extperm->msts) )
   {
      SCIPdebugMessage("msts not empty! \n");
      return FALSE;
   }

   if( !graph_csrdepo_isEmpty(extperm->msts_reduced) )
   {
      SCIPdebugMessage("msts_reduced not empty! \n");
      return FALSE;
   }

   if( !extreduce_mldistsIsEmpty(extperm->sds_vertical) )
   {
      SCIPdebugMessage("sds_vertical not empty! size=%d \n", extreduce_mldistsNlevels(extperm->sds_vertical));
      return FALSE;
   }

   if( !extreduce_mldistsIsEmpty(extperm->sds_horizontal) )
   {
      SCIPdebugMessage("sds_horizontal not empty! \n");
      return FALSE;
   }

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
void extreduce_extPermaFreeMembers(
   SCIP*                 scip,               /**< SCIP */
   EXTPERMA*             extperm             /**< extension data */
)
{
   assert(scip && extperm);

   extreduce_mldistsFree(scip, &(extperm->sds_horizontal));
   extreduce_mldistsFree(scip, &(extperm->sds_vertical));
   graph_csrdepo_free(scip, &(extperm->msts_reduced));
   graph_csrdepo_free(scip, &(extperm->msts));
   reduce_dcmstFree(scip, &(extperm->dcmst));

   SCIPfreeMemoryArrayNull(scip, &(extperm->pcSdToNode));
   SCIPfreeMemoryArray(scip, &(extperm->bottleneckDistNode));
   SCIPfreeMemoryArray(scip, &(extperm->tree_deg));
   SCIPfreeMemoryArray(scip, &(extperm->isterm));

   extperm->nnodes = -1;
}


/** cleans extension data */
void extreduce_extdataClean(
   EXTDATA*              extdata             /**< extension data */
)
{
   assert(extdata);

   extdata->extstack_ncomponents = 0;
   extdata->tree_nDelUpArcs = 0;
   extdata->tree_nleaves = 0;
   extdata->tree_nedges = 0;
   extdata->tree_depth = 0;
   extdata->tree_root = -1;
   extdata->tree_redcost = 0.0;
}


/** cleans reduction data */
void extreduce_reddataClean(
   REDDATA*              reddata             /**< reduction data */
)
{
   assert(reddata);


   // todo maybe we need to clean up some MSTs here?
}


/** is the extension data clean? */
SCIP_Bool extreduce_extdataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   if( extdata->extstack_ncomponents != 0 )
   {
      printf("extdata->extstack_ncomponents %d \n", extdata->extstack_ncomponents);
      return FALSE;
   }

   if( extdata->tree_nDelUpArcs != 0 )
   {
      printf("extdata->tree_nDelUpArcs %d \n", extdata->tree_nDelUpArcs);
      return FALSE;
   }

   if( !EQ(extdata->tree_redcost, 0.0) )
   {
      printf("extdata->tree_redcost %f \n", extdata->tree_redcost);
      return FALSE;
   }

   if( extdata->tree_nleaves != 0 )
   {
      printf("extdata->tree_nleaves %d \n", extdata->tree_nleaves);
      return FALSE;
   }

   if( extdata->tree_root != -1 )
   {
      printf("extdata->tree_root %d \n", extdata->tree_root);
      return FALSE;
   }

   if( extdata->tree_nedges != 0 )
   {
      printf("extdata->tree_nedges %d \n", extdata->tree_nedges);
      return FALSE;
   }

   if( extdata->tree_depth != 0 )
   {
      printf("extdata->tree_depth %d \n", extdata->tree_depth);
      return FALSE;
   }

   for( int i = 0; i < graph->knots; i++ )
   {
      if( !(extdata->tree_deg[i] == 0 || extdata->tree_deg[i] == -1) )
      {
         printf("extdata->tree_deg[i] %d \n", extdata->tree_deg[i]);
         return FALSE;
      }

      if( !EQ(extdata->tree_bottleneckDistNode[i], -1.0) )
      {
         printf("extdata->bottleneckDistNode[i] %f \n", extdata->tree_bottleneckDistNode[i]);
         return FALSE;
      }
   }

   return TRUE;
}


/** is the reduction data clean? */
SCIP_Bool extreduce_reddataIsClean(
   const GRAPH*          graph,              /**< graph data structure */
   const REDDATA*        reddata             /**< reduction data */
)
{
   assert(graph && reddata);

   for( int i = 0; i < graph->knots; i++ )
   {
      if( reddata->pseudoancestor_mark[i] != 0 )
      {
         printf("pseudoancestor_mark %d \n", reddata->pseudoancestor_mark[i]);
         return FALSE;
      }
   }

   return TRUE;
}


/** initializes multi-level distances structure */
SCIP_RETCODE extreduce_mldistsInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnlevels,         /**< maximum number of levels that can be handled */
   int                   maxnslots,          /**< maximum number of of slots (per level) that can be handled */
   int                   maxntargets,        /**< maximum number of of targets (per slot) that can be handled */
   MLDISTS**             mldistances         /**< to be initialized */
)
{
   MLDISTS* mldists;

   assert(scip && mldistances);
   assert(maxnlevels >= 1 && maxnslots >= 1 && maxntargets >= 1);

   SCIP_CALL( SCIPallocMemory(scip, mldistances) );

   mldists = *mldistances;
   mldists->nlevels = 0;
   mldists->maxnlevels = maxnlevels;
   mldists->level_maxntargets = maxntargets;
   mldists->level_maxnslots = maxnslots;
   mldists->maxnslots = maxnlevels * maxnslots;
   mldists->maxntargets = maxnlevels * maxnslots * maxntargets;
   mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;

#ifndef NDEBUG
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->target_ids), mldists->maxntargets) );
#endif
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->target_dists), mldists->maxntargets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->base_ids), mldists->maxnslots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->level_basestart), maxnlevels + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->level_targetstart), maxnlevels + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->level_ntargets), maxnlevels) );

   mldists->level_basestart[0] = 0;
   mldists->level_targetstart[0] = 0;

   return SCIP_OKAY;
}


/** frees  multi-level distances structure */
void extreduce_mldistsFree(
   SCIP*                 scip,               /**< SCIP */
   MLDISTS**             mldistances         /**< to be freed */
)
{
   MLDISTS* mldists;

   assert(scip && mldistances);

   mldists = *mldistances;

   SCIPfreeMemoryArray(scip, &(mldists->level_ntargets));
   SCIPfreeMemoryArray(scip, &(mldists->level_targetstart));
   SCIPfreeMemoryArray(scip, &(mldists->level_basestart));
   SCIPfreeMemoryArray(scip, &(mldists->base_ids));
   SCIPfreeMemoryArray(scip, &(mldists->target_dists));
#ifndef NDEBUG
   SCIPfreeMemoryArray(scip, &(mldists->target_ids));
#endif


   SCIPfreeMemory(scip, mldistances);
}

/** empty? */
SCIP_Bool extreduce_mldistsIsEmpty(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);

   return (mldists->nlevels == 0);
}


/** does an empty slot exits? (on current level) */
SCIP_Bool extreduce_mldistsEmptySlotExists(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(MLDISTS_EMPTYSLOT_NONE == mldists->emptyslot_number || mldists->emptyslot_number >= 0);

   return (mldists->emptyslot_number != MLDISTS_EMPTYSLOT_NONE);
}


/** gets targets IDs memory from clean slot (to be filled in) */
int* extreduce_mldistsEmptySlotTargetIds(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
#ifndef NDEBUG
   const int start = mldistsGetPosEmptyTargetsStart(mldists);
   const int level = mldistsGetTopLevel(mldists);
   const int end = start + mldists->level_ntargets[level];

   for( int i = start; i < end; ++i )
   {
      assert(mldists->target_ids[i] == MLDISTS_ID_UNSET);
   }

   return &(mldists->target_ids[start]);

#else
   return NULL;
#endif
}


/** gets targets distances memory from clean slot (to be filled in) */
SCIP_Real* extreduce_mldistsEmptySlotTargetDists(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

#ifndef NDEBUG
   const int level = mldistsGetTopLevel(mldists);
   const int end = start + mldists->level_ntargets[level];

   for( int i = start; i < end; ++i )
   {
      assert(EQ(mldists->target_dists[i], MLDISTS_DIST_UNSET));
   }
#endif

   return &(mldists->target_dists[start]);
}


/** gets level of current empty slot */
int extreduce_mldistsEmptySlotLevel(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(extreduce_mldistsEmptySlotExists(mldists));

   return mldistsGetTopLevel(mldists);
}


/** sets base of empty slot */
void extreduce_mldistsEmtpySlotSetBase(
   int                   baseid,             /**< base */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(position >= 0 && position < mldists->maxnslots);
   assert(baseid >= 0);
   assert(mldists->base_ids[position] == MLDISTS_ID_UNSET);

   mldists->base_ids[position] = baseid;
}


/** Resets all changes (especially bases) of empty slot.
 *  NOTE: Assumes that at least the basis of the slot has been set */
void extreduce_mldistsEmtpySlotReset(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(position >= 0 && position < mldists->maxnslots);
   assert(mldists->base_ids[position] != MLDISTS_ID_UNSET);

   mldists->base_ids[position] = MLDISTS_ID_UNSET;

#ifndef NDEBUG
   {
      const int target_start = mldistsGetPosEmptyTargetsStart(mldists);
      const int target_end = target_start + mldists->level_ntargets[level];

      for( int i = target_start; i != target_end; i++ )
      {
         mldists->target_dists[i] = MLDISTS_DIST_UNSET;
         mldists->target_ids[i] = MLDISTS_ID_UNSET;
      }
   }
#endif
}


/** marks current empty slot as filled */
void extreduce_mldistsEmptySlotSetFilled(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->emptyslot_number < extreduce_mldistsLevelNSlots(mldists, level));


#ifndef NDEBUG
   {
      const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

      assert(mldists->base_ids[position] != MLDISTS_ID_UNSET);
   }
#endif


   mldists->emptyslot_number++;

   /* all slots of current level used? */
   if( mldists->emptyslot_number >= extreduce_mldistsLevelNSlots(mldists, level) )
   {
      mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;
   }
}


/** adds another level of slots at top */
void extreduce_mldistsLevelAddTop(
   int                   nslots,             /**< number of slots per this level */
   int                   nslottargets,       /**< number of targets per slot */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   int nlevels;

   assert(!extreduce_mldistsEmptySlotExists(mldists));
   assert(nslots > 0 && nslots <= mldists->level_maxnslots);
   assert(nslottargets > 0 && nslottargets <= mldists->level_maxntargets);
   assert(mldists->nlevels < mldists->maxnlevels);

   mldists->emptyslot_number = 0;
   mldists->nlevels++;

   nlevels = mldists->nlevels;

   assert(nlevels > 0);

   mldists->level_basestart[nlevels] = mldists->level_basestart[nlevels - 1] + nslots;
   mldists->level_targetstart[nlevels] = mldists->level_targetstart[nlevels - 1] + nslots * nslottargets;

   mldists->level_ntargets[nlevels - 1] = nslottargets;

   assert(extreduce_mldistsEmptySlotExists(mldists));

   mldistsTopLevelUnset(mldists);
}


/** closes the top level for further extensions */
void extreduce_mldistsLevelCloseTop(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(!extreduce_mldistsIsEmpty(mldists));

   if( extreduce_mldistsEmptySlotExists(mldists) )
   {
      const int nlevels = mldists->nlevels;
      const int emptyslot = mldists->emptyslot_number;
      int nslottargets;

      assert(nlevels >= 1);
      assert(emptyslot >= 0);

      nslottargets = mldists->level_ntargets[nlevels - 1];

      assert(nslottargets >= 1);

      mldists->level_basestart[nlevels] = mldists->level_basestart[nlevels - 1] + emptyslot;
      mldists->level_targetstart[nlevels] = mldists->level_targetstart[nlevels - 1] + emptyslot * nslottargets;

      mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;
   }

   assert(!extreduce_mldistsEmptySlotExists(mldists));
}


/** removes top level of slots */
void extreduce_mldistsLevelRemoveTop(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(!extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->nlevels > 0);

   mldistsTopLevelUnset(mldists);

   mldists->nlevels--;

   assert(!extreduce_mldistsEmptySlotExists(mldists));
}


/** removes top level of slots, which is not yet closed */
void extreduce_mldistsLevelRemoveTopNonClosed(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->nlevels > 0);

   mldistsTopLevelUnset(mldists);

   mldists->nlevels--;
   mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;

   assert(!extreduce_mldistsEmptySlotExists(mldists));
}


/** gets number of targets (per slots) for given level */
int extreduce_mldistsLevelNTargets(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< level */
)
{
   assert(mldists);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_ntargets[level] >= 0 && mldists->level_ntargets[level] <= mldists->level_maxntargets);

   return mldists->level_ntargets[level];
}

/** gets number of targets (per slots) for top level */
int extreduce_mldistsLevelNTopTargets(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);

   return extreduce_mldistsLevelNTargets(mldists, mldistsGetTopLevel(mldists));
}


/** gets number of slots for given level */
int extreduce_mldistsLevelNSlots(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< level */
)
{
   assert(mldists);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_basestart[level + 1] - mldists->level_basestart[level] >= 0);
   assert(mldists->level_basestart[level + 1] - mldists->level_basestart[level] <= mldists->level_maxnslots);

   return (mldists->level_basestart[level + 1] - mldists->level_basestart[level]);
}


/** gets number of levels */
int extreduce_mldistsNlevels(
   const MLDISTS*        mldists            /**< multi-level distances */
)
{
   assert(mldists);
   assert(mldists->nlevels >= 0);

   return (mldists->nlevels);
}


/** is the base contained in a slot of the given level? */
SCIP_Bool extreduce_mldistsLevelContainsBase(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base id */
)
{
   assert(mldists);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(baseid >= 0);

   return (mldistsGetPosBase(mldists, level, baseid) != -1);
}


/** gets targets ids */
const int* extreduce_mldistsTargetIds(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base */
)
{
#ifndef NDEBUG
   const int targetpos = mldistsGetPosTargetsStart(mldists, level, baseid);

   return &(mldists->target_ids[targetpos]);
#else
   return NULL;
#endif
}


/** gets targets distances */
const SCIP_Real* extreduce_mldistsTargetDists(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base */
)
{
   const int targetpos = mldistsGetPosTargetsStart(mldists, level, baseid);

   return &(mldists->target_dists[targetpos]);
}


/** Gets targets ids */
const int* extreduce_mldistsTopTargetIds(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   baseid              /**< the base */
)
{
   return extreduce_mldistsTargetIds(mldists, mldistsGetTopLevel(mldists), baseid);
}


/** gets targets distances */
const SCIP_Real* extreduce_mldistsTopTargetDists(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   baseid              /**< the base */
)
{
   const int targetpos = mldistsGetPosTargetsStart(mldists, mldistsGetTopLevel(mldists), baseid);

   return &(mldists->target_dists[targetpos]);
}
