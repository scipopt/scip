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

// #define SCIP_DEBUG
#include "extreduce.h"
#include "misc_stp.h"
#include "portab.h"


#ifdef RED_UTIL_TIME
#include <time.h>
#endif

#define MLDISTS_EMPTYSLOT_NONE -1
#define EDGE_FORBIDDEN_NONE -2

/** Structure for storing distances in the extension tree.
 *  Organized in slots that can be filled by the user.
 *  On each level there are a number of slots available (specified by the user).
 *  Each slots consists of a base (id) and a number of targets. Each target has a distance and an ID.
 *  Each slot on a level has the same number of targets, namely level_ntargets[level].  */
struct multi_level_distances_storage
{
   int*                  target_ids;        /**< target ids only in DEBUG mode! */
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
   SCIP_Bool             target_withids;    /**< use ids? */
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
   assert(ntargets >= 0);

   targetpos = mldists->level_targetstart[level] + offset * ntargets;

   assert(targetpos >= 0 && targetpos < mldists->maxntargets);

   return targetpos;
}


/** gets (internal) position of targets for given base id and target id */
static inline
int mldistsGetPosTargets(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid,             /**< the base id */
   int                   targetid            /**< the target id */
)
{
   const int start = mldistsGetPosTargetsStart(mldists, level, baseid);
   const int end = start + mldists->level_ntargets[level];
   const int* const target_ids = mldists->target_ids;
   int i;

   assert(mldists->target_withids);

   for( i = start; i != end; i++ )
   {
      if( target_ids[i] == targetid )
      {
         break;
      }
   }

   assert(i != end);

   return i;
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
   assert(ntargets >= 0 && ntargets <= mldists->level_maxntargets);

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
      mldists->base_ids[i] = STP_MLDISTS_ID_UNSET;
   }

   assert(start_target >= 0);
   assert(start_target <= end_target);
   assert(end_target <= mldists->maxntargets);

   for( int i = start_target; i < end_target; ++i )
   {
      mldists->target_dists[i] = STP_MLDISTS_DIST_UNSET;
   }

   if( mldists->target_withids )
   {
      assert(mldists->target_ids);

      for( int i = start_target; i < end_target; ++i )
      {
         mldists->target_ids[i] = STP_MLDISTS_ID_UNSET;
      }
   }
#endif
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

   assert(size > 0);

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
   const DISTDATA*       distdata,           /**< data */
   int                   startvertex,        /**< start vertex */
   DIJK*                 dijkdata,           /**< limited Dijkstra data */
   int**                 prededge,           /**< predecessors */
   SCIP_Bool**           edgemark            /**< debug only */
)
{
   SCIP_Real* const dist = dijkdata->node_distance;
   DHEAP* const dheap = dijkdata->dheap;
   const int nnodes = g->knots;

   assert(dheap->size == 0 && distdata->closenodes_maxnumber > 0);
   assert(startvertex >= 0 && startvertex < g->knots);
   assert(distdata->closenodes_range[startvertex].start == distdata->closenodes_range[startvertex].end);

   SCIP_CALL( SCIPallocBufferArray(scip, prededge, nnodes) );

#ifndef NDEBUG
   SCIP_CALL( SCIPallocMemoryArray(scip, edgemark, g->edges / 2) );

   for( int e = 0; e < g->edges / 2; e++ )
      (*edgemark)[e] = FALSE;

   for( int k = 0; k < nnodes; k++ )
   {
      (*prededge)[k] = -1;
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
   int                   startvertex,        /**< start vertex */
   SCIP_Bool             is_buildphase,      /**< in build up phase? */
   DIJK*                 dijkdata,           /**< limited Dijkstra data */
   DISTDATA*             distdata,           /**< to be initialized */
   int*                  prededge,           /**< predecessors */
   SCIP_Bool*            edgemark            /**< debug only */
   )
{
   int* const visitlist = dijkdata->visitlist;
   SCIP_Real* const dist = dijkdata->node_distance;
   DHEAP* const dheap = dijkdata->dheap;
   STP_Bool* const visited = dijkdata->node_visited;
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
   const int closenodes_limit = distdata->closenodes_maxnumber;
   const SCIP_Bool isPc = graph_pc_isPc(g);

   assert(dcsr && dist && visitlist && visited && dheap && prededge && edgemark);
   assert(range_closenodes && closenodes_indices && closenodes_prededges);
   assert(dheap->size == 1);
   assert(!isPc || pc_costshifts);

   /* main loop */
   while( dheap->size > 0 )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;
      const SCIP_Real k_dist = dist[k];

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
      }

      /* correct adjacent nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];
         assert(g->mark[m]);

         if( state[m] != CONNECT )
         {
            const SCIP_Real distnew = isPc ?
               k_dist + cost_csr[e] - pc_costshifts[m]
             : k_dist + cost_csr[e];

            assert(GE(distnew, 0.0));

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
   int**                 prededge,           /**< predecessors */
   SCIP_Bool**           edgemark            /**< debug only */
)
{
   SCIPfreeMemoryArrayNull(scip, edgemark);
   SCIPfreeBufferArray(scip, prededge);
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
static inline
SCIP_RETCODE distDataComputeCloseNodesSD(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   startvertex,        /**< start vertex */
   SCIP_Bool             is_buildphase,      /**< in build up phase? */
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
   DIJK*                 dijkdata,           /**< limited Dijkstra data */
   DISTDATA*             distdata            /**< to be initialized */
   )
{
   int* prededge = NULL;
   SCIP_Bool* edgemark = NULL;

   assert(scip && g && dijkdata && distdata);

   SCIP_CALL( closeNodesRunInit(scip, g, distdata, startvertex, dijkdata, &prededge, &edgemark) );

   SCIP_CALL( closeNodesRunCompute( scip, g, startvertex, is_buildphase,
         dijkdata, distdata, prededge, edgemark) );

   closeNodesRunExit(scip, &prededge, &edgemark);

   /* sort close nodes according to their index */
   closeNodesRunSort(g, startvertex, distdata);

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &(distdata->closenodes_prededges), closenodes_totalsize) );

   return SCIP_OKAY;
}



static inline
void distDataRecomputeNormalDist(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   DISTDATA*             distdata            /**< distance data */
)
{
   assert(distdata->pathroot_isdirty[vertex1]);

   distdata->pathroot_nrecomps[vertex1]++;
   SCIP_CALL_ABORT( distDataComputeCloseNodes(scip, g, vertex1, FALSE, distdata->dijkdata, distdata) );

   graph_dijkLimited_reset(g, distdata->dijkdata);

   SCIPdebugMessage("vertex %d is dirty, recompute \n", vertex1);

   distdata->pathroot_isdirty[vertex1] = FALSE;
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



/*
 * Interface methods
 */


/** reverts extension component (makes root the only extension leaf) */
void extreduce_extCompRevert(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   EXTCOMP*              extcomp             /**< component to be cleaned for */
)
{
   int* const extleaves = extcomp->extleaves;
   int* const compedges = extcomp->compedges;
   const int comproot_org = extcomp->comproot;
   const SCIP_Bool compIsSingleEdge = (extcomp->ncompedges == 1);

   assert(extcomp->ncompedges >= 1);
   assert(extcomp->comproot == graph->tail[extcomp->compedges[0]]);

   if( compIsSingleEdge )
   {
      assert(extcomp->nextleaves == 1);
   }
   else
   {
      const int firstedge = compedges[0];

      assert(extcomp->ncompedges >= 3);
      assert(extcomp->nextleaves >= 2);
      assert(graph->head[compedges[1]] == extleaves[0]);

      compedges[0] = compedges[1];
      compedges[1] = flipedge(firstedge);
   }

   compedges[0] = flipedge(compedges[0]);

   assert(extleaves[0] != extcomp->comproot);

   extcomp->comproot = extleaves[0];
   extleaves[0] = comproot_org;
   extcomp->nextleaves = 1;
}


/** is extension component promising candidate for extension? */
SCIP_Bool extreduce_extCompIsPromising(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   const EXTCOMP*        extcomp             /**< component to be cleaned for */)
{
   const int* const extleaves = extcomp->extleaves;
   const SCIP_Bool* const isterm = extperma->isterm;
   const int nextleaves = extcomp->nextleaves;

   assert(extcomp->ncompedges == 1 || extcomp->ncompedges >= 3);

   /* we want to at least check each star! */
   if( extcomp->ncompedges >= 3 )
   {
      return TRUE;
   }

   /* go over all possible extensions and see whether any of them are promising */
   for( int i = 0; i < nextleaves; ++i )
   {
      const int leaf = extleaves[i];

      if( extLeafIsExtendable(graph, isterm, leaf) )
         return TRUE;
   }

   return FALSE;
}


/** is extension component or its reversed version promising candidate for extension? */
SCIP_Bool extreduce_extCompFullIsPromising(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTPERMA*       extperma,           /**< extension data */
   const EXTCOMP*        extcomp             /**< component to be cleaned for */)
{
   assert(graph && extperma && extcomp);

   if( extreduce_extCompIsPromising(graph, extperma, extcomp) )
      return TRUE;

   if( extLeafIsExtendable(graph, extperma->isterm, extcomp->comproot) )
      return TRUE;

   return FALSE;
}


/** initializes distance data */
SCIP_RETCODE extreduce_distDataInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   int                   maxnclosenodes,     /**< maximum number of close nodes to each node */
   SCIP_Bool             computeSD,          /**< also compute special distances? */
   DISTDATA*             distdata            /**< to be initialized */
)
{
   const int nnodes = g->knots;
   RANGE* range_closenodes;
   DIJK* dijkdata;
   DHEAP* dheap;

   assert(distdata && g && scip && g->dcsr_storage);
   assert(graph_valid_dcsr(g, FALSE));
   assert(!graph_pc_isPcMw(g) || !g->extended);

   distDataInitSizes(g, maxnclosenodes, distdata);
   SCIP_CALL( distDataAllocateNodesArrays(scip, g, computeSD, distdata) );

   SCIP_CALL( graph_dijkLimited_init(scip, g, &(distdata->dijkdata)) );
   dijkdata = distdata->dijkdata;

   if( graph_pc_isPc(g) )
   {
      SCIP_CALL( graph_dijkLimited_initPcShifts(scip, g, dijkdata) );
   }

   if( computeSD )
   {
      SCIP_CALL( reduce_sdInit(scip, g, &(distdata->sdistdata)) );
   }
   else
   {
      distdata->sdistdata = NULL;
   }

   range_closenodes = distdata->closenodes_range;

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

#if 0
      if( computeSD )
         SCIP_CALL( distDataComputeCloseNodesSD(scip, g, k, TRUE, dijkdata, distdata) );
      else
#endif
      SCIP_CALL( distDataComputeCloseNodes(scip, g, k, TRUE, dijkdata, distdata) );

      graph_dijkLimited_reset(g, dijkdata);
   }

   /* store for each edge the roots of all paths it is used for */
   SCIP_CALL( distDataPathRootsInitialize(scip, g, distdata) );

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &dheap) );
   distdata->dheap = dheap;

   assert(distdata->closenodes_prededges);
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
   SCIP_Real dist;

   assert(distdata);

   dist = distDataGetNormalDist(scip, g, vertex1, vertex2, distdata);

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
      {
       //  printf("%f->%f \n", dist, dist_sd);

         dist = dist_sd;
      }

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
   SCIPfreeMemoryArrayNull(scip, &(distdata->closenodes_prededges));

   if( distdata->sdistdata )
      reduce_sdFree(scip, &(distdata->sdistdata));

   distDataPathRootsFree(scip, graph, distdata);
   graph_dijkLimited_free(scip, &(distdata->dijkdata));
}


/** initializes multi-level distances structure */
SCIP_RETCODE extreduce_mldistsInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnlevels,         /**< maximum number of levels that can be handled */
   int                   maxnslots,          /**< maximum number of of slots (per level) that can be handled */
   int                   maxntargets,        /**< maximum number of of targets (per slot) that can be handled */
   int                   emptyslot_nbuffers, /**< buffer entries for empty slot (dists and IDs array is that much longer) */
   SCIP_Bool             use_targetids,      /**< use target IDs? */
   MLDISTS**             mldistances         /**< to be initialized */
)
{
   MLDISTS* mldists;

   assert(scip && mldistances);
   assert(maxnlevels >= 1 && maxnslots >= 1 && maxntargets >= 1);
   assert(emptyslot_nbuffers >= 0);

   SCIP_CALL( SCIPallocMemory(scip, mldistances) );

   mldists = *mldistances;
   mldists->nlevels = 0;
   mldists->maxnlevels = maxnlevels;
   mldists->level_maxntargets = maxntargets;
   mldists->level_maxnslots = maxnslots;
   mldists->target_withids = use_targetids;
   mldists->maxnslots = maxnlevels * maxnslots;
   mldists->maxntargets = maxnlevels * maxnslots * maxntargets + emptyslot_nbuffers;
   mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;

   if( use_targetids )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->target_ids), mldists->maxntargets) );
   }
   else
   {
      mldists->target_ids = NULL;
   }

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

   if( mldists->target_withids )
   {
      assert(mldists->target_ids);

      SCIPfreeMemoryArray(scip, &(mldists->target_ids));
   }

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
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

#ifndef NDEBUG
   const int level = mldistsGetTopLevel(mldists);
   const int end = start + mldists->level_ntargets[level];

   assert(mldists->target_withids);
   assert(mldists->target_ids);

   for( int i = start; i < end; ++i )
   {
      assert(mldists->target_ids[i] == STP_MLDISTS_ID_UNSET);
   }
#endif

   return &(mldists->target_ids[start]);
}

/** gets targets IDs memory from clean slot (to be filled in) */
int* extreduce_mldistsEmptySlotTargetIdsDirty(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

   assert(mldists->target_ids);

   return &(mldists->target_ids[start]);
}


/** Gets targets distances memory from clean slot (to be filled in).
 *  NOTE: Can only be called as long as this empty slots' distances have not not modified! */
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
      assert(EQ(mldists->target_dists[i], STP_MLDISTS_DIST_UNSET));
   }
#endif

   return &(mldists->target_dists[start]);
}


/** Gets targets distances memory from empty slot.
 *  NOTE: This method does not make sure that the distances are clean! (i.e. not already set) */
SCIP_Real* extreduce_mldistsEmptySlotTargetDistsDirty(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

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
void extreduce_mldistsEmptySlotSetBase(
   int                   baseid,             /**< base */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(position >= 0 && position < mldists->maxnslots);
   assert(baseid >= 0);
   assert(mldists->base_ids[position] == STP_MLDISTS_ID_UNSET);

   mldists->base_ids[position] = baseid;
}


/** Resets all changes (especially bases) of empty slot.
 *  NOTE: Assumes that at least the basis of the slot has been set */
void extreduce_mldistsEmptySlotReset(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(position >= 0 && position < mldists->maxnslots);
   assert(mldists->base_ids[position] != STP_MLDISTS_ID_UNSET);

   mldists->base_ids[position] = STP_MLDISTS_ID_UNSET;

#ifndef NDEBUG
   {
      const int target_start = mldistsGetPosEmptyTargetsStart(mldists);
      const int target_end = target_start + mldists->level_ntargets[level];

      for( int i = target_start; i != target_end; i++ )
      {
         mldists->target_dists[i] = STP_MLDISTS_DIST_UNSET;
      }

      if( mldists->target_withids )
      {
         assert(mldists->target_ids);

         for( int i = target_start; i != target_end; i++ )
         {
            mldists->target_ids[i] = STP_MLDISTS_ID_UNSET;
         }
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

      assert(mldists->base_ids[position] != STP_MLDISTS_ID_UNSET);
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
   assert(nslottargets >= 0 && nslottargets <= mldists->level_maxntargets);
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


/** gets number of slots for top level */
int extreduce_mldistsTopLevelNSlots(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   return (extreduce_mldistsLevelNSlots(mldists, extreduce_mldistsTopLevel(mldists)));
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


/** gets top level */
int extreduce_mldistsTopLevel(
   const MLDISTS*        mldists            /**< multi-level distances */
)
{
   assert(mldists);
   assert(mldists->nlevels >= 1);

   return (mldists->nlevels - 1);
}


/** gets top level bases*/
const int* extreduce_mldistsTopLevelBases(
   const MLDISTS*        mldists            /**< multi-level distances */
)
{
   const int toplevel = mldistsGetTopLevel(mldists);
   const int basestart_pos = mldistsGetPosBasesStart(mldists, toplevel);

   return &(mldists->base_ids[basestart_pos]);
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
   const int targetpos = mldistsGetPosTargetsStart(mldists, level, baseid);

   assert(mldists->target_withids);
   assert(mldists->target_ids);

   return &(mldists->target_ids[targetpos]);
}


/** gets targets distances */
const SCIP_Real* extreduce_mldistsTargetDists(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base */
)
{
   int targetpos;

   assert(mldists);
   assert(level >= 0 && baseid >= 0);

   targetpos = mldistsGetPosTargetsStart(mldists, level, baseid);

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


/** gets (one!) target distance for given target ID and base ID */
SCIP_Real extreduce_mldistsTopTargetDist(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   baseid,             /**< the base */
   int                   targetid            /**< the identifier */
)
{
   int targetpos;

   assert(mldists);
   assert(baseid >= 0 && targetid >= 0);

   targetpos = mldistsGetPosTargets(mldists, mldistsGetTopLevel(mldists), baseid, targetid);

   assert(GE(mldists->target_dists[targetpos], 0.0));

   return (mldists->target_dists[targetpos]);
}

/** gets (one!) target distance for given level, target ID and base ID */
SCIP_Real extreduce_mldistsTargetDist(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid,             /**< the base */
   int                   targetid            /**< the identifier */
)
{
   int targetpos;

   assert(mldists);
   assert(baseid >= 0 && targetid >= 0);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));

   targetpos = mldistsGetPosTargets(mldists, level, baseid, targetid);

   assert(GE(mldists->target_dists[targetpos], 0.0));

   return (mldists->target_dists[targetpos]);
}
