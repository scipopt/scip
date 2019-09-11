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

/**@file   graph_ancestors.c
 * @brief  includes graph ancestor methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * Graph ancestor methods for Steiner problems. Usually used to reconstruct a graph after reductions
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h) */

#include <assert.h>
#include "graph.h"
#include "portab.h"


/** fixed graph components, typedef FIXED */
struct fixed_graph_components
{
   IDX*                  fixedges;           /**< fixed edges */
   int*                  fixpseudonodes;     /**< fixed psuedo eliminated nodes  */
   int                   nfixnodes;          /**< number of fixed nodes  */
};

/** blocked pseudo ancestors */
typedef struct blocked_pseudo_ancestor
{
   int**                 blocks;             /**< blocks of ancestors (of size  nblocks) */
   int*                  sizes;              /**< current number of ancestors (of size  nblocks)  */
   int*                  capacities;         /**< current capacities of ancestors (of size nblocks) */
   int                   nblocks;            /**< number of ancestor blocks */
} BLOCKANS;

/** node ancestors resulting from pseudo-elimination, typedef PSEUDOANS */
struct pseudo_ancestors
{
   BLOCKANS*             ans_halfedges;      /**< (half) edge ancestors  */
   BLOCKANS*             ans_nodes;          /**< (pc/mw) node ancestors  */
   int                   halfnedges;         /**< half number of edges */
   int                   nnodes;             /**< number of nodes */
};


/** get next power of 2 number; assumes number >= 1 and number < UINT_MAX */
static inline
int getNextPow2(int number)
{
   uint32_t n = (uint32_t) number;
   assert(number >= 1);

   n--;
   n |= n >> 1;
   n |= n >> 2;
   n |= n >> 4;
   n |= n >> 8;
   n |= n >> 16;
   n++;

   assert(n < INT_MAX );
   assert(number >= (int) n && number <= 2 * (int) n);
   assert(number >= 2);

   return (int) n;
}


/** initialize */
static
SCIP_RETCODE blockedAncestors_init(
   int                   nblocks,            /**< number of ancestor blocks */
   BLOCKANS**            blockedans          /**< blocked pseudo-ancestors */
)
{
   int** blocks;
   int* sizes;
   int* capacities;

   assert(nblocks > 0);

   SCIP_CALL( SCIPallocMemory(scip, blockedans) );

   (*blockedans)->nblocks = nblocks;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(blocks), nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sizes), nblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(capacities), nblocks) );

   for( int k = 0; k < nblocks; k++ )
   {
      blocks[k] = NULL;
      sizes[k] = 0;
      capacities[k] = 0;
   }

   (*blockedans)->blocks = blocks;
   (*blockedans)->sizes = sizes;
   (*blockedans)->capacities = capacities;

   return SCIP_OKAY;
}


/** free */
static inline
void blockedAncestors_freeBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   block_id,           /**< id for which to free pseudo ancestors */
   BLOCKANS*             blockedans         /**< blocked pseudo-ancestors */
)
{
   const int capacity = blockedans->capacities[block_id];

   assert(scip && blockedans && blockedans->sizes && blockedans->capacities);
   assert(block_id >= 0 && block_id < blockedans->nblocks);
   assert(capacity >= 0);

   if( capacity > 0 )
   {
      assert(blockedans->blocks[block_id]);
      SCIPfreeBlockMemoryArray(scip, &(blockedans->blocks[block_id]), capacity);

      blockedans->sizes[block_id] = 0;
      blockedans->capacities[block_id] = 0;
   }

   assert(blockedans->sizes[block_id] == 0);
   assert(blockedans->capacities[block_id] == 0);
   assert(blockedans->blocks[block_id] == NULL);
}


/** free */
static
void blockedAncestors_free(
   SCIP*                 scip,               /**< SCIP data structure */
   BLOCKANS**            blockedans          /**< blocked pseudo-ancestors */
)
{
   BLOCKANS* const ancestors = (*blockedans);
   int* const * blocks = ancestors->blocks;
   const int* sizes = ancestors->sizes;
   const int* capacities = ancestors->capacities;
   const int nblocks = ancestors->nblocks;

   assert(blockedans && ancestors);
   assert(nblocks >= 1);
   assert(capacities && sizes && blocks);

   for( int e = nblocks - 1; e >= 0; --e )
      blockedAncestors_freeBlock(scip, e, ancestors);

   SCIPfreeMemoryArray(scip, &sizes);
   SCIPfreeMemoryArray(scip, &capacities);
   SCIPfreeMemoryArray(scip, &blocks);

   SCIPfreeMemoryArray(scip, blockedans);
}


/** hash ancestors of given entry */
static inline
void blockedAncestors_hash(
   const BLOCKANS*       blockedans,         /**< blocked pseudo-ancestors */
   int                   block_id,           /**< entry for which to reallocate */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int nancestors = blockedans->sizes[block_id];
   const int* const ancestorsblock = blockedans->blocks[block_id];

   assert(blockedans && hasharr);
   assert(block_id >= 0 && block_id < blockedans->nblocks);
   assert(nancestors >= 0);

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestorsblock[k];
      assert(a >= 0 && hasharr[a] == 0);

      hasharr[a] = 1;
   }
}


/** unhash nAncestorsToClean many ancestors  */
static inline
void blockedAncestors_unhashPartial(
   const BLOCKANS*       blockedans,         /**< blocked pseudo-ancestors */
   int                   block_id,           /**< entry for which to reallocate */
   int                   nAncestorsToClean,  /**< number of (first) ancestors to use for clean, or -1 to clean all */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int* const ancestorsblock = blockedans->blocks[block_id];
   const int n = (nAncestorsToClean >= 0) ? nAncestorsToClean : blockedans->sizes[block_id];

   assert(blockedans && hasharr);
   assert(block_id >= 0 && block_id < blockedans->nblocks);
   assert(nAncestorsToClean >= 0 || nAncestorsToClean == -1);

   for( int k = 0; k < n; k++ )
   {
      const int a = ancestorsblock[k];
      assert(a >= 0);
      assert(hasharr[a] == 1);

      hasharr[a] = 0;
   }
}


/** hash ancestors of given entry with possible conflicts */
static inline
void blockedAncestors_hashDirty(
   const BLOCKANS*       blockedans,         /**< blocked pseudo-ancestors */
   int                   block_id,           /**< entry for which to reallocate */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   SCIP_Bool*            conflict,           /**< conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int nancestors = blockedans->sizes[block_id];
   const int* const ancestorsblock = blockedans->blocks[block_id];

   assert(blockedans && hasharr);
   assert(block_id >= 0 && block_id < blockedans->nblocks);
   assert(nancestors >= 0);

   *conflict = FALSE;

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestorsblock[k];
      assert(a >= 0 );
      assert(hasharr[a] == 0 || hasharr[a] == 1);

      if( hasharr[a] == 1 )
      {
         *conflict = TRUE;

         if( breakOnConflict )
            return;
      }

      hasharr[a] = 1;
   }
}


/** unhash ancestors of given entry with possible conflicts */
static inline
void blockedAncestors_unhashDirty(
   const BLOCKANS*       blockedans,         /**< blocked pseudo-ancestors */
   int                   block_id,           /**< entry for which to reallocate */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int nancestors = blockedans->sizes[block_id];
   const int* const ancestorsblock = blockedans->blocks[block_id];

   assert(blockedans && hasharr);
   assert(block_id >= 0 && block_id < blockedans->nblocks);
   assert(nancestors >= 0);

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestorsblock[k];
      assert(a >= 0 );
      assert(hasharr[a] == 0 || hasharr[a] == 1);

      if( breakOnConflict && hasharr[a] == 0 )
         return;

      hasharr[a] = 0;
   }
}


/** ancestor already hashed? */
static inline
SCIP_Bool blockedAncestors_hashIsHit(
   const BLOCKANS*       blockedans,         /**< blocked pseudo-ancestors */
   int                   ancestor,           /**< ancestor to check */
   const int*            hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   assert(blockedans && hasharr);
   assert(ancestor >= 0 && ancestor < blockedans->nblocks);
   assert(hasharr[ancestor] == 0 || hasharr[ancestor] == 1);

   return (hasharr[ancestor] == 1);
}


/** reallocate ancestor array of given entry */
static inline
SCIP_RETCODE blockedAncestors_realloc(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   block_id,           /**< entry for which to reallocate */
   int                   maxsize_new,        /**< the new size */
   BLOCKANS*             blockedans          /**< blocked pseudo-ancestors */
)
{
   const int maxsize = blockedans->capacities[block_id];

   assert(scip && blockedans);
   assert(block_id >= 0 && block_id < blockedans->nblocks);
   assert(maxsize >= 0 && maxsize_new > maxsize);
   assert(maxsize_new >= 2);

   if( maxsize == 0 )
   {
      assert(blockedans->blocks[block_id] == NULL);
      assert(blockedans->sizes[block_id] == 0);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(blockedans->blocks[block_id]), maxsize_new) );
   }
   else
   {
      assert(blockedans->blocks[block_id] != NULL);
      assert(maxsize >= 2);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(blockedans->blocks[block_id]), maxsize, maxsize_new) );
   }

   blockedans->capacities[block_id] = maxsize_new;

   return SCIP_OKAY;
}


/** pseudo ancestors of given block without conflicts? */
static inline
SCIP_Bool blockedAncestors_blockIsValid(
   int                   block_id,           /**< entry for which to reallocate */
   const BLOCKANS*       blockedans,         /**< blocked pseudo-ancestors */
   int*                  hasharr             /**< clean hash array of size nnodes (wrt pseudo ancestors) */
)
{
   SCIP_Bool conflict;

   assert(blockedans);
   assert(block_id >= 0 && block_id < blockedans->nblocks);

   blockedAncestors_hashDirty(blockedans, block_id, FALSE, &conflict, hasharr);
   blockedAncestors_unhashDirty(blockedans, block_id, FALSE, hasharr);

   return !conflict;
}


/** appends copy of pseudo ancestors of source to target */
static
SCIP_RETCODE blockedAncestors_appendCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   block_target,       /**< target block */
   int                   block_source,       /**< source block */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   int                   nnodes,             /**< number of nodes of underlying graph */
   BLOCKANS*             blockedans,         /**< blocked pseudo-ancestors */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   int position_target;
   const int size_source = blockedans->sizes[block_source];

   assert(scip && blockedans && conflict);
   assert(nnodes >= 0 && size_source >= 0);
   assert(block_target >= 0 && block_target < blockedans->nblocks);
   assert(block_source >= 0 && block_source < blockedans->nblocks);

   *conflict = FALSE;

   /* anything to append? */
   if( size_source > 0 )
   {
      const int* const ancestors_source = blockedans->blocks[block_source];
      int* ancestors_target;
      int* hasharr;
      const int size_target_org = blockedans->sizes[block_target];
      const int size_targetPlusSource = size_target_org + size_source;

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, nnodes) );

      /* need to realloc target ancestor array? */
      if( size_targetPlusSource > blockedans->capacities[block_target] )
      {
         const int maxsize_up = getNextPow2(size_targetPlusSource);
         SCIP_CALL( blockedAncestors_realloc(scip, block_target, maxsize_up, blockedans) );
      }

      ancestors_target = blockedans->blocks[block_target];

      /* mark ancestors of target */
      blockedAncestors_hash(blockedans, block_target, hasharr);

      position_target = size_target_org;

      /* add source to target */
      for( int e = 0; e < size_source; e++ )
      {
         const int ancestor = ancestors_source[e];

         if( blockedAncestors_hashIsHit(blockedans, ancestor, hasharr) )
         {
            *conflict = TRUE;

            if( breakOnConflict )
               break;

            continue;
         }

         assert(position_target < blockedans->capacities[block_target]);

         ancestors_target[position_target++] = ancestor;
      }

      assert(position_target <= size_targetPlusSource);
      assert(position_target >= blockedans->sizes[block_target]);

      blockedans->sizes[block_target] = position_target;

      /* unmark ancestors of target */
      blockedAncestors_unhashPartial(blockedans, block_target, size_target_org, hasharr);

      SCIPfreeCleanBufferArray(scip, &hasharr);
   }

   return SCIP_OKAY;
}


/** adds pseudo ancestor to ancestor list of given block */
static
SCIP_RETCODE blockedAncestors_addAncestor(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   block_target,       /**< block for which to add pseudo ancestor */
   int                   ancestor,           /**< (index of) pseudo ancestor */
   int                   nnodes,             /**< number of nodes of underlying graph */
   BLOCKANS*             blockedans          /**< blocked pseudo-ancestors */
)
{
   int* const sizes = blockedans->sizes;
   int* const capacities = blockedans->capacities;

   assert(scip && blockedans);
   assert(block_target >= 0 && block_target < blockedans->nblocks);

   /* need to reallocate? */
   if( sizes[block_target] == capacities[block_target] )
   {
      const int maxsize_up = getNextPow2(sizes[block_target] + 1);
      SCIP_CALL( blockedAncestors_realloc(scip, block_target, maxsize_up, blockedans) );
   }

   assert(sizes[block_target] < capacities[block_target]);

   blockedans->blocks[block_target][sizes[block_target]++] = ancestor;

#ifndef NDEBUG
   {
      int* hasharr;
      SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &hasharr, nnodes) );
      assert(blockedAncestors_blockIsValid(block_target, blockedans, hasharr));
      SCIPfreeCleanBufferArray(scip, &hasharr);
   }
#endif

   return SCIP_OKAY;
}


/** valid pseudo-ancestors (no conflicts)? */
static
SCIP_Bool blockedAncestors_isValid(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nnodes,             /**< number of nodes of underlying graph */
   const BLOCKANS*       blockedans          /**< blocked pseudo-ancestors */
)
{
   int* hasharr;
   SCIP_Bool isValid = TRUE;

   assert(scip && scip && blockedans);
   assert(nnodes >= 0);

   /* check whether sizes/capacities are correct */
   for( int e = 0; e < blockedans->nblocks && isValid; e++ )
   {
      if( blockedans->sizes[e] > blockedans->capacities[e] )
         isValid = FALSE;

      if( blockedans->sizes[e] > 0 && blockedans->blocks[e] == NULL )
         isValid = FALSE;

      if( blockedans->capacities[e] > 0 && blockedans->blocks[e] == NULL )
         isValid = FALSE;
   }

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &hasharr, nnodes) );

   /* check whether there are conflict within the ancestor blocks */
   for( int e = 0; e < blockedans->nblocks && isValid; e++ )
   {
      if( blockedans->blocks[e] != NULL && !blockedAncestors_blockIsValid(e, blockedans, hasharr) )
         isValid = FALSE;
   }

   SCIPfreeCleanBufferArray(scip, &hasharr);

   return isValid;
}

/** hash ancestors of given edge */
void graph_pseudoAncestors_hashEdge(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to hash */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   blockedAncestors_hash(pseudoancestors->ans_halfedges, halfedge, hasharr);
}

/** hash ancestors of given node */
void graph_pseudoAncestors_hashNode(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   node,               /**< node for which to hash */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   assert(pseudoancestors && hasharr);
   assert(node >= 0 && node < pseudoancestors->nnodes);

   blockedAncestors_hash(pseudoancestors->ans_nodes, node, hasharr);
}


/** unhash ancestors of given edge */
void graph_pseudoAncestors_unhashEdge(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to unhash */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   blockedAncestors_unhashPartial(pseudoancestors->ans_halfedges, halfedge, -1, hasharr);
}


/** hash ancestors of given node */
void graph_pseudoAncestors_unhashNode(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   node,               /**< node for which to unhash */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   assert(pseudoancestors && hasharr);
   assert(node >= 0 && node < pseudoancestors->nnodes);

   blockedAncestors_unhashPartial(pseudoancestors->ans_nodes, node, -1, hasharr);
}


/** hash ancestors of given edge with possible conflicts */
void graph_pseudoAncestors_hashEdgeDirty(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to hash */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   SCIP_Bool*            conflict,           /**< conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   blockedAncestors_hashDirty(pseudoancestors->ans_halfedges, halfedge, breakOnConflict, conflict, hasharr);
}


/** hash ancestors of given node with possible conflicts */
void graph_pseudoAncestors_hashNodeDirty(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   node,               /**< node for which to hash */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   SCIP_Bool*            conflict,           /**< conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   assert(pseudoancestors && hasharr);
   assert(node >= 0 && node < pseudoancestors->nnodes);

   blockedAncestors_hashDirty(pseudoancestors->ans_nodes, node, breakOnConflict, conflict, hasharr);
}


/** unhash ancestors of given edge with possible conflicts */
void graph_pseudoAncestors_unhashEdgeDirty(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to unhash */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   blockedAncestors_unhashDirty(pseudoancestors->ans_halfedges, halfedge, breakOnConflict, hasharr);
}


/** hash ancestors of given node with possible conflicts */
void graph_pseudoAncestors_unhashNodeDirty(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   node,               /**< node for which to unhash */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   assert(pseudoancestors && hasharr);
   assert(node >= 0 && node < pseudoancestors->nnodes);

   blockedAncestors_unhashDirty(pseudoancestors->ans_nodes, node, breakOnConflict, hasharr);
}


/** initializes pseudo ancestors */
SCIP_RETCODE graph_init_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
)
{
   PSEUDOANS* pseudoancestors;

   assert(scip && g);
   assert(g->pseudoancestors == NULL);

   SCIP_CALL( SCIPallocMemory(scip, &pseudoancestors) );

   g->pseudoancestors = pseudoancestors;

   pseudoancestors->nnodes = g->knots;
   pseudoancestors->halfnedges = g->edges / 2;

   blockedAncestors_init(pseudoancestors->halfnedges, &(pseudoancestors->ans_halfedges));

   if( graph_pc_isPcMw(g) )
   {
      blockedAncestors_init(g->knots, &(pseudoancestors->ans_nodes));
   }
   else
   {
      pseudoancestors->ans_nodes = NULL;
   }

   return SCIP_OKAY;
}


/** frees pseudo ancestors */
void graph_free_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   PSEUDOANS* pseudoancestors;

   assert(scip && g && g->pseudoancestors);
   assert(g->pseudoancestors->nnodes >= 1);

   pseudoancestors = g->pseudoancestors;

   blockedAncestors_free(scip, &(pseudoancestors->ans_halfedges));
   blockedAncestors_free(scip, &(pseudoancestors->ans_nodes));

   assert(pseudoancestors->ans_halfedges == NULL);
   assert(pseudoancestors->ans_nodes == NULL);

   SCIPfreeMemoryArray(scip, &(g->pseudoancestors));

   assert(g->pseudoancestors == NULL);
}


/** frees pseudo ancestor block for given edge */
void graph_free_pseudoAncestorsEdgeBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_free,          /**< edge for which to free pseudo ancestors */
   GRAPH*                g                   /**< the graph */
)
{
   const int block_id = edge_free / 2;

   assert(scip && g && g->pseudoancestors && g->pseudoancestors->ans_halfedges);
   assert(block_id >= 0 && block_id < g->pseudoancestors->halfnedges);

   blockedAncestors_freeBlock(scip, block_id, g->pseudoancestors->ans_halfedges);
}


/** frees pseudo ancestor block for given node */
void graph_free_pseudoAncestorsNodeBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node_free,          /**< node for which to free pseudo ancestors */
   GRAPH*                g                   /**< the graph */
)
{
   assert(scip && g && g->pseudoancestors && g->pseudoancestors->ans_nodes);
   assert(graph_pc_isPcMw(g));
   assert(node_free >= 0 && node_free < g->pseudoancestors->nnodes);

   blockedAncestors_freeBlock(scip, node_free, g->pseudoancestors->ans_nodes);
}


/** returns number of pseudo ancestors for given edge */
int graph_get_nPseudoAncestorsEdge(
   const GRAPH*          g,            /**< the graph */
   int                   edge          /**< edge for which to return number of pseudo ancestors */
   )
{
   const int halfedge = edge / 2;

   assert(g && g->pseudoancestors && g->pseudoancestors->ans_halfedges);
   assert(halfedge >= 0 && halfedge < g->pseudoancestors->halfnedges);

   return g->pseudoancestors->ans_halfedges->sizes[halfedge];
}


/** returns number of pseudo ancestors for given node */
int graph_get_nPseudoAncestorsNode(
   const GRAPH*          g,            /**< the graph */
   int                   node          /**< node for which to return number of pseudo ancestors */
   )
{
   assert(g && g->pseudoancestors && g->pseudoancestors->ans_nodes);
   assert(node >= 0 && node < g->pseudoancestors->nnodes);

   return g->pseudoancestors->ans_nodes->sizes[node];
}


/** returns pseudo ancestors for given edge */
const int* graph_get_pseudoAncestorsEdge(
   const GRAPH*          g,            /**< the graph */
   int                   edge          /**< edge for which to return pseudo ancestors */
   )
{
   const int halfedge = edge / 2;

   assert(g && g->pseudoancestors && g->pseudoancestors->ans_halfedges);
   assert(halfedge >= 0 && halfedge < g->pseudoancestors->halfnedges);

   return g->pseudoancestors->ans_halfedges->blocks[halfedge];
}


/** returns pseudo ancestors for given node */
const int* graph_get_pseudoAncestorsNode(
   const GRAPH*          g,            /**< the graph */
   int                   node          /**< node for which to return pseudo ancestors */
   )
{
   assert(g && g->pseudoancestors && g->pseudoancestors->ans_nodes);
   assert(node >= 0 && node < g->pseudoancestors->nnodes);

   return g->pseudoancestors->ans_nodes->blocks[node];
}


/** returns array number of nodes */
int graph_pseudoAncestors_getNnodes(
   const GRAPH*          g            /**< the graph */
   )
{
   assert(g && g->pseudoancestors);

   return g->pseudoancestors->nnodes;
}


/** appends copy of pseudo ancestors of edge_source to edge_target */
SCIP_RETCODE graph_pseudoAncestors_appendCopyEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_target,        /**< edge target */
   int                   edge_source,        /**< edge source */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   PSEUDOANS* const pseudoancestors = g->pseudoancestors;
   const int target = edge_target / 2;
   const int source = edge_source / 2;

   assert(scip && g && pseudoancestors && conflict);

   SCIP_CALL( blockedAncestors_appendCopy(scip, target, source, breakOnConflict,
         pseudoancestors->nnodes, pseudoancestors->ans_halfedges, conflict) );

   return SCIP_OKAY;
}


/** appends copy of pseudo ancestors of node source to node target */
SCIP_RETCODE graph_pseudoAncestors_appendCopyNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node_target,        /**< node target */
   int                   node_source,        /**< node source */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   PSEUDOANS* const pseudoancestors = g->pseudoancestors;

   assert(scip && g && pseudoancestors && conflict);

   SCIP_CALL( blockedAncestors_appendCopy(scip, node_target, node_source, breakOnConflict,
         pseudoancestors->nnodes, pseudoancestors->ans_nodes, conflict) );

   return SCIP_OKAY;
}


/** appends pseudo ancestors of edge_source to edge_target, ancestors for edge_source are deleted */
SCIP_RETCODE graph_pseudoAncestors_appendMoveEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_target,        /**< target edge */
   int                   edge_source,        /**< source edge  */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   SCIP_CALL( graph_pseudoAncestors_appendCopyEdge(scip, edge_target, edge_source, breakOnConflict, g, conflict) );
   graph_free_pseudoAncestorsEdgeBlock(scip, edge_source, g);

   return SCIP_OKAY;
}


/** appends pseudo ancestors of node source to node target, ancestors for source are deleted */
SCIP_RETCODE graph_pseudoAncestors_appendMoveNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node_target,        /**< target node */
   int                   node_source,        /**< source node */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   SCIP_CALL( graph_pseudoAncestors_appendCopyNode(scip, node_target, node_source, breakOnConflict, g, conflict) );
   graph_free_pseudoAncestorsNodeBlock(scip, node_source, g);

   return SCIP_OKAY;
}


/** adds pseudo ancestor to ancestor list of given edge */
SCIP_RETCODE graph_pseudoAncestors_addToEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_target,        /**< edge for which to add pseudo ancestor */
   int                   ancestor,           /**< (index of) pseudo ancestor */
   GRAPH*                g                   /**< the graph */
)
{
   PSEUDOANS* const pseudoancestors = g->pseudoancestors;
   const int target = edge_target / 2;

   assert(scip && g && pseudoancestors);
   assert(target >= 0 && target < g->pseudoancestors->halfnedges);
   assert(ancestor >= 0 && ancestor < g->pseudoancestors->nnodes);

   SCIP_CALL( blockedAncestors_addAncestor(scip, target, ancestor, pseudoancestors->nnodes, pseudoancestors->ans_halfedges) );

   return SCIP_OKAY;
}


/** adds pseudo ancestor to ancestor list of given node */
SCIP_RETCODE graph_pseudoAncestors_addToNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node_target,        /**< node for which to add pseudo ancestor */
   int                   ancestor,           /**< (index of) pseudo ancestor */
   GRAPH*                g                   /**< the graph */
)
{
   PSEUDOANS* const pseudoancestors = g->pseudoancestors;

   assert(scip && g && pseudoancestors);
   assert(node_target >= 0 && node_target < g->pseudoancestors->nnodes);
   assert(ancestor >= 0 && ancestor < g->pseudoancestors->nnodes);

   SCIP_CALL( blockedAncestors_addAncestor(scip, node_target, ancestor, pseudoancestors->nnodes, pseudoancestors->ans_nodes) );

   return SCIP_OKAY;
}


/** valid pseudo-ancestors (no conflicts)? */
SCIP_Bool graph_valid_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   const PSEUDOANS* const pseudoancestors = g->pseudoancestors;
   const int nnodes = pseudoancestors->nnodes;
   SCIP_Bool isValid;

   assert(scip && g && pseudoancestors);

   isValid = blockedAncestors_isValid(scip, nnodes, pseudoancestors->ans_halfedges);

   if( !isValid )
      return FALSE;

   if( graph_pc_isPcMw(g) )
   {
      isValid = blockedAncestors_isValid(scip, nnodes, pseudoancestors->ans_nodes);

      if( !isValid )
         return FALSE;
   }

   return TRUE;
}

#if 0
/** check whether conflict for one edge */
SCIP_RETCODE graph_checkConflict1_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< the graph */
   int                   edge1,              /**< first edge */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   const int block_id1 = edge1 / 2;


   assert(scip && g && g->pseudoancestors);
   assert(block_id1 >= 0 && block_id1 < g->pseudoancestors->halfnedges);


   return SCIP_OKAY;
}
#endif


/** initializes fixed */
SCIP_RETCODE graph_init_fixed(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
)
{
   FIXED* fixedcomponents;

   assert(scip && g);
   assert(g->fixedcomponents == NULL);

   SCIP_CALL( SCIPallocMemory(scip, &(g->fixedcomponents)) );

   fixedcomponents = g->fixedcomponents;

   fixedcomponents->fixedges = NULL;
   fixedcomponents->fixpseudonodes = NULL;
   fixedcomponents->nfixnodes = 0;

   return SCIP_OKAY;
}


/** frees fixed */
void graph_free_fixed(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   IDX* curr;
   FIXED* fixedcomponents = g->fixedcomponents;

   assert(scip && g);
   assert(fixedcomponents != NULL);
   assert(fixedcomponents->nfixnodes >= 0);

   if( fixedcomponents->nfixnodes > 0 )
      SCIPfreeBlockMemoryArray(scip, &(fixedcomponents->fixpseudonodes), fixedcomponents->nfixnodes);

   curr = fixedcomponents->fixedges;
   while( curr != NULL )
   {
      fixedcomponents->fixedges = curr->parent;
      SCIPfreeBlockMemory(scip, &(curr));

      curr = fixedcomponents->fixedges;
   }

   SCIPfreeMemoryArray(scip, &(g->fixedcomponents));
}


/** adds new fixed components */
SCIP_RETCODE graph_fixed_add(
   SCIP*                 scip,               /**< SCIP data structure */
   IDX*                  edges,              /**< edge to add or NULL */
   const int*            pseudonodes,        /**< nodes to add */
   int                   npseudonodes,       /**< number of nodes to add  */
   GRAPH*                g                   /**< the graph */
)
{
   FIXED* fixedcomponents;

   assert(scip && g);
   assert(npseudonodes >= 0);

   if( g->fixedcomponents == NULL )
      SCIP_CALL( graph_init_fixed(scip, g) );

   fixedcomponents = g->fixedcomponents;

   if( edges )
   {
      SCIP_Bool conflict = FALSE;

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(fixedcomponents->fixedges), edges, &conflict) );
      assert(!conflict);
   }

   if( npseudonodes > 0 )
   {
      int fixnnodes = fixedcomponents->nfixnodes;
      const int fixnnodes_new = fixnnodes + npseudonodes;

      assert(pseudonodes);

      if( fixnnodes == 0 )
      {
         assert(!fixedcomponents->fixpseudonodes);
         assert(fixnnodes_new == npseudonodes);

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(fixedcomponents->fixpseudonodes), fixnnodes_new) );
      }
      else
      {
         assert(fixedcomponents->fixpseudonodes);

         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(fixedcomponents->fixpseudonodes), fixnnodes, fixnnodes_new) );
      }

      for( int i = 0; i < npseudonodes; i++ )
         fixedcomponents->fixpseudonodes[fixnnodes++] = pseudonodes[i];

      assert(fixnnodes == fixnnodes_new);

      fixedcomponents->nfixnodes = fixnnodes_new;
   }

   return SCIP_OKAY;
}

/** adds ancestors from given edges */
SCIP_RETCODE graph_fixed_addEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge,               /**< edge */
   GRAPH*                g                   /**< the graph */
)
{
   assert(scip && g);
   assert(edge >= 0 && edge < g->edges);
   assert(g->ancestors);

   SCIP_CALL( graph_fixed_add(scip, g->ancestors[edge], NULL, 0, g) );

#if 0
   SCIP_CALL( graph_fixed_add(scip, g->ancestors[edge], graph_get_pseudoAncestors(g, edge),
         graph_get_nPseudoAncestors(g, edge), g) );
#endif
   return SCIP_OKAY;
}

/** adds ancestors from given edges */
SCIP_RETCODE graph_fixed_addNodePc(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node,               /**< node */
   GRAPH*                g                   /**< the graph */
)
{
   assert(scip && g);
   assert(node >= 0 && node < g->knots);
   assert(graph_pc_isPcMw(g));
   assert(g->ancestors && g->pcancestors);

   SCIP_CALL( graph_fixed_add(scip, g->pcancestors[node], NULL, 0, g) );

   return SCIP_OKAY;
}

/** gets fixed edges */
IDX* graph_get_fixedges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   assert(g && scip);

   if( !g->fixedcomponents )
      return NULL;

   return g->fixedcomponents->fixedges;
}

/** gets fixed pseudo eliminated nodes */
const int* graph_get_fixpseudonodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   assert(g && scip);

   if( !g->fixedcomponents )
      return NULL;

   return g->fixedcomponents->fixpseudonodes;
}

/** gets number of fixed pseudo eliminated nodes */
int graph_get_nFixpseudonodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   assert(g && scip);

   if( !g->fixedcomponents )
      return 0;

   return g->fixedcomponents->nfixnodes;
}
