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

#include "graph.h"
#include "portab.h"


/** fixed graph components, typedef FIXED */
struct fixed_graph_components
{
   IDX*                  fixedges;           /**< fixed edges */
   int*                  fixpseudonodes;     /**< fixed psuedo eliminated nodes  */
   int                   nfixnodes;          /**< number of fixed nodes  */
};

/** node ancestors resulting from pseudo-elimination, typedef FIXED */
struct pseudo_ancestors
{
   int**                 edgeblocks;         /**< blocks of ancestors for each halfedge */
   int**                 pcnodeblocks;       /**< blocks of ancestors for each node (only used for PC/MW) */
   int*                  sizes;              /**< number of ancestors for each halfedge */
   int*                  maxsizes;           /**< current maximum number of ancestors for each halfedge */
   int                   halfnedges;         /**< half number of edges */
   int                   nnodes;             /**< number of nodes */
};


/** get next power of 2 number */
static inline
uint32_t getNextPow2(uint32_t n)
{
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
    n++;
    return n;
}


/** unhash some ancestors of given edge */
static inline
void pseudoAncestors_hashCleanLimited(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   halfedge,           /**< edge for which to hash */
   int                   nAncestorsToClean,  /**< number of (first) ancestors to use for clean */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int* ancestors = pseudoancestors->edgeblocks[halfedge];

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);
   assert(nAncestorsToClean >= 0 && nAncestorsToClean <= pseudoancestors->sizes[halfedge]);

   for( int k = 0; k < nAncestorsToClean; k++ )
   {
      const int a = ancestors[k];
      assert(a >= 0 && a < pseudoancestors->nnodes);
      assert(hasharr[a] == 1);

      hasharr[a] = 0;
   }
}


/** ancestor already hashed? */
static inline
SCIP_Bool pseudoAncestors_hashIsHit(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   ancestor,           /**< ancestor to check */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   assert(pseudoancestors && hasharr);
   assert(ancestor >= 0 && ancestor < pseudoancestors->nnodes);
   assert(hasharr[ancestor] == 0 || hasharr[ancestor] == 1);

   return (hasharr[ancestor] == 1);
}


/** reallocate ancestor array of given edge */
static inline
SCIP_RETCODE pseudoAncestors_realloc(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   halfedge,           /**< edge for which to reallocate */
   int                   min_maxsize_new,    /**< the minimum new maximum size, next power of 2 value will be taken */
   PSEUDOANS*            pseudoancestors     /**< pseudo-ancestors */
)
{
   const int maxsize_up = (int) getNextPow2(min_maxsize_new);
   const int maxsize = pseudoancestors->maxsizes[halfedge];

   assert(pseudoancestors);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);
   assert(min_maxsize_new >= 1 && min_maxsize_new > maxsize);
   assert(maxsize_up >= min_maxsize_new && maxsize_up <= 2 * min_maxsize_new);
   assert(maxsize_up >= 2);

   if( maxsize == 0 )
   {
      assert(pseudoancestors->edgeblocks[halfedge] == NULL);
      assert(pseudoancestors->sizes[halfedge] == 0);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pseudoancestors->edgeblocks[halfedge]), maxsize_up) );
   }
   else
   {
      assert(pseudoancestors->edgeblocks[halfedge] != NULL);
      assert(maxsize >= 2);

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(pseudoancestors->edgeblocks[halfedge]), maxsize, maxsize_up) );
   }

   printf("new sizes: %d %d \n", min_maxsize_new, maxsize_up);

   pseudoancestors->maxsizes[halfedge] = maxsize_up;

   return SCIP_OKAY;
}


/** pseudo ancestors of given edge without conflicts? */
static inline
SCIP_Bool pseudoAncestors_isValid(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   halfedge,           /**< edge for which to check */
   int*                  hasharr             /**< clean hash array of size nnodes (wrt pseudo ancestors) */
)
{
   SCIP_Bool conflict;

   assert(pseudoancestors);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   graph_pseudoAncestors_hashConflicting(pseudoancestors, 2 * halfedge, FALSE, &conflict, hasharr);
   graph_pseudoAncestors_hashCleanConflicting(pseudoancestors, 2 * halfedge, FALSE, hasharr);

   return !conflict;
}


/** hash ancestors of given edge */
void graph_pseudoAncestors_hash(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to hash */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;
   const int* ancestors = pseudoancestors->edgeblocks[halfedge];
   const int nancestors = pseudoancestors->sizes[halfedge];

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestors[k];
      assert(a >= 0 && a < pseudoancestors->nnodes);
      assert(hasharr[a] == 0);

      hasharr[a] = 1;
   }
}


/** hash ancestors of given edge, check for conflicts */
void graph_pseudoAncestors_hashConflicting(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to hash */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   SCIP_Bool*            conflict,           /**< conflict? */
   int*                  hasharr             /**< clean hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;
   const int* ancestors = pseudoancestors->edgeblocks[halfedge];
   const int nancestors = pseudoancestors->sizes[halfedge];

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   *conflict = FALSE;

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestors[k];
      assert(a >= 0 && a < pseudoancestors->nnodes);
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


/** unhash ancestors of given edge */
void graph_pseudoAncestors_hashClean(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to hash */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;
   const int* ancestors = pseudoancestors->edgeblocks[halfedge];
   const int nancestors = pseudoancestors->sizes[halfedge];

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestors[k];
      assert(a >= 0 && a < pseudoancestors->nnodes);
      assert(hasharr[a] == 1);

      hasharr[a] = 0;
   }
}


/** unhash ancestors of given edge */
void graph_pseudoAncestors_hashCleanConflicting(
   const PSEUDOANS*      pseudoancestors,    /**< pseudo-ancestors */
   int                   edge,               /**< edge for which to hash */
   SCIP_Bool             breakOnConflict,    /**< break on conflict? */
   int*                  hasharr             /**< hash array of size nnodes (wrt pseudo ancestors) */
)
{
   const int halfedge = edge / 2;
   const int* ancestors = pseudoancestors->edgeblocks[halfedge];
   const int nancestors = pseudoancestors->sizes[halfedge];

   assert(pseudoancestors && hasharr);
   assert(halfedge >= 0 && halfedge < pseudoancestors->halfnedges);

   for( int k = 0; k < nancestors; k++ )
   {
      const int a = ancestors[k];
      assert(a >= 0 && a < pseudoancestors->nnodes);
      assert(hasharr[a] == 0 || hasharr[a] == 1);

      if( breakOnConflict && hasharr[a] == 0 )
         return;

      hasharr[a] = 0;
   }
}


/** initializes pseudo ancestors */
SCIP_RETCODE graph_init_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
)
{
   int** blocks;
   int* sizes;
   int* maxsizes;
   const int halfnedges = g->edges / 2;
   PSEUDOANS* pseudoancestors;

   assert(scip && g);
   assert(g->pseudoancestors == NULL);

   SCIP_CALL( SCIPallocMemory(scip, &pseudoancestors) );

   g->pseudoancestors = pseudoancestors;

   pseudoancestors->nnodes = g->knots;
   pseudoancestors->halfnedges = halfnedges;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(blocks), halfnedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sizes), halfnedges) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(maxsizes), halfnedges) );

   for( int e = 0; e < halfnedges; e++ )
   {
      blocks[e] = NULL;
      sizes[e] = 0;
      maxsizes[e] = 0;
   }

   pseudoancestors->edgeblocks = blocks;
   pseudoancestors->sizes = sizes;
   pseudoancestors->maxsizes = maxsizes;

   return SCIP_OKAY;
}


/** frees pseudo ancestors */
void graph_free_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< the graph */
   )
{
   PSEUDOANS* pseudoancestors;

   assert(scip && g);
   assert(g->pseudoancestors != NULL && g->pseudoancestors->nnodes >= 1);

   pseudoancestors = g->pseudoancestors;

   for( int e = pseudoancestors->halfnedges - 1; e >= 0; --e )
   {
      const int size = pseudoancestors->sizes[e];

      assert(size >= 0);

      if( size > 0 )
      {
         assert(pseudoancestors->edgeblocks[e]);
         SCIPfreeBlockMemoryArray(scip, &(pseudoancestors->edgeblocks[e]), size);
      }
   }

   SCIPfreeMemoryArray(scip, &(pseudoancestors->maxsizes));
   SCIPfreeMemoryArray(scip, &(pseudoancestors->sizes));
   SCIPfreeMemoryArray(scip, &(pseudoancestors->edgeblocks));

   SCIPfreeMemoryArray(scip, &(g->pseudoancestors));
   assert(g->pseudoancestors == NULL);
}


/** frees pseudo ancestor block for given edge */
SCIP_RETCODE graph_free_pseudoAncestorsBlock(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_free,          /**< edge for which to free pseudo ancestors */
   GRAPH*                g                   /**< the graph */
)
{
   PSEUDOANS* pseudoancestors = g->pseudoancestors;
   const int block_id = edge_free / 2;
   const int size = pseudoancestors->sizes[block_id];

   assert(scip && g && g->pseudoancestors);
   assert(block_id >= 0 && block_id < g->pseudoancestors->halfnedges);
   assert(size >= 0);

   if( size > 0 )
   {
      assert(pseudoancestors->edgeblocks[block_id]);
      SCIPfreeBlockMemoryArray(scip, &(pseudoancestors->edgeblocks[block_id]), size);
   }

   pseudoancestors->sizes[block_id] = 0;
   pseudoancestors->maxsizes[block_id] = 0;

   assert(pseudoancestors->edgeblocks[block_id] == NULL);

   return SCIP_OKAY;
}


/** returns number of pseudo ancestors for given edge */
int graph_get_nPseudoAncestors(
   const GRAPH*          g,            /**< the graph */
   int                   edge          /**< edge for which to return number of pseudo ancestors */
   )
{
   const int halfedge = edge / 2;

   assert(g && g->pseudoancestors);
   assert(halfedge >= 0 && halfedge < g->pseudoancestors->halfnedges);

   return g->pseudoancestors->sizes[halfedge];
}

/** returns array of pseudo ancestors for given edge (possibly NULL) */
const int* graph_get_pseudoAncestors(
   const GRAPH*          g,            /**< the graph */
   int                   edge          /**< edge for which to return array of pseudo ancestors */
   )
{
   const int halfedge = edge / 2;

   assert(g && g->pseudoancestors);
   assert(halfedge >= 0 && halfedge < g->pseudoancestors->halfnedges);

   return g->pseudoancestors->edgeblocks[halfedge];
}

/** appends copy of pseudo ancestors of edge_source to edge_target */
SCIP_RETCODE graph_appendCopy_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_target,        /**< edge target */
   int                   edge_source,        /**< edge source */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   PSEUDOANS* const pseudoancestors = g->pseudoancestors;
   int position_target;
   const int source = edge_source / 2;
   const int size_source = pseudoancestors->sizes[source];

   assert(scip && g && pseudoancestors && conflict);
   assert(source >= 0 && source < pseudoancestors->halfnedges);
   assert(size_source >= 0);

   *conflict = FALSE;

   /* anything to append? */
   if( size_source > 0 )
   {
      const int* const ancestors_source = pseudoancestors->edgeblocks[source];
      int* ancestors_target;
      int* hasharr;
      const int target = edge_target / 2;
      const int size_target_org = pseudoancestors->sizes[target];
      const int size_targetPlusSource = size_target_org + size_source;

      assert(target >= 0 && target < pseudoancestors->halfnedges);

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, pseudoancestors->nnodes) );

      /* need to realloc target ancestor array? */
      if( size_targetPlusSource > pseudoancestors->maxsizes[target] )
      {
         SCIP_CALL( pseudoAncestors_realloc(scip, target, size_targetPlusSource, pseudoancestors) );
      }

      ancestors_target = pseudoancestors->edgeblocks[target];

      /* mark ancestors of target */
      graph_pseudoAncestors_hash(pseudoancestors, edge_target, hasharr);

      position_target = size_target_org;

      /* add source to target */
      for( int e = 0; e < size_source; e++ )
      {
         const int ancestor = ancestors_source[e];

         if( pseudoAncestors_hashIsHit(pseudoancestors, ancestor, hasharr) )
         {
            *conflict = TRUE;

            if( breakOnConflict )
               break;

            continue;
         }

         assert(position_target < pseudoancestors->maxsizes[target]);

         ancestors_target[position_target++] = ancestor;
      }

      assert(position_target <= size_targetPlusSource);
      assert(position_target >= pseudoancestors->sizes[target]);

      pseudoancestors->sizes[target] = position_target;
      pseudoAncestors_hashCleanLimited(pseudoancestors, target, size_target_org, hasharr);

      SCIPfreeCleanBufferArray(scip, &hasharr);
   }

   return SCIP_OKAY;
}

/** appends pseudo ancestors of edge_source to edge_target, ancestors for edge_source are deleted */
SCIP_RETCODE graph_appendMove_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_target,        /**< edge target */
   int                   edge_source,        /**< edge source */
   SCIP_Bool             breakOnConflict,    /**< break on conflict */
   GRAPH*                g,                  /**< the graph */
   SCIP_Bool*            conflict            /**< conflict? */
)
{
   SCIP_CALL( graph_appendCopy_pseudoAncestors(scip, edge_target, edge_source, breakOnConflict, g, conflict) );
   SCIP_CALL( graph_free_pseudoAncestorsBlock(scip, edge_source, g) );

   return SCIP_OKAY;
}


/** adds pseudo ancestor to ancestor list of given edge */
SCIP_RETCODE graph_add_pseudoAncestor(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edge_target,        /**< edge for which to add pseudo ancestor */
   int                   ancestor,           /**< (index of) pseudo ancestors */
   GRAPH*                g                   /**< the graph */
)
{
   PSEUDOANS* const pseudoancestors = g->pseudoancestors;
   const int target = edge_target / 2;
   int* const sizes = pseudoancestors->sizes;
   int* const maxsizes = pseudoancestors->maxsizes;

   assert(scip && g && pseudoancestors);
   assert(target >= 0 && target < g->pseudoancestors->halfnedges);
   assert(ancestor >= 0 && ancestor < g->pseudoancestors->nnodes);
   assert(sizes[target] <= maxsizes[target]);

   /* need to reallocate? */
   if( sizes[target] == maxsizes[target] )
   {
      SCIP_CALL( pseudoAncestors_realloc(scip, target, sizes[target] + 1, pseudoancestors) );
   }

   assert(sizes[target] < maxsizes[target]);

   pseudoancestors->edgeblocks[edge_target][sizes[target]++] = ancestor;

#ifndef NDEBUG
   {
      int* hasharr;
      SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &hasharr, pseudoancestors->nnodes) );
      assert(pseudoAncestors_isValid(pseudoancestors, target, hasharr));
      SCIPfreeCleanBufferArray(scip, &hasharr);
   }
#endif

   return SCIP_OKAY;
}

/** valid pseudo-ancestors (no conflicts)? */
SCIP_Bool graph_valid_pseudoAncestors(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g                   /**< the graph */
)
{
   const PSEUDOANS* const pseudoancestors = g->pseudoancestors;
   const int halfnedges = pseudoancestors->halfnedges;
   int* hasharr;
   int** blocks = pseudoancestors->edgeblocks;
   SCIP_Bool isValid = TRUE;

   assert(scip && g && pseudoancestors);

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &hasharr, pseudoancestors->nnodes) );

   for( int e = 0; e < halfnedges; e++ )
   {
      if( blocks[e] != NULL && !pseudoAncestors_isValid(pseudoancestors, e, hasharr) )
      {
         isValid = FALSE;
         break;
      }
   }

   SCIPfreeCleanBufferArray(scip, &hasharr);

   return isValid;
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
