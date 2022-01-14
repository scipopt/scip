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

/**@file   extreduce_core.c
 * @brief  extended reduction techniques for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements the core algorithms for the extended reduction techniques, namely for the tree search.
 * Starting from a given component (e.g. a single edge), a number of possible extensions are checked to be able
 * to verify that this component is not part of at least one optimal solution.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
// #define STP_DEBUG_EXT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"
#include "extreducedefs.h"


#define EXT_COSTS_RECOMPBOUND 10

/*
 * Local methods
 */


/** gets next subset of given size, returns FALSE if no further subset possible, otherwise TRUE */
static inline
SCIP_Bool ksubsetGetNext(
   int                   k,                  /**< size of subset */
   int                   n,                  /**< size of underlying set */
   SCIP_Bool             isFirstSubset,      /**< first call? */
   int*                  ksubset             /**< the k-subset IN/OUT */
)
{
   assert(0 < k && k <= n);

   if( isFirstSubset )
   {
      for( int i = 0; i < k; i++ )
         ksubset[i] = i;

      return TRUE;
   }

   for( int i = k - 1; i >= 0; i-- )
   {
      /* can index be incremented? */
      if( ksubset[i] < n - k + i )
      {
         ksubset[i]++;

         /* set right entries in incremental order */
         while( ++i < k )
            ksubset[i] = ksubset[i - 1] + 1;

         return TRUE;
      }
   }

   return FALSE;
}


/** are the given extension vertices in conflict with the extension conditions? */
static inline
SCIP_Bool extensionHasImplicationConflict(
   const GRAPH*          graph,              /**< graph data structure */
   const STP_Vectype(int)  implications,     /**< implications for extroot */
   const int*            tree_deg,           /**< tree degree or NULL */
   int                   tree_root,          /**< tree root */
   const int*            extedges,           /**< extension edges */
   int                   nextedges           /**< number of extension edges */
)
{
   const int nimplications = StpVecGetSize(implications);

   assert(nextedges > 0);
   assert(graph_knot_isInRange(graph, tree_root));

   for( int i = 0; i < nimplications; i++ )
   {
      int j;
      const int impnode = implications[i];
      assert(graph_knot_isInRange(graph, impnode));
      assert(Is_term(graph->term[impnode]));

      if( tree_root == impnode )
         continue;

      if( tree_deg && tree_deg[impnode] > 0 )
         continue;

      for( j = 0; j < nextedges; j++ )
      {
         const int extnode = graph->head[extedges[j]];
         assert(graph_knot_isInRange(graph, extnode));

         if( impnode == extnode )
            break;
      }

      /* implication node not contained in extension nodes? */
      if( j == nextedges )
      {
         return TRUE;
      }
   }

   return FALSE;
}

/** returns position of last marked component before the current one */
static inline
int extStackGetPrevMarked(
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const extstack_state = extdata->extstack_state;
   int stackpos = extStackGetPosition(extdata) - 1;

   assert(stackpos >= 0);

   while( extstack_state[stackpos] != EXT_STATE_MARKED )
   {
      stackpos--;
      assert(stackpos >= 0);
   }

   return stackpos;
}


/** can the extension stack hold new components? */
static inline
SCIP_Bool extStackIsExtendable(
   int                   nextedges,          /**< number of edges for extension */
   int                   nnewcomps,          /**< number of new components */
   int                   stack_datasize,     /**< datasize of stack */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int newstacksize_upper = (stack_datasize + nnewcomps * (nextedges + 1) / 2);
   const int newncomponents_upper = extdata->extstack_ncomponents + nnewcomps;

   if( newstacksize_upper > extdata->extstack_maxsize || newncomponents_upper >= extdata->extstack_maxncomponents )
   {
      assert(extdata->extstack_state[extStackGetPosition(extdata)] == EXT_STATE_NONE);

      SCIPdebugMessage("stack too full, cannot expand \n");

      return FALSE;
   }

   return TRUE;
}


/** Adds non-expanded components (i.e. sets of edges extending at a certain leaf) to the stack.
 *  Components are ordered such that smallest component is added last. */
static inline
void extStackAddCompsNonExpanded(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            extedgesstart,      /**< starts of extension edges for one components */
   const int*            extedges,           /**< extensions edges */
   int                   nextendable_leaves, /**< number of leaves that can be extended */
   EXTDATA*              extdata             /**< extension data */
   )
{
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extStackGetPosition(extdata);
   int datasize = extstack_start[stackpos + 1];
   int extsize[STP_EXT_MAXGRAD];
   int extindex[STP_EXT_MAXGRAD];

   assert(nextendable_leaves > 0);
   assert(extedgesstart[nextendable_leaves] - extedgesstart[0] > 0);
   assert(datasize + (extedgesstart[nextendable_leaves] - extedgesstart[0]) <= extdata->extstack_maxsize);

   for( int i = 0; i < nextendable_leaves; i++ )
   {
      assert(extedgesstart[i + 1] >= 0);

      extsize[i] = extedgesstart[i + 1] - extedgesstart[i];
      extindex[i] = i;
      assert(extsize[i] >= 0);
   }

   SCIPsortDownIntInt(extsize, extindex, nextendable_leaves);

   /* put the non-empty extensions on the stack, with smallest last */
   for( int i = 0; i < nextendable_leaves; i++ )
   {
      const int index = extindex[i];

      if( extsize[i] == 0 )
      {
#ifndef NDEBUG
         assert(i > 0);
         assert(extedgesstart[index + 1] - extedgesstart[index] == 0);

         for( int j = i; j < nextendable_leaves; j++ )
            assert(extsize[j] == 0);
#endif

         break;
      }

      for( int j = extedgesstart[index]; j < extedgesstart[index + 1]; j++ )
      {
         assert(extedges[j] >= 0);
         extstack_data[datasize++] = extedges[j];
      }

      assert(stackpos < extdata->extstack_maxsize - 2);

      extstack_state[++stackpos] = EXT_STATE_NONE;
      extstack_start[stackpos + 1] = datasize;
   }

#ifdef SCIP_DEBUG
   printf("added extending edges:  \n");

   for( int i = extstack_start[extdata->extstack_ncomponents]; i < extstack_start[stackpos + 1]; i++ )
      graph_edge_printInfo(graph, extstack_data[i]);
#endif

   extdata->extstack_ncomponents = stackpos + 1;

   assert(extdata->extstack_ncomponents <= extdata->extstack_maxncomponents);
}


/** adds a new leaf */
static inline
void extLeafAdd(
   int                   leaf,              /**< leaf to add */
   EXTDATA*              extdata            /**< extension data */
)
{
   assert(extdata && extdata->tree_leaves);
   assert(leaf >= 0 && extdata->tree_deg[leaf] == 0);

   extdata->tree_leaves[(extdata->tree_nleaves)++] = leaf;
}


/** remove entry from leaves list */
static inline
void extLeafRemove(
   int                   leaf,              /**< leaf to remove */
   EXTDATA*              extdata            /**< extension data */
)
{
   int* const tree_leaves = extdata->tree_leaves;
   int nleaves;
   int position;

   assert(extdata->tree_deg[leaf] == 1);

   /* switch last leaf and leaf to be removed */
   assert(extdata->tree_nleaves > 1);

   position = extLeafFindPos(extdata, leaf);
   assert(position > 0);

   extdata->tree_nleaves--;

   nleaves = extdata->tree_nleaves;

   for( int p = position; p < nleaves; p++ )
   {
      tree_leaves[p] = tree_leaves[p + 1];
   }
}


/** Remove top component from leaves, and restore the root of the component as a leaf
 *  NOTE: assumes that the tree_deg has already been adapted */
static inline
void extLeafRemoveTop(
   const GRAPH*          graph,             /**< graph data structure */
   int                   topsize,           /**< size of top to remove */
   int                   toproot,           /**< root of the top component */
   EXTDATA*              extdata            /**< extension data */
)
{
   const int* const extstack_data = extdata->extstack_data;
   int* const tree_leaves = extdata->tree_leaves;
   const int stackpos_prev = extStackGetPrevMarked(extdata);
   const int prevedges_start = extStackGetOutEdgesStart(extdata, stackpos_prev);
   const int prevedges_end = extStackGetOutEdgesEnd(extdata, stackpos_prev);
   const int compsize = prevedges_end - prevedges_start;
   int compleaves_start;

   assert(extdata);
   assert(topsize >= 1 && topsize < extdata->tree_nleaves);
   assert(toproot >= 0);
   assert(extdata->tree_deg[toproot] == 1);
   assert(compsize > 0);

#ifndef NDEBUG
   {
      const int* const tree_deg = extdata->tree_deg;
      const int nleaves = extdata->tree_nleaves;

      for( int i = 1; i <= topsize; i++ )
      {
         const int leaf = tree_leaves[nleaves - i];
         assert(tree_deg[leaf] == 0);
      }
   }
#endif

   extdata->tree_nleaves -= (topsize - 1);

#ifndef NDEBUG
   tree_leaves[(extdata->tree_nleaves) - 1] = toproot;

   for( int i = prevedges_start; i != prevedges_end; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      assert(extLeafFindPos(extdata, head) >= extdata->tree_nleaves - compsize);
   }
#endif

   compleaves_start = extdata->tree_nleaves - compsize;

   for( int i = prevedges_start; i < prevedges_end; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      tree_leaves[compleaves_start++] = head;
   }

   assert(compleaves_start == extdata->tree_nleaves);
   assert(extdata->tree_nleaves >= 1);
}


/** adds a new inner node */
static inline
void extInnerNodeAdd(
   const GRAPH*          graph,             /**< graph data structure */
   int                   innernode,         /**< node to add */
   EXTDATA*              extdata            /**< extension data */
)
{
   const SCIP_Bool isPc = (graph->prize != NULL);

   assert(extdata->tree_innerNodes);
   assert(0 <= innernode && innernode < graph->knots);
   assert(extdata->tree_deg[innernode] > 0);
   assert(isPc == graph_pc_isPc(graph));

   if( isPc )
   {
      extdata->pcdata->tree_innerPrize += graph->prize[innernode];
   }

   extdata->tree_innerNodes[(extdata->tree_ninnerNodes)++] = innernode;
}


/** removes a new inner node */
static inline
void extInnerNodeRemoveTop(
   const GRAPH*          graph,             /**< graph data structure */
   int                   innernode,         /**< (top) node to remove */
   EXTDATA*              extdata            /**< extension data */
)
{
   assert(extdata && extdata->tree_innerNodes);
   assert(innernode >= 0);

   /* update tree leaves array */
   if( innernode != extdata->tree_root )
   {
      const SCIP_Bool isPc = (graph->prize != NULL);

      assert(extdata->tree_ninnerNodes >= 1);
      assert(innernode == extdata->tree_innerNodes[extdata->tree_ninnerNodes - 1]);
      assert(isPc == graph_pc_isPc(graph));

      extdata->tree_ninnerNodes--;

      if( isPc )
      {
         extdata->pcdata->tree_innerPrize -= graph->prize[innernode];

         assert(GE(extdata->pcdata->tree_innerPrize, 0.0));
      }
   }
   else
   {
      assert(!extInitialCompIsEdge(extdata) || extdata->tree_ninnerNodes == 0);
      assert(!extInitialCompIsStar(extdata) || extdata->tree_ninnerNodes == 1);
      assert(!extInitialCompIsGenStar(extdata) || extdata->tree_ninnerNodes == 2);
   }
}


/** adds edge to tree */
static inline
void extTreeAddEdge(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to be added */
   EXTDATA*              extdata             /**< extension data */
)
{
   int* const tree_deg = extdata->tree_deg;
   int* const tree_edges = extdata->tree_edges;
   int* const tree_parentNode = extdata->tree_parentNode;
   SCIP_Real* const tree_parentEdgeCost = extdata->tree_parentEdgeCost;
   const SCIP_Real edgecost = graph->cost[edge];
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   const int genstar_centerhead = (extdata->genstar_centeredge == -1) ? -1 : graph->head[extdata->genstar_centeredge];

   assert(tree_deg[head] == 0);
   assert(tree_deg[tail] > 0 || tail == extdata->tree_root);

   if( extdata->tree_starcenter != head && genstar_centerhead != head )
   {
      extLeafAdd(head, extdata);
   }
   else
   {
      assert(extInitialCompIsStar(extdata) || extInitialCompIsGenStar(extdata));
      assert(extIsAtInitialComp(extdata));
   }

   /* NOTE: a bit hacky, but also works for general stars, because of the order
    * in which the initial edges are processed */
   extdata->tree_cost += edgecost;
   tree_deg[head] = 1;
   tree_edges[(extdata->tree_nedges)++] = edge;
   tree_parentNode[head] = tail;
   tree_parentEdgeCost[head] = edgecost;
   tree_deg[tail]++;
}


/** removes root of stack top component from tree */
static inline
void extTreeStackTopRootRemove(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int comproot = extStackGetTopRoot(graph, extdata);

   /* update tree leaves array */
   if( !extIsAtInitialComp(extdata) )
   {
      assert(comproot != extdata->tree_root);

      extLeafRemove(comproot, extdata);
      extInnerNodeAdd(graph, comproot, extdata);
   }
   else
   {
      assert(extdata->tree_nleaves == 1);
      assert(extdata->tree_deg[comproot] == 1);
      assert(comproot == extdata->tree_root);

      /* will be reset afterwards! */
      extdata->tree_deg[comproot] = 0;

      if( extInitialCompIsStar(extdata) )
      {
         const int starcenter = extdata->tree_starcenter;
         assert(graph_knot_isInRange(graph, starcenter));

         extdata->tree_deg[starcenter] = 0;
      }
      else if( extInitialCompIsGenStar(extdata) )
      {
         const int centeredge = extdata->genstar_centeredge;
         const int starcenter = extdata->tree_starcenter;

         assert(graph_edge_isInRange(graph, centeredge));
         assert(graph_knot_isInRange(graph, starcenter));
         assert(graph->tail[centeredge] == starcenter);

         extdata->tree_deg[starcenter] = 0;
         extdata->tree_deg[graph->head[centeredge]] = 0;
      }
   }
}


/** adds top component of stack to tree */
static
void extTreeStackTopAdd(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            conflict            /**< conflict found? */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   REDDATA* const reddata = extdata->reddata;
   int* const pseudoancestor_mark = reddata->pseudoancestor_mark;
   const int stackpos = extStackGetPosition(extdata);
   int conflictIteration = -1;

   assert(!(*conflict));
   assert(stackpos >= 0);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);
   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   extreduce_redcostInitExpansion(graph, extdata);
   extTreeStackTopRootRemove(graph, extdata);

   /* add top expanded component to tree data */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];

      assert(extdata->tree_nedges < extdata->extstack_maxsize);
      assert(edge >= 0 && edge < graph->edges);

      extreduce_redcostAddEdge(graph, edge, reddata, extdata);
      extTreeAddEdge(graph, edge, extdata);

      /* no conflict found yet? */
      if( conflictIteration == -1 )
      {
         assert(*conflict == FALSE);

         graph_pseudoAncestors_hashEdgeDirty(graph->pseudoancestors, edge, TRUE, conflict, pseudoancestor_mark);

         if( *conflict )
         {
            SCIPdebugMessage("pseudoancestor conflict for edge %d \n", edge);
            conflictIteration = i;
            assert(conflictIteration >= 0);
         }
      }
   }

   /* conflict found? */
   if( conflictIteration != -1 )
   {
      assert(*conflict && conflictIteration >= 0);

      for( int i = extstack_start[stackpos]; i < conflictIteration; i++ )
      {
         const int edge = extstack_data[i];
         graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, edge, pseudoancestor_mark);
      }
   }

   extdata->tree_depth++;

   assert(!extreduce_treeIsFlawed(scip, graph, extdata));
}


/** removes top component of stack from tree */
static inline
void extTreeStackTopRemove(
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Bool             ancestor_conflict,  /**< with ancestor conflict? */
   EXTDATA*              extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   int* const pseudoancestor_mark = reddata->pseudoancestor_mark;
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const tree_deg = extdata->tree_deg;
   const int stackpos = extStackGetPosition(extdata);
   const int stackstart = extstack_start[stackpos];
   const int comproot = extStackGetTopRoot(graph, extdata);
   const int compsize = extstack_start[stackpos + 1] - extstack_start[stackpos];

   assert(compsize > 0);
   assert(stackstart > 0);
   assert(tree_deg[comproot] > 1);
   assert(extdata->extstack_state[stackpos] != EXT_STATE_NONE);

   /* remove top component */
   for( int i = stackstart; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      assert(edge >= 0 && edge < graph->edges);
      assert(comproot == graph->tail[edge]);
      assert(tree_deg[head] == 1 && tree_deg[comproot] > 1);

      extreduce_redcostRemoveEdge(edge, reddata, extdata);

      extdata->tree_cost -= graph->cost[edge];
      tree_deg[head] = 0;
      tree_deg[comproot]--;

      if( !ancestor_conflict ) /* in case of a conflict, edge is unhashed already */
         graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, edge, pseudoancestor_mark);
   }

   assert(tree_deg[comproot] == 1);

   (extdata->tree_nedges) -= compsize;
   (extdata->tree_depth)--;

   /* finally, remove top component from leaves and MST storages and restore the component root */
   extLeafRemoveTop(graph, compsize, comproot, extdata);
   extInnerNodeRemoveTop(graph, comproot, extdata);
   extreduce_mstCompRemove(graph, extdata);

   assert(extdata->tree_nedges >= 0 && extdata->tree_depth >= 0);
}


/** can any extension via edge be ruled out by using simple test?? */
static inline
SCIP_Bool extTreeRuleOutEdgeSimple(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   int                   edge                /**< edge to be tested */
)
{
   const REDDATA* const reddata = extdata->reddata;
   const int extvert = graph->head[edge];
   const int* const tree_deg = extdata->tree_deg;

   if( tree_deg[extvert] != 0 )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMessage("simple rule-out (edge-to-tree)  ");
      graph_edge_printInfo(graph, edge);
#endif

      return TRUE;
   }

   if( graph_pseudoAncestors_edgeIsHashed(graph->pseudoancestors, edge, reddata->pseudoancestor_mark) )
   {
#ifdef SCIP_DEBUG

      SCIPdebugMessage("simple rule-out (ancestor)  ");
      graph_edge_printInfo(graph, edge);
#endif

      return TRUE;
   }

   return FALSE;
}


#ifdef SCIP_DISABLED
/** can any extension via edge except for only the edge itself be ruled out? */
static
SCIP_Bool extRuleOutEdgeCombinations(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   extedge             /**< edge to be tested */
)
{
}
#endif



/** Can current singleton extension be ruled out?
 *  NOTE: Also stores vertical SDs for this singleton if not ruled out! */
static inline
SCIP_Bool extTreeRuleOutSingletonFull(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Bool ruledOutFull = FALSE;
   const int edge = extdata->extstack_data[extdata->extstack_start[extStackGetPosition(extdata)]];
   const int base = graph->tail[edge];
   const int leaf = graph->head[edge];

   SCIPdebugMessage("adding SDs for edge %d->%d \n", graph->tail[edge], graph->head[edge]);

   extreduce_mstLevelVerticalReopen(extdata);

   assert(extdata->tree_deg[base] == 2);
   assert(extdata->tree_deg[leaf] == 1);

   /* NOTE the following is needed to keep invariants */
   extdata->tree_deg[base] = 1;
   extdata->tree_deg[leaf] = 0;
   extLeafRemoveTop(graph, 1, base, extdata);
   extInnerNodeRemoveTop(graph, base, extdata);

   extreduce_mstLevelVerticalAddLeaf(scip, graph, edge, extdata, &ruledOutFull);

   extLeafRemove(base, extdata);
   extLeafAdd(leaf, extdata);
   extInnerNodeAdd(graph, base, extdata);
   extdata->tree_deg[base] = 2;
   extdata->tree_deg[leaf] = 1;

   extreduce_mstLevelVerticalClose(extdata->reddata);

   return ruledOutFull;
}


/** Can current singleton extension be ruled out by implication argument? */
static inline
SCIP_Bool extTreeRuleOutSingletonImplied(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int edge = extdata->extstack_data[extdata->extstack_start[extStackGetPosition(extdata)]];
   const int base = graph->tail[edge];
   const int *tree_deg = graph_pc_isPc(graph) ? extdata->tree_deg : NULL;
   const STP_Vectype(int) implications = extdata->reddata->nodes_implications[base];

   if( extensionHasImplicationConflict(graph, implications, tree_deg,
         extdata->tree_root, &edge, 1) )
   {
      SCIPdebugMessage("implication singleton rule-out \n");
      return TRUE;
   }

   return FALSE;
}


/** Can current tree be peripherally ruled out?
 *  NOTE: If tree cannot be ruled-out, the current component will be put into the MST storage 'reddata->msts' */
static inline
SCIP_Bool extTreeRuleOutPeriph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
//#define EXT_PRINT_STATS
#ifdef EXT_PRINT_STATS
   static SCIP_Longint contracts = 0;
   static SCIP_Longint mst = 0;
   static SCIP_Longint red = 0;
#endif

   /* if we have a singleton edge, we first check whether the edge is redundant in all extension components
    * and concomitantly compute the SDs to the current tree leafs */
   if( extStackTopIsSingleton(extdata) && !extIsAtInitialComp(extdata) )
   {
      /* NOTE SDs will also be computed if not ruled out */
      SCIP_Bool ruledOutFull = extTreeRuleOutSingletonFull(scip, graph, extdata);


      if( ruledOutFull )
         return TRUE;

      if( extTreeRuleOutSingletonImplied(scip, graph, extdata) )
         return TRUE;
   }

   if( extreduce_mstRuleOutPeriph(scip, graph, extdata) )
   {
#ifdef EXT_PRINT_STATS
      mst++;
      if( mst % 10000 == 0 )
      {
         printf("rule-out-red=%lld \n", red);
         printf("rule-out-mst=%lld \n", mst);
         printf("rule-out-contracts=%lld \n", contracts);
      }
#endif

      return TRUE;
   }

   if( extreduce_redcostRuleOutPeriph(graph, extdata) )
   {
#ifdef EXT_PRINT_STATS
      red++;
      if( red % 10000 == 0 )
      {
         printf("rule-out-red=%lld \n", red);
         printf("rule-out-mst=%lld \n", mst);
         printf("rule-out-contracts=%lld \n", contracts);
      }
#endif

      return TRUE;
   }

   if( extreduce_contractionRuleOutPeriph(scip, graph, extdata) )
   {
#ifdef EXT_PRINT_STATS
      contracts++;
      if( contracts % 10000 == 0 )
      {
         printf("rule-out-red=%lld \n", red);
         printf("rule-out-mst=%lld \n", mst);
         printf("rule-out-contracts=%lld \n", contracts);
      }
#endif

      return TRUE;
   }

   return FALSE;
}


/** Stores extensions of tree from current (expanded and marked) stack top that cannot be ruled-out. */
static inline
void extTreeFindExtensions(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   int*                  extedgesstart,      /**< starts of extension edges for one components */
   int*                  extedges,           /**< extensions edges */
   int*                  nextensions,        /**< number of all extensions */
   int*                  nextendableleaves,  /**< number of leaves that can be extended */
   SCIP_Bool*            with_ruledout_leaf  /**< one leaf could already be ruled out? */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const SCIP_Bool* const isterm = extdata->node_isterm;
   const int stackpos = extStackGetPosition(extdata);
   const int topedges_start = extStackGetTopOutEdgesStart(extdata, stackpos);
   const int topedges_end = extStackGetTopOutEdgesEnd(extdata, stackpos);

#ifndef NDEBUG
   assert(graph && extdata && extedgesstart && extedges && nextensions && nextendableleaves && with_ruledout_leaf);
   assert(stackpos >= 0);
   assert(EXT_STATE_MARKED == extdata->extstack_state[stackpos]);
   assert(!(*with_ruledout_leaf));
   assert(*nextensions == 0 && *nextendableleaves == 0);
   assert(topedges_start < topedges_end);

   extreduce_extendInitDebug(extedgesstart, extedges);
#endif

   extedgesstart[0] = 0;

   /* loop over all leaves of top component */
   for( int i = topedges_start; i < topedges_end; i++ )
   {
      int nleafextensions = 0;
      const int leaf = graph->head[extstack_data[i]];

      assert(extstack_data[i] >= 0 && extstack_data[i] < graph->edges);

      /* extensions from leaf not possible? */
      if( !extLeafIsExtendable(graph, isterm, leaf) )
         continue;

      // todo extra method...skipping all edges except for the reverted one
      if( extIsAtInitialStar(extdata) && extdata->extcomp->nextleaves == 1 && extdata->extcomp->extleaves[0] != leaf )
         continue;

      /* assemble feasible single edge extensions from leaf */
      for( int e = graph->outbeg[leaf]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( extTreeRuleOutEdgeSimple(graph, extdata, e) )
         {
            continue;
         }

         assert(*nextensions < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD);
         extedges[(*nextensions)++] = e;
         nleafextensions++;
      }

      if( nleafextensions == 0 )
      {
         *with_ruledout_leaf = TRUE;
         return;
      }

      extedgesstart[++(*nextendableleaves)] = *nextensions;
   }

   assert(*nextensions >= *nextendableleaves);
   assert(*nextendableleaves <= STP_EXT_MAXGRAD);
}


/** synchronize tree with the stack */
static inline
void extTreeSyncWithStack(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            conflict            /**< conflict found? */
)
{
   const int stackposition = extStackGetPosition(extdata);

   assert(scip && graph && extdata && conflict);
   assert(!(*conflict));

#ifdef SCIP_DEBUG
   extreduce_printStack(graph, extdata);
#endif

   /* is current component expanded? */
   if( extdata->extstack_state[stackposition] == EXT_STATE_EXPANDED && !extStackTopIsWrapped(extdata) )
   {
      /* add top component to tree */
      extTreeStackTopAdd(scip, graph, extdata, conflict);
   }

   /* recompute edge costs and reduced costs? */
   if( ++(extdata->ncostupdatestalls) > EXT_COSTS_RECOMPBOUND )
   {
      extreduce_treeRecompCosts(scip, graph, extdata);
      extdata->ncostupdatestalls = 0;
   }

#ifndef NDEBUG
   if( !(*conflict) )
   {
      assert(extreduce_treeIsHashed(graph, extdata));
   }
#endif
}


/** should we truncate from current component? */
static
SCIP_Bool extTruncate(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const SCIP_Bool* const  isterm = extdata->node_isterm;
   const int stackpos = extStackGetPosition(extdata);

   assert(extdata->extstack_state[stackpos] == EXT_STATE_MARKED);
   assert(extstack_start[stackpos] < extdata->extstack_maxsize);

   if( extdata->tree_depth >= extdata->tree_maxdepth )
   {
      SCIPdebugMessage("truncate (depth too high) \n");
      return TRUE;
   }

   if( extdata->tree_nedges >= extdata->tree_maxnedges )
   {
      SCIPdebugMessage("truncate (too many tree edges) \n");
      return TRUE;
   }

   if( extdata->tree_nleaves >= extdata->tree_maxnleaves )
   {
      SCIPdebugMessage("truncate (too many leaves) \n");
      return TRUE;
   }

   if( extdata->extstack_ncomponents >= extdata->extstack_maxncomponents - 1 )
   {
      SCIPdebugMessage("truncate (too many stack components) \n");
      return TRUE;
   }

   /* check whether at least one leaf is extendable */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int leaf = graph->head[edge];

      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[leaf] > 0);

      if( extLeafIsExtendable(graph, isterm, leaf) )
         return FALSE;
   }

   SCIPdebugMessage("truncate (non-promising) \n");
   return TRUE;
}


/** top component is rebuilt, and
 *  if success == TRUE: goes back to first marked component
 *  if success == FALSE: goes back to first marked or non-expanded component
 *   */
static
void extBacktrack(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Bool             success,            /**< backtrack from success? */
   SCIP_Bool             ancestor_conflict,  /**< backtrack triggered by ancestor conflict? */
   EXTDATA*              extdata             /**< extension data */
)
{
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extStackGetPosition(extdata);
   const int extstate_top = extstack_state[stackpos];

   assert(graph && extdata);
   assert(extdata->extstack_start[stackpos + 1] - extdata->extstack_start[stackpos] > 0);

   /* top component already expanded (including marked)? */
   if( extstate_top != EXT_STATE_NONE )
   {
      extTreeStackTopRemove(graph, ancestor_conflict, extdata);
   }

   stackpos--;

   /* backtrack */
   if( success )
   {
      if( extstack_state[stackpos] != EXT_STATE_EXPANDED  )
      {
         /* the MST level associated top component cannot be used anymore, because the next component is not a sibling */
         extreduce_mstLevelRemove(extdata->reddata);
      }

      while( extstack_state[stackpos] == EXT_STATE_NONE )
      {
         stackpos--;
         assert(stackpos >= 0);
      }

      SCIPdebugMessage("backtrack SUCCESS \n");
      assert(extstack_state[stackpos] == EXT_STATE_EXPANDED || extstack_state[stackpos] == EXT_STATE_MARKED);
   }
   else
   {
      /* the MST level associated with top component cannot be used anymore, because siblings will be removed */
      extreduce_mstLevelRemove(extdata->reddata);

      while( extstack_state[stackpos] == EXT_STATE_EXPANDED )
      {
         stackpos--;
         assert(stackpos >= 0);
      }

      SCIPdebugMessage("backtrack FAILURE \n");
      assert(extstack_state[stackpos] == EXT_STATE_NONE || extstack_state[stackpos] == EXT_STATE_MARKED);
   }

   extdata->extstack_ncomponents = stackpos + 1;

   assert(extdata->extstack_ncomponents <= extdata->extstack_maxncomponents);
   assert(!extreduce_treeIsFlawed(scip, graph, extdata));
}

/** Builds components from top edges and adds them.
 *  Backtracks if stack is too full.
 *  Helper method for 'extStackTopExpand' */
static inline
void extStackAddCompsExpanded(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   nextedges,          /**< number of edges for extension */
   const int*            extedges,           /**< array of edges for extension */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success,            /**< success pointer */
   SCIP_Bool*            ruledOut            /**< all ruled out? */
)
{
   int ksubset[STP_EXT_MAXGRAD];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   const int extroot = graph->tail[extedges[0]];
   const STP_Vectype(int) implications = extdata->reddata->nodes_implications[extroot];
   int stackpos = extStackGetPosition(extdata);
   int datasize = extstack_start[stackpos];
   const uint32_t powsize = (uint32_t) pow(2.0, nextedges);
   const int* tree_deg = graph_pc_isPc(graph) ? extdata->tree_deg : NULL;

   assert(0 < nextedges && nextedges < STP_EXT_MAXGRAD);
   assert(powsize >= 2);
   assert(extstack_data[datasize] == EXT_EDGE_WRAPPED);

   *ruledOut = FALSE;

   /* stack too full? */
   if( !extStackIsExtendable(nextedges, (int) powsize, datasize, extdata) )
   {
      *success = FALSE;
      extBacktrack(scip, graph, *success, FALSE, extdata);

      return;
   }

   /* remove the dummy level */
   extreduce_mstLevelHorizontalRemove(extdata->reddata);

   extreduce_mstLevelHorizontalAdd(scip, graph, nextedges, extedges, extdata);

   /* compute and add components (overwrite previous, non-expanded component) */
   // todo since single extensions are used first, might make sense to also have a special bottleneck test for this case!
   // -also good if we have the method that excludes everything except for the singletons!
   for( int setsize = nextedges; setsize >= 2; setsize--  )
   {
      SCIP_Bool isFirst = TRUE;

      while( ksubsetGetNext(setsize, nextedges, isFirst, ksubset) )
      {
         const int datasize_prev = datasize;
         isFirst = FALSE;

         for( int j = 0; j < setsize; j++ )
         {
            const int pos = ksubset[j];

            assert(datasize < extdata->extstack_maxsize);
            assert(graph->tail[extedges[pos]] == extroot);

            extstack_data[datasize++] = extedges[pos];
            SCIPdebugMessage(" head %d \n", graph->head[extedges[pos]]);
         }

         SCIPdebugMessage("... added \n");
         assert(stackpos < extdata->extstack_maxsize - 1);

         if( extensionHasImplicationConflict(graph, implications, tree_deg, extdata->tree_root,
            &(extstack_data[datasize_prev]), datasize - datasize_prev) )
         {
            SCIPdebugMessage("implication conflict found for root %d \n", extroot);
            datasize = datasize_prev;
         }
         else
         {
            extstack_state[stackpos] = EXT_STATE_EXPANDED;
            extstack_start[++stackpos] = datasize;
         }

         assert(extstack_start[stackpos] - extstack_start[stackpos - 1] > 0);
      }
   }

   /* nothing added? */
   if( stackpos == extStackGetPosition(extdata) )
   {
      assert(datasize == extstack_start[stackpos]);
      *ruledOut = TRUE;
      extstack_data[datasize] = EXT_EDGE_WRAPPED;
   }
   else
   {
      assert(stackpos > extStackGetPosition(extdata));
      assert(stackpos >= extdata->extstack_ncomponents);
      assert(stackpos <= extdata->extstack_maxncomponents);

      extdata->extstack_ncomponents = stackpos;
   }
}


/** Builds singleton components from top edges and adds them.
 *  Backtracks if stack is too full */
static inline
void extStackAddCompsExpandedSing(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   nextedges,          /**< number of edges for extension */
   const int*            extedges,           /**< array of edges for extension */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success,            /**< success pointer */
   SCIP_Bool*            ruledOut            /**< all ruled out? */
)
{
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   const int extroot = graph->tail[extedges[0]];
   int stackpos = extStackGetPosition(extdata);
   int datasize = extstack_start[stackpos];

   assert(nextedges > 0 && nextedges < STP_EXT_MAXGRAD);
   *ruledOut = FALSE;

   /* stack too full? */
   if( !extStackIsExtendable(nextedges, nextedges + 1, datasize, extdata) )
   {
      *success = FALSE;
      extBacktrack(scip, graph, *success, FALSE, extdata);
      return;
   }

   extreduce_mstLevelHorizontalAddEmpty(graph, extdata);
   extreduce_mstLevelClose(scip, graph, extroot, extdata);

   /* first we add the wrapped component */
   extstack_data[datasize++] = EXT_EDGE_WRAPPED;
   extstack_state[stackpos] = EXT_STATE_EXPANDED;
   extstack_start[++stackpos] = datasize;

   /* compute and add singleton components (overwrite previous, non-expanded component) */
   for( int pos = 0; pos < nextedges; pos++  )
   {
      assert(datasize < extdata->extstack_maxsize && stackpos < extdata->extstack_maxsize - 1);
      assert(graph->tail[extedges[pos]] == extroot);

      extstack_data[datasize++] = extedges[pos];
      SCIPdebugMessage("add singleton head %d \n", graph->head[extedges[pos]]);

      extstack_state[stackpos] = EXT_STATE_EXPANDED;
      extstack_start[++stackpos] = datasize;

      assert(extstack_start[stackpos] - extstack_start[stackpos - 1] > 0);
   }

   assert(stackpos > extStackGetPosition(extdata));

   /* only wrapped component added? */
   if( stackpos == extStackGetPosition(extdata) + 1 )
   {
      *ruledOut = TRUE;
   }
   else
   {
      assert(stackpos > extStackGetPosition(extdata));
      assert(stackpos >= extdata->extstack_ncomponents);
      assert(stackpos <= extdata->extstack_maxncomponents);

      extdata->extstack_ncomponents = stackpos;
   }
}


/** adds initial component as 'expanded' to stack (including computation of horizontal SDs) */
static inline
void extStackAddCompInitialExpanded(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int* extedges;
   int nextedges;

   /* NOTE: in case of star component we want to skip the root edge for the horizontal SD computation */
   if( extInitialCompIsStar(extdata) )
   {
      extedges = &(extstack_data[1]);
      nextedges = (extstack_start[1] - 1);
   }
   else if( extInitialCompIsGenStar(extdata) )
   {
      extedges = &(extstack_data[2]);
      nextedges = (extstack_start[1] - 2);
   }
   else
   {
      extedges = extstack_data;
      nextedges = extstack_start[1];
   }

   assert(nextedges > 0);
   assert(graph->tail[extstack_data[0]] == extdata->tree_root);

   extreduce_mstLevelHorizontalAdd(scip, graph, nextedges, extedges, extdata);
   extreduce_mstLevelClose(scip, graph, extdata->tree_root, extdata);

   extdata->extstack_state[0] = EXT_STATE_EXPANDED;
   extdata->extstack_ncomponents = 1;
}


/** Collects edges top component of stack that we need to consider for extension
 *  (i.e. which cannot be ruled out).
 *  Helper method for 'extStackTopExpand' */
static inline
void extStackTopCollectExtEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   int*                  extedges,           /**< array of collected edges */
   int*                  nextedges           /**< number of edges */
)
{
   const MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const int stackpos = extStackGetPosition(extdata);
   const int toplevel = extreduce_mldistsTopLevel(sds_vertical);
   const int comproot = graph->tail[extdata->extstack_data[extdata->extstack_start[stackpos] + 1]];
   // todo make that less hacky

   assert(*nextedges == 0);
   assert(extreduce_mstTopCompInSync(scip, graph, extdata));

   for( int e = graph->outbeg[comproot]; e != EAT_LAST; e = graph->oeat[e] )
   {
      if( extdata->tree_deg[graph->head[e]] != 0 )
         continue;

      if( extreduce_mldistsLevelContainsBase(sds_vertical, toplevel, graph->head[e]) )
      {
         SCIPdebugMessage("collecting edge %d->%d \n", graph->tail[e], graph->head[e]);
         extedges[(*nextedges)++] = e;
      }
   }
}



/** Collects edges top component of stack that we need to consider for singletons extension */
static inline
void extStackTopCollectExtEdgesSing(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   int*                  extedges,           /**< array of collected edges */
   int*                  nextedges           /**< number of edges */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extStackGetPosition(extdata);

   assert(*nextedges == 0);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] >= 1);
   assert(extreduce_mstTopCompInSync(scip, graph, extdata));

   /* collect edges for components (and try to rule each of them out) */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];

      assert(*nextedges < STP_EXT_MAXGRAD);
      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[graph->head[edge]] == 0);
      assert(!extTreeRuleOutEdgeSimple(graph, extdata, edge));

      extedges[(*nextedges)++] = edge;
   }
}


/** Gets start of data for initial component */
static inline
int extStackTopGetInitalDataStart(
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const extstack_start = extdata->extstack_start;
   int start;

   assert(extStackGetPosition(extdata) == 0);

   if( extInitialCompIsEdge(extdata) )
      start = extstack_start[0];
   else if( extInitialCompIsGenStar(extdata) )
      start = extstack_start[0] + 2;
   else
      start = extstack_start[0] + 1;

   return start;
}


/** Computes ancestor SDs for leaves of initial component.
 *  Also checks for possible rule-out. */
static inline
void extStackTopProcessInitialEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            initialRuleOut      /**< pointer to mark early rule-out */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extStackGetPosition(extdata);
   const SCIP_Bool compIsEdge = extInitialCompIsEdge(extdata);
   const int data_start = extStackTopGetInitalDataStart(extdata);
   const int data_end = extstack_start[stackpos + 1];

   assert(*initialRuleOut == FALSE);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] >= 1);
   assert(extreduce_mstTopCompInSync(scip, graph, extdata));
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] == extstack_start[1]);
   assert(data_start < data_end);

   /* compute SDs for all leaves (and try to rule each of them out) */
   for( int i = data_start; i < data_end; i++ )
   {
      const int edge = extstack_data[i];

      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[graph->head[edge]] == 0);
      assert(!extTreeRuleOutEdgeSimple(graph, extdata, edge));

      /* computes the SDs from 'leaf' to all tree leaves in 'sds_vertical', unless the edge is ruled out */
      extreduce_mstLevelVerticalAddLeafInitial(scip, graph, edge, extdata, initialRuleOut);

      if( *initialRuleOut )
      {
         SCIPdebugMessage("initial component ruled-out! \n");
         return;
      }

#if 0
     if( extRuleOutEdgeCombinations(graph, extdata, e) )
     {
        // todo: need some marker here to say that we have a single edge!
        continue;
     }
#endif
   }

   if( !compIsEdge )
   {
      const STP_Vectype(int) implications = extdata->reddata->nodes_implications[extdata->tree_starcenter];

      assert(extdata->tree_deg[graph->tail[extstack_data[0]]] == 1);
      assert(graph_knot_isInRange(graph, extdata->tree_starcenter));

      /* NOTE: need to be careful here, because the fist edge is an in-edge!  */
      if( extensionHasImplicationConflict(graph, implications, extdata->tree_deg, extdata->tree_root,
          &(extstack_data[data_start]), data_end - data_start) )
       {
          SCIPdebugMessage("implication conflict found for initial star component \n");

          *initialRuleOut = TRUE;
       }
   }
}


/** Expands top component of stack to singletons */
static
void extStackTopExpandSingletons(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD];
   int nextedges = 0;
   SCIP_Bool ruledOut = FALSE;

#ifndef NDEBUG
   const int stackpos = extStackGetPosition(extdata);
   const int* const extstack_state = extdata->extstack_state;
   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;

   assert(scip && graph && success);
   assert(EXT_STATE_NONE == extstack_state[stackpos]);
#endif
   SCIPdebugMessage("expanding singletons \n");

   extStackTopCollectExtEdgesSing(scip, graph, extdata, extedges, &nextedges);
   assert(nextedges > 0);

   /* add a placeholder level to keep invariants */
   extreduce_mstLevelVerticalAddEmpty(graph, extdata);

   /* use the just collected edges 'extedges' to build singleton components and add them to the stack */
   extStackAddCompsExpandedSing(scip, graph, nextedges, extedges, extdata, success, &ruledOut);

#ifndef NDEBUG
   {
      const int stackpos_new = extStackGetPosition(extdata);
      assert(extstack_state[stackpos_new] == EXT_STATE_EXPANDED || (stackpos_new < stackpos) );
   }
#endif
}


/** Expands wrapped component of stack
 *  by adding all possible subsets of the component that cannot be ruled-out.
 */
static
void extStackTopExpandWrapped(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD];
   int* const extstack_state = extdata->extstack_state;
   const int stackpos = extStackGetPosition(extdata);
   int nextedges = 0;
   SCIP_Bool ruledOut = FALSE;

#ifndef NDEBUG
   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;

   assert(scip && graph && success);
   assert(EXT_STATE_EXPANDED == extstack_state[stackpos]);
#endif
   SCIPdebugMessage("expanding wrapped components \n");

   /* collect surviving singleton edges */
   extStackTopCollectExtEdges(scip, graph, extdata, extedges, &nextedges);

   //extreduce_mstLevelVerticalClose(reddata);

   /* everything ruled out already?
    * NOTE in the case of 1 extension edge, this edge has already been ruled out! */
   if( nextedges == 0 || nextedges == 1 )
   {
      ruledOut = TRUE;
   }
   else
   {
      /* use the just collected edges 'extedges' to build components and add them to the stack */
      extStackAddCompsExpanded(scip, graph, nextedges, extedges, extdata, success, &ruledOut);

#ifndef NDEBUG
      if( !ruledOut )
      {
         const int stackpos_new = extStackGetPosition(extdata);
         assert(extstack_state[stackpos_new] == EXT_STATE_EXPANDED || (stackpos_new < stackpos) );
      }
#endif
   }

   if( ruledOut )
   {
      *success = TRUE;

      if( extStackTopIsWrapped(extdata) )
         extstack_state[stackpos] = EXT_STATE_NONE;

      assert(extstack_state[stackpos] == EXT_STATE_NONE);
      assert(stackpos != 0);  /* not the initial component! */

      extBacktrack(scip, graph, *success, FALSE, extdata);
   }
}


/** Expands top component of stack.
 *  Note: This method can backtrack:
 *        1. If stack is full (with success set to FALSE),
 *        2. If all edges of the component can be ruled-out (with success set to TRUE).
 */
static
void extStackTopExpand(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   if( extStackTopIsWrapped(extdata) )
      extStackTopExpandWrapped(scip, graph, extdata, success);
   else
      extStackTopExpandSingletons(scip, graph, extdata, success);
}


/** same as 'extStackTopExpand', but for initial component only */
static inline
void extStackTopExpandInitial(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            initialRuleOut      /**< pointer to mark early rule-out */
)
{
   REDDATA* const reddata = extdata->reddata;

#ifndef NDEBUG
   const int nextedges = extdata->extstack_start[1];
   const int* const extstack_state = extdata->extstack_state;

   assert(EXT_STATE_NONE == extstack_state[0]);
   assert(extStackGetPosition(extdata) == 0);
   assert(extIsAtInitialComp(extdata));
   assert(nextedges == 1 || nextedges >= 3);
#endif

   extreduce_mstLevelInitialInit(reddata, extdata);
   extStackTopProcessInitialEdges(scip, graph, extdata, initialRuleOut);
   extreduce_mstLevelVerticalClose(reddata);

   assert(extstack_state[0] == EXT_STATE_NONE);

   /* initial component ruled out already? */
   if( *initialRuleOut )
   {
      extreduce_mstLevelVerticalRemove(reddata);
   }
   else
   {
      extStackAddCompInitialExpanded(scip, graph, extdata);

      assert(extstack_state[0] != EXT_STATE_NONE);
   }

   assert(*initialRuleOut || extdata->extstack_state[0] != EXT_STATE_NONE );
}


/** Extends top component of stack.
 *  Backtracks if stack is full.
 *  Will not add anything in case of rule-out of at least one extension node. */
static
void extExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer, FALSE iff no extension was possible */
)
{
   int* extedges = NULL;
   int* extedgesstart = NULL;
   int* const extstack_start = extdata->extstack_start;
   const int stackpos = extStackGetPosition(extdata);
   int nextendable_leaves = 0;
   int nextensions = 0;
   SCIP_Bool has_ruledout_leaf = FALSE;

   assert(stackpos >= 0);
   assert(EXT_STATE_MARKED == extdata->extstack_state[stackpos]);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] <= STP_EXT_MAXGRAD);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &extedges, STP_EXT_MAXGRAD * STP_EXT_MAXGRAD) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &extedgesstart, STP_EXT_MAXGRAD + 1) );

   extTreeFindExtensions(scip, graph, extdata, extedgesstart, extedges, &nextensions,
         &nextendable_leaves, &has_ruledout_leaf);

   if( has_ruledout_leaf )
   {
      *success = TRUE;

      SCIPdebugMessage("ruled-out one leaf \n");
   }
   else if( nextendable_leaves == 0 )  /* found no valid extensions? */
   {
      *success = FALSE;

      assert(nextensions == 0);
      SCIPdebugMessage("no valid extensions found \n");
   }
   else  /* found non-empty valid extensions */
   {
      const int stacksize_new = extstack_start[stackpos + 1] + (extedgesstart[nextendable_leaves] - extedgesstart[0]);
      assert(nextendable_leaves > 0 && nextensions > 0);

      /* stack too small? */
      if( stacksize_new > extdata->extstack_maxsize )
      {
         *success = FALSE;
         assert(EXT_STATE_MARKED == extdata->extstack_state[extStackGetPosition(extdata)]);

         SCIPdebugMessage("stack is full! need to backtrack \n");

         extBacktrack(scip, graph, *success, FALSE, extdata);
      }
      else
      {
         *success = TRUE;

         extStackAddCompsNonExpanded(graph, extedgesstart, extedges, nextendable_leaves, extdata);

         /* try to expand last (and smallest!) component, which currently is just a set of edges */
         extStackTopExpand(scip, graph, extdata, success);
      }
   }

   SCIPfreeBufferArray(scip, &extedgesstart);
   SCIPfreeBufferArray(scip, &extedges);
}


/** adds initial single edge to stack */
static inline
void extPreprocessInitialEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int edge = extcomp->compedges[0];

   assert(0 <= edge && edge < graph->edges);
   assert(extcomp->ncompedges == 1);
   assert(graph->tail[edge] == extcomp->comproot);

#ifdef SCIP_DEBUG
   printf("\n --- ADD initial EDGE component --- \n\n");
   SCIPdebugMessage("...initial edge %d: %d->%d \n\n", edge, graph->tail[edge], graph->head[edge]);
#endif

   extdata->extstack_data[0] = edge;
}


/** adds initial star component edges to stack */
/*  NOTE: it is vital the the first edge of the star component comes from the root! */
static inline
void extPreprocessInitialStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int* const compedges = extcomp->compedges;
   const int ncompedges = extcomp->ncompedges;
   const int rootedge = compedges[0];
   const int starcenter = graph->head[rootedge];
   const int comproot = extcomp->comproot;

   assert(ncompedges >= 3);
   assert(graph->tail[rootedge] == comproot);

#ifdef SCIP_DEBUG
   printf("\n --- ADD initial STAR component --- \n\n");
   SCIPdebugMessage("...root star edge %d: %d->%d \n", rootedge, comproot, starcenter);
#endif

   extdata->extstack_data[0] = rootedge;
   extdata->tree_deg[starcenter] = ncompedges;
   extdata->tree_parentNode[starcenter] = comproot;
   extdata->tree_parentEdgeCost[starcenter] = graph->cost[rootedge];
   extdata->tree_starcenter = starcenter;

   for( int i = 1; i < ncompedges; i++ )
   {
      const int e = compedges[i];

      SCIPdebugMessage("...star edge %d: %d->%d \n", e, graph->tail[e], graph->head[e]);
      assert(graph->tail[e] == starcenter);
      extdata->extstack_data[i] = e;
   }

   extInnerNodeAdd(graph, starcenter, extdata);

#ifdef SCIP_DEBUG
   printf(" \n");
#endif
}



/** adds initial general star component edges to stack */
/*  NOTE: it is vital the the first edge of the star component comes from the root! */
static inline
void extPreprocessInitialGenStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int* const compedges = extcomp->compedges;
   const int ncompedges = extcomp->ncompedges;
   const int rootedge = compedges[0];
   const int starcenter = graph->head[rootedge];
   const int comproot = extcomp->comproot;
   const int centeredge = extdata->genstar_centeredge;
   const int centerhead = graph->head[centeredge];

   assert(extInitialCompIsGenStar(extdata));
   assert(ncompedges >= 3);
   assert(graph->tail[rootedge] == comproot);
   assert(graph_edge_isInRange(graph, extdata->genstar_centeredge));
   assert(graph->tail[centeredge] == starcenter);

#ifdef SCIP_DEBUG
   printf("\n --- ADD initial GENERAL star component --- \n\n");
   SCIPdebugMessage("...root star edge %d: %d->%d \n", rootedge, comproot, starcenter);
   SCIPdebugMessage("...center star edge %d: %d->%d \n", centeredge, starcenter, centerhead);
#endif

   extdata->extstack_data[0] = rootedge;
   /* NOTE: one time for root edge, one time for center edge */
   extdata->tree_deg[starcenter] = 2;
   extdata->tree_parentNode[starcenter] = comproot;
   extdata->tree_parentEdgeCost[starcenter] = graph->cost[rootedge];
   extdata->tree_starcenter = starcenter;

   extdata->extstack_data[1] = centeredge;
   extdata->tree_deg[centerhead] = 1;
   extdata->tree_parentNode[centerhead] = starcenter;
   extdata->tree_parentEdgeCost[centerhead] = graph->cost[centeredge];

   for( int i = 1; i < ncompedges; i++ )
   {
      const int e = compedges[i];

      SCIPdebugMessage("...star edge %d: %d->%d \n", e, graph->tail[e], graph->head[e]);

      extdata->tree_deg[graph->tail[e]]++;
      extdata->extstack_data[i + 1] = e;
   }

   assert(extdata->tree_deg[starcenter] >= 2);
   assert(extdata->tree_deg[centerhead] >= 2);

   extInnerNodeAdd(graph, starcenter, extdata);
   extInnerNodeAdd(graph, centerhead, extdata);

#ifdef SCIP_DEBUG
   printf(" \n");
#endif
}


/** Initial component preprocessing:
 *  The component root is added to the tree and the stack,
 *  and the remainder is added to the stack to allow for further expansion. */
static inline
void extPreprocessInitialComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Real* redcost_treenodeswap = extdata->reddata->redcost_treenodeswaps;
   const int ncompedges = extcomp->ncompedges;
   const int comproot = extcomp->comproot;
   const SCIP_Bool compIsEdge = (ncompedges == 1);
   const SCIP_Bool compIsGenStar = extInitialCompIsGenStar(extdata);
   const int redcost_nlevels = extdata->reddata->redcost_nlevels;
   const int nnodes = graph->knots;

   assert(ncompedges >= 1 && ncompedges < STP_EXT_MAXGRAD);
   assert(comproot >= 0 && comproot < graph->knots);
   assert(redcost_nlevels >= 1);

   if( compIsEdge )
      extPreprocessInitialEdge(scip, graph, extcomp, extdata);
   else if( compIsGenStar )
      extPreprocessInitialGenStar(scip, graph, extcomp, extdata);
   else
      extPreprocessInitialStar(scip, graph, extcomp, extdata);

   /* for now, only the component root is marked as a leaf */
   extreduce_mstAddRootLevel(scip, comproot, extdata);
   extLeafAdd(comproot, extdata);

   extdata->tree_root = comproot;
   extdata->extstack_ncomponents = 1;
   extdata->extstack_state[0] = EXT_STATE_NONE;
   extdata->extstack_start[0] = 0;
   extdata->extstack_start[1] = ncompedges;
   extdata->tree_parentNode[comproot] = -1;
   extdata->tree_parentEdgeCost[comproot] = -1.0;

   for( int i = 0; i < redcost_nlevels; i++ )
      redcost_treenodeswap[comproot + i * nnodes] = 0.0;

   if( compIsGenStar )
   {
      /* center edge is also included... */
      extdata->extstack_start[1]++;
   }

   assert(extdata->tree_leaves[0] == comproot);
   assert(extdata->tree_deg[comproot] == 0);
   assert(extdata->tree_nleaves == 1);

   extdata->tree_deg[comproot] = 1;
}


/** helper for rule-out during initial component processing */
static inline
void extUnhashInitialComponent(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata             /**< extension data */
)
{
   const PSEUDOANS* const pseudoancestors = graph->pseudoancestors;
   const int* const compedges = extcomp->compedges;
   int* const pseudoancestor_mark = extdata->reddata->pseudoancestor_mark;
   const int ncompedges = extcomp->ncompedges;

   /* necessary because these edges will be deleted in clean-up otherwise */
   for( int i = 0; i < ncompedges; i++ )
   {
      graph_pseudoAncestors_unhashEdge(pseudoancestors, compedges[i], pseudoancestor_mark);
   }

   if( extInitialCompIsGenStar(extdata) )
   {
      const int centeredge = extcomp->genstar_centeredge;
      assert(graph_edge_isInRange(graph, centeredge));

      graph_pseudoAncestors_unhashEdge(pseudoancestors, centeredge, pseudoancestor_mark);
   }
}


/** Adds extensions initial component to stack (needs to be star component rooted in root).
  * If no extensions are added, then the component has been ruled-out. */
static inline
void extProcessInitialComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            ruledOut,           /**< initial component ruled out? */
   SCIP_Bool*            isExtendable        /**< extension possible? */
)
{
   SCIP_Bool conflict = FALSE;

   assert(FALSE == (*ruledOut) && TRUE == (*isExtendable));

   extPreprocessInitialComponent(scip, graph, extcomp, extdata);
   extStackTopExpandInitial(scip, graph, extdata, ruledOut);

   assert(extStackGetPosition(extdata) == 0);

   /* early rule-out? */
   if( *ruledOut )
   {
      return;
   }

   assert(extdata->extstack_state[0] == EXT_STATE_EXPANDED);

   extTreeStackTopAdd(scip, graph, extdata, &conflict);

   if( conflict )
   {
      assert(extInitialCompIsStar(extdata) || extInitialCompIsGenStar(extdata));
      *ruledOut = TRUE;
      return;
   }

   assert(extdata->tree_deg[extdata->tree_root] == 1);
   assert(extdata->tree_deg[graph->head[extdata->extstack_data[0]]] == extcomp->ncompedges
         || extInitialCompIsGenStar(extdata));

   /* NOTE: anyway necessary to keep the MST graph up-to-date! */
   if( extTreeRuleOutPeriph(scip, graph, extdata) )
   {
      *ruledOut = TRUE;
      extUnhashInitialComponent(graph, extcomp, extdata);
      return;
   }

   /* the initial component could not be ruled-out, so set its stage to 'marked' */
   extdata->extstack_state[0] = EXT_STATE_MARKED;
   assert(extStackGetPosition(extdata) == 0);
   extExtend(scip, graph, extdata, isExtendable);

   if( !(*isExtendable) )
   {
      extUnhashInitialComponent(graph, extcomp, extdata);
   }

   assert(*isExtendable || extcomp->ncompedges >= 3);
}


/** Checks whether component 'extcomp' (star or single edge) can be ruled out. */
static
void extProcessComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< initial component to be checked */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            deletable           /**< is component deletable? */
)
{
   int* const extstack_state = extdata->extstack_state;
   SCIP_Bool success = TRUE;
   SCIP_Bool conflict = FALSE;

   assert(extreduce_extdataIsClean(graph, extdata));
   assert(extreduce_reddataIsClean(graph, extdata->reddata));
   assert(extreduce_pcdataIsClean(graph, extdata->pcdata));
   assert(FALSE == (*deletable));

   /* put initial component on the stack */
   extProcessInitialComponent(scip, graph, extcomp, extdata, deletable, &success);

   /* early rule-out? or no extension possible? */
   if( *deletable || !success )
   {
      extreduce_extCompClean(scip, graph, extcomp, FALSE, extdata);
      return;
   }

   assert(extstack_state[0] == EXT_STATE_MARKED);

   /* limited DFS backtracking; stops once back at initial component */
   while( extdata->extstack_ncomponents > 1 )
   {
      const int stackposition = extStackGetPosition(extdata);
      conflict = FALSE;

      extTreeSyncWithStack(scip, graph, extdata, &conflict);

      /* has current component already been extended? */
      if( extstack_state[stackposition] == EXT_STATE_MARKED )
      {
         extBacktrack(scip, graph, success, FALSE, extdata);
         continue;
      }

      /* component not expanded yet or is wrapped? */
      if( extstack_state[stackposition] != EXT_STATE_EXPANDED || extStackTopIsWrapped(extdata) )
      {
         assert(extstack_state[stackposition] == EXT_STATE_NONE || extStackTopIsWrapped(extdata));

         extStackTopExpand(scip, graph, extdata, &success);
         continue;
      }

      assert(extstack_state[stackposition] == EXT_STATE_EXPANDED);

      if( conflict || extTreeRuleOutPeriph(scip, graph, extdata) )
      {
         success = TRUE;
         extBacktrack(scip, graph, success, conflict, extdata);
         continue;
      }

      /* the component could not be ruled-out, so set its stage to 'marked' */
      extstack_state[stackposition] = EXT_STATE_MARKED;

      if( extTruncate(graph, extdata) )
      {
         success = FALSE;
         extBacktrack(scip, graph, success, FALSE, extdata);
         continue;
      }

      /* neither ruled out nor truncated, so extend */
      extExtend(scip, graph, extdata, &success);

   } /* DFS loop */

   *deletable = success;

   extreduce_extCompClean(scip, graph, extcomp, TRUE, extdata);
}

/*
 * Interface methods
 */


/** check component for possible elimination */
SCIP_RETCODE extreduce_checkComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   EXTCOMP*              extcomp,            /**< component to be checked (might be reverted) */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            compIsDeletable     /**< is component deletable? */
)
{
   int* extstack_data;
   int* extstack_start;
   int* extstack_state;
   int* tree_edges;
   int* tree_leaves;
   int* tree_innerNodes;
   int* tree_parentNode;
   int* pcSdCands = NULL;
   SCIP_Real* tree_parentEdgeCost;
   SCIP_Real* redcost_treenodeswaps;
   SCIP_Real* redcost_treecosts;
   int* pseudoancestor_mark = NULL;
   SCIP_Bool* redcost_noReversedTree;
   SCIP_Bool* sdeq_edgesIsForbidden;
   const int nnodes = graph->knots;
   const int redcosts_nlevels = redcosts_getNlevels(redcostdata);
   const int maxstacksize = extreduce_getMaxStackSize();
   const int maxncomponents = extreduce_getMaxStackNcomponents(graph);

   assert(!(*compIsDeletable));
   assert(extreduce_extCompFullIsPromising(graph, extpermanent, extcomp) || (extcomp->ncompedges > 1));
   assert(redcosts_nlevels > 0 && maxstacksize > 0 && maxncomponents > 0 && nnodes > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &extstack_data, maxstacksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &extstack_start, maxncomponents + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &extstack_state, maxncomponents + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_edges, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_leaves, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_innerNodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentNode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentEdgeCost, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost_treenodeswaps, nnodes * redcosts_nlevels) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost_treecosts, redcosts_nlevels) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost_noReversedTree, redcosts_nlevels) );

   if( graph_pc_isPc(graph) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pcSdCands, nnodes) );
   }

   if( graph_pseudoAncestorsGetHashArraySize(graph->pseudoancestors) > 0 )
   {
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &pseudoancestor_mark, graph_pseudoAncestorsGetHashArraySize(graph->pseudoancestors)) );
   }
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &sdeq_edgesIsForbidden, graph->edges / 2) );

   {
      const SCIP_Bool* isterm = extpermanent->isterm;

      PCDATA pcdata = { .pcSdToNode = extpermanent->pcSdToNode, .pcSdCands = pcSdCands, .nPcSdCands = -1,
         .pcSdStart = -1, .tree_innerPrize = 0.0 };
      REDDATA reddata = { .contration = extpermanent->contration,
         .dcmst = extpermanent->dcmst, .msts_comp = extpermanent->msts_comp,
         .msts_levelbase = extpermanent->msts_levelbase,
         .sds_horizontal = extpermanent->sds_horizontal, .sds_vertical = extpermanent->sds_vertical,
         .sdsbias_horizontal = extpermanent->sdsbias_horizontal, .sdsbias_vertical = extpermanent->sdsbias_vertical,
         .edgedeleted = extpermanent->edgedeleted, .pseudoancestor_mark = pseudoancestor_mark,
         .nodes_implications = extpermanent->nodes_implications,
         .redcost_treenodeswaps = redcost_treenodeswaps, .redcost_treecosts = redcost_treecosts,
         .redcost_noReversedTree = redcost_noReversedTree,
         .redcost_nlevels = redcosts_nlevels, .redcost_allowEquality = extpermanent->redcostEqualAllow,
         .sdsbias_use = extpermanent->useSdBias };
      EXTDATA extdata = { .extstack_data = extstack_data, .extstack_start = extstack_start,
         .extstack_state = extstack_state, .extstack_ncomponents = 0, .tree_leaves = tree_leaves,
         .tree_edges = tree_edges, .tree_deg = extpermanent->tree_deg, .tree_nleaves = 0,
         .tree_bottleneckDistNode = extpermanent->bottleneckDistNode, .tree_parentNode = tree_parentNode,
         .tree_parentEdgeCost = tree_parentEdgeCost, .tree_cost = 0.0, .ncostupdatestalls = 0,
         .tree_nDelUpArcs = 0, .tree_root = -1, .tree_starcenter = -1, .tree_nedges = 0, .tree_depth = 0,
         .extstack_maxsize = maxstacksize, .extstack_maxncomponents = maxncomponents,
         .pcdata = &pcdata, .redcostdata = redcostdata,
         .sdeq_resetStack = NULL, .sdeq_edgesIsForbidden = sdeq_edgesIsForbidden, .sdeq_hasForbiddenEdges = FALSE,
         .genstar_centeredge = extcomp->genstar_centeredge,
         .tree_innerNodes = tree_innerNodes, .tree_ninnerNodes = 0,
         .tree_maxdepth = extpermanent->tree_maxdepth,
         .tree_maxnleaves = extpermanent->tree_maxnleaves,
         .tree_maxnedges = extpermanent->tree_maxnedges, .node_isterm = isterm, .reddata = &reddata,
         .distdata = extpermanent->distdata_default, .distdata_biased = extpermanent->distdata_biased,
         .mode = extpermanent->mode, .extcomp = extcomp };

#ifdef STP_DEBUG_EXT
      extreduce_extdataCleanArraysDbg(graph, &extdata);
#endif

      extreduce_reddataClean(&reddata);

      if( extreduce_extCompIsPromising(graph, extpermanent, extcomp) )
         extProcessComponent(scip, graph, extcomp, &extdata, compIsDeletable);

      /* also try the other way? */
      if( !(*compIsDeletable) && extcomp->allowReversion )
      {
         extreduce_extCompRevert(graph, extpermanent, extcomp);

         if( extreduce_extCompIsPromising(graph, extpermanent, extcomp)  )
         {
            extProcessComponent(scip, graph, extcomp, &extdata, compIsDeletable);
         }
      }
   }

   assert(extreduce_extPermaIsClean(graph, extpermanent));

   SCIPfreeCleanBufferArray(scip, &sdeq_edgesIsForbidden);
   if( pseudoancestor_mark )
      SCIPfreeCleanBufferArray(scip, &pseudoancestor_mark);

   SCIPfreeBufferArrayNull(scip, &pcSdCands);

   SCIPfreeBufferArray(scip, &redcost_noReversedTree);
   SCIPfreeBufferArray(scip, &redcost_treecosts);
   SCIPfreeBufferArray(scip, &redcost_treenodeswaps);
   SCIPfreeBufferArray(scip, &tree_parentEdgeCost);
   SCIPfreeBufferArray(scip, &tree_parentNode);
   SCIPfreeBufferArray(scip, &tree_innerNodes);
   SCIPfreeBufferArray(scip, &tree_leaves);
   SCIPfreeBufferArray(scip, &tree_edges);
   SCIPfreeBufferArray(scip, &extstack_state);
   SCIPfreeBufferArray(scip, &extstack_start);
   SCIPfreeBufferArray(scip, &extstack_data);

   return SCIP_OKAY;
}
