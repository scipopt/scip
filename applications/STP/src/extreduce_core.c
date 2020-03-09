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
//  #define SCIP_DEBUG
// #define STP_DEBUG_EXT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"

#define EXT_COSTS_RECOMPBOUND 10

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
   const int* const extstack_start = extdata->extstack_start;
   const int* const extstack_data = extdata->extstack_data;
   int* const tree_leaves = extdata->tree_leaves;
   const int stackpos_prev = extStackGetPrevMarked(extdata);
   const int stackstart = extstack_start[stackpos_prev];
   const int stackend = extstack_start[stackpos_prev + 1];
   const int compsize = stackend - stackstart;
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

   for( int i = stackstart; i < stackend; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      assert(extLeafFindPos(extdata, head) >= extdata->tree_nleaves - compsize);
   }
#endif

   compleaves_start = extdata->tree_nleaves - compsize;

   for( int i = stackstart; i < stackend; i++ )
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
   int                   comproot,          /**< root of component that is removed */
   EXTDATA*              extdata            /**< extension data */
)
{
   assert(extdata && extdata->tree_innerNodes);
   assert(comproot >= 0);

   /* update tree leaves array todo might need to be changed for pseudo-elimination */
   if( comproot != extdata->tree_root )
   {
      const SCIP_Bool isPc = (graph->prize != NULL);

      assert(extdata->tree_ninnerNodes >= 1);
      assert(comproot == extdata->tree_innerNodes[extdata->tree_ninnerNodes - 1]);
      assert(isPc == graph_pc_isPc(graph));

      extdata->tree_ninnerNodes--;

      if( isPc )
      {
         extdata->pcdata->tree_innerPrize -= graph->prize[comproot];

         assert(GE(extdata->pcdata->tree_innerPrize, 0.0));
      }
   }
   else
   {
      assert(extdata->tree_ninnerNodes == 0);
   }
}


/** adds edge to tree */
static inline
void extTreeAddEdge(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to be added */
   int                   comproot,           /**< component root */
   EXTDATA*              extdata             /**< extension data */
)
{
   int* const tree_deg = extdata->tree_deg;
   int* const tree_edges = extdata->tree_edges;
   int* const tree_parentNode = extdata->tree_parentNode;
   SCIP_Real* const tree_parentEdgeCost = extdata->tree_parentEdgeCost;
   const SCIP_Real edgecost = graph->cost[edge];
   const int head = graph->head[edge];

   assert(comproot == graph->tail[edge]);
   assert(tree_deg[head] == 0);
   assert(tree_deg[comproot] > 0 || comproot == extdata->tree_root);

   extLeafAdd(head, extdata);

   extdata->tree_cost += edgecost;
   tree_deg[head] = 1;
   tree_edges[(extdata->tree_nedges)++] = edge;
   tree_parentNode[head] = comproot;
   tree_parentEdgeCost[head] = edgecost;
   tree_deg[comproot]++;
}


/** removes root of stack top component from tree */
static inline
void extTreeStackTopRootRemove(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int comproot = extStackGetTopRoot(graph, extdata);

   /* update tree leaves array todo might need to be changed for pseudo-elimination */
   if( comproot != extdata->tree_root )
   {
      extLeafRemove(comproot, extdata);
      extInnerNodeAdd(graph, comproot, extdata);
   }
   else
   {
      assert(extdata->tree_nleaves == 1);
      assert(extdata->tree_deg[comproot] == 1);

      extdata->tree_deg[comproot] = 0;
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
   const int comproot = extStackGetTopRoot(graph, extdata);
   int conflictIteration = -1;
   const SCIP_Bool noReversedRedCostTree = extreduce_redcostReverseTreeRuledOut(graph, extdata);

   assert(!(*conflict));
   assert(stackpos >= 0);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);
   assert(comproot >= 0 && comproot < graph->knots);
   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   extTreeStackTopRootRemove(graph, extdata);

   /* add top expanded component to tree data */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];

      assert(extdata->tree_nedges < extdata->extstack_maxsize);
      assert(edge >= 0 && edge < graph->edges);

      extreduce_redcostAddEdge(graph, edge, noReversedRedCostTree, reddata, extdata);
      extTreeAddEdge(graph, edge, comproot, extdata);

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

  // for( int k = extdata->tree_nleaves - 1; k >= extdata->tree_nleaves - compsize; k-- ) printf("bt remove leaf %d \n", extdata->tree_leaves[k]);

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


#if 0
/** can any extension via edge except for only the edge itself be ruled out? */
static
SCIP_Bool extRuleOutEdgeCombinations(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   extedge             /**< edge to be tested */
)
{
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real tree_redcost = extTreeGetRedcostBound(scip, graph, extdata);
   const SCIP_Real cutoff = reddata->cutoff;
   const int base = graph->tail[extedge];
   const int extvert = graph->head[extedge];
   const SCIP_Real* const redcost = reddata->redCosts;

   assert(extedge >= 0 && extedge < graph->edges);
   assert(extdata->tree_deg[base] == 1 && base != extdata->reddata->redCostRoot);

   // add to leaves?

#ifdef SCIP_DEBUG
        printf("edge combinations ruled out: ");
        graph_edge_printInfo(graph, e);
#endif

   if( extvert == extdata->reddata->redCostRoot )
   {
      assert(Is_term(graph->term[extvert]));

      tree_redcost = FARAWAY;
   }
   else
      tree_redcost += redcost[extedge] + nodeTo3TermsPaths[extvert].dist - nodeTo3TermsPaths[base].dist;

   if( reddata->edgedeleted != NULL && reddata->edgedeleted[flipedge(extedge)] )
      checkReverseTrees = FALSE;

}
#endif


/** Can current tree be peripherally ruled out?
 *  NOTE: If tree cannot be ruled-out, the current component will be put into the MST storage 'reddata->msts' */
static inline
SCIP_Bool extTreeRuleOutPeriph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   if( extreduce_redcostRuleOutPeriph(graph, extdata) )
      return TRUE;

   if( extreduce_mstRuleOutPeriph(scip, graph, extdata) )
      return TRUE;

   return FALSE;
}


/** Stores extensions of tree from current (expanded and marked) stack top that cannot be ruled-out.
 *  Also computes SDs from leaves of any extension to the tree leaves and other extension leaves. */
static
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
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extStackGetPosition(extdata);
   const SCIP_Bool* const isterm = extdata->node_isterm;

#ifndef NDEBUG
   assert(graph && extdata && extedgesstart && extedges && nextensions && nextendableleaves && with_ruledout_leaf);
   assert(stackpos >= 0);
   assert(EXT_STATE_MARKED == extdata->extstack_state[stackpos]);
   assert(!(*with_ruledout_leaf));
   assert(*nextensions == 0 && *nextendableleaves == 0);

   extreduce_extendInitDebug(extedgesstart, extedges);
#endif

   extedgesstart[0] = 0;

   /* loop over all leaves of top component */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      int nleafextensions = 0;
      const int leaf = graph->head[extstack_data[i]];

      assert(extstack_data[i] >= 0 && extstack_data[i] < graph->edges);

      /* extensions from leaf not possible? */
      if( !extLeafIsExtendable(graph, isterm, leaf) )
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
   if( extdata->extstack_state[stackposition] == EXT_STATE_EXPANDED )
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
   SCIP_Bool*            success             /**< success pointer */
)
{
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extStackGetPosition(extdata);
   int datasize = extstack_start[stackpos];
   const uint32_t powsize = (uint32_t) pow(2.0, nextedges);

   assert(nextedges > 0 && nextedges < 32);
   assert(powsize >= 2);

   /* stack too full? */
   if( !extStackIsExtendable(nextedges, (int) powsize, datasize, extdata) )
   {
      *success = FALSE;
      extBacktrack(scip, graph, *success, FALSE, extdata);

      return;
   }

   /* todo try to rule out pairs of edges with simple, 2-edge bottleneck */
   // int* antipairs_starts;
   // int* antipairs_edges;
   // int* anitpairs_hasharr; size nedges, clean

   extreduce_mstLevelHorizontalAdd(scip, graph, nextedges, extedges, extdata);

   extreduce_mstLevelClose(scip, graph, graph->tail[extedges[0]], extdata);

   /* compute and add components (overwrite previous, non-expanded component) */
   // todo we probably want to order so that the smallest components are put last!
   // todo we might just ignore the single edge extensions if counter & (counter - 1) == 0, then add them at the end
   // todo extra method!
   // todo if single extensions are used first, might make sense to also have a special bottleneck test for this case!
   // also good if we have the method that excludes everything except for the singletons!
   for( uint32_t counter = powsize - 1; counter >= 1; counter-- )
   {
      for( uint32_t j = 0; j < (uint32_t) nextedges; j++ )
      {
         /* Check if jth bit in counter is set */
         if( counter & ((uint32_t) 1 << j) )
         {
            assert(datasize < extdata->extstack_maxsize);
            assert(extedges[j] >= 0);

            extstack_data[datasize++] = extedges[j];
            SCIPdebugMessage(" head %d \n", graph->head[extedges[j]]);
         }
      }

      SCIPdebugMessage("... added \n");
      assert(stackpos < extdata->extstack_maxsize - 1);

      // todo check with hashing whether antipairs extis. If so, remove again

      extstack_state[stackpos] = EXT_STATE_EXPANDED;
      extstack_start[++stackpos] = datasize;

      assert(extstack_start[stackpos] - extstack_start[stackpos - 1] > 0);
   }

   assert(stackpos > extStackGetPosition(extdata));
   assert(stackpos >= extdata->extstack_ncomponents);
   assert(stackpos <= extdata->extstack_maxncomponents);

   extdata->extstack_ncomponents = stackpos;
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
      SCIP_Bool ruledOut = FALSE;

      assert(*nextedges < STP_EXT_MAXGRAD);
      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[graph->head[edge]] == 0);
      assert(!extTreeRuleOutEdgeSimple(graph, extdata, edge));

      /* computes the SDs from 'leaf' to all tree leaves in 'sds_vertical', unless the edge is ruled out */
      extreduce_mstLevelVerticalAddLeaf(scip, graph, edge, extdata, &ruledOut);

      if( ruledOut )
      {
         continue;
      }

#if 0
     if( extRuleOutEdgeCombinations(graph, extdata, e) )
     {
        // todo: need some marker here to say that we have a single edge!
        continue;
     }
#endif

      extedges[(*nextedges)++] = edge;
   }
}


/** Expands top component of stack.
 *  I.e.. adds all possible subsets of the top component that cannot be ruled-out.
 *  Note: This method can backtrack:
 *        1. If stack is full (with success set to FALSE),
 *        2. If all edges of the component can be ruled-out (with success set to TRUE).
 *  Note: This method computes SDs for newly added leaves! */
static
void extStackTopExpand(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD];
   REDDATA* const reddata = extdata->reddata;
   int nextedges = 0;
   const int stackpos = extStackGetPosition(extdata);

#ifndef NDEBUG
   const int* const extstack_state = extdata->extstack_state;
   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;
#endif

   assert(scip && graph && success);
   assert(EXT_STATE_NONE == extstack_state[stackpos]);

   extreduce_mstLevelInit(reddata, extdata);

   /* Note: Also computes ancestor SDs for leaves that are not ruled-out
    * and adds them to vertical level! */
   extStackTopCollectExtEdges(scip, graph, extdata, extedges, &nextedges);

   extreduce_mstLevelVerticalClose(reddata);

   /* everything ruled out already? */
   if( nextedges == 0 )
   {
      *success = TRUE;

      assert(extstack_state[stackpos] == EXT_STATE_NONE);

      /* not the initial component? */
      if( stackpos != 0 )
      {
         extBacktrack(scip, graph, *success, FALSE, extdata);
      }
      else
      {
         extreduce_mstLevelVerticalRemove(reddata);
      }
   }
   else
   {
      /* use the just collected edges 'extedges' to build components and add them to the stack */
      extStackAddCompsExpanded(scip, graph, nextedges, extedges, extdata, success);

#ifndef NDEBUG
      {
         const int stackpos_new = extStackGetPosition(extdata);
         assert(extstack_state[stackpos_new] == EXT_STATE_EXPANDED || (stackpos_new < stackpos) );
      }
#endif
   }
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


/** Adds extensions initial component to stack (needs to be star component rooted in root).
  * If no extensions are added, then the component has been ruled-out. */
static
void extProcessInitialComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            ruledOut,           /**< initial component ruled out? */
   SCIP_Bool*            success             /**< extension successful? */
)
{
   const int* const compedges = extcomp->compedges;
   const int ncompedges = extcomp->ncompedges;
   const int comproot = extcomp->comproot;
   SCIP_Bool conflict;

   assert(compedges);
   assert(ncompedges >= 1 && ncompedges < STP_EXT_MAXGRAD);
   assert(ncompedges < extdata->extstack_maxsize);
   assert(comproot >= 0 && comproot < graph->knots);
   assert(FALSE == (*ruledOut) && TRUE == (*success));

#ifdef SCIP_DEBUG
   printf("\n --- ADD initial component --- \n\n");
#endif

   for( int i = 0; i < ncompedges; i++ )
   {
      const int e = compedges[i];
      const int tail = graph->tail[e];

      assert(e >= 0 && e < graph->edges);
      assert(tail == comproot);

      SCIPdebugMessage("edge %d: %d->%d \n", e, graph->tail[e], graph->head[e]);

      extdata->extstack_data[i] = e;
      extLeafAdd(tail, extdata);
      extreduce_mstAddRootLevel(scip, tail, extdata);
   }

   extdata->tree_root = comproot;
   extdata->extstack_ncomponents = 1;
   extdata->extstack_state[0] = EXT_STATE_NONE;
   extdata->extstack_start[0] = 0;
   extdata->extstack_start[1] = ncompedges;
   extdata->tree_parentNode[comproot] = -1;
   extdata->tree_redcostSwap[comproot] = 0.0;
   extdata->tree_parentEdgeCost[comproot] = -1.0;

   assert(ncompedges > 1 || extdata->tree_leaves[0] == comproot);
   assert(extdata->tree_deg[comproot] == 0);
   assert(extdata->tree_nleaves == ncompedges);

   extdata->tree_deg[comproot] = 1;

   /* expand the single edge */
   extStackTopExpand(scip, graph, extdata, success);

   assert(*success);
   assert(extStackGetPosition(extdata) == 0);

   /* early rule-out? */
   if( extdata->extstack_state[0] == EXT_STATE_NONE )
   {
      *ruledOut = TRUE;

      /* necessary because this edge will be deleted in clean-up otherwise */
      graph_pseudoAncestors_hashEdge(graph->pseudoancestors, compedges[0], extdata->reddata->pseudoancestor_mark);

      return;
   }

   assert(extdata->extstack_state[0] == EXT_STATE_EXPANDED);

   conflict = FALSE;
   extTreeStackTopAdd(scip, graph, extdata, &conflict);

   assert(!conflict);

   /* NOTE: necessary to keep the MST graph up-to-date */
   if( extTreeRuleOutPeriph(scip, graph, extdata) )
   {
      *ruledOut = TRUE;
      return;
   }

   /* the single edge component could not be ruled-out, so set its stage to 'marked' */
   extdata->extstack_state[extStackGetPosition(extdata)] = EXT_STATE_MARKED;

   extExtend(scip, graph, extdata, success);

   assert(success || ncompedges >= 3);
}


/** Checks whether component can be ruled out. */
static
void extProcessComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            deletable           /**< is arc deletable? */
)
{
   int* const extstack_state = extdata->extstack_state;
   SCIP_Bool success = TRUE;
   SCIP_Bool conflict = FALSE;

   assert(extreduce_extdataIsClean(graph, extdata));
   assert(extreduce_reddataIsClean(graph, extdata->reddata));
   assert(extreduce_pcdataIsClean(graph, extdata->pcdata));

   assert(!(*deletable));

   /* put 'extcomp' on the stack */
   extProcessInitialComponent(scip, graph, extcomp, extdata, deletable, &success);

   /* early rule-out? or no success? */
   if( *deletable || !success )
   {
      extreduce_extCompClean(scip, graph, extcomp, extdata);
      return;
   }

   assert(extstack_state[0] == EXT_STATE_MARKED);

   /* limited DFS backtracking; stops once back at 'edge' */
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

      /* component not expanded yet? */
      if( extstack_state[stackposition] != EXT_STATE_EXPANDED )
      {
         assert(extstack_state[stackposition] == EXT_STATE_NONE);

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

   extreduce_extCompClean(scip, graph, extcomp, extdata);
}


/** check (directed) arc */
SCIP_RETCODE extreduce_checkComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   EXTCOMP*              extcomp,            /**< component to be checked (might be reverted) */
   DISTDATA*             distdata,           /**< data for distance computations */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            compIsDeletable     /**< is component deletable? */
)
{
   const SCIP_Bool* isterm = extpermanent->isterm;
   const int root = redcostdata->redCostRoot;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real* rootdist = redcostdata->rootToNodeDist;
   const PATH* nodeToTermpaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real cutoff = redcostdata->cutoff;
   int* extstack_data;
   int* extstack_start;
   int* extstack_state;
   int* tree_edges;
   int* tree_leaves;
   int* tree_innerNodes;
   int* tree_parentNode;
   int* pcSdCands = NULL;
   SCIP_Real* tree_parentEdgeCost;
   SCIP_Real* tree_redcostSwap;
   int* pseudoancestor_mark;
   const int nnodes = graph->knots;
   const int maxstacksize = extreduce_getMaxStackSize();
   const int maxncomponents = extreduce_getMaxStackNcomponents(graph);
   const SCIP_Bool isPseudoElim = (extcomp->ncompedges > 1);

   SCIP_CALL( SCIPallocBufferArray(scip, &extstack_data, maxstacksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &extstack_start, maxncomponents + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &extstack_state, maxncomponents + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_edges, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_leaves, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_innerNodes, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentNode, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentEdgeCost, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_redcostSwap, nnodes) );
   if( graph_pc_isPc(graph) )
      SCIP_CALL( SCIPallocBufferArray(scip, &pcSdCands, nnodes) );

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &pseudoancestor_mark, nnodes) );

   /* is the component (or the reverted one) promising? Or are we trying pseudo elimination? */
   if( extreduce_extCompFullIsPromising(graph, extpermanent, extcomp) || isPseudoElim )
   {
      PCDATA pcdata = { .pcSdToNode = extpermanent->pcSdToNode, .pcSdCands = pcSdCands, .nPcSdCands = -1, .pcSdStart = -1,
         .tree_innerPrize = 0.0 };
      REDDATA reddata = { .dcmst = extpermanent->dcmst, .msts_comp = extpermanent->msts_comp,
         .msts_levelbase = extpermanent->msts_levelbase,
         .sds_horizontal = extpermanent->sds_horizontal, .sds_vertical = extpermanent->sds_vertical,
         .redCosts = redcost, .rootToNodeDist = rootdist, .nodeTo3TermsPaths = nodeToTermpaths,
         .nodeTo3TermsBases = redcostdata->nodeTo3TermsBases, .edgedeleted = extpermanent->edgedeleted,
         .pseudoancestor_mark = pseudoancestor_mark, .cutoff = cutoff, .equality = extpermanent->redcostEqualAllow, .redCostRoot = root };
      EXTDATA extdata = { .extstack_data = extstack_data, .extstack_start = extstack_start,
         .extstack_state = extstack_state, .extstack_ncomponents = 0, .tree_leaves = tree_leaves,
         .tree_edges = tree_edges, .tree_deg = extpermanent->tree_deg, .tree_nleaves = 0,
         .tree_bottleneckDistNode = extpermanent->bottleneckDistNode, .tree_parentNode = tree_parentNode,
         .tree_parentEdgeCost = tree_parentEdgeCost, .tree_redcostSwap = tree_redcostSwap,
         .tree_cost = 0.0, .tree_redcost = 0.0, .ncostupdatestalls = 0,
         .tree_nDelUpArcs = 0, .tree_root = -1, .tree_nedges = 0, .tree_depth = 0,
         .extstack_maxsize = maxstacksize, .extstack_maxncomponents = maxncomponents,
         .pcdata = &pcdata,
         .tree_innerNodes = tree_innerNodes, .tree_ninnerNodes = 0,
         .tree_maxdepth = extreduce_getMaxTreeDepth(graph),
         .tree_maxnleaves = STP_EXTTREE_MAXNLEAVES,
         .tree_maxnedges = STP_EXTTREE_MAXNEDGES, .node_isterm = isterm, .reddata = &reddata, .distdata = distdata };

      extProcessComponent(scip, graph, extcomp, &extdata, compIsDeletable);

      if( !(*compIsDeletable) )
      {
         extreduce_extCompRevert(graph, extpermanent, extcomp);

         if( extreduce_extCompIsPromising(graph, extpermanent, extcomp)  )
         {
            extProcessComponent(scip, graph, extcomp, &extdata, compIsDeletable);
         }
      }
   }

   assert(extreduce_extPermaIsClean(graph, extpermanent));

   SCIPfreeCleanBufferArray(scip, &pseudoancestor_mark);
   SCIPfreeBufferArrayNull(scip, &pcSdCands);
   SCIPfreeBufferArray(scip, &tree_redcostSwap);
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
