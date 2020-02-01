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

/**@file   extreduce_extmst.c
 * @brief  extended-reduction specific MST algorithms for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements MST algorithms for extended reduction techniques for Steiner problems.
 * Allows to efficiently compute and store special distance MSTs between the leaves of extension tree.
 * Furthermore, can check for tree bottlenecks.
 *
 * A 'level' of the extension tree consists of all possible extension edges from the leaf used for extension.
 * For each level there are a number of 'components': all the subsets that were not already ruled-out.
 * Once a level is initiated, all SDs to the other leaves of the tree are computed ('vertical'),
 * as well as the SDs among the level ('horizontal').
 * These SDs are kept until the level has been removed again.
 * Furthermore, for each level we store the MST corresponding to the extension tree without the level.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


//#define SCIP_DEBUG
//#define STP_DEBUG_EXT

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"

#define EXT_PC_SDMAXVISITS 10  /**< maximum visits for PC specific SD computation */

/** returns special distance */
static inline
SCIP_Real extGetSD(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Real sd = extreduce_distDataGetSD(scip, g, vertex1, vertex2, extdata->distdata);

   assert((extdata->pcSdToNode != NULL) == graph_pc_isPcMw(g));

   if( extdata->pcSdToNode )
   {
      const SCIP_Real sdpc = extdata->pcSdToNode[vertex2];

      assert(SCIPisEQ(scip, sdpc, -1.0) || SCIPisGE(scip, sdpc, 0.0));

      if( sdpc > -0.5 && (sdpc < sd || sd < -0.5) )
      {
         SCIPdebugMessage("special distance update for pc: %f to %f \n", sd, sdpc);
         sd = sdpc;
      }
   }

   assert(SCIPisEQ(scip, sd, -1.0) || SCIPisGE(scip, sd, 0.0));

   return sd;
}


/** returns size of top component on the stack */
static inline
int extStackGetTopSize(
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int stackpos = extStackGetPosition(extdata);
   const int* const stack_start = extdata->extstack_start;
   const int size = stack_start[stackpos + 1] - stack_start[stackpos];

   assert(extdata->extstack_state[stackpos] != EXT_STATE_NONE);
   assert(size > 0 && size < STP_EXT_MAXGRAD);

   return size;
}


/** returns number of ancestor leaves (i.e. number of leaves below current level) */
static inline
int extGetNancestorLeaves(
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int compsize = extStackGetTopSize(extdata);
   const int nleaves = extdata->tree_nleaves;
   const int nleaves_ancestors = nleaves - compsize;

   assert(nleaves_ancestors > 0 && nleaves_ancestors < nleaves);

   return nleaves_ancestors;
}


/** is given SD non-trivial? */
static inline
SCIP_Real sdIsNonTrivial(
   SCIP_Real             specialDist         /**< SD */
  )
{
   assert(specialDist >= 0 || EQ(specialDist, -1.0));
   assert(LT(specialDist, FARAWAY));

   return (specialDist >= -0.5);
}


/** marks single PcSd array entry */
static inline
void pcSdMarkSingle(
   const GRAPH*          graph,              /**< graph data structure */
   int                   entry,              /**< entry to mark */
   SCIP_Real             value,              /**< value to mark with */
   SCIP_Real*            pcSdToNode,         /**< node mark array */
   int*                  pcSdCands,          /**< marked candidates list */
   int*                  nPcSdCands          /**< pointer to store number of candidates */
)
{
   /* entry not marked yet? */
   if( pcSdToNode[entry] < -0.5 )
   {
      assert(EQ(pcSdToNode[entry], -1.0));
      assert(*nPcSdCands < graph->knots);

      pcSdCands[(*nPcSdCands)++] = entry;
      pcSdToNode[entry] = value;
   }
   else if( value < pcSdToNode[entry] )
   {
      pcSdToNode[entry] = value;
   }

   assert(GE(pcSdToNode[entry], 0.0));
}


/** marks PcSd array */
static
void pcSdToNodeMark(
   const GRAPH*          graph,              /**< graph data structure */
   int                   startvertex,        /**< vertex to start from */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const pcSdToNode = extdata->pcSdToNode;
   int* const pcSdCands = extdata->pcSdCands;
   const DCSR* const dcsr = graph->dcsr_storage;
   const RANGE* const range_csr = dcsr->range;
   const int* const head_csr = dcsr->head;
   const SCIP_Real* const cost_csr = dcsr->cost;
   const SCIP_Real* const prize = graph->prize;
   const int* const tree_deg = extdata->tree_deg;
   const int start = range_csr[startvertex].start;
   const int end = range_csr[startvertex].end;
   int count1 = 0;
   int count2 = 0;

   assert(graph_pc_isPcMw(graph));
   assert(pcSdCands && pcSdToNode && prize);
   assert(extdata->nPcSdCands == -1);

   extdata->nPcSdCands = 0;

   for( int i = start; i != end; i++ )
   {
      const SCIP_Real edgecost = cost_csr[i];
      const int head = head_csr[i];
      assert(tree_deg[head] >= 0);

      if( tree_deg[head] == 0 )
      {
         const int start2 = range_csr[head].start;
         const int end2 = range_csr[head].end;

         for( int i2 = start2; i2 != end2; i2++ )
         {
            const int head2 = head_csr[i2];
            assert(tree_deg[head2] >= 0);

            /* tree reached? */
            if( tree_deg[head2] > 0 && head2 != startvertex )
            {
               const SCIP_Real edgecost2 = cost_csr[i2];
               const SCIP_Real maxedgecost = MAX(edgecost, edgecost2);
               SCIP_Real dist2 = MAX(maxedgecost, edgecost + edgecost2 - prize[head]);

               assert(0.0 == prize[head] || Is_term(graph->term[head]));

               pcSdMarkSingle(graph, head2, dist2, pcSdToNode, pcSdCands, &(extdata->nPcSdCands));
            }

            if( count2++ > EXT_PC_SDMAXVISITS )
               break;
         }
      }
      else
      {
         assert(head != startvertex);
         pcSdMarkSingle(graph, head, edgecost, pcSdToNode, pcSdCands, &(extdata->nPcSdCands));
      }

      if( count1++ > EXT_PC_SDMAXVISITS )
         break;
   }
}


/** unmarks PcSd array */
static inline
void pcSdToNodeUnmark(
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const pcSdToNode = extdata->pcSdToNode;
   const int* const pcSdCands = extdata->pcSdCands;
   const int nPcSdCands = extdata->nPcSdCands;

   assert(graph_pc_isPcMw(graph));
   assert(pcSdCands && pcSdToNode);
   assert(nPcSdCands >= 0);

   for( int i = 0; i < nPcSdCands; i++ )
   {
      const int cand = pcSdCands[i];

      assert(pcSdToNode[cand] >= 0.0);

      pcSdToNode[cand] = -1.0;
   }

#ifndef NDEBUG
   extdata->nPcSdCands = -1;
#endif
}


/** marks bottleneck array on path to tree root */
static
void bottleneckMarkRootPath(
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex,             /**< vertex to start from */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const SCIP_Real* const parentEdgeCost = extdata->tree_parentEdgeCost;
   const int* const parentNode = extdata->tree_parentNode;
   const int* const tree_deg = extdata->tree_deg;
   const int tree_root = extdata->tree_root;

   assert(bottleneckDist_node && parentEdgeCost && parentNode && tree_deg);
   assert(vertex >= 0 && vertex < graph->knots);
   assert(bottleneckDist_node[vertex] == -1.0);
   assert(bottleneckDist_node[tree_root] == -1.0);


   if( vertex == tree_root )
   {
      bottleneckDist_node[vertex] = 0.0;
   }
   else
   {
      /* go down from vertex */

      SCIP_Real bottleneck = 0.0;
      SCIP_Real bottleneck_local = 0.0;
      int childNode = vertex;
      int currentNode = parentNode[vertex];
      const SCIP_Bool isPc = graph_pc_isPc(graph);

      assert(currentNode != -1);
      assert(tree_deg[childNode] == 1);

      while( currentNode != -1 )
      {
         assert(currentNode >= 0 && tree_deg[currentNode] >= 0);
         assert(parentEdgeCost[childNode] >= 0.0 && bottleneckDist_node[currentNode] == -1.0);
         assert(currentNode != vertex);
         assert(!isPc || !graph_pc_knotIsDummyTerm(graph, currentNode));

         if( tree_deg[childNode] == 2 )
         {
            bottleneck_local += parentEdgeCost[childNode];
            if( isPc && Is_term(graph->term[childNode]) )
            {
               assert(graph_pc_termIsNonLeafTerm(graph, childNode) && graph->prize[childNode] > 0.0);
               bottleneck_local -= graph->prize[childNode];
            }
         }
         else
            bottleneck_local = parentEdgeCost[childNode];

         if( bottleneck < bottleneck_local )
            bottleneck = bottleneck_local;

         bottleneckDist_node[currentNode] = bottleneck;
         childNode = currentNode;
         currentNode = parentNode[currentNode];
      }

      assert(childNode == tree_root);
   }
}

/** unmarks bottleneck array on path to tree root */
static
void bottleneckUnmarkRootPath(
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex,             /**< vertex to start from */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const int* const parentNode = extdata->tree_parentNode;
   const int tree_root = extdata->tree_root;

   assert(extdata && bottleneckDist_node && parentNode);
   assert(bottleneckDist_node[vertex] == -1.0 || vertex == tree_root);
   assert(bottleneckDist_node[tree_root] >= 0.0);

   if( vertex == tree_root )
   {
      bottleneckDist_node[vertex] = -1.0;
      assert(parentNode[vertex] == -1);
   }
   else
   {
      assert(parentNode[vertex] >= 0);
   }

   /* go down from vertex and reset bottleneckDist_node */
   for( int currentNode = parentNode[vertex]; currentNode != -1; currentNode = parentNode[currentNode]  )
   {
      assert(currentNode >= 0);
      assert(extdata->tree_deg[currentNode] >= 0);
      assert(bottleneckDist_node[currentNode] >= 0.0);

      bottleneckDist_node[currentNode] = -1.0;
   }

   assert(bottleneckDist_node[tree_root] == -1.0);
}


/** computes the tree bottleneck between vertices in the current tree,
 * for which vertex_pathmarked root path has been marked already */
static
SCIP_Real bottleneckGetDist(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
#ifndef NDEBUG
   int                   vertex_pathmarked,  /**< vertex with marked rootpath */
#endif
   int                   vertex_unmarked     /**< second vertex */
   )
{
   const SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const SCIP_Real* const parentEdgeCost = extdata->tree_parentEdgeCost;
   const int* const parentNode = extdata->tree_parentNode;
   int* const tree_deg = extdata->tree_deg;
   SCIP_Real bottleneck;
   const int tree_root = extdata->tree_root;
   int currentNode;

   assert(bottleneckDist_node && parentEdgeCost && parentNode);
   assert(bottleneckDist_node[vertex_pathmarked] == -1.0 || vertex_pathmarked == tree_root);
   assert(bottleneckDist_node[vertex_unmarked] == -1.0 || vertex_unmarked == tree_root || tree_deg[vertex_unmarked] > 1);
   assert(bottleneckDist_node[tree_root] >= 0.0);
   assert(vertex_pathmarked != vertex_unmarked);

   /* go down from vertex_unmarked up to lowest common ancestor with vertex_pathmarked  */
   bottleneck = 0.0;

   if( vertex_unmarked == tree_root )
   {
      currentNode = vertex_unmarked;
   }
   else
   {
      SCIP_Real bottleneck_local = 0.0;
      const SCIP_Bool isPc = graph_pc_isPc(graph);

      assert(parentNode[vertex_unmarked] >= 0);

      for( currentNode = vertex_unmarked; bottleneckDist_node[currentNode] < -0.5; currentNode = parentNode[currentNode] )
      {
         assert(tree_deg[currentNode] >= 0 && parentEdgeCost[currentNode] >= 0.0);
         assert(bottleneckDist_node[currentNode] == -1.0);
         assert(currentNode != vertex_pathmarked);

         if( tree_deg[currentNode] == 2 )
         {
            bottleneck_local += parentEdgeCost[currentNode];
            if( isPc && Is_term(graph->term[currentNode]) )
            {
               assert(graph_pc_termIsNonLeafTerm(graph, currentNode) && graph->prize[currentNode] > 0.0);
               bottleneck_local -= graph->prize[currentNode];
            }
         }
         else
            bottleneck_local = parentEdgeCost[currentNode];

         if( bottleneck < bottleneck_local )
            bottleneck = bottleneck_local;

         assert(parentNode[currentNode] >= 0 && parentNode[currentNode] != vertex_unmarked);
      }
   }

   bottleneck = MAX(bottleneck, bottleneckDist_node[currentNode]);

   return bottleneck;
}


/** does a special distance approximation dominate the tree bottleneck distance between
 *  vertex_pathmarked and vertex_unmarked in the current tree? */
static inline
SCIP_Bool bottleneckIsDominated(
   const GRAPH*          graph,              /**< graph data structure */
   int                   extedge,            /**< edge along which we want to extend the tree, or -1 */
   int                   vertex_pathmarked,  /**< vertex for which bottleneck path to root has been marked */
   int                   vertex_unmarked,    /**< second vertex */
   SCIP_Real             specialDist,        /**< best computed special distance approximation (-1.0 if unknown) */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real bottleneckDist;
   const SCIP_Bool hasSpecialDist = sdIsNonTrivial(specialDist);

   assert(vertex_pathmarked >= 0 && vertex_pathmarked < graph->knots);
   assert(vertex_unmarked >= 0 && vertex_unmarked < graph->knots);
   assert(extedge == -1 || vertex_pathmarked == graph->tail[extedge]);

   if( !hasSpecialDist )
   {
      return FALSE;
   }

   assert(GE(specialDist, 0.0));

   if( extedge >= 0 )
   {
      if( LT(specialDist, graph->cost[extedge]) )
         return TRUE;
   }

   if( vertex_pathmarked == vertex_unmarked )
   {
      return FALSE;
   }

#ifndef NDEBUG
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_pathmarked, vertex_unmarked);
#else
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_unmarked);
#endif

   SCIPdebugMessage("%d->%d: sd=%f bottleneck=%f \n", vertex_pathmarked, vertex_unmarked, specialDist, bottleneckDist);

   if( LT(specialDist, bottleneckDist) )
      return TRUE;
   else if( LE(specialDist, bottleneckDist) && 0 ) /* todo cover equality */
      return TRUE;

   return FALSE;
}


/** does a special distance approximation dominate the tree bottleneck distance between
 *  vertex_pathmarked and vertex_unmarked in the current tree? */
static inline
SCIP_Bool bottleneckToSiblingIsDominated(
   const GRAPH*          graph,              /**< graph data structure */
   int                   extedge,            /**< edge for extension */
   int                   edge2sibling,       /**< edge to sibling of extedge head */
   SCIP_Real             specialDist         /**< best computed special distance approximation (-1.0 if unknown) */
   )
{
   const SCIP_Bool hasSpecialDist = sdIsNonTrivial(specialDist);

   assert(extedge >= 0 && edge2sibling >= 0);
   assert(extedge != edge2sibling);
   assert(graph->tail[extedge] == graph->tail[edge2sibling]);

   if( !hasSpecialDist )
   {
      return FALSE;
   }
   else
   {
      const SCIP_Real* const edgecost = graph->cost;

      assert(GE(specialDist, 0.0));

      if( LT(specialDist, edgecost[edge2sibling]) )
         return TRUE;

      if( LT(specialDist, edgecost[extedge]) )
         return TRUE;

      if( 0 && LE(specialDist, edgecost[edge2sibling]) ) // todo cover equality!
         return TRUE;

      if( 0 && LE(specialDist, edgecost[extedge]) ) // todo cover equality!
         return TRUE;
   }

   return FALSE;
}


// todo check here always the bottleneck distances to some or all internal tree
// nodes. Apart from the base! Need to keep them in a list similar to the leaves!

/** checks tree bottleneck distances to non-leaves of the tree */
static inline
void bottleneckCheckNonLeaves(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2neighbor,      /**< the edge from the tree to the neighbor */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            ruledOut            /**< could the extension be ruled out */
)
{
   const int* const pcSdCands = extdata->pcSdCands;
   const int* const tree_deg = extdata->tree_deg;
   const int nPcSdCands = extdata->nPcSdCands;
   const int neighbor = graph->head[edge2neighbor];
   const int neighbor_base = graph->tail[edge2neighbor];

   assert(pcSdCands);
   assert(ruledOut);
   assert(!(*ruledOut));
   assert(nPcSdCands >= 0);

   /* also check non-leaves */
   for( int c = 0; c < nPcSdCands; c++ )
   {
      SCIP_Real specialDist;
      const int cand = pcSdCands[c];

      assert(cand >= 0 && cand < graph->knots);

      /* leaf or not contained? */
      if( tree_deg[cand] <= 1 )
         continue;

      specialDist = extGetSD(scip, graph, neighbor, cand, extdata);

      if( bottleneckIsDominated(graph, edge2neighbor, neighbor_base, cand, specialDist, extdata) )
      {
         SCIPdebugMessage("---non-leaf bottleneck rule-out---\n");
         *ruledOut = TRUE;

         return;
      }
   }
}


#ifndef NDEBUG
/** has the leaf a dominated bottleneck with other leaves? */
static
SCIP_Bool dbgBottleneckFromLeafIsDominated(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   topleaf,            /**< component leaf to check for */
   EXTDATA*              extdata             /**< extension data */
   )
{
   const int* const leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Bool ruleOut = FALSE;

   bottleneckMarkRootPath(graph, topleaf, extdata);

   for( int j = 0; j < nleaves; j++ )
   {
      const int leaf = leaves[j];

      if( leaf != topleaf )
      {
         const SCIP_Real specialDist = extreduce_extGetSD(scip, graph, topleaf, leaf, extdata);

         if( bottleneckIsDominated(graph, -1, topleaf, leaf, specialDist, extdata) )
         {
            ruleOut = TRUE;
            break;
         }
      }
   }

   bottleneckUnmarkRootPath(graph, topleaf, extdata);

   return ruleOut;
}

#endif


/** Gets SDs from leaf of top tree component to siblings for MST calculation.
 *  Returns early (with leafRuledOut == TRUE) if extension via 'edge2leaf' can be ruled out already.
 *  NOTE: Only restricted bottleneck tests are performed! */
static inline
void mstCompLeafGetSDsToSiblings(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2top,           /**< edge to the top component leaf */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Real             sds[],              /**< array to store the SDs */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
   )
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int* const ghead = graph->head;
   const int stackpos = extStackGetPosition(extdata);
   const int topleaf = ghead[edge2top];

   SCIP_Bool hitTopLeaf = FALSE;

   assert(!(*leafRuledOut));

   for( int i = extstack_start[stackpos], j = 0; i < extstack_start[stackpos + 1]; i++, j++ )
   {
      const int edge2sibling = extstack_data[i];
      const int sibling = ghead[edge2sibling];
      SCIP_Real specialDist;

      assert(extreduce_nodeIsInStackTop(graph, extdata, sibling));
      assert(extdata->tree_deg[sibling] == 1);
      assert(graph->tail[edge2top] == graph->tail[edge2sibling]);
      assert(EQ(sds[j], -1.0));

      if( sibling == topleaf )
      {
         hitTopLeaf = TRUE;
         sds[j] = FARAWAY;

         continue;
      }

      /* only make bottleneck test for 'left' siblings to avoid double checks */
      if( hitTopLeaf )
      {
         //continue; // todo
      }

      /* todo here get the adjcosts from horizontal storage */

      specialDist = extGetSD(scip, graph, topleaf, graph->head[edge2sibling], extdata);

      sds[j] = specialDist; // todo

      if( bottleneckToSiblingIsDominated(graph, edge2top, edge2sibling, specialDist) )
      {
         SCIPdebugMessage("---bottleneck rule-out component (siblings test)---\n");
         *leafRuledOut = TRUE;
         break;
      }
   }

   assert(hitTopLeaf || *leafRuledOut);
}


/** Gets SDs from leaf of top tree component to ancestors for MST calculation.
 *  Returns early (with leafRuledOut == TRUE) if extension via 'edge2leaf' can be ruled out already.
 *  NOTE: Only restricted bottleneck tests are performed, UNLESS the leaf has no siblings! */
static inline
void mstCompLeafGetSDsToAncestors(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2leaf,          /**< edge to the top component leaf */
   int                   nleaves_ancestors,  /**< number of leaves to ancestors */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Real             sds[],              /**< array to store the SDs */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
   )
{
   const MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const int* const leaves = extdata->tree_leaves;
   const int topleaf = graph->head[edge2leaf];
   const SCIP_Real* const adjedgecosts = extreduce_mldistsTopTargetDists(sds_vertical, topleaf);
   const SCIP_Bool hasSiblings = (extStackGetTopSize(extdata) > 1);
   const SCIP_Bool isPc = graph_pc_isPc(graph);
#ifndef NDEBUG
   const int* const adjids = extreduce_mldistsTopTargetIds(sds_vertical, topleaf);
#endif

   assert(adjedgecosts);
   assert(!(*leafRuledOut));
   assert(extreduce_mldistsLevelNTopTargets(sds_vertical) == nleaves_ancestors);

   /* expensive check, maybe only do if STP_DEBUG_EXT is set? */
   assert(extreduce_sdsverticalInSync(scip, graph, extStackGetTopSize(extdata), topleaf, extdata));

   /* if there are no siblings, then there is a chance to find a non-trivial bottleneck rule-out */
   if( !hasSiblings )
   {
      bottleneckMarkRootPath(graph, topleaf, extdata);

      if( isPc )
         pcSdToNodeMark(graph, topleaf, extdata);
   }

   /* get the SDs to the ancestor (lower) leafs (and try bottleneck rule out if there are no siblings) */
   for( int j = 0; j < nleaves_ancestors; j++ )
   {
      const int leaf = leaves[j];
      const SCIP_Real sd = adjedgecosts[j];
      const SCIP_Real specialDist = EQ(sd, FARAWAY) ? -1.0 : sd;

      assert(!extreduce_nodeIsInStackTop(graph, extdata, leaf));
      assert(extdata->tree_deg[leaf] == 1);
      assert(leaf != topleaf);
      assert(adjids[j] == leaf);
      assert(EQ(specialDist, extreduce_extGetSD(scip, graph, topleaf, leaf, extdata)));

      /* any chance for rule-out? */
      if( !hasSiblings )
      {
         if( bottleneckIsDominated(graph, -1, topleaf, leaf, specialDist, extdata) )
         {
            SCIPdebugMessage("---bottleneck rule-out component (standard test)---\n");
            *leafRuledOut = TRUE;
            break;
         }
      }
   }

   /* clean up */
   if( !hasSiblings )
   {
      bottleneckUnmarkRootPath(graph, topleaf, extdata);

      if( isPc )
         pcSdToNodeUnmark(graph, extdata);
   }
}


/** Gets SDs from leaf (head of 'edge2leaf') to all other leaves of the tree.
 *  Returns early (with leafRuledOut == TRUE) if extension via 'edge2leaf' can be ruled out already.
 *  NOTE: Only restricted bottleneck tests are performed! */
static inline
void mstCompLeafGetSDs(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2leaf,          /**< edge to the top component leaf */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Real             sds[],              /**< array to store the SDs */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
   )
{
   const int nleaves_ancestors = extGetNancestorLeaves(extdata);
#ifndef NDEBUG
   const int compleaf = graph->head[edge2leaf];
#endif

   assert(leafRuledOut && !(*leafRuledOut));

   /* fill in the first part of the sds array */
   mstCompLeafGetSDsToAncestors(scip, graph, edge2leaf, nleaves_ancestors, extdata, sds, leafRuledOut);

   if( *leafRuledOut )
   {
      assert(dbgBottleneckFromLeafIsDominated(scip, graph, compleaf, extdata));
      return;
   }

   /* fill in the second part of the sds array */
   mstCompLeafGetSDsToSiblings(scip, graph, edge2leaf, extdata, &(sds[nleaves_ancestors]), leafRuledOut);

   if( *leafRuledOut )
   {
      assert(dbgBottleneckFromLeafIsDominated(scip, graph, compleaf, extdata));
      return;
   }

   //    if( !graph_pc_isPc(graph) ) //  remove once all is udated from storages
   assert(!dbgBottleneckFromLeafIsDominated(scip, graph, compleaf, extdata));

   //assert(extreduce_sdsTopInSync(scip, graph, sds, compleaf, extdata)); // todo wont work right now

}


/** computes SDs from head of extension edge to all leaves of the tree */
static inline
void mstLevelLeafSetVerticalSDs(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2neighbor,      /**< the edge from the tree to the neighbor */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            ruledOut            /**< early rule out? */
)
{
   MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   SCIP_Real* const adjedgecosts = extreduce_mldistsEmptySlotTargetDists(sds_vertical);
   const int* const leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const int neighbor = graph->head[edge2neighbor];
   const int neighbor_base = graph->tail[edge2neighbor];
   const int neighbor_base_proper = (neighbor_base == extdata->tree_root)? -1 : neighbor_base;

#ifndef NDEBUG
   int* const adjids = extreduce_mldistsEmptySlotTargetIds(sds_vertical);
   SCIP_Bool basehit = FALSE;
#endif

   assert(adjedgecosts && leaves && ruledOut);
   assert(*ruledOut == FALSE);

   for( int j = 0, k = 0; j < nleaves; j++ )
   {
      SCIP_Real specialDist;
      const int leaf = leaves[j];

      assert(extdata->tree_deg[leaf] == 1 && leaf != neighbor);

      specialDist = extGetSD(scip, graph, neighbor, leaf, extdata);

      /* save the SD? */
      if( leaf != neighbor_base_proper )
      {
         adjedgecosts[k] = (specialDist >= -0.5) ? specialDist : FARAWAY;
#ifndef NDEBUG
         adjids[k] = leaf;
#endif
         k++;
      }
      else
      {
#ifndef NDEBUG
         assert(!basehit);
         basehit = TRUE;
#endif
      }

      if( bottleneckIsDominated(graph, edge2neighbor, neighbor_base, leaf, specialDist, extdata) )
      {
         SCIPdebugMessage("---bottleneck rule-out---\n");
         assert(*ruledOut == FALSE);

         *ruledOut = TRUE;

         break;
      }
   }

#ifndef NDEBUG
   if( !(*ruledOut) )
   {
      assert(basehit || neighbor_base_proper != neighbor_base);
   }
#endif

}


/** initialization for adding a leaf to a level */
static inline
void mstLevelLeafInit(
   const GRAPH*          graph,              /**< graph data structure */
   int                   neighbor_base,      /**< neighbor base */
   int                   neighbor,           /**< neighbor */
   EXTDATA*              extdata             /**< extension data */
)
{
   MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   extreduce_mldistsEmtpySlotSetBase(neighbor, sds_vertical);

   /* Initialization for bottleneck. We start from the base of the neighbor! */
   bottleneckMarkRootPath(graph, neighbor_base, extdata);

   if( isPc )
   {
      pcSdToNodeMark(graph, neighbor, extdata);
   }
}


/** finalization for adding a leaf to a level */
static inline
void mstLevelLeafExit(
   const GRAPH*          graph,              /**< graph data structure */
   int                   neighbor_base,      /**< neighbor base */
   SCIP_Bool             ruledOut,           /**< early rule out? */
   EXTDATA*              extdata             /**< extension data */
)
{
   MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   if( ruledOut )
      extreduce_mldistsEmtpySlotReset(sds_vertical);
   else
      extreduce_mldistsEmptySlotSetFilled(sds_vertical);

   bottleneckUnmarkRootPath(graph, neighbor_base, extdata);

   if( isPc )
      pcSdToNodeUnmark(graph, extdata);
}


/** adds the initial vertex */
void extreduce_mstAddRoot(
   SCIP*                 scip,               /**< SCIP */
   int                   root,               /**< the root */
   REDDATA*              reddata             /**< reduction data */
)
{
   int todo; // move mst and bottleneck computation to extra file
   // maybe extreduce_extmst.c . Maybe need to add interface method for some static methods here
   // also these weird pcsd methods could be put there
   CSR* mst1;
   CSRDEPO* const mstdepo = reddata->msts;

   assert(root >= 0);
   assert(graph_csrdepo_isEmpty(mstdepo));

   // todo
#if 0
   graph_csrdepo_addEmptyTop(mstdepo, 1, 0);
   graph_csrdepo_getEmptyTop(mstdepo, mst1);

   assert(!graph_csrdepo_hasEmptyTop(mstdepo));

   reduce_dcmstGet1NodeMst(scip, mst1);
#endif
}




/** Adds leaf from top component of current tree to MST. I.e., adds SD adjacency costs updates MST.
 * 'edge2leaf' must be in top component of the stack.
 *  Returns early (with leafRuledOut == TRUE) if extension via 'edge2leaf' can be ruled out already.
 *  NOTE: SDs are not computed but taken from storage! */
void extreduce_mstCompAddLeaf(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2leaf,          /**< edge to the top component leaf */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
)
{
   DCMST* const dcmst = extdata->reddata->dcmst;
   SCIP_Real* const adjcosts = reduce_dcmstGetAdjcostBuffer(dcmst);

   assert(leafRuledOut && !(*leafRuledOut));

   mstCompLeafGetSDs(scip, graph, edge2leaf, extdata, adjcosts, leafRuledOut);

   // todo update mst from mst (in-place)
   // reduce_dcmstAddNodeInplace


}


/** adds current component (subset of the top level) */
void extreduce_mstCompInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2leaf,          /**< edge to the top component leaf */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
)
{
   DCMST* const dcmst = extdata->reddata->dcmst;
   SCIP_Real* const adjcosts = reduce_dcmstGetAdjcostBuffer(dcmst);

   assert(leafRuledOut && !(*leafRuledOut));
   assert(reduce_dcmstGetMaxnnodes(dcmst) >= extdata->tree_nleaves);

   mstCompLeafGetSDs(scip, graph, edge2leaf, extdata, adjcosts, leafRuledOut);

   // todo update mst from mst_reduce of previous!:
   // reduce_dcmstAddNode


}


/** Removes current component (subset of the top level) from MST storages */
void extreduce_mstCompRemove(
   const GRAPH*          graph,             /**< graph data structure */
   EXTDATA*              extdata            /**< extension data */
   )
{
   REDDATA* const reddata = extdata->reddata;
   CSRDEPO* const msts = reddata->msts;

//   assert(extreduce_stackTopMstDepoInSync(graph, extdata));
//   assert(graph_csrdepo_getNcsrs(msts) == extdata->tree_depth);

//   graph_csrdepo_removeTop(msts); todo
}


/** adds a full new level at the top */
void extreduce_mstLevelInit(
   REDDATA*              reddata,            /**< reduction data */
   EXTDATA*              extdata             /**< extension data */
)
{
   MLDISTS* const sds_vertical = reddata->sds_vertical;

   /* Reserve space for the SDs from each potential vertex of the new level to all leaves
    * of the tree except for the extending vertex.
    * But for the initial component we need to keep the root! */
   if( extStackGetPosition(extdata) == 0 )
      extreduce_mldistsLevelAddTop(STP_EXT_MAXGRAD, extdata->tree_nleaves, sds_vertical);
   else
      extreduce_mldistsLevelAddTop(STP_EXT_MAXGRAD, extdata->tree_nleaves - 1, sds_vertical);

   SCIPdebugMessage("init MST level %d \n", extreduce_mldistsNlevels(sds_vertical) - 1);

   assert(extdata->tree_depth == extreduce_mldistsNlevels(sds_vertical) - 1);
}


/** Adds neighbor of tree for MST calculation.
 *  Basically, SDs to all leafs are computed and stored in 'reddata->sds_vertical'.
 *  Neighbor is given by head of edge 'edge2neighbor'.
 *  Returns early (with leafRuledOut == TRUE, and without adding the neighbor)
 *  if extension via this edge can be ruled out already by using a bottleneck argument or MST. */
void extreduce_mstLevelAddLeaf(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2neighbor,      /**< the edge from the tree to the neighbor */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
)
{
   const int neighbor = graph->head[edge2neighbor];
   const int neighbor_base = graph->tail[edge2neighbor];
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   assert(leafRuledOut);
   assert(extdata->tree_deg[neighbor_base] == 1);
   assert(extdata->tree_deg[neighbor] == 0);
   assert(*leafRuledOut == FALSE);

   mstLevelLeafInit(graph, neighbor_base, neighbor, extdata);

   /* compute and store SDs to all leaves */
   mstLevelLeafSetVerticalSDs(scip, graph, edge2neighbor, extdata, leafRuledOut);

   /* if not yet ruled out, build the MST */
   if( !(*leafRuledOut) )
   {
       // todo compute MST (call add node method)

       // check for rule out

#ifdef STP_DEBUG_EXT
      {
         // todo compute weight!
         const SCIP_Real mstweight = extreduce_treeGetSdMstExtWeight(scip, graph, extvert, extdata);

         assert(GE(mstweight, 0.0));

        // printf("ext. mstobj=%f \n", mstweight);
      }
#endif
   }

   /* if not yet ruled out, try bottleneck distances to non-leaves of the tree */
   // todo do this also for STP!
   if( isPc && !(*leafRuledOut) )
   {
      bottleneckCheckNonLeaves(scip, graph, edge2neighbor, extdata, leafRuledOut);
   }

   mstLevelLeafExit(graph, neighbor_base, *leafRuledOut, extdata);
}

/** closes top MST level for further additions */
void extreduce_mstLevelClose(
   REDDATA*              reddata             /**< reduction data */
)
{
   MLDISTS* const sds_vertical = reddata->sds_vertical;
   int nslots;

   extreduce_mldistsLevelCloseTop(sds_vertical);

   nslots = extreduce_mldistsLevelNSlots(sds_vertical, extreduce_mldistsNlevels(sds_vertical) - 1);

   SCIPdebugMessage("close MST level %d, nslots=%d\n", extreduce_mldistsNlevels(sds_vertical) - 1,
         nslots);

   if( nslots > 0 )
   {
      // todo initialize mst_reduced for this level! If there is a positive number of slots!

   }
   else
   {
      // todo just add dummy entry?
   }

}


/** Removes top MST level.
 *  NOTE: SDs from level vertices to all leafs will be discarded! */
void extreduce_mstLevelRemove(
   REDDATA*              reddata             /**< reduction data */
)
{
   MLDISTS* const sds_vertical = reddata->sds_vertical;

   SCIPdebugMessage("remove MST level %d \n", extreduce_mldistsNlevels(sds_vertical) - 1);

   extreduce_mldistsLevelRemoveTop(sds_vertical);

   // todo also remove msts_reduced here?
}


/** Returns special distance.
 *  NOTE: might lead different result if 'vertex1' and 'vertex2' are swapped.  */
SCIP_Real extreduce_extGetSD(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   return extGetSD(scip, g, vertex1, vertex2, extdata);
}
