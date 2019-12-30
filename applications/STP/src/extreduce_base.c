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

/**@file   reduce_ext2.c
 * @brief  extended reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements extended reduction techniques for several Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"

#define EXT_STATE_NONE     0
#define EXT_STATE_EXPANDED 1
#define EXT_STATE_MARKED   2
#define EXT_REDCOST_NRECOMP 10
#define EXT_SDMAXVISITS 10


#ifndef NDEBUG
/** is current tree flawed? */
static
SCIP_Bool extTreeIsFlawed(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   int* edgecount;
   int* degreecount;
   const int* const tree_edges = extdata->tree_edges;
   const int* const tree_deg = extdata->tree_deg;
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const int tree_nedges = extdata->tree_nedges;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   int leavescount;
   const SCIP_Bool isPc = graph_pc_isPcMw(graph);

   SCIP_Bool flawed = FALSE;

   assert(nleaves >= 1);

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &edgecount, nedges) );
   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &degreecount, nnodes) );

   for( int i = 0; i < tree_nedges; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];
      const int tail = graph->tail[e];

      assert(e >= 0 && e < nedges);

      if( edgecount[e] > 0 || edgecount[flipedge(e)] > 0  )
      {
         printf("tree_nedges %d \n", tree_nedges);
         printf("FLAW: double edge \n");
         graph_edge_printInfo(graph, e);
         flawed = TRUE;
      }

      if( tree_deg[tail] <= 0 )
      {
         printf("FLAW: non-positive degree for %d (%d) \n", tail, tree_deg[tail]);
         flawed = TRUE;
      }

      if( tree_deg[head] <= 0 )
      {
         printf("FLAW: non-positive degree for %d (%d) \n", head, tree_deg[head]);
         flawed = TRUE;
      }

      degreecount[tail]++;
      degreecount[head]++;

      edgecount[e]++;
   }

   leavescount = 1; /* for tail of initial edge */

   /* degree check */
   for( int i = 0; i < tree_nedges && !flawed; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];

      if( degreecount[head] == 1 )
         leavescount++;

      if( degreecount[head] != tree_deg[head] )
      {
         printf("FLAW: wrong degree  \n");
         flawed = TRUE;
      }
   }

   /* leaves check */
   if( !flawed && leavescount != nleaves )
   {
      printf("FLAW wrong leaves count %d != %d \n", leavescount, nleaves);
      flawed = TRUE;
   }

   for( int i = 0; i < nleaves && !flawed; i++ )
   {
      const int leaf = tree_leaves[i];
      if( degreecount[leaf] != 1 )
      {
         printf("FLAW wrong leaf %d degree %d != %d \n", leaf, degreecount[leaf], 1);
         flawed = TRUE;
      }
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( isPc && extdata->pcSdToNode[i] >= -0.5 )
      {
         printf("FLAW wrong pcSdToNode entry[%d]=%f \n", i, extdata->pcSdToNode[i]);
         flawed = TRUE;
      }
      if( extdata->tree_bottleneckDistNode[i] >= -0.5 )
      {
         printf("FLAW wrong tree_bottleneckDistNode entry[%d]=%f \n", i, extdata->tree_bottleneckDistNode[i]);
         flawed = TRUE;
      }
   }

   /* clean-up */
   for( int i = 0; i < tree_nedges; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];
      const int tail = graph->tail[e];

      edgecount[e] = 0;
      degreecount[tail] = 0;
      degreecount[head] = 0;
   }

   SCIPfreeCleanBufferArray(scip, &degreecount);
   SCIPfreeCleanBufferArray(scip, &edgecount);

   return flawed;
}
#endif

/** insertion sort; note: could be speed-up by use of sentinel value at position 0 */
static inline
void sortDescendingIntRealReal(
   int*                  keyArr,             /**< key array of size 'nentries' */
   SCIP_Real*            dataArr1,           /**< array of size 'nentries' */
   SCIP_Real*            dataArr2,           /**< array of size 'nentries' */
   int                   nentries            /**< number of entries */
)
{
   assert(keyArr && dataArr1 && dataArr2);
   assert(nentries >= 1);

   for( int i = 1; i < nentries; i++ )
   {
      int j;
      const int currKey = keyArr[i];
      const SCIP_Real currData1 = dataArr1[i];
      const SCIP_Real currData2 = dataArr2[i];

      for( j = i - 1; j >= 0 && currKey > keyArr[j]; j-- )
      {
         keyArr[j + 1] = keyArr[j];
         dataArr1[j + 1] = dataArr1[j];
         dataArr2[j + 1] = dataArr2[j];
      }

      keyArr[j + 1] = currKey;
      dataArr1[j + 1] = currData1;
      dataArr2[j + 1] = currData2;
   }

#ifndef NDEBUG
   for( int i = 1; i < nentries; i++ )
      assert(keyArr[i - 1] >= keyArr[i] );
#endif
}


/** helper for rooted tree reduced cost computation */
static inline
SCIP_Real getMinDistCombination(
   const SCIP_Real*      firstTermDist,      /**< array of size 'nentries' */
   const SCIP_Real*      secondTermDist,     /**< array of size 'nentries' */
   int                   nentries            /**< number of entries to check */
)
{
   SCIP_Real min;

   assert(firstTermDist && secondTermDist);
   assert(nentries >= 1);

   if( nentries == 1 )
   {
      min = firstTermDist[0];
   }
   else
   {
      int i;
      SCIP_Real secondSum = 0.0;
      min = FARAWAY;

      for( i = 0; i < nentries; i++ )
      {
         assert(LE(firstTermDist[i], secondTermDist[i]));
         assert(LE(secondTermDist[i], FARAWAY));

         if( EQ(secondTermDist[i], FARAWAY) )
            break;

         secondSum += secondTermDist[i];
      }

      assert(LT(secondSum, FARAWAY));

      /* is there an index i with secondTermDist[i] == FARAWAY? */
      if( i < nentries )
      {
         assert(EQ(secondTermDist[i], FARAWAY));

         min = firstTermDist[i];
         for( int j = 0; j < nentries; j++ )
         {
            if( j == i )
               continue;

            min += secondTermDist[j];
         }
      }
      else
      {
         for( i = 0; i < nentries; i++ )
         {
            const SCIP_Real distCombination = secondSum + firstTermDist[i] - secondTermDist[i];

            assert(LT(secondTermDist[i], FARAWAY));

            if( distCombination < min )
               min = distCombination;
         }

         assert(LE(min, secondSum));
      }
   }

   return min;
}


/** is the edge valid? */
static
SCIP_Bool edgeIsValid(
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


/** deletes an edge and makes corresponding adaptations */
static
void removeEdge(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   DISTDATA*             distdata            /**< distance data (in/out) */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("remove edge ");
   graph_edge_printInfo(graph, edge);
#endif

   graph_edge_delFull(scip, graph, edge, TRUE);
   extreduce_distDataDeleteEdge(scip, graph, edge, distdata);

   if( graph->grad[tail] == 0 )
   {
      if( Is_term(graph->term[tail])  )
      {
         assert(graph_pc_isPcMw(graph) || tail == graph->source);
      }
      else
      {
         graph->mark[tail] = FALSE;
      }
   }

   if( graph->grad[head] == 0 )
   {
      if( Is_term(graph->term[head]) || head == graph->source )
      {
         assert(graph_pc_isPcMw(graph));
      }
      else
      {
         graph->mark[head] = FALSE;
      }
   }
}


/** get maximum allow depth for extended tree in given graph */
static
int getMaxTreeDepth(
   const GRAPH*          graph               /**< graph data structure */
)
{
   const int maxdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;

   assert(maxdepth > 0);

   return maxdepth;
}


/** returns current position in the stack */
static inline
int extStackGetPosition(
   const EXTDATA*        extdata             /**< extension data */
)
{
   assert(extdata->extstack_ncomponents > 0);
   return (extdata->extstack_ncomponents - 1);
}


/** prints the current stack */
static
void extStackPrint(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
#ifdef SCIP_DEBUG
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extdata->extstack_ncomponents - 1;

   for( int j = 0; j <= stackpos; j++ )
   {
      if( extdata->extstack_state[j] == EXT_STATE_NONE )
         printf("pos=%d state=NONE \n", j);
      else if( extdata->extstack_state[j] == EXT_STATE_EXPANDED )
         printf("pos=%d state=EXPANDED \n", j);
      else
      {
         assert(extdata->extstack_state[j] == EXT_STATE_MARKED);

         printf("pos=%d state=MARKED \n", j);
      }

      /* check all leaves of current component */
      for( int i = extstack_start[j]; i < extstack_start[j + 1]; i++ )
      {
         const int edge = extstack_data[i];
         assert(edge >= 0 && edge < graph->edges);

         printf("  ");
         graph_edge_printInfo(graph, edge);
      }
   }
#endif
}


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


/** can we extend the tree from given leaf? */
inline static
SCIP_Bool extLeafIsExtendable(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   int                   leaf                /**< the leaf */
)
{
   assert(graph && isterm);
   assert(leaf >= 0 && leaf < graph->knots);

   // todo if not a terminal, check whether number of neigbhors not contained in current tree < STP_EXT_MAXGRAD

   return (!isterm[leaf] && graph->grad[leaf] <= STP_EXT_MAXGRAD);
}


/** Finds position of given leaf in leaves data.
 *  Returns -1 if leaf could not be found. */
static inline
int extLeafFindPos(
   const EXTDATA*        extdata,            /**< extension data */
   int                   leaf,               /**< leaf to find */
   int                   startpos            /**< position to start from (going backwards) */
)
{
   int i;
   const int* const tree_leaves = extdata->tree_leaves;

   assert(extdata && tree_leaves);
   assert(startpos > 0);
   assert(leaf >= 0 && extdata->tree_deg[leaf] >= 1);

   for( i = startpos; i >= 0; i-- )
   {
      const int currleaf = tree_leaves[i];

      if( currleaf == leaf )
         break;
   }

   return i;
}


/** adds a new leaf */
static inline
void extLeafAdd(
   EXTDATA*              extdata,            /**< extension data */
   CGRAPH*               cgraph,             /**< complete graph data */
   int                   leaf                /**< leaf to add */
)
{
   assert(extdata && cgraph && extdata->tree_leaves);
   assert(leaf >= 0 && extdata->tree_deg[leaf] == 0);

   extdata->tree_leaves[(extdata->tree_nleaves)++] = leaf;

   cgraph_node_append(cgraph, leaf);
}


/** restores a previous leaf */
static inline
void extLeafRestore(
   EXTDATA*              extdata,            /**< extension data */
   CGRAPH*               cgraph,             /**< complete graph data */
   int                   leaf                /**< leaf to add */
)
{
    // maybe copy saved edge costs here if possible, check whether it is valid and if so, reinsert.
     // somehow need to remember whether leaf has been reinstated or not!


   assert(extdata && cgraph && extdata->tree_leaves);
   assert(leaf >= 0 && extdata->tree_deg[leaf] == 1);

   extdata->tree_leaves[(extdata->tree_nleaves)++] = leaf;

   cgraph_node_append(cgraph, leaf);
}


/** remove entry from leaves list */
static inline
void extLeafRemove(
   int                   leaf,              /**< leaf to remove */
   CGRAPH*               cgraph,            /**< complete graph data structure */
   EXTDATA*              extdata            /**< extension data */
)
{
   int* const tree_leaves = extdata->tree_leaves;
   int position;

   assert(cgraph);
   assert(extdata->tree_deg[leaf] == 1);

   /* switch last leaf and leaf to be removed */
   extdata->tree_nleaves--;
   assert(extdata->tree_nleaves > 0);

   position = extLeafFindPos(extdata, leaf, extdata->tree_nleaves);
   assert(position > 0);

   assert(cgraph->nodeids[position] == leaf);

   /* is leaf the last entry? */
   if( position == extdata->tree_nleaves)
      cgraph_node_deleteTop(cgraph);
   else
      cgraph_node_repositionTop(cgraph, position);

   tree_leaves[position] = tree_leaves[extdata->tree_nleaves];

   assert(cgraph_idsInSync(cgraph, tree_leaves, extdata->tree_nleaves));
}

/** remove top entry from leaves list */
static inline
void extLeafRemoveTop(
   int                   topsize,           /**< size of top to remove */
   CGRAPH*               cgraph,            /**< complete graph data structure */
   EXTDATA*              extdata            /**< extension data */
)
{
   assert(cgraph && extdata);
   assert(topsize >= 0 && topsize < extdata->tree_nleaves);

   for( int i = 0; i < topsize; i++ )
      cgraph_node_deleteTop(cgraph);

   extdata->tree_nleaves -= topsize;

   assert(cgraph_idsInSync(cgraph, extdata->tree_leaves, extdata->tree_nleaves));
}


/** marks single PcSd array entry */
static inline
void extPcSdMarkSingle(
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
      assert(*nPcSdCands < graph->knots);
      pcSdCands[(*nPcSdCands)++] = entry;
      pcSdToNode[entry] = value;
   }
   else if( value < pcSdToNode[entry] )
   {
      pcSdToNode[entry] = value;
   }
}


/** marks PcSd array */
static
void extPcSdToNodeMark(
   const GRAPH*          graph,              /**< graph data structure */
   int                   startvertex,        /**< vertex to start from */
   EXTDATA*              extdata,            /**< extension data */
   int*                  pcSdCands,          /**< marked candidates list */
   int*                  nPcSdCands          /**< pointer to store number of candidates */
   )
{
   SCIP_Real* const pcSdToNode = extdata->pcSdToNode;
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
   assert(pcSdToNode && prize && nPcSdCands);
   assert(*nPcSdCands == 0);

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
            if( tree_deg[head2] > 0 )
            {
               const SCIP_Real edgecost2 = cost_csr[i2];
               const SCIP_Real maxedgecost = MAX(edgecost, edgecost2);
               SCIP_Real dist2 = MAX(maxedgecost, edgecost + edgecost2 - prize[head]);

               assert(0.0 == prize[head] || Is_term(graph->term[head]));

               extPcSdMarkSingle(graph, head2, dist2, pcSdToNode, pcSdCands, nPcSdCands);
            }

            if( count2++ > EXT_SDMAXVISITS )
               break;
         }
      }
      else
      {
         extPcSdMarkSingle(graph, head, edgecost, pcSdToNode, pcSdCands, nPcSdCands);
      }

      if( count1++ > EXT_SDMAXVISITS )
         break;
   }
}


/** unmarks PcSd array */
static
void extPcSdToNodeUnmark(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            pcSdCands,          /**< marked candidates list */
   int                   nPcSdCands,         /**< number of candidates */
   EXTDATA*              extdata             /**< extension data */

   )
{
   SCIP_Real* const pcSdToNode = extdata->pcSdToNode;

   assert(graph_pc_isPcMw(graph));
   assert(pcSdToNode);

   for( int i = 0; i < nPcSdCands; i++ )
   {
      const int cand = pcSdCands[i];
      assert(pcSdToNode[cand] >= 0.0);
      pcSdToNode[cand] = -1.0;
   }
}


/** marks bottleneck array on path to tree root */
static
void extTreeBottleneckMarkRootPath(
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
void extTreeBottleneckUnmarkRootPath(
   int                   vertex,             /**< vertex to start from */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const int* const parentNode = extdata->tree_parentNode;
   const int tree_root = extdata->tree_root;

   assert(extdata && bottleneckDist_node && parentNode);
   assert(bottleneckDist_node[vertex] == -1.0);
   assert(bottleneckDist_node[tree_root] >= 0.0);

   if( vertex == tree_root )
   {
      bottleneckDist_node[vertex] = 0.0;
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
SCIP_Real extTreeBottleneckGetDist(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   vertex_pathmarked,  /**< vertex with marked rootpath */
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
static
SCIP_Bool extTreeBottleneckIsDominated(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex_pathmarked,  /**< vertex for which bottleneck path to root has been marked */
   int                   vertex_unmarked,    /**< second vertex */
   SCIP_Real             specialDist,        /**< best computed special distance approximation (-1.0 if unknown) */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Bool success = FALSE;

   assert(vertex_pathmarked >= 0 && vertex_pathmarked < graph->knots && vertex_unmarked >= 0 && vertex_unmarked < graph->knots);
   assert(SCIPisEQ(scip, specialDist, -1.0) || SCIPisGE(scip, specialDist, 0.0));

   if( specialDist >= -0.5 )
   {
      const SCIP_Real bottleneckDist = extTreeBottleneckGetDist(graph, extdata, vertex_pathmarked, vertex_unmarked);

      assert(SCIPisGE(scip, specialDist, 0.0));

      SCIPdebugMessage("%d->%d: sd=%f bottleneck=%f \n", vertex_pathmarked, vertex_unmarked, specialDist, bottleneckDist);

      if( SCIPisLT(scip, specialDist, bottleneckDist) )
         success = TRUE;
      else if( SCIPisLE(scip, specialDist, bottleneckDist) && 0 ) /* todo cover equality */
         success = TRUE;
   }

   return success;
}


/** gets reduced cost of current tree rooted at leave 'root', called direct if tree cannot */
static
SCIP_Real extTreeGetDirectedRedcostProper(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   root                /**< the root for the orientation */
)
{
   int nearestTerms[STP_EXTTREE_MAXNLEAVES_GUARD];
   SCIP_Real firstTermDist[STP_EXTTREE_MAXNLEAVES_GUARD];
   SCIP_Real secondTermDist[STP_EXTTREE_MAXNLEAVES_GUARD];
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const SCIP_Real swapcost = extdata->tree_redcostSwap[root];
   const REDDATA* const reddata = extdata->reddata;
   const PATH* const nodeTo3TermsPaths = reddata->nodeTo3TermsPaths;
   const int* const next3Terms = reddata->nodeTo3TermsBases;
   const SCIP_Bool* const isterm = extdata->node_isterm;
   const int* tree_deg = extdata->tree_deg;

   SCIP_Real redcost_directed = extdata->tree_redcost + reddata->rootToNodeDist[root] + swapcost;
   const int nnodes = graph->knots;
   int leavescount = 0;

#ifndef NDEBUG
   SCIP_Real redcost_debug = redcost_directed;
   for( int i = 0; i < STP_EXTTREE_MAXNLEAVES_GUARD; i++ )
   {
      nearestTerms[i] = -1;
      firstTermDist[i] = -1.0;
      secondTermDist[i] = -1.0;
   }

   assert(SCIPisLT(scip, redcost_directed, FARAWAY));
#endif

   for( int j = 0; j < nleaves; j++ )
   {
      int i;
      int term;
      const int leaf = tree_leaves[j];

      if( leaf == root || isterm[leaf] )
      {
         assert(leaf == root || EQ(0.0, nodeTo3TermsPaths[leaf].dist));
         continue;
      }

      /* find closest valid terminal for extension */
      for( i = 0; i < 3; i++ )
      {
         term = next3Terms[leaf + i * nnodes];

         if( term == UNKNOWN )
            break;

         assert(graph_pc_isPcMw(graph) || Is_term(graph->term[term]));
         assert(term >= 0 && term < graph->knots);
         assert(term != leaf);

         /* terminal not in current tree?*/
         if( tree_deg[term] == 0 )
            break;

         if( tree_deg[term] < 0 )
         {
            assert(graph_pc_isPcMw(graph) && tree_deg[term] == -1);
            assert(Is_pseudoTerm(graph->term[term]));
            break;
         }
      }

      /* no terminal reachable? */
      if( term == UNKNOWN )
      {
         assert(i < 3 && SCIPisGE(scip, nodeTo3TermsPaths[leaf + i * nnodes].dist, FARAWAY));
         return FARAWAY;
      }

      /* all terminals in current tree? */
      if( i == 3 )
         i = 2;

      nearestTerms[leavescount] = term;
      firstTermDist[leavescount] = nodeTo3TermsPaths[leaf + i * nnodes].dist;
      secondTermDist[leavescount] = (i < 2) ? nodeTo3TermsPaths[leaf + (i + 1) * nnodes].dist : firstTermDist[leavescount];

      assert(LE(nodeTo3TermsPaths[leaf].dist, firstTermDist[leavescount]));
      assert(LE(firstTermDist[leavescount], secondTermDist[leavescount]));

      //   printf("i %d \n", i);
 // printf("term=%d, first=%f second=%f def=%f \n", term, firstTermDist[leavescount], secondTermDist[leavescount], nodeTo3TermsPaths[leaf].dist);
      leavescount++;

#ifndef NDEBUG
      redcost_debug += nodeTo3TermsPaths[leaf].dist;
      assert(leavescount <= STP_EXTTREE_MAXNLEAVES_GUARD);
#endif
   }

   if( leavescount > 0 )
   {
      int first = 0;

      sortDescendingIntRealReal(nearestTerms, firstTermDist, secondTermDist, leavescount);

      for( int i = 1; i < leavescount; i++ )
      {
         assert(nearestTerms[i] >= 0 && firstTermDist[i] >= 0.0 && secondTermDist[i] >= 0.0);
         if( nearestTerms[i] != nearestTerms[i - 1] )
         {
            const int n = i - first;
            redcost_directed += getMinDistCombination(firstTermDist + first, secondTermDist + first, n);
            first = i;
         }
      }

      redcost_directed += getMinDistCombination(firstTermDist + first, secondTermDist + first, leavescount - first);
   }

 //  printf("redcost_directed=%f redcost_debug=%f  \n", redcost_directed, redcost_debug );

   assert(SCIPisGE(scip, redcost_directed, redcost_debug));

   return redcost_directed;
}


/** gets reduced cost of current tree rooted at leave 'root' */
static inline
SCIP_Real extTreeGetDirectedRedcost(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   root                /**< the root for the orientation */
)
{
   const SCIP_Real* const tree_redcostSwap = extdata->tree_redcostSwap;

   assert(extdata->tree_nleaves > 1 && extdata->tree_nleaves < STP_EXTTREE_MAXNLEAVES_GUARD);
   assert(extdata->tree_leaves[0] == extdata->tree_root);
   assert(root >= 0 && root < graph->knots);

   /* are there any deleted arcs in the directed tree? */
   if( root == extdata->tree_root && extdata->tree_nDelUpArcs > 0 )
   {
      return FARAWAY;
   }

   /* is the rooting possible? */
   if( SCIPisLT(scip, tree_redcostSwap[root], FARAWAY) )
   {
      return extTreeGetDirectedRedcostProper(scip, graph, extdata, root);
   }

   return FARAWAY;
}


/** gets reduced cost bound of current tree */
static
SCIP_Real extTreeGetRedcostBound(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const tree_leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Real tree_redcost;

   assert(graph && extdata);
   assert(nleaves > 1 && tree_leaves[0] == extdata->tree_root);

   tree_redcost = FARAWAY;

   // todo perhaps needs to be adapted for pseudo elimination tests...

   /* take each leaf as root of the tree */
   for( int i = 0; i < nleaves; i++ )
   {
      const int leaf = tree_leaves[i];
      const SCIP_Real tree_redcost_new = extTreeGetDirectedRedcost(scip, graph, extdata, leaf);

      tree_redcost = MIN(tree_redcost, tree_redcost_new);
   }

   return tree_redcost;
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
   int* const tree_edges = extdata->tree_edges;
   int* const tree_deg = extdata->tree_deg;
   int* const tree_parentNode = extdata->tree_parentNode;
   SCIP_Real* const tree_parentEdgeCost = extdata->tree_parentEdgeCost;
   SCIP_Real* const tree_redcostSwap = extdata->tree_redcostSwap;
   REDDATA* const reddata = extdata->reddata;
   CGRAPH* const cgraph = reddata->cgraph;
   const SCIP_Real* const redcost = reddata->redCosts;
   int* const pseudoancestor_mark = reddata->pseudoancestor_mark;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   const int stackpos = extStackGetPosition(extdata);
   const int comproot = graph->tail[extstack_data[extstack_start[stackpos]]];
   int conflictIteration = -1;
   const SCIP_Bool noReversedRedCostTree = (reddata->redCostRoot == extdata->tree_root || SCIPisGE(scip, tree_redcostSwap[comproot], FARAWAY));

   assert(!(*conflict));
   assert(stackpos >= 0);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);
   assert(comproot >= 0 && comproot < graph->knots);
   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   /* update tree leaves array todo might need to be changed for pseudo-elimination */
   if( comproot != extdata->tree_root )
      extLeafRemove(comproot, reddata->cgraph, extdata);
   else
   {
      assert(extdata->tree_nleaves == 1);
   }

   /* add top expanded component to tree data */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];
      const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);

      assert(extdata->tree_nedges < extdata->extstack_maxsize);
      assert(edge >= 0 && edge < graph->edges);
      assert(tree_deg[head] == 0);
      assert(tree_deg[comproot] > 0 || comproot == extdata->tree_root);
      assert(comproot == graph->tail[edge]);

      if( noReversedRedCostTree || (edgedeleted && edgedeleted[flipedge(edge)]) )
      {
         tree_redcostSwap[head] = FARAWAY;
      }
      else
      {
         tree_redcostSwap[head] = tree_redcostSwap[comproot] + redcost[flipedge(edge)];

         if( !edgeIsDeleted )
            tree_redcostSwap[head] -= redcost[edge];

         assert(LT(tree_redcostSwap[head], FARAWAY));
      }

      if( !edgeIsDeleted )
      {
         extdata->tree_redcost += redcost[edge];
         assert(LT(extdata->tree_redcost, FARAWAY));
      }
      else
      {
         extdata->tree_nDelUpArcs++;
      }

      extLeafAdd(extdata, cgraph, head);
      tree_deg[head] = 1;
      tree_edges[(extdata->tree_nedges)++] = edge;
      tree_parentNode[head] = comproot;
      tree_parentEdgeCost[head] = graph->cost[edge];
      tree_deg[comproot]++;

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

   assert(!extTreeIsFlawed(scip, graph, extdata));
}


/** removes top component of stack from tree */
static
void extTreeStackTopRemove(
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Bool             ancestor_conflict,  /**< with ancestor conflict? */
   EXTDATA*              extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real* const redcost = reddata->redCosts;
   int* const pseudoancestor_mark = reddata->pseudoancestor_mark;
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const tree_deg = extdata->tree_deg;
   const int stackpos = extStackGetPosition(extdata);
   const int stackstart = extstack_start[stackpos];
   const int comproot = graph->tail[extstack_data[extstack_start[stackpos]]];
   const int compsize = extstack_start[stackpos + 1] - extstack_start[stackpos];
   const STP_Bool* const edgedeleted = reddata->edgedeleted;

   assert(compsize > 0);
   assert(tree_deg[comproot] > 1);
   assert(extdata->extstack_state[stackpos] != EXT_STATE_NONE);

   /* remove top component */
   for( int i = stackstart; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];
      const int tail = graph->tail[edge];
      const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);

      assert(edge >= 0 && edge < graph->edges);
      assert(tree_deg[head] == 1 && tree_deg[tail] > 1);

      if( !edgeIsDeleted )
      {
         extdata->tree_redcost -= redcost[edge];
         assert(LT(extdata->tree_redcost, FARAWAY));
      }
      else
      {
         extdata->tree_nDelUpArcs--;
         assert(extdata->tree_nDelUpArcs >= 0);
      }

      tree_deg[head] = 0;
      tree_deg[tail]--;

      if( !ancestor_conflict ) /* in case of a conflict, edge is unhashed already */
         graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, edge, pseudoancestor_mark);
   }

   (extdata->tree_nedges) -= compsize;

  // for( int k = extdata->tree_nleaves - 1; k >= extdata->tree_nleaves - compsize; k-- ) printf("bt remove leaf %d \n", extdata->tree_leaves[k]);

   extLeafRemoveTop(compsize, reddata->cgraph, extdata);

   (extdata->tree_depth)--;

   /* restore component root as a leaf */
   assert(tree_deg[comproot] == 1);
   extLeafRestore(extdata, reddata->cgraph, comproot);

   assert(extdata->tree_nedges >= 0 && extdata->tree_depth >= 0);
}


/** recompute reduced costs */
static
void extTreeRecompRedCosts(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
#ifndef NDEBUG
   const int tree_nDelUpArcs = extdata->tree_nDelUpArcs;
#endif
   REDDATA* const reddata = extdata->reddata;
   const STP_Bool* const edgedeleted = reddata->edgedeleted;
   SCIP_Real treecost = 0.0;
   const SCIP_Real* const redcost = reddata->redCosts;
   const int* const tree_edges = extdata->tree_edges;
   const int tree_nedges = extdata->tree_nedges;

   extdata->tree_nDelUpArcs = 0;

   assert(!extTreeIsFlawed(scip, graph, extdata));

   for( int i = 0; i < tree_nedges; i++ )
   {
      const int edge = tree_edges[i];
      const SCIP_Bool edgeIsDeleted = (edgedeleted && edgedeleted[edge]);

      assert(edge >= 0 && edge < graph->edges);

      if( !edgeIsDeleted )
      {
         treecost += redcost[edge];
         assert(LT(treecost, FARAWAY));
      }
      else
      {
         extdata->tree_nDelUpArcs++;
      }
   }

   assert(SCIPisEQ(scip, treecost, extdata->tree_redcost));
   assert(tree_nDelUpArcs == extdata->tree_nDelUpArcs);

   extdata->tree_redcost = treecost;
}


/** some updates */
static
void extTreeSyncWithStack(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   int*                  nupdatestalls,      /**< update stalls counter */
   SCIP_Bool*            conflict            /**< conflict found? */
)
{
   const int stackposition = extStackGetPosition(extdata);

   assert(scip && graph && extdata && nupdatestalls && conflict);
   assert(!(*conflict));

   extStackPrint(graph, extdata);

   /* is current component expanded? */
   if( extdata->extstack_state[stackposition] == EXT_STATE_EXPANDED )
      extTreeStackTopAdd(scip, graph, extdata, conflict); /* add component to tree */

   /* recompute reduced costs? */
   if( ++(*nupdatestalls) > EXT_REDCOST_NRECOMP )
   {
      extTreeRecompRedCosts(scip, graph, extdata);
      *nupdatestalls = 0;
   }

   /* assert that the entire tree is hashed */
#ifndef NDEBUG
   if( *conflict == FALSE )
   {
      const REDDATA* const reddata = extdata->reddata;

      for( int i = 0; i < extdata->tree_nedges; i++ )
      {
         const int edge = extdata->tree_edges[i];
         assert(graph_edge_nPseudoAncestors(graph, edge) == 0 ||
            graph_pseudoAncestors_edgeIsHashed(graph->pseudoancestors, edge, reddata->pseudoancestor_mark));
      }
   }
#endif
}


#ifndef NDEBUG
/** gets reduced cost bound of current tree */
static
SCIP_Bool extMSTisInSync(
   const CGRAPH*         cgraph              /**< complete graph data */
   )
{
// assert that the costs of the tree all coincide with the actual SD etc distances!
// might be good to have this and reddata extdata stuff in extra method reduce_ext_util.c or just reduce_util.c
// need some flag (in cgraph?) to see whether a leaf in the cgraph does not actually have valid costs (or any)
// and should be recomputed!
   return TRUE;
}


#endif


/** gets reduced cost bound of current tree */
static
void extMSTaddLeaf(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   extleaf,            /**< the leaf */
   EXTDATA*              extdata,            /**< extension data */
   CGRAPH*               cgraph,             /**< complete graph data */
   int*                  pcSdCands,          /**< != NULL iff PC/RPC */
   SCIP_Bool*            leafRuledOut        /**< could the extension already by ruled out */
)
{
   SCIP_Real* const adjedgecosts = cgraph->adjedgecosts;
   const int* const leaves = extdata->tree_leaves;
   int leafpos = -1;
   int nPcSdCands = 0;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Bool ruledOut = FALSE;
   const SCIP_Bool isPc = (pcSdCands != NULL);

   assert(adjedgecosts && leaves);
   assert(isPc == graph_pc_isPcMw(graph));

   extTreeBottleneckMarkRootPath(graph, extleaf, extdata);

   /* for PC/RPC initialize pcSdToNode array */
   if( isPc )
      extPcSdToNodeMark(graph, extleaf, extdata, pcSdCands, &nPcSdCands);

   for( int j = 0; j < nleaves; j++ )
   {
      SCIP_Real specialDist;
      const int leaf = leaves[j];

      assert(extdata->tree_deg[leaf] == 1);

      if( leaf == extleaf )
      {
         leafpos = j;
         adjedgecosts[j] = FARAWAY;
         continue;
      }

      specialDist = extGetSD(scip, graph, extleaf, leaf, extdata);
      adjedgecosts[j] = specialDist >= -0.5 ? specialDist : FARAWAY;

      if( extTreeBottleneckIsDominated(scip, graph, extleaf, leaf, specialDist, extdata) )
      {
         SCIPdebugMessage("---bottleneck rule-out---\n");
         ruledOut = TRUE;
         break;
      }
   }


   if( isPc && !ruledOut )
   {
      const int* const tree_deg = extdata->tree_deg;
      // todo if not successful so far, perhaps try bottleneck distances for inner vertices of STP as well!

      /* also check non-leaves */
      for( int c = 0; c < nPcSdCands; c++ )
      {
         SCIP_Real specialDist;
         const int cand = pcSdCands[c];

         assert(cand >= 0 && cand < graph->knots);

         /* leaf or not contained? */
         if( tree_deg[cand] <= 1 )
            continue;

         specialDist = extGetSD(scip, graph, extleaf, cand, extdata);

         if( extTreeBottleneckIsDominated(scip, graph, extleaf, cand, specialDist, extdata) )
         {
            SCIPdebugMessage("---non-leaf bottleneck rule-out---\n");
            ruledOut = TRUE;
            break;
         }
      }
   }

   if( !ruledOut && 0 )
   {
      assert(leafpos >= 0);

      int todo;
      cgraph_node_applyMinAdjCosts(cgraph, leafpos, extleaf);
   }

   extTreeBottleneckUnmarkRootPath(extleaf, extdata);

   if( isPc )
      extPcSdToNodeUnmark(graph, pcSdCands, nPcSdCands, extdata);

   *leafRuledOut = ruledOut;
}


/** adds initial component to stack (needs to be star component rooted in root) */
static
void extAddInitialComponent(
   const GRAPH*          graph,              /**< graph data structure */
   const int*            compedges,          /**< component edges */
   int                   ncompedges,         /**< number of component edges */
   int                   root,               /**< root of the component */
   EXTDATA*              extdata             /**< extension data */
)
{
   CGRAPH* const cgraph = extdata->reddata->cgraph;

   assert(compedges);
   assert(ncompedges >= 1 && ncompedges < STP_EXT_MAXGRAD);
   assert(ncompedges < extdata->extstack_maxedges);
   assert(root >= 0 && root < graph->knots);

#ifdef SCIP_DEBUG
   printf("\n --- ADD initial component --- \n\n");
#endif

   // todo method needs to be adapted for pseudo-elimination!

   for( int i = 0; i < ncompedges; i++ )
   {
      const int e = compedges[i];
      const int tail = graph->tail[e];

      assert(e >= 0 && e < graph->edges);
      assert(tail == root);

      SCIPdebugMessage("edge %d: %d->%d \n", e, graph->tail[e], graph->head[e]);

      extdata->extstack_data[i] = e;
      extLeafAdd(extdata, cgraph, tail);
   }

   extdata->tree_root = root;
   extdata->extstack_ncomponents = 1;
   extdata->extstack_state[0] = EXT_STATE_EXPANDED;
   extdata->extstack_start[0] = 0;
   extdata->extstack_start[1] = ncompedges;
   extdata->tree_parentNode[root] = -1;
   extdata->tree_redcostSwap[root] = 0.0;
   extdata->tree_parentEdgeCost[root] = -1.0;
   assert(ncompedges > 1 || extdata->tree_leaves[0] == root);
   assert(extdata->tree_deg[root] == 0);
   assert(extdata->tree_nleaves == ncompedges);
}


/** can any extension via edge be ruled out? */
static
SCIP_Bool extRuleOutEdge(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   edge                /**< edge to be tested */
)
{
   REDDATA* const reddata = extdata->reddata;
   const int* const tree_deg = extdata->tree_deg;
   const int extvert = graph->head[edge];

   assert(scip && graph && reddata && extdata);
   assert(edge >= 0 && edge < graph->edges);

   if( tree_deg[extvert] != 0 )
   {
      return TRUE;
   }
   else
   {
      if( graph_pseudoAncestors_edgeIsHashed(graph->pseudoancestors, edge, reddata->pseudoancestor_mark) )
      {
         SCIPdebugMessage("ancestor simple rule-out  ");
         return TRUE;
      }
   }

   return FALSE;
}

#if 0
/** can any extension via edge except for only the edge itself be ruled out? */
static
SCIP_Bool extRuleOutEdgeCombinations(
   SCIP*                 scip,               /**< SCIP */
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


/** can the extension with only the edge itself be ruled out? */
static
SCIP_Bool extRuleOutEdgeSingle(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   edge                /**< edge to be tested */
)
{
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real tree_redcost = extTreeGetRedcostBound(scip, graph, extdata, edge);
   const SCIP_Real cutoff = reddata->cutoff;

   if( reddata->equality ? (SCIPisGE(scip, tree_redcost, cutoff)) : SCIPisGT(scip, tree_redcost, cutoff) )
   {
      SCIPdebugMessage("redcost simple rule-out (%f >= %f)  ", tree_redcost, cutoff);
      return TRUE;
   }

}
#endif


/** can subtree be peripherally ruled out? */
static
SCIP_Bool extRuleOutPeriph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real tree_redcost = extTreeGetRedcostBound(scip, graph, extdata);
   const SCIP_Real cutoff = reddata->cutoff;

   if( reddata->equality ? (SCIPisGE(scip, tree_redcost, cutoff)) : SCIPisGT(scip, tree_redcost, cutoff) )
   {
      SCIPdebugMessage("Rule-out periph (with red.cost=%f) \n", tree_redcost);
      return TRUE;
   }
   else
   {
      CGRAPH* const cgraph = reddata->cgraph;
      int* pcSdCands = NULL;
      const int* const extstack_data = extdata->extstack_data;
      const int* const extstack_start = extdata->extstack_start;
      const int stackpos = extStackGetPosition(extdata);
      const SCIP_Bool isPc = (STP_PCSPG == graph->stp_type || STP_RPCSPG == graph->stp_type);
      SCIP_Bool ruledOut = FALSE;

      assert(EXT_STATE_EXPANDED == extdata->extstack_state[stackpos]);
      assert(cgraph_idsInSync(cgraph, extdata->tree_leaves, extdata->tree_nleaves));

      if( isPc )
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &pcSdCands, graph->knots) );

      /* compute special distances for MST test (and compare with tree bottleneck distances for early rule-out) */
      for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
      {
         const int extleaf = graph->head[extstack_data[i]];

         extMSTaddLeaf(scip, graph, extleaf, extdata, cgraph, pcSdCands, &ruledOut);

         /* early rule out? */
         if( ruledOut )
            break;
      }

      SCIPfreeBufferArrayNull(scip, &pcSdCands);

      if( ruledOut )
         return TRUE;


      /* assert the that cgraph is clean now... */

      /* now compute the MST */


#ifndef NDEBUG
      for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
         assert( graph_edge_nPseudoAncestors(graph, extstack_data[i]) == 0 || graph_pseudoAncestors_edgeIsHashed(graph->pseudoancestors, extstack_data[i], reddata->pseudoancestor_mark));
#endif
   }

   return FALSE;
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

   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

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

   if( extstack_start[stackpos] >= extdata->extstack_maxedges )
   {
      SCIPdebugMessage("truncate (too many edges on stack) \n");
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

   assert(graph && extdata);
   assert(extdata->extstack_start[stackpos + 1] - extdata->extstack_start[stackpos] > 0);

   /* top component already expanded? */
   if( extstack_state[stackpos] != EXT_STATE_NONE )
   {
      extTreeStackTopRemove(graph, ancestor_conflict, extdata);
   }

   stackpos--;

   /* backtrack */
   if( success )
   {
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
      while( extstack_state[stackpos] == EXT_STATE_EXPANDED )
      {
         stackpos--;
         assert(stackpos >= 0);
      }

      SCIPdebugMessage("backtrack FAILURE \n");
      assert(extstack_state[stackpos] == EXT_STATE_NONE || extstack_state[stackpos] == EXT_STATE_MARKED);
   }

   extdata->extstack_ncomponents = stackpos + 1;

   assert(!extTreeIsFlawed(scip, graph, extdata));
}

/** expands top component of stack (backtracks if stack is full) */
static
void extStackExpand(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extStackGetPosition(extdata);
   int datasize = extstack_start[stackpos];
   const int setsize = extstack_start[stackpos + 1] - extstack_start[stackpos];
   const uint32_t powsize = (uint32_t) pow(2.0, setsize);

#ifndef NDEBUG
   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;
#endif

   assert(extdata && scip && graph && success);
   assert(setsize <= STP_EXT_MAXGRAD);
   assert(setsize > 0 && setsize <= 32);
   assert(stackpos >= 1);
   assert(extstack_state[stackpos] == EXT_STATE_NONE);

   /* stack too full? */
   if( (datasize + (int) powsize * (setsize + 1) / 2) > extdata->extstack_maxsize )
   {
      *success = FALSE;
      extBacktrack(scip, graph, *success, FALSE, extdata);

      return;
   }

   /* collect edges for new component and find conflicts */
   for( int i = extstack_start[stackpos], j = 0; i < extstack_start[stackpos + 1]; i++, j++ )
   {
      const int edge = extstack_data[i];
      assert(j < STP_EXT_MAXGRAD);
      assert(edge >= 0 && edge < graph->edges);
      assert(extdata->tree_deg[graph->head[edge]] == 0);

      // todo find excluding pairs, use ancestormark and bottleneck

      extedges[j] = edge;
   }

   /* compute and add components (overwrite previous, non-expanded component) */
   for( uint32_t counter = 1; counter < powsize; counter++ )
   {
      for( unsigned int j = 0; j < (unsigned) setsize; j++ )
      {
         /* Check if jth bit in counter is set */
         if( counter & (1 << j) )
         {
            assert(datasize < extdata->extstack_maxsize);
            assert(extedges[j] >= 0);

            extstack_data[datasize++] = extedges[j];
            SCIPdebugMessage(" head %d \n", graph->head[extedges[j]]);
         }
      }

      SCIPdebugMessage("... added \n");
      assert(stackpos < extdata->extstack_maxsize - 1);

      extstack_state[stackpos] = EXT_STATE_EXPANDED;
      extstack_start[++stackpos] = datasize;

      assert(extstack_start[stackpos] - extstack_start[stackpos - 1] > 0);
   }

   assert(stackpos >= extdata->extstack_ncomponents);

   extdata->extstack_ncomponents = stackpos;
}


/** extends top component of stack (or backtracks if stack is full) */
static
void extExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD * STP_EXT_MAXGRAD];
   int extedgesstart[STP_EXT_MAXGRAD + 1];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   const SCIP_Bool* const isterm = extdata->node_isterm;
   int stackpos = extStackGetPosition(extdata);
   int nfullextensions;
   int nsingleextensions;
   SCIP_Bool with_ruledout_leaf;

#ifndef NDEBUG
   assert(stackpos >= 0);
   assert(extstack_state[stackpos] == EXT_STATE_EXPANDED );
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] <= STP_EXT_MAXGRAD);

   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedgesstart[i] = -1;

   for( int i = 0; i < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;
#endif

   extstack_state[stackpos] = EXT_STATE_MARKED;

   nfullextensions = 0;
   nsingleextensions = 0;
   extedgesstart[0] = 0;

   with_ruledout_leaf = FALSE;

   /* loop over all leaves of extension */
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
         if( !extRuleOutEdge(scip, graph, extdata, e) )
         {
            int todo;
            // check whether we can kill extension with Combinations test

            // yes?? then check whether we can rule out single

            assert(nsingleextensions < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD);
            extedges[nsingleextensions++] = e;
            nleafextensions++;
         }
#ifdef SCIP_DEBUG
         else
         {
            printf("simple rule out: ");
            graph_edge_printInfo(graph, e);
         }
#endif
      }

      if( nleafextensions == 0 )
      {
         with_ruledout_leaf = TRUE;
         break;
      }

      extedgesstart[++nfullextensions] = nsingleextensions;
   }

   assert(nfullextensions <= STP_EXT_MAXGRAD);

   if( with_ruledout_leaf )
   {
      SCIPdebugMessage("ruled-out one leaf \n");
      *success = TRUE;
   }
   else if( nfullextensions == 0 )  /* found no valid extensions? */
   {
      assert(nsingleextensions == 0);
      SCIPdebugMessage("no valid extensions found \n");
      *success = FALSE;
   }
   else  /* found non-empty valid extensions */
   {
      int datasize = extstack_start[stackpos + 1];
      int extsize[STP_EXT_MAXGRAD];
      int extindex[STP_EXT_MAXGRAD];

      assert(nsingleextensions > 0);
      assert(extedgesstart[nfullextensions] - extedgesstart[0] > 0);

      /* stack too small? */
      if( datasize + (extedgesstart[nfullextensions] - extedgesstart[0]) > extdata->extstack_maxsize )
      {
         *success = FALSE;
         extBacktrack(scip, graph, *success, FALSE, extdata);

         return;
      }

      for( int i = 0; i < nfullextensions; i++ )
      {
         assert(extedgesstart[i + 1] >= 0);

         extsize[i] = extedgesstart[i + 1] - extedgesstart[i];
         extindex[i] = i;
         assert(extsize[i] >= 0);
      }

      SCIPsortDownIntInt(extsize, extindex, nfullextensions);

      /* put the non-empty extensions on the stack, with smallest last */
      for( int i = 0; i < nfullextensions; i++ )
      {
         const int index = extindex[i];

         if( extsize[i] == 0 )
         {
            assert(i > 0);
            assert(extedgesstart[index + 1] - extedgesstart[index] == 0);

            for( int j = i; j < nfullextensions; j++ )
               assert(extsize[j] == 0);

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

      *success = TRUE;

      /* try to expand last (smallest) component */
      extStackExpand(scip, graph, extdata, success);
   }
}


/** check (directed) arc (internal method) */
static
void extCheckArc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< directed edge to be checked */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            deletable           /**< is arc deletable? */
)
{
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   int* const tree_deg = extdata->tree_deg;
   int* const extstack_state = extdata->extstack_state;
   int nupdatestalls = 0;
   SCIP_Bool success = TRUE;
   SCIP_Bool conflict = FALSE;

   assert(extreduce_reddataIsClean(graph, extdata->reddata) && extreduce_extdataIsClean(graph, extdata));

   /* put 'edge' on the stack */
   extAddInitialComponent(graph, &edge, 1, tail, extdata);

   extTreeSyncWithStack(scip, graph, extdata, &nupdatestalls, &conflict);

   assert(!conflict);
   assert(extdata->tree_parentNode[head] == tail && extdata->tree_parentEdgeCost[head] == graph->cost[edge]);

   extExtend(scip, graph, extdata, &success);

   assert(extstack_state[0] == EXT_STATE_MARKED);
   assert(success || 1 == extdata->extstack_ncomponents);

   /* limited DFS backtracking; stops once back at 'edge' */
   while( extdata->extstack_ncomponents > 1 )
   {
      const int stackposition = extStackGetPosition(extdata);
      conflict = FALSE;

      extTreeSyncWithStack(scip, graph, extdata, &nupdatestalls, &conflict);

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

         extStackExpand(scip, graph, extdata, &success);
         continue;
      }

      assert(extstack_state[stackposition] == EXT_STATE_EXPANDED);

      if( conflict || extRuleOutPeriph(scip, graph, extdata) )
      {
         success = TRUE;
         extBacktrack(scip, graph, success, conflict, extdata);
         continue;
      }

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
   assert(tree_deg[head] == 1 && tree_deg[tail] == 1 && extdata->tree_nedges == 1);

   tree_deg[head] = 0;
   tree_deg[tail] = 0;

   graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, edge, extdata->reddata->pseudoancestor_mark);
   extreduce_extdataClean(extdata);
   extreduce_reddataClean(extdata->reddata);

   assert(extreduce_reddataIsClean(graph, extdata->reddata) && extreduce_extdataIsClean(graph, extdata));
}


/** check (directed) arc */
SCIP_RETCODE extreduce_checkArc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   int                   edge,               /**< edge to be checked */
   SCIP_Bool             equality,           /**< delete edge also in case of reduced cost equality? */
   DISTDATA*             distdata,           /**< data for distance computations */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            edgeIsDeletable     /**< is edge deletable? */
)
{
   const SCIP_Bool* isterm = extpermanent->isterm;
   const int root = redcostdata->redCostRoot;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real* rootdist = redcostdata->rootToNodeDist;
   const PATH* nodeToTermpaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real cutoff = redcostdata->cutoff;
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   const SCIP_Real edgebound = redcost[edge] + rootdist[tail] + nodeToTermpaths[head].dist;
   SCIP_Bool restoreAntiArcDeleted = FALSE;
   STP_Bool* const edgedeleted = extpermanent->edgedeleted;

   assert(scip && graph && redcost && rootdist && nodeToTermpaths && distdata);
   assert(edge >= 0 && edge < graph->edges);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph->mark[tail] && graph->mark[head]);
   assert(graph_isMarked(graph));
   assert(extreduce_extPermaIsClean(graph, extpermanent));

   /* trivial rule-out? */
   if( SCIPisGT(scip, edgebound, cutoff) || (equality && SCIPisEQ(scip, edgebound, cutoff)) || head == root )
   {
      *edgeIsDeletable = TRUE;
      return SCIP_OKAY;
   }

   if( edgedeleted && !edgedeleted[flipedge(edge)] )
   {
      edgedeleted[flipedge(edge)] = TRUE;
      restoreAntiArcDeleted = TRUE;
   }

   *edgeIsDeletable = FALSE;

   /* can we extend from 'edge'? */
   if( extLeafIsExtendable(graph, isterm, head) )
   {
      int* extstack_data;
      int* extstack_start;
      int* extstack_state;
      int* tree_edges;
      int* tree_leaves;
      int* tree_parentNode;
      SCIP_Real* tree_parentEdgeCost;
      SCIP_Real* tree_redcostSwap;
      int* pseudoancestor_mark;
      const int nnodes = graph->knots;
      const int maxstackedges = MIN(nnodes / 2, STP_EXT_MAXEDGES);

      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_data, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_state, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_edges, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_leaves, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentNode, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentEdgeCost, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_redcostSwap, nnodes) );

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &pseudoancestor_mark, nnodes) );

      {
         REDDATA reddata = { .cmst = extpermanent->cmst, .cgraph = extpermanent->cgraph,
            .cgraphEdgebuffer = extpermanent->cgraphEdgebuffer,
            .redCosts = redcost, .rootToNodeDist = rootdist, .nodeTo3TermsPaths = nodeToTermpaths,
            .nodeTo3TermsBases = redcostdata->nodeTo3TermsBases, .edgedeleted = edgedeleted,
            .pseudoancestor_mark = pseudoancestor_mark, .cutoff = cutoff, .equality = equality, .redCostRoot = root };
         EXTDATA extdata = { .extstack_data = extstack_data, .extstack_start = extstack_start,
            .extstack_state = extstack_state, .extstack_ncomponents = 0, .tree_leaves = tree_leaves,
            .tree_edges = tree_edges, .tree_deg = extpermanent->tree_deg, .tree_nleaves = 0,
            .tree_bottleneckDistNode = extpermanent->bottleneckDistNode, .tree_parentNode = tree_parentNode,
            .tree_parentEdgeCost = tree_parentEdgeCost, .tree_redcostSwap = tree_redcostSwap, .tree_redcost = 0.0,
            .tree_nDelUpArcs = 0, .tree_root = -1, .tree_nedges = 0, .tree_depth = 0, .extstack_maxsize = nnodes - 1,
            .pcSdToNode = extpermanent->pcSdToNode, .extstack_maxedges = maxstackedges,
            .tree_maxnleaves = STP_EXTTREE_MAXNLEAVES, .tree_maxdepth = getMaxTreeDepth(graph),
            .tree_maxnedges = STP_EXTTREE_MAXNEDGES, .node_isterm = isterm, .reddata = &reddata, .distdata = distdata };

         extCheckArc(scip, graph, edge, &extdata, edgeIsDeletable);
      }

      assert(extreduce_extPermaIsClean(graph, extpermanent));

      SCIPfreeCleanBufferArray(scip, &pseudoancestor_mark);
      SCIPfreeBufferArray(scip, &tree_redcostSwap);
      SCIPfreeBufferArray(scip, &tree_parentEdgeCost);
      SCIPfreeBufferArray(scip, &tree_parentNode);
      SCIPfreeBufferArray(scip, &tree_leaves);
      SCIPfreeBufferArray(scip, &tree_edges);
      SCIPfreeBufferArray(scip, &extstack_state);
      SCIPfreeBufferArray(scip, &extstack_start);
      SCIPfreeBufferArray(scip, &extstack_data);
   }

   if( restoreAntiArcDeleted )
   {
      assert(edgedeleted);
      edgedeleted[flipedge(edge)] = FALSE;
   }

   return SCIP_OKAY;
}


/** check edge */
SCIP_RETCODE extreduce_checkEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   int                   edge,               /**< edge to be checked */
   SCIP_Bool             equality,           /**< delete edge also in case of reduced cost equality? */
   DISTDATA*             distdata,           /**< data for distance computations */
   EXTPERMA*             extpermanent,       /**< extension data */
   SCIP_Bool*            edgeIsDeletable     /**< is edge deletable? */
)
{
   const SCIP_Bool* isterm = extpermanent->isterm;
   const int root = redcostdata->redCostRoot;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real* rootdist = redcostdata->rootToNodeDist;
   const PATH* nodeToTermpaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real cutoff = redcostdata->cutoff;
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];

   assert(scip && graph && redcost && rootdist && nodeToTermpaths && edgeIsDeletable && distdata && extpermanent);
   assert(edge >= 0 && edge < graph->edges);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph->mark[tail] && graph->mark[head]);
   assert(graph_isMarked(graph));
   assert(extreduce_extPermaIsClean(graph, extpermanent));

   *edgeIsDeletable = FALSE;

   /* can we extend from 'edge'? */
   if( extLeafIsExtendable(graph, isterm, tail) || extLeafIsExtendable(graph, isterm, head) )
   {
      int* extstack_data;
      int* extstack_start;
      int* extstack_state;
      int* tree_edges;
      int* tree_leaves;
      int* tree_parentNode;
      SCIP_Real* tree_parentEdgeCost;
      SCIP_Real* tree_redcostSwap;
      int* pseudoancestor_mark;
      const int nnodes = graph->knots;
      const int maxdfsdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;
      const int maxstackedges = MIN(nnodes / 2, STP_EXT_MAXEDGES);

      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_data, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_start, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &extstack_state, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_edges, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_leaves, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentNode, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_parentEdgeCost, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &tree_redcostSwap, nnodes) );

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &pseudoancestor_mark, nnodes) );

      {
         REDDATA reddata = { .cmst = extpermanent->cmst, .cgraph = extpermanent->cgraph,
            .cgraphEdgebuffer = extpermanent->cgraphEdgebuffer,
            .redCosts = redcost, .rootToNodeDist = rootdist, .nodeTo3TermsPaths = nodeToTermpaths,
            .nodeTo3TermsBases = redcostdata->nodeTo3TermsBases, .edgedeleted = extpermanent->edgedeleted,
            .pseudoancestor_mark = pseudoancestor_mark, .cutoff = cutoff, .equality = equality, .redCostRoot = root };
         EXTDATA extdata = { .extstack_data = extstack_data, .extstack_start = extstack_start,
            .extstack_state = extstack_state, .extstack_ncomponents = 0, .tree_leaves = tree_leaves,
            .tree_edges = tree_edges, .tree_deg = extpermanent->tree_deg, .tree_nleaves = 0,
            .tree_bottleneckDistNode = extpermanent->bottleneckDistNode, .tree_parentNode = tree_parentNode,
            .tree_parentEdgeCost = tree_parentEdgeCost, .tree_redcostSwap = tree_redcostSwap, .tree_redcost = 0.0,
            .tree_nDelUpArcs = 0, .tree_root = -1, .tree_nedges = 0, .tree_depth = 0, .extstack_maxsize = nnodes - 1,
            .pcSdToNode = extpermanent->pcSdToNode,
            .extstack_maxedges = maxstackedges, .tree_maxnleaves = STP_EXTTREE_MAXNLEAVES, .tree_maxdepth = maxdfsdepth,
            .tree_maxnedges = STP_EXTTREE_MAXNEDGES, .node_isterm = isterm, .reddata = &reddata, .distdata = distdata };

         /* can we extend from head? */
         if( extLeafIsExtendable(graph, isterm, head) )
            extCheckArc(scip, graph, edge, &extdata, edgeIsDeletable);

         /* try to extend from tail? */
         if( !(*edgeIsDeletable) && extLeafIsExtendable(graph, isterm, tail) )
            extCheckArc(scip, graph, flipedge(edge), &extdata, edgeIsDeletable);
      }

      assert(extreduce_extPermaIsClean(graph, extpermanent));

      SCIPfreeCleanBufferArray(scip, &pseudoancestor_mark);
      SCIPfreeBufferArray(scip, &tree_redcostSwap);
      SCIPfreeBufferArray(scip, &tree_parentEdgeCost);
      SCIPfreeBufferArray(scip, &tree_parentNode);
      SCIPfreeBufferArray(scip, &tree_leaves);
      SCIPfreeBufferArray(scip, &tree_edges);
      SCIPfreeBufferArray(scip, &extstack_state);
      SCIPfreeBufferArray(scip, &extstack_start);
      SCIPfreeBufferArray(scip, &extstack_data);
   }

   return SCIP_OKAY;
}


/** Extended reduction test for arcs.
 * This method will also set edgedeletable[a] to TRUE if arc 'a' can be deleted, but its anti-parallel arc not. */
SCIP_RETCODE extreduce_deleteArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   int*                  nelims              /**< number of eliminations */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const redcost = redcostdata->redEdgeCost;

   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   graph_mark(graph);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, graph, edgedeletable, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( edgeIsValid(graph, e) )
      {
         const int erev = e + 1;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) && SCIPisZero(scip, redcost[erev]) )
            continue;

         if( !edgedeletable[e] )
         {
            SCIP_Bool deletable;
            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, e, allowequality, &distdata, &extpermanent,
                  &deletable) );

            if( deletable )
               edgedeletable[e] = TRUE;
         }

         if( !edgedeletable[erev] )
         {
            SCIP_Bool erevdeletable = FALSE;

            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, erev, allowequality, &distdata, &extpermanent,
                  &erevdeletable) );

            if( erevdeletable )
               edgedeletable[erev] = TRUE;
         }

         if( edgedeletable[e] && edgedeletable[erev] )
         {
            assert(edgedeletable[e] && edgedeletable[erev]);

            removeEdge(scip, e, graph, &distdata);

            (*nelims)++;
         }
      }
   }

   extreduce_extPermaFreeMembers(scip, &extpermanent);
   extreduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

#ifndef NDEBUG
   for( int k = 0; k < graph->knots; k++ )
      if( graph->grad[k] == 0 && k != redcostdata->redCostRoot && !Is_term(graph->term[k]) )
         assert(!graph->mark[k]);
#endif

   return SCIP_OKAY;
}


/** extended reduction test for edges */
SCIP_RETCODE extreduce_deleteEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const redcost = redcostdata->redEdgeCost;

   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   graph_mark(graph);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, graph, edgedeletable, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( edgeIsValid(graph, e) )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) && SCIPisZero(scip, redcost[erev]) )
            continue;

         SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, e, allowequality, &distdata, &extpermanent, &deletable) );

         if( deletable )
         {
            removeEdge(scip, e, graph, &distdata);

            (*nelims)++;
         }
      }
   }

   extreduce_extPermaFreeMembers(scip, &extpermanent);
   extreduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

#ifndef NDEBUG
   for( int k = 0; k < graph->knots; k++ )
      if( graph->grad[k] == 0 && k != redcostdata->redCostRoot && !Is_term(graph->term[k]) )
         assert(!graph->mark[k]);
#endif

   return SCIP_OKAY;
}
