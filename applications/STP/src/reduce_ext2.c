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
#include "reduce.h"

#define EXT_ANCESTORS_MAX  16
#define EXT_STATE_NONE     0
#define EXT_STATE_EXPANDED 1
#define EXT_STATE_MARKED   2
#define EXT_REDCOST_NRECOMP 10

#define STP_EXT_MAXDFSDEPTH 6
#define STP_EXT_MINDFSDEPTH 4
#define STP_EXT_MAXGRAD 8
#define STP_EXT_MAXEDGES 500
#define STP_EXT_MAXTREESIZE 20
#define STP_EXT_MAXNLEAVES 20
#define STP_EXT_EDGELIMIT 50000

#define EXEDGE_FREE 0
#define EXEDGE_FIXED 1
#define EXEDGE_KILLED 2

/** reduction data */
typedef struct reduction_data
{
   const SCIP_Real* const reducedcosts;
   const SCIP_Real* const rootdist;
   const PATH* const termpaths;
   const STP_Bool* const edgedeleted;
   int*  const ancestormark;
   const SCIP_Real cutoff;
   const SCIP_Real treeredcostoffset;
   const SCIP_Bool equality;
} REDDATA;

/** extension data */
typedef struct extension_data
{
   int* const extstack_data;
   int* const extstack_start;
   int* const extstack_state;
   int* const tree_leaves;
   int* const tree_edges;
   int* const tree_deg;                      /**< array of size nnodes */
   SCIP_Real* const tree_bottleneckDistNode; /**< needs to be set to -1.0 (for all nodes) */
   int* const tree_parentNode;
   SCIP_Real* const tree_parentEdgeCost;
   SCIP_Real tree_redcost;
   int tree_nleaves;
   int tree_size;
   int tree_depth;
   int extstack_size;
   const int extstack_maxsize;
   const int extstack_maxedges;
   const int tree_maxnleaves;
   const int tree_maxdepth;
   const int tree_maxsize;
   REDDATA* const reddata;
   DISTDATA* const distdata;
} EXTDATA;

/** mark ancestors of given edge that might already be in conflict */
static
SCIP_Bool markEdgeAncestorsConflicting(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      const unsigned idx = ((unsigned) curr->index) / 2;

      assert(curr->index >= 0 && idx < (unsigned) (MAX(graph->edges, graph->orgedges) / 2));

      if( ancestormark[idx] )
         return TRUE;

      ancestormark[idx] = 1;
   }

   return FALSE;
}

/** unmark ancestors of given edge that might already be in conflict */
static
void unmarkEdgeAncestorsConflicting(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      ancestormark[((unsigned) curr->index) / 2] = 0;
   }
}

/** unmark ancestors of given edge */
static
void unmarkEdgeAncestors(
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge,               /**< edge to use */
   int*                  ancestormark        /**< ancestor mark array */
   )
{
   int count = 0;
   assert(edge >= 0);

   for( IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      const unsigned idx = ((unsigned) curr->index) / 2;

      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      assert(ancestormark[idx] == 1);

      ancestormark[idx] = 0;
   }
}


#ifndef NDEBUG
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
   const int treesize = extdata->tree_size;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   int leavescount;

   SCIP_Bool flawed = FALSE;

   assert(nleaves >= 1);

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &edgecount, nedges) );
   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &degreecount, nnodes) );

   for( int i = 0; i < treesize; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];
      const int tail = graph->tail[e];

      assert(e >= 0 && e < nedges);

      if( edgecount[e] > 0 || tree_deg[tail] <= 0 || tree_deg[head] <= 0 )
         flawed = TRUE;

      degreecount[tail]++;
      degreecount[head]++;

      edgecount[e]++;
   }


   leavescount = 1; /* for tail of initial edge */

   /* degree check */
   for( int i = 0; i < treesize && !flawed; i++ )
   {
      const int e = tree_edges[i];
      const int head = graph->head[e];

      if( degreecount[head] == 1 )
         leavescount++;

      if( degreecount[head] != extdata->tree_deg[head] )
         flawed = TRUE;
   }

   /* leaves check */
   if( !flawed && leavescount != nleaves )
      flawed = TRUE;

   for( int i = 0; i < nleaves && !flawed; i++ )
   {
      const int leaf = tree_leaves[i];
      if( degreecount[leaf] != 1 )
         flawed = TRUE;
   }

   /* clean-up */
   for( int i = 0; i < treesize; i++ )
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


/** should we truncate from current component? */
static
void printStack(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
#ifdef SCIP_DEBUG
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extdata->extstack_size - 1;

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


inline static
SCIP_Bool extTreeEdgeAncestorConflict(
   const REDDATA*        reddata,            /**< reduction data */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge                /**< edge to check for conflict */
)
{
   const int* const ancestormark = reddata->ancestormark;
   int count = 0;
   for( const IDX* curr = graph->ancestors[edge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
   {
      assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
      if( ancestormark[((unsigned) curr->index) / 2] )
      {
#ifdef SCIP_DEBUG
         printf("conflict found for ");
         graph_edge_printInfo(graph, edge);
#endif
         return TRUE;
      }
   }
   return FALSE;
}

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

/** finds position of given leaf in leaves data */
static
int extFindLeafPos(
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


/** computes the tree bottleneck between vertex1 and vertex2 in the current tree */
static
SCIP_Real extTreeGetBottleneckDist(
   const EXTDATA*        extdata,            /**< extension data */
   int                   vertex1,            /**< first vertex */
   int                   vertex2             /**< second vertex */
   )
{
   SCIP_Real* bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const SCIP_Real* parentEdgeCost = extdata->tree_parentEdgeCost;
   const int* parentNode = extdata->tree_parentNode;
   int* const tree_deg = extdata->tree_deg;
   SCIP_Real bottleneck;
   int currentNode;
   int childNode;

   assert(bottleneckDist_node && parentEdgeCost && parentNode);

   /* go down from vertex1 */
   bottleneck = 0.0;
   childNode = vertex1;
   currentNode = parentNode[vertex1];

   while( currentNode != - 1 )
   {
      assert(currentNode >= 0);
      assert(tree_deg[currentNode] >= 0);
      assert(parentEdgeCost[childNode] >= 0.0);
      assert(bottleneckDist_node[currentNode] == -1.0);

      if( tree_deg[childNode] == 1 )
         bottleneck += parentEdgeCost[childNode];
      else
         bottleneck = parentEdgeCost[childNode];

      bottleneckDist_node[currentNode] = bottleneck;

      childNode = currentNode;
      currentNode = parentNode[currentNode];
   }

   /* go down from vertex2 up to least common ancestor with vertex1 */
   bottleneck = 0.0;
   childNode = vertex2;
   currentNode = parentNode[vertex2];

   while( bottleneckDist_node[currentNode] < -0.5 )
   {
      assert(tree_deg[currentNode] >= 0);
      assert(parentEdgeCost[childNode] >= 0.0);
      assert(bottleneckDist_node[currentNode] == -1.0);

      if( tree_deg[childNode] == 1 )
         bottleneck += parentEdgeCost[childNode];
      else
         bottleneck = parentEdgeCost[childNode];

      childNode = currentNode;
      currentNode = parentNode[currentNode];

      assert(currentNode >= 0);
   }

   assert(bottleneckDist_node[childNode] == -1.0);
   assert(bottleneckDist_node[currentNode] >= 0.0);

   if( tree_deg[childNode] == 1 )
      bottleneck += parentEdgeCost[childNode];
   else
      bottleneck = parentEdgeCost[childNode];

   bottleneck = MAX(bottleneck, bottleneckDist_node[currentNode]);

   /* go down from vertex1 and reset bottleneckDist_node */
   for( currentNode = parentNode[vertex1]; currentNode != -1; currentNode = parentNode[currentNode]  )
   {
      assert(currentNode >= 0);
      assert(tree_deg[currentNode] >= 0);
      assert(bottleneckDist_node[currentNode] >= -0.0);

      bottleneckDist_node[currentNode] = -1.0;
   }

   int todo;
   bottleneck = 0.0;

   return bottleneck;
}

/** adds top component of stack to tree */
static
void extTreeAddStackTop(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            conflict            /**< conflict found? */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   int* const tree_edges = extdata->tree_edges;
   int* const tree_leaves = extdata->tree_leaves;
   int* const tree_deg = extdata->tree_deg;
   int* const tree_parentNode = extdata->tree_parentNode;
   SCIP_Real* const tree_parentEdgeCost = extdata->tree_parentEdgeCost;
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real* const redcost = reddata->reducedcosts;
   const int stackpos = extdata->extstack_size - 1;
   const int comproot = graph->tail[extstack_data[extstack_start[stackpos]]];

   assert(!(*conflict));
   assert(stackpos >= 0);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);
   assert(comproot >= 0 && comproot < graph->knots);
   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   /* not the initial edge? */
   if( tree_deg[comproot] != graph->knots )
   {
      /* update tree leaves array */

      int comprootpos;

      assert(tree_deg[comproot] == 1);

      /* switch last leaf and root component */
      extdata->tree_nleaves--;
      assert(extdata->tree_nleaves > 0);

      comprootpos = extFindLeafPos(extdata, comproot, extdata->tree_nleaves);
      assert(comprootpos > 0);

      tree_leaves[comprootpos] = tree_leaves[extdata->tree_nleaves];
   }

   /* add top expanded component to tree data */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      assert(extdata->tree_size < extdata->extstack_maxsize);
      assert(edge >= 0 && edge < graph->edges);
      assert(tree_deg[graph->head[edge]] == 0 && tree_deg[graph->tail[edge]] > 0);

      extdata->tree_redcost += redcost[edge];

      tree_edges[(extdata->tree_size)++] = edge;
      tree_leaves[(extdata->tree_nleaves)++] = head;
      tree_deg[head] = 1;
      tree_parentNode[head] = comproot;
      tree_parentEdgeCost[head] = graph->cost[edge];

      /* not the initial edge? */
      if( stackpos > 0 )
         tree_deg[graph->tail[edge]]++;

      // todo break and use unmark conflict in extBacktrack
      if( markEdgeAncestorsConflicting(graph, edge, reddata->ancestormark) )
         *conflict = TRUE;
   }

   extdata->tree_depth++;

   assert(!extTreeIsFlawed(scip, graph, extdata));
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
   const int stackposition = extdata->extstack_size - 1;
   REDDATA* const reddata = extdata->reddata;

   assert(scip && graph && extdata && reddata && nupdatestalls && conflict);
   assert(!(*conflict));

   printStack(graph, extdata);

   /* is current component expanded? */
   if( extdata->extstack_state[stackposition] == EXT_STATE_EXPANDED )
      extTreeAddStackTop(scip, graph, extdata, conflict); /* add component to tree */

   /* recompute reduced costs? */
   if( ++(*nupdatestalls) > EXT_REDCOST_NRECOMP )
   {
      SCIP_Real treecost = reddata->treeredcostoffset;
      const SCIP_Real* const redcost = reddata->reducedcosts;
      const int* const tree_edges = extdata->tree_edges;
      const int tree_size = extdata->tree_size;

      *nupdatestalls = 0;

      assert(!extTreeIsFlawed(scip, graph, extdata));

      for( int i = 0; i < tree_size; i++ )
      {
         const int edge = tree_edges[i];
         assert(edge >= 0 && edge < graph->edges);

         treecost += redcost[edge];
      }

      assert(SCIPisEQ(scip, treecost, extdata->tree_redcost));

      extdata->tree_redcost = treecost;
   }
}

/** should we truncate from current component? */
static
SCIP_Bool extTruncate(
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extdata->extstack_size - 1;

   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED);

   if( extdata->tree_depth >= extdata->tree_maxdepth )
   {
      SCIPdebugMessage("truncate (depth) \n");
      return TRUE;
   }

   if( extdata->tree_size >= extdata->tree_maxsize )
   {
      SCIPdebugMessage("truncate (tree size) \n");
      return TRUE;
   }

   if( extdata->tree_nleaves >= extdata->tree_maxnleaves )
   {
      SCIPdebugMessage("truncate (number of leaves) \n");
      return TRUE;
   }

   if( extstack_start[stackpos] >= extdata->extstack_maxedges )
   {
      SCIPdebugMessage("truncate (edges on stack) \n");
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

/** get peripheral reduced cost of current tree including  */
static
SCIP_Real extTreeGetRedcosts(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   extedge             /**< edge for extension or -1 */
)
{
   REDDATA* const reddata = extdata->reddata;
   const int* const tree_leaves = extdata->tree_leaves;
   const PATH* const termpaths = reddata->termpaths;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Real tree_redcost = extdata->tree_redcost; /* includes reduced costs to initial tail */

   assert(graph && reddata && extdata);
   assert(nleaves > 1);

   for( int i = 1; i < nleaves; i++ )
   {
      const int leaf = tree_leaves[i];
      assert(extdata->tree_deg[leaf] == 1);

      tree_redcost += termpaths[leaf].dist;
   }

   // todo for the general case we also must consider the case that tail[edge] is the root!
   if( extedge != -1 )
   {
      const int base = graph->tail[extedge];
      const int extvert = graph->head[extedge];
      const SCIP_Real* const redcost = reddata->reducedcosts;

      assert(extedge >= 0 && extedge < graph->edges);
      assert(extdata->tree_deg[base] == 1);

      tree_redcost += redcost[extedge] + termpaths[extvert].dist - termpaths[base].dist;
   }

   return tree_redcost;
}

/** can any extension via edge be ruled out? */
static
SCIP_Bool extRuleOutSimple(
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
      return TRUE;

   if( reddata->edgedeleted && reddata->edgedeleted[edge] )
      return TRUE;
   else
   {
      const SCIP_Real tree_redcost = extTreeGetRedcosts(graph, extdata, edge);
      const SCIP_Real cutoff = reddata->cutoff;

      if( reddata->equality ? (SCIPisGE(scip, tree_redcost, cutoff)) : SCIPisGT(scip, tree_redcost, cutoff) )
         return TRUE;

      if( extTreeEdgeAncestorConflict(reddata, graph, edge) )
         return TRUE;
   }

   return FALSE;
}

/** can subtree be peripherally ruled out? */
static
SCIP_Bool extRuleOutPeriph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   REDDATA* const reddata = extdata->reddata;
   const SCIP_Real tree_redcost = extTreeGetRedcosts(graph, extdata, -1);
   const SCIP_Real cutoff = reddata->cutoff;

   if( reddata->equality ? (SCIPisGE(scip, tree_redcost, cutoff)) : SCIPisGT(scip, tree_redcost, cutoff) )
   {
      SCIPdebugMessage("Rule-out periph (red. cost) \n");
      return TRUE;
   }
   else
   {
      // todo extra method?
      /* tree bottleneck test */
      {
         DISTDATA* const distdata = extdata->distdata;
         const int stackpos = extdata->extstack_size - 1;
         const int* const extstack_data = extdata->extstack_data;
         const int* const extstack_start = extdata->extstack_start;
         const int* const leaves = extdata->tree_leaves;
         const int nleaves = extdata->tree_nleaves;

         for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
         {
            const int extleaf = graph->head[extstack_data[i]];

            assert(extleaf >= 0 && extleaf < graph->knots);

            for( int j = 0; j < nleaves; j++ )
            {
               SCIP_Real specialDist;
               const int leaf = leaves[j];

               assert(leaf != extleaf && leaf >= 0 && leaf < graph->knots);

               specialDist = reduce_distDataGetSD(distdata, extleaf, leaf);

               /* could a valid special distance be found?  */
               if( specialDist >= -0.5  )
               {
                  const SCIP_Real treeBottleneckDist = extTreeGetBottleneckDist(extdata, extleaf, leaf);

                  if( SCIPisLT(scip, specialDist, treeBottleneckDist) ) /* todo cover equality */
                     return TRUE;
               }
            }
         }
      }

      // todo do MST test

#ifndef NDEBUG
      {
         const int stackpos = extdata->extstack_size - 1;
         const int* const extstack_data = extdata->extstack_data;
         const int* const extstack_start = extdata->extstack_start;
         const int* const ancestormark = reddata->ancestormark;

         for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
         {
            const int curredge = extstack_data[i];
            int count = 0;
            for( IDX* curr = graph->ancestors[curredge]; curr != NULL && count <= EXT_ANCESTORS_MAX; curr = curr->parent, count++ )
            {
               assert(curr->index >= 0 && curr->index / 2 < (MAX(graph->edges, graph->orgedges) / 2));
               assert(ancestormark[((unsigned) curr->index) / 2] == 1);
            }
         }
      }
#endif
   }

   return FALSE;
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
   REDDATA* const reddata = extdata->reddata;
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int* const tree_deg = extdata->tree_deg;
   int stackpos = extdata->extstack_size - 1;
   const int stackstart = extstack_start[stackpos];

   assert(graph && reddata && extdata);
   assert(extstack_start[stackpos + 1] - extstack_start[stackpos] > 0);

   /* top component already expanded? */
   if( extstack_state[stackpos] != EXT_STATE_NONE )
   {
      const SCIP_Real* const redcost = reddata->reducedcosts;
      int* const tree_leaves = extdata->tree_leaves;
      const int comproot = graph->tail[extstack_data[extstack_start[stackpos]]];
      const int compsize = extstack_start[stackpos + 1] - extstack_start[stackpos];

      assert(compsize > 0);
      assert(tree_deg[comproot] > 1);

      /* remove top component */
      for( int i = stackstart; i < extstack_start[stackpos + 1]; i++ )
      {
         const int edge = extstack_data[i];
         const int head = graph->head[edge];
         const int tail = graph->tail[edge];

         assert(edge >= 0 && edge < graph->edges);
         assert(tree_deg[head] == 1 && tree_deg[tail] > 1);

         (extdata->tree_redcost) -= redcost[edge];
         tree_deg[head] = 0;
         tree_deg[tail]--;

         if( ancestor_conflict )
            unmarkEdgeAncestorsConflicting(graph, edge, reddata->ancestormark);
         else
            unmarkEdgeAncestors(graph, edge, reddata->ancestormark);
      }

      (extdata->tree_size) -= compsize;
#if 0
      for( int k = extdata->tree_nleaves - 1; k >= extdata->tree_nleaves - compsize; k-- )
         printf("bt remove leaf %d \n", tree_leaves[k]);
#endif
      (extdata->tree_nleaves) -= compsize;
      (extdata->tree_depth)--;

      /* add component root to leaves array and remove current one */
      assert(tree_deg[comproot] == 1);
      tree_leaves[extdata->tree_nleaves++] = comproot;

      assert(extdata->tree_size >= 0 && extdata->tree_depth >= 0);
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

   extdata->extstack_size = stackpos + 1;

   assert(!extTreeIsFlawed(scip, graph, extdata));
}

/** expands top component of stack (backtracks if stack is full) */
static
void extStackExpand(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extdata->extstack_size - 1;
   int datasize = extstack_start[stackpos];
   const int setsize = extstack_start[stackpos + 1] - extstack_start[stackpos];
   const uint32_t powsize = pow(2, setsize);

   assert(extdata && scip && graph && isterm && success);
   assert(setsize <= STP_EXT_MAXGRAD);
   assert(setsize > 0 && setsize <= 32);
   assert(stackpos >= 1);
   assert(extstack_state[stackpos] == EXT_STATE_NONE);

   /* stack too full? */
   if( (datasize + (int) pow(2, setsize)) > extdata->extstack_maxsize )
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
      for( int j = 0; j < setsize; j++ )
      {
         /* Check if jth bit in counter is set */
         if( counter & (1 << j) )
         {
            extstack_data[datasize++] = extedges[j];
            SCIPdebugMessage("  %d \n", graph->head[extedges[j]]);
         }
      }

      SCIPdebugMessage("... added \n");
      assert(stackpos < extdata->extstack_maxsize - 1);

      extstack_state[stackpos] = EXT_STATE_EXPANDED;
      extstack_start[++stackpos] = datasize;

      assert(extstack_start[stackpos] - extstack_start[stackpos - 1] > 0);
   }

   assert(stackpos >= extdata->extstack_size);

   extdata->extstack_size = stackpos;
}

/** extend top component of stack (backtracks if stack is full) */
static
void extExtend(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            success             /**< success pointer */
)
{
   int extedges[STP_EXT_MAXGRAD * STP_EXT_MAXGRAD];
   int extedgesstart[STP_EXT_MAXGRAD + 1];
   int* const extstack_data = extdata->extstack_data;
   int* const extstack_start = extdata->extstack_start;
   int* const extstack_state = extdata->extstack_state;
   int stackpos = extdata->extstack_size - 1;
   int nfullextensions;
   int nsingleextensions;

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

   /* loop over all leaves of extension */
   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int leaf = graph->head[extstack_data[i]];

      assert(extstack_data[i] >= 0 && extstack_data[i] < graph->edges);

      /* extensions from leaf not possible? */
      if( !extLeafIsExtendable(graph, isterm, leaf) )
         continue;

      /* assemble feasible single edge extensions from leaf */
      for( int e = graph->outbeg[leaf]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( !extRuleOutSimple(scip, graph, extdata, e) )
         {
            assert(nsingleextensions < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD);
            extedges[nsingleextensions++] = e;
         }
#ifdef SCIP_DEBUG
         else
         {
            printf("simple rule out: ");
            graph_edge_printInfo(graph, e);
         }
#endif
      }

      extedgesstart[++nfullextensions] = nsingleextensions;
   }

   assert(nfullextensions <= STP_EXT_MAXGRAD);

   /* found no valid extensions? */
   if( nfullextensions == 0 )
   {
      *success = FALSE;
   }
   /* found valid extensions, but all ruled out already? */
   else if( nsingleextensions == 0 )
   {
      *success = TRUE;
   }
   /* found non-empty valid extensions */
   else
   {
      int datasize = extstack_start[stackpos + 1];
      int extsize[STP_EXT_MAXGRAD];
      int extindex[STP_EXT_MAXGRAD];

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

      for( int i = extstack_start[extdata->extstack_size]; i < extstack_start[stackpos + 1]; i++ )
         graph_edge_printInfo(graph, extstack_data[i]);
#endif

      extdata->extstack_size = stackpos + 1;

      *success = TRUE;

      /* try to expand last (smallest) component */
      extStackExpand(scip, graph, isterm, extdata, success);
   }
}


/** check (directed) arc */
SCIP_RETCODE reduce_extendedCheckArc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data structures */
   const STP_Bool*       edgedeleted,        /**< edge array to mark which directed edge can be removed or NULL */
   const SCIP_Bool*      isterm,             /**< marks whether node is a terminal (or proper terminal for PC) */
   int                   edge,               /**< directed edge to be checked */
   SCIP_Bool             equality,           /**< delete edge also in case of reduced cost equality? */
   DISTDATA*             distdata,           /**< data for distance computations */
   SCIP_Real*            bottleneckDistNode, /**< needs to be set to -1.0 (size nnodes) */
   int*                  tree_deg,           /**< -1 for forbidden nodes (e.g. PC terminals), 0 otherwise; in method: position ( > 0) for nodes in tree */
   SCIP_Bool*            deletable           /**< is edge deletable? */
)
{
   const int root = redcostdata->redCostRoot;
   const SCIP_Real* redcost = redcostdata->redEdgeCost;
   const SCIP_Real* rootdist = redcostdata->rootToNodeDist;
   const PATH* node2termpaths = redcostdata->nodeTo3TermsPaths;
   const SCIP_Real cutoff = redcostdata->cutoff;
   const int head = graph->head[edge];
   const int tail = graph->tail[edge];
   SCIP_Real edgebound = redcost[edge] + rootdist[tail] + node2termpaths[head].dist;

#ifndef NDEBUG
   assert(scip && graph && redcost && rootdist && node2termpaths && deletable && distdata && tree_deg && bottleneckDistNode);
   assert(edge >= 0 && edge < graph->edges);
   assert(!graph_pc_isPcMw(graph) || !graph->extended);
   assert(graph->mark[tail] && graph->mark[head]);

   if( !graph_pc_isPcMw(graph) )
      for( int k = 0; k < graph->knots; k++ )
         assert(graph->mark[k] == (graph->grad[k] > 0));
#endif

   /* trivial rule-out? */
   if( SCIPisGT(scip, edgebound, cutoff) || (equality && SCIPisEQ(scip, edgebound, cutoff)) || head == root )
   {
      *deletable = TRUE;
      return SCIP_OKAY;
   }

   *deletable = FALSE;

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
      int* ancestormark;
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

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &ancestormark, (MAX(graph->edges, graph->orgedges) / 2)) );

      tree_deg[root] = -1;

      /* can we extend from head? */
      if( extLeafIsExtendable(graph, isterm, head) )
      {
         const SCIP_Real treeredcostoffset = rootdist[tail];
         int nupdatestalls = 0;
         SCIP_Bool success = TRUE;
         SCIP_Bool conflict = FALSE;
         REDDATA reddata = {redcost, rootdist, node2termpaths, edgedeleted, ancestormark, cutoff, treeredcostoffset, equality};
         EXTDATA extdata = {extstack_data, extstack_start, extstack_state, tree_leaves, tree_edges,
            tree_deg, bottleneckDistNode, tree_parentNode, tree_parentEdgeCost, 0.0, 0, 0, 0, 0, nnodes - 1, maxstackedges,
            STP_EXT_MAXNLEAVES,  maxdfsdepth, STP_EXT_MAXTREESIZE, &reddata, distdata};

         extdata.tree_redcost = treeredcostoffset;
         extdata.tree_depth = 0;
         tree_deg[tail] = nnodes;

         /* put 'edge' on the stack */
         extdata.extstack_size = 1;
         extstack_start[0] = 0;
         extstack_start[1] = 1;
         extstack_data[0] = edge;
         extstack_state[0] = EXT_STATE_EXPANDED;
         tree_leaves[0] = tail;
         tree_parentNode[tail] = -1;
         extdata.tree_nleaves = 1;

         int todo; // internal method here that takes pointer to REDATA and EXTDATA
         extTreeSyncWithStack(scip, graph, &extdata, &nupdatestalls, &conflict);

         assert(!conflict);
         assert(tree_parentNode[head] == tail && tree_parentEdgeCost[head] == graph->cost[edge]);

         extExtend(scip, graph, isterm, &extdata, &success);

         assert(extstack_state[0] == EXT_STATE_MARKED);
         assert(success || 1 == extdata.extstack_size);

         /* limited DFS backtracking; stops once back at 'edge' */
         while( extdata.extstack_size > 1 )
         {
            const int stackposition = extdata.extstack_size - 1;
            conflict = FALSE;

            extTreeSyncWithStack(scip, graph, &extdata, &nupdatestalls, &conflict);

            /* has current component already been extended? */
            if( extstack_state[stackposition] == EXT_STATE_MARKED )
            {
               extBacktrack(scip, graph, success, FALSE, &extdata);
               continue;
            }

            /* component not expanded yet? */
            if( extstack_state[stackposition] != EXT_STATE_EXPANDED )
            {
               assert(extstack_state[stackposition] == EXT_STATE_NONE);

               extStackExpand(scip, graph, isterm, &extdata, &success);
               continue;
            }

            assert(extstack_state[stackposition] == EXT_STATE_EXPANDED);

            if( conflict || extRuleOutPeriph(scip, graph, &extdata) )
            {
               success = TRUE;
               extBacktrack(scip, graph, success, conflict, &extdata);
               continue;
            }

            if( extTruncate(graph, isterm, &extdata) )
            {
               success = FALSE;
               extBacktrack(scip, graph, success, FALSE, &extdata);
               continue;
            }

            /* neither ruled out nor truncated, so extend */
            extExtend(scip, graph, isterm, &extdata, &success);

         } /* DFS loop */

         *deletable = success;
         assert(tree_deg[head] == 1 && tree_deg[tail] == nnodes);

         tree_deg[head] = 0;
         tree_deg[tail] = 0;

         unmarkEdgeAncestors(graph, edge, ancestormark);
      } /* extend from head */


      /* finalize arrays */

      tree_deg[root] = 0;

#ifndef NDEBUG
      for( int i = 0; i < MAX(graph->edges, graph->orgedges) / 2; i++ )
         assert(ancestormark[i] == 0);

      for( int i = 0; i < nnodes; i++ )
      {
         assert(tree_deg[i] == 0 || tree_deg[i] == -1);
         assert(bottleneckDistNode[i] == -1.0);
      }
#endif
      SCIPfreeCleanBufferArray(scip, &ancestormark);
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

/** extended reduction test for edges */
SCIP_RETCODE reduce_extendedEdge2(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   const PATH*           termpaths,          /**< paths to nearest terminals */
   const SCIP_Real*      redcost,            /**< dual ascent costs */
   const SCIP_Real*      rootdist,           /**< shortest path distance (w.r.t reduced costs) from root to any node */
   const int*            result,             /**< solution array */
   SCIP_Real             minpathcost,        /**< the required reduced path cost to be surpassed */
   int                   root,               /**< the root */
   SCIP_Bool             markdirected,       /**< try to also mark edge if anti-parallel is not marked */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   int*                  nelims              /**< number of eliminations */
)
{
   SCIP_Bool* isterm;
   int* tree_deg;
   SCIP_Real* bottleneckDist;
   const int nedges = graph->edges;
   const int nnodes = graph->knots;
   const SCIP_Bool pcmw = graph_pc_isPcMw(graph);
   REDCOST redcostdata = {redcost, rootdist, termpaths, minpathcost, root};
   DISTDATA distdata;

   assert(scip && graph && redcost);
   assert(!pcmw || !graph->extended);
   assert(root >= 0 && root < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, minpathcost) )
      return SCIP_OKAY;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( reduce_distDataInit(scip, graph, EXT_CLOSENODES_MAXN, FALSE, &distdata) );

   SCIP_CALL( SCIPallocBufferArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_deg, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bottleneckDist, nnodes) );

   graph_get_isTerm(graph, isterm);

   if( !pcmw )
      for( int k = 0; k < nnodes; k++ )
         graph->mark[k] = (graph->grad[k] > 0);

   for( int k = 0; k < nnodes; k++ )
   {
      bottleneckDist[k] = -1.0;
      if( graph->mark[k] )
         tree_deg[k] = 0;
      else
         tree_deg[k] = -1;
   }

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      const int tail = graph->tail[e];
      const int head = graph->head[e];

      if( pcmw && (!graph->mark[tail] || !graph->mark[head]) )
         continue;

      if( graph->oeat[e] != EAT_FREE )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev);
         assert(SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) || SCIPisZero(scip, redcost[erev]) )
            continue;

         if( !edgedeletable[e] )
         {
            SCIP_CALL( reduce_extendedCheckArc(scip, graph, &redcostdata, edgedeletable,
                  isterm, e, allowequality, &distdata, bottleneckDist, tree_deg, &deletable) );

            if( deletable )
               edgedeletable[e] = TRUE;
         }

         if( !edgedeletable[erev] && (deletable || markdirected) )
         {
            SCIP_Bool erevdeletable = TRUE;

            SCIP_CALL( reduce_extendedCheckArc(scip, graph, &redcostdata, edgedeletable,
                  isterm, erev, allowequality, &distdata, bottleneckDist, tree_deg, &erevdeletable) );

            if( erevdeletable )
               edgedeletable[erev] = TRUE;

            deletable = (deletable && erevdeletable);
         }

         if( deletable )
         {
            graph_edge_delFull(scip, graph, e, TRUE);

            if( graph->grad[tail] == 0 )
               graph->mark[tail] = FALSE;

            if( graph->grad[head] == 0 )
               graph->mark[head] = FALSE;

            (*nelims)++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &bottleneckDist);
   SCIPfreeBufferArray(scip, &tree_deg);
   SCIPfreeBufferArray(scip, &isterm);

   reduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
      if( graph->grad[k] == 0 && k != root )
         assert(!graph->mark[k]);
#endif

   return SCIP_OKAY;
}
