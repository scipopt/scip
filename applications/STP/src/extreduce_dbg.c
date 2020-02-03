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

/**@file   extreduce_dbg.c
 * @brief  extended reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements extended reduction debugging routines for several Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
// #define STP_DEBUG_EXTPC // use if special sds for PC are deactivated

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"


/** get SD MST weight
 *  NOTE: might deviate because only getSd is used...maybe only use double in the code? */
static
SCIP_Real sdmstGetWeight(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const int*            nodes,              /**< nodes (from graph) for MST computation */
   int                   nnodes,             /**< number of nodes for MST computation*/
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real mstweight;
   CGRAPH* cgraph;
   CMST* cmst;
   SCIP_Real* adjcosts;

   SCIP_CALL_ABORT( cgraph_init(scip, &cgraph, STP_EXTTREE_MAXNLEAVES_GUARD) );
   SCIP_CALL_ABORT( cmst_init(scip, &cmst, STP_EXTTREE_MAXNLEAVES_GUARD) );

   adjcosts = cgraph->adjedgecosts;

   assert(adjcosts);
   assert(nnodes > 0 && nnodes <= STP_EXTTREE_MAXNLEAVES_GUARD);

   /* build the MST graph */
   for( int i = 0; i < nnodes; i++ )
      cgraph_node_append(cgraph, i);

   for( int i = 0; i < nnodes; i++ )
   {
      const int startnode = nodes[i];

      for( int j = 0; j < nnodes; j++ )
      {
         SCIP_Real specialDist;

         if( i == j )
         {
            specialDist = FARAWAY;
         }
         else
         {
            const int endnode = nodes[j];
            specialDist = extreduce_extGetSd(scip, graph, startnode, endnode, extdata);

            if( specialDist <= 0.0 )
            {
               assert(EQ(specialDist, -1.0));
               specialDist = BLOCKED;
            }
         }

         adjcosts[j] = specialDist;
      }

      cgraph_node_applyMinAdjCosts(cgraph, i, i);
   }

   /* compute the MST */
   cmst_computeMst(cgraph, 0, cmst);

   mstweight = cmst->mstobj;

   cmst_free(scip, &cmst);
   cgraph_free(scip, &cgraph);

   assert(GE(mstweight, 0.0));

   return mstweight;
}


/** Helper.
 *  Gives maximum cost among all close nodes. */
static inline
SCIP_Real distCloseNodesGetMaxCost(
   int                   vertex,             /**< vertex for which to get the maximum cost*/
   const DISTDATA*       distdata            /**< distance data */
   )
{
   const SCIP_Real* const distances = distdata->closenodes_distances;
   const RANGE* const range = distdata->closenodes_range;
   const int range_start = range[vertex].start;
   const int range_end = range[vertex].end;
   SCIP_Real maxcost = -FARAWAY;

   for( int i = range_start; i < range_end; ++i )
   {
      const SCIP_Real dist = distances[i];

      if( dist > maxcost )
         maxcost = dist;
   }

   assert(GE(maxcost, 0.0) || range_start == range_end);

   return maxcost;
}


/** Gets close nodes and corresponding distances.
 *  NOTE: needs to correspond to 'distDataComputeCloseNodes' in 'extreduce_dbg.c' */
static
SCIP_RETCODE distCloseNodesCompute(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< distance data */
   int                   startvertex,        /**< start vertex */
   int*                  closenodes_indices, /**< indices of close nodes */
   SCIP_Real*            closenodes_dists,   /**< distances of close nodes */
   int*                  nclosenodes         /**< number of added close nodes */
   )
{
   SCIP_Real* dist;
   DHEAP* dheap = NULL;
   int* state;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const int nnodes = g->knots;
   const SCIP_Real closenodes_maxcost = distCloseNodesGetMaxCost(startvertex, distdata);
   int clodenode_count;

   assert(dcsr && g && distdata);
   assert(startvertex >= 0 && startvertex < g->knots);

   SCIP_CALL( graph_heap_create(scip, nnodes, NULL, NULL, &dheap) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &dist, nnodes) );

   state = dheap->position;

   for( int k = 0; k < nnodes; k++ )
   {
      dist[k] = FARAWAY;
      assert(state[k] == UNKNOWN);
   }

   clodenode_count = 0;
   dist[startvertex] = 0.0;
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
         if( GT(dist[k], closenodes_maxcost) )
            break;

         closenodes_indices[clodenode_count] = k;
         closenodes_dists[clodenode_count] = dist[k];

         clodenode_count++;
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
               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
      }
   }

   SCIPfreeMemoryArray(scip, &dist);
   graph_heap_free(scip, TRUE, TRUE, &dheap);

   *nclosenodes = clodenode_count;

   return SCIP_OKAY;
}

/* Is each original close node to 'vertex' included in the 'closenodes_indices' array?
 * And if so, with correct costs? */
static
SCIP_Bool distCloseNodesIncluded(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata,           /**< distance data */
   int                   vertex,             /**< vertex for check */
   const int*            closenodes_indices, /**< indices of newly found close nodes */
   const SCIP_Real*      closenodes_dists,   /**< distances of newly found close nodes */
   int                   nclosenodes         /**< number of newly found close nodes */
)
{
   RANGE* const org_range = distdata->closenodes_range;
   int* const org_indices = distdata->closenodes_indices;
   SCIP_Real* const org_dists = distdata->closenodes_distances;
   int* newnodeindex;
   const int nnodes = graph_get_nNodes(g);
   const int org_start = org_range[vertex].start;
   const int org_end = org_range[vertex].end;
   const int org_nclosenodes = org_end - org_start;
   SCIP_Bool isIncluded = TRUE;

   if( nclosenodes < org_nclosenodes )
   {
      SCIPdebugMessage("too few new closenodes! %d < %d \n", nclosenodes, org_nclosenodes);
      return FALSE;
   }

   SCIP_CALL_ABORT( SCIPallocCleanBufferArray(scip, &newnodeindex, nnodes) );

   /* mark the indices of new close nodes */
   for( int j = 0; j < nclosenodes; j++ )
   {
      const int cnode = closenodes_indices[j];

      assert(newnodeindex[cnode] == 0);

      newnodeindex[cnode] = j + 1;
   }

   /* the actual checks */
   for( int i = org_start; i < org_end; ++i )
   {
      const SCIP_Real org_dist = org_dists[i];
      const int org_node = org_indices[i];
      const int new_index = newnodeindex[org_node] - 1;

      assert(new_index >= -1);

      if( new_index == -1 )
      {
         SCIPdebugMessage("could not find vertex %d in new close nodes \n", org_node);
#ifdef SCIP_DEBUG
         printf("new nodes: \n");

         for( int k = 0; k < nclosenodes; ++k )
         {
            printf("(idx=%d dist=%f) ", k, closenodes_dists[k]);
            graph_knot_printInfo(g, closenodes_indices[k]);
         }

         printf("original nodes: \n");

         for( int k = org_start; k < org_end; ++k )
         {
            printf("(dist=%f) ", org_dists[k]);
            graph_knot_printInfo(g, org_indices[k]);
         }
#endif

         isIncluded = FALSE;
         break;
      }

      assert(org_node == closenodes_indices[new_index]);

      if( !EQ(closenodes_dists[new_index], org_dist) )
      {
         SCIPdebugMessage("wrong distances: %f != %f \n", org_dists[i], closenodes_dists[new_index]);

         isIncluded = FALSE;
         break;
      }
   }

   /* unmark the new close nodes */
   for( int j = 0; j < nclosenodes; j++ )
   {
      const int cnode = closenodes_indices[j];

      assert(newnodeindex[cnode] == j + 1);

      newnodeindex[cnode] = 0;
   }

   SCIPfreeCleanBufferArray(scip, &newnodeindex);

   return isIncluded;
}


/** is current tree flawed? */
SCIP_Bool extreduce_treeIsFlawed(
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

   for( int i = 0; i < nnodes; i++ )
   {
      assert(degreecount[i] == 0);
   }

   for( int i = 0; i < nedges; i++ )
   {
      assert(edgecount[i] == 0);
   }

   SCIPfreeCleanBufferArray(scip, &degreecount);
   SCIPfreeCleanBufferArray(scip, &edgecount);

   return flawed;
}


/** is current tree completely hashed? */
SCIP_Bool extreduce_treeIsHashed(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
   const REDDATA* const reddata = extdata->reddata;

   for( int i = 0; i < extdata->tree_nedges; i++ )
   {
      const int edge = extdata->tree_edges[i];
      const int nAncestors = graph_edge_nPseudoAncestors(graph, edge);

      assert(nAncestors >= 0);

      if( nAncestors == 0 )
         continue;

      if( !graph_pseudoAncestors_edgeIsHashed(graph->pseudoancestors, edge, reddata->pseudoancestor_mark) )
         return FALSE;
   }

   return TRUE;
}


/** gets MST weight for SD MST spanning all leaves
 *  NOTE: only for debugging! very slow!
 * */
SCIP_Real extreduce_treeGetSdMstWeight(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata             /**< extension data */
)
{
   const int* const leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;

   return sdmstGetWeight(scip, graph, leaves, nleaves, extdata);
}


/** gets MST weight for SD MST spanning all leaves and extension vertex
 *  NOTE: only for debugging! very slow!
 * */
SCIP_Real extreduce_treeGetSdMstExtWeight(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   extvert,            /**< extended vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Real mstweight;
   int* extleaves;
   const int* const leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;

   assert(extvert >= 0 && extvert < graph->knots);

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &extleaves, STP_EXTTREE_MAXNLEAVES_GUARD + 1) );

   assert(nleaves <= STP_EXTTREE_MAXNLEAVES_GUARD);

   for( int i = 0; i < nleaves; ++i )
   {
      const int leaf = leaves[i];
      assert(extvert != leaf);

      extleaves[i] = leaf;
   }

   extleaves[nleaves] = extvert;

   mstweight = sdmstGetWeight(scip, graph, extleaves, nleaves + 1, extdata);

   SCIPfreeBufferArray(scip, &extleaves);

   return mstweight;
}


/** prints the current stack */
void extreduce_printStack(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata             /**< extension data */
)
{
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
}


/** is the node in the current top component of the stack? */
SCIP_Bool extreduce_nodeIsInStackTop(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
   int                   node                /**< the node */
   )
{
   const int stackpos = extStackGetPosition(extdata);
   const int* const stack_start = extdata->extstack_start;
   const int* const stack_data = extdata->extstack_data;

   assert(graph);
   assert(stack_start && stack_data);
   assert(node >= 0 && node < graph->knots);
   assert(extdata->extstack_state[stackpos] == EXT_STATE_EXPANDED || extdata->extstack_state[stackpos] == EXT_STATE_MARKED);

   for( int i = stack_start[stackpos]; i < stack_start[stackpos + 1]; i++ )
   {
      const int topnode = graph->head[stack_data[i]];

      assert(topnode >= 0 && topnode < graph->knots);

      if( topnode == node )
         return TRUE;
   }

   return FALSE;
}


/** Are the close-nodes still valid?
 *  NOTE: expensive method, just designed for debugging! */
SCIP_Bool extreduce_distCloseNodesAreValid(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   const DISTDATA*       distdata            /**< distance data */
)
{
   SCIP_Real* closenodes_dists;
   int* closenodes_indices;
   int nclosenodes;
   const int nnodes = graph_get_nNodes(g);
   SCIP_Bool isValid = TRUE;

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &closenodes_indices, nnodes) );
   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &closenodes_dists, nnodes) );

   assert(scip && distdata);
   assert(distdata->pathroot_isdirty);

   for( int i = 0; i < nnodes; i++ )
   {
      if( distdata->pathroot_isdirty[i] )
         continue;

      SCIP_CALL_ABORT( distCloseNodesCompute(scip, g, distdata, i,
            closenodes_indices, closenodes_dists, &nclosenodes) );

      if( !distCloseNodesIncluded(scip, g, distdata, i, closenodes_indices, closenodes_dists, nclosenodes) )
      {
#ifdef SCIP_DEBUG
         SCIPdebugMessage("corrupted node: ");
         graph_knot_printInfo(g, i);
#endif

         isValid = FALSE;
         break;
      }
   }

   SCIPfreeMemoryArray(scip, &closenodes_dists);
   SCIPfreeMemoryArray(scip, &closenodes_indices);

   return isValid;
}


/** debug initialization */
void extreduce_extendInitDebug(
   int*                  extedgesstart,      /**< array */
   int*                  extedges            /**< array */
)
{
   assert(extedgesstart && extedges);

   for( int i = 0; i < STP_EXT_MAXGRAD; i++ )
      extedgesstart[i] = -1;

   for( int i = 0; i < STP_EXT_MAXGRAD * STP_EXT_MAXGRAD; i++ )
      extedges[i] = -1;
}


/** check whether vertical SDs are up to date for given leaf of component */
SCIP_Bool extreduce_sdsverticalInSync(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   compsize,           /**< size of component */
   int                   nleaves_ancestors,  /**< number of leaves to ancestors */
   int                   topleaf,            /**< component leaf to check for */
   EXTDATA*              extdata             /**< extension data */
   )
{
   const MLDISTS* const sds_vertical = extdata->reddata->sds_vertical;
   const int* const leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   const SCIP_Real* const adjedgecosts = extreduce_mldistsTopTargetDists(sds_vertical, topleaf);
   const int nleaves_old = nleaves - compsize;
   SCIP_Bool isInSync = TRUE;
#ifndef NDEBUG
   const int* const adjids = extreduce_mldistsTopTargetIds(sds_vertical, topleaf);
#endif

   assert(adjedgecosts);
   assert(nleaves_old == nleaves_ancestors);
   assert(nleaves_old > 0 && nleaves_old < nleaves);

#ifndef STP_DEBUG_EXTPC
   if( graph_pc_isPc(graph) )
      return TRUE;
#endif

   /* get the SDs to the ancestor (lower) leafs and compare */
   for( int j = 0; j < nleaves_old; j++ )
   {
      const int leaf = leaves[j];
      const SCIP_Real sd_old = adjedgecosts[j];
      const SCIP_Real specialDist_new = extreduce_extGetSd(scip, graph, topleaf, leaf, extdata);
      const SCIP_Real sd_new = (specialDist_new >= -0.5) ? specialDist_new : FARAWAY;

      assert(!extreduce_nodeIsInStackTop(graph, extdata, leaf));
      assert(extdata->tree_deg[leaf] == 1);
      assert(leaf != topleaf);
      assert(adjids[j] == leaf);

       if( !EQ(sd_old, sd_new) )
       {
          SCIPdebugMessage("vertical SDs are wrong! %f!=%f \n", sd_old, sd_new);

          isInSync = FALSE;
          break;
       }
   }


   return isInSync;
}


/** check whether horizontal SDs are up to date for given leaf of component */
SCIP_Bool extreduce_sdshorizontalInSync(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   topleaf,            /**< component leaf to check for */
   EXTDATA*              extdata             /**< extension data */
   )
{
   const MLDISTS* const sds_horizontal = extdata->reddata->sds_horizontal;
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int* const ghead = graph->head;
   const int stackpos = extStackGetPosition(extdata);
   SCIP_Bool isInSync = TRUE;
#ifndef NDEBUG
   SCIP_Bool hitTopLeaf = FALSE;
#endif

#ifndef STP_DEBUG_EXTPC
   if( graph_pc_isPc(graph) )
      return TRUE;
#endif

   for( int i = extstack_start[stackpos], j = 0; i < extstack_start[stackpos + 1]; i++, j++ )
   {
      const int edge2sibling = extstack_data[i];
      const int sibling = ghead[edge2sibling];

      assert(extreduce_nodeIsInStackTop(graph, extdata, sibling));
      assert(extdata->tree_deg[sibling] == 1);

      if( sibling == topleaf )
      {
#ifndef NDEBUG
         hitTopLeaf = TRUE;
#endif
         continue;
      }
      else
      {
         const SCIP_Real sd_old = extreduce_mldistsTopTargetDist(sds_horizontal, topleaf, sibling);
         const SCIP_Real specialDist_new = extreduce_extGetSdDouble(scip, graph, topleaf, sibling, extdata);
         const SCIP_Real sd_new = (specialDist_new >= -0.5) ? specialDist_new : FARAWAY;

         if( !EQ(sd_old, sd_new) )
         {
            SCIPdebugMessage("vertical SDs are wrong! %f!=%f \n", sd_old, sd_new);

            isInSync = FALSE;
            break;
         }
      }
   }

   assert(hitTopLeaf || !isInSync);

   return isInSync;
}


/** are sds from top component leaf corresponding to current tree? */
SCIP_Bool extreduce_sdsTopInSync(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real       sds[],              /**< SDs from top leaf */
   int                   topleaf,            /**< component leaf to check for */
   EXTDATA*              extdata             /**< extension data */
   )
{
   const int* const leaves = extdata->tree_leaves;
   const int nleaves = extdata->tree_nleaves;
   SCIP_Bool isInSync = TRUE;

#ifndef STP_DEBUG_EXTPC
   if( graph_pc_isPc(graph) )
      return TRUE;
#endif

   for( int j = 0; j < nleaves; j++ )
   {
      const int leaf = leaves[j];

      if( leaf != topleaf )
      {
         const SCIP_Real specialDist_double = extreduce_extGetSdDouble(scip, graph, topleaf, leaf, extdata);
         const SCIP_Real sd_double = (specialDist_double > -0.5)? specialDist_double : FARAWAY;

         const SCIP_Real specialDist = extreduce_extGetSd(scip, graph, topleaf, leaf, extdata);
         const SCIP_Real sd = (specialDist > -0.5)? specialDist : FARAWAY;

         if( !EQ(sds[j], sd_double) && !EQ(sds[j], sd) )
         {
            SCIPdebugMessage("SD from %d to %d not correct! \n", topleaf, leaf);
            SCIPdebugMessage("new sds (double, single): %f, %f ... old sd: %f  \n", sd_double, sd, sds[j]);

            isInSync = FALSE;
            break;
         }
      }
      else
      {
         if( !EQ(sds[j], FARAWAY) )
         {
            SCIPdebugMessage("SD to topleaf not FARAWAY! (but %f) \n", sds[j]);

            isInSync = FALSE;
            break;
         }
      }
   }

   return isInSync;
}


#if 0
/** does the stack top correspond to MST depository top? */
SCIP_Bool extreduce_stackTopMstDepoInSync(
   const GRAPH*          graph,             /**< graph data structure */
   const EXTDATA*        extdata            /**< extension data */
)
{
   const REDDATA* const reddata = extdata->reddata;
   const CSRDEPO* const msts = reddata->msts;
   CSR topmst;
   const int* const extstack_data = extdata->extstack_data;
   const int* const extstack_start = extdata->extstack_start;
   const int stackpos = extStackGetPosition(extdata);
   const int topsize = (extstack_start[stackpos + 1] - extstack_start[stackpos]);

   graph_csrdepo_getTop(msts, &topmst);

   assert(topsize > 0);

   for( int i = extstack_start[stackpos]; i < extstack_start[stackpos + 1]; i++ )
   {
      const int edge = extstack_data[i];
      const int head = graph->head[edge];

      if( reduce )
      {

      }

   }

   assert(extdata->extstack_state[stackpos] != EXT_STATE_NONE);

   return TRUE;
}
#endif


// assert that the costs of the tree all coincide with the actual SD etc distances!
// might be good to have this and reddata extdata stuff in extra method reduce_ext_util.c or just reduce_util.c
// need some flag (in cgraph?) to see whether a leaf in the cgraph does not actually have valid costs (or any)
// and should be recomputed!
