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
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"


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



/** does the cgraph correspond to the current tree? */
SCIP_Bool extreduce_cgraphInSyncWithTree(
   const EXTDATA*        extdata             /**< extension data */
   )
{
   CGRAPH* cgraph;

   assert(extdata);
   assert(extdata->reddata);

   cgraph = extdata->reddata->cgraph;

   assert(cgraph);


// assert that the costs of the tree all coincide with the actual SD etc distances!
// might be good to have this and reddata extdata stuff in extra method reduce_ext_util.c or just reduce_util.c
// need some flag (in cgraph?) to see whether a leaf in the cgraph does not actually have valid costs (or any)
// and should be recomputed!

   if( !cgraph_idsInSync(cgraph, extdata->tree_leaves, extdata->tree_nleaves) )
      return FALSE;

   return TRUE;
}
