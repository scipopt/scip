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

/**@file   stptest_extutils.c
 * @brief  tests for Steiner tree extended reduction utilites
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problem reductions.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "portab.h"
#include "graph.h"
#include "reduce.h"
#include "extreduce.h"


/** tests building and un-building */
static
SCIP_RETCODE testMldistsBuilding(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   const int maxnlevels = 3;
   const int maxnslots = 12;
   const int maxntargets = 7;
#ifndef NDEBUG
   const SCIP_Real* dists;
#endif

   MLDISTS* mldists;

   SCIP_CALL( extreduce_mldistsInit(scip, maxnlevels, maxnslots, maxntargets, 0, TRUE, &mldists) );

   /* add */

   extreduce_mldistsLevelAddTop(2, 7, mldists);

#ifndef NDEBUG
   dists = extreduce_mldistsEmptySlotTargetDists(mldists);
   assert(dists);
#endif

   for( int i = 0; i < 2; i++ )
   {
      extreduce_mldistsEmptySlotSetBase(i + 10, mldists);
      extreduce_mldistsEmptySlotSetFilled(mldists);
   }

   extreduce_mldistsLevelAddTop(12, 5, mldists);

   for( int i = 0; i < 12; i++ )
   {
      extreduce_mldistsEmptySlotSetBase(i, mldists);
      extreduce_mldistsEmptySlotSetFilled(mldists);
   }

   extreduce_mldistsLevelAddTop(2, 7, mldists);

   for( int i = 0; i < 2; i++ )
   {
      extreduce_mldistsEmptySlotSetBase(i, mldists);
      extreduce_mldistsEmptySlotSetFilled(mldists);
   }

   /* remove */

   extreduce_mldistsLevelRemoveTop(mldists);
   extreduce_mldistsLevelRemoveTop(mldists);
   extreduce_mldistsLevelRemoveTop(mldists);

   STPTEST_ASSERT_MSG(extreduce_mldistsIsEmpty(mldists), "MLDISTS not empty \n");

   extreduce_mldistsFree(scip, &mldists);

   return SCIP_OKAY;
}



/** test close nodes are computed correctly */
static
SCIP_RETCODE testDistCloseNodesPcAreValid1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DISTDATA* distdata;
   GRAPH* graph;
   const int nnodes = 4;
   const int nedges = 8;
   const int nclosenodes = 3;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 2.1);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   graph->prize[0] = 5.0;
   graph->prize[1] = 0.6;
   graph->prize[2] = 0.5;
   graph->prize[3] = 0.5;

   SCIP_CALL( stptest_graphSetUpPcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );

   /* 2. do the actual test */

   /* test close nodes */
   {
      const RANGE* node_range = distdata->closenodes_range;
      const int* node_idx = distdata->closenodes_indices;
      const SCIP_Real* node_dist = distdata->closenodes_distances;

      STPTEST_ASSERT(node_idx[node_range[0].start] == 1);
      STPTEST_ASSERT(node_idx[node_range[0].start + 1] == 2);
      STPTEST_ASSERT(node_idx[node_range[0].start + 2] == 3);

      STPTEST_ASSERT(EQ(node_dist[node_range[0].start], 1.0));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 1], 1.0));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 2], 1.5));
   }

   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);
   stptest_extreduceTearDown(scip, graph, NULL);

   return SCIP_OKAY;
}


/** test close nodes are computed correctly */
static
SCIP_RETCODE testDistCloseNodesPcAreValid2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DISTDATA* distdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 10;
   const int nclosenodes = 4;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.5);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 2.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 2, 0.5);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   graph->prize[0] = 5.0;
   graph->prize[1] = 2.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 1.0;
   graph->prize[4] = 0.0;

   SCIP_CALL( stptest_graphSetUpPcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );

   /* 2. do the actual test */

   /* test close nodes */
   {
      const RANGE* node_range = distdata->closenodes_range;
      const SCIP_Real* node_dist = distdata->closenodes_distances;

      STPTEST_ASSERT(EQ(node_dist[node_range[0].start], 1.5));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 1], 1.5));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 2], 1.5));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 3], 1.5));
   }

   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);
   stptest_extreduceTearDown(scip, graph, NULL);

   return SCIP_OKAY;
}


/** test close nodes are computed correctly */
static
SCIP_RETCODE testDistCloseNodesPcAreValidAfterDeletion(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DISTDATA* distdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 12;
   const int nclosenodes = 4;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.5);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 0.2);
   graph_edge_addBi(scip, graph, 3, 4, 0.7);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   graph->prize[0] = 0.0;
   graph->prize[1] = 1.2;
   graph->prize[2] = 0.4;
   graph->prize[3] = 1.0;
   graph->prize[4] = 2.0;

   SCIP_CALL( stptest_graphSetUpPcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );

   /* 2. do the actual test */

   /* test close nodes */
   {
      const RANGE* node_range = distdata->closenodes_range;
      const int* node_idx = distdata->closenodes_indices;
      const SCIP_Real* node_dist = distdata->closenodes_distances;

      STPTEST_ASSERT(EQ(node_dist[node_range[0].start], 1.0));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 1], 1.0));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 2], 1.5));
      STPTEST_ASSERT(EQ(node_dist[node_range[0].start + 3], 1.0));

      STPTEST_ASSERT(node_idx[node_range[0].start + 3] == 4);

      extreduce_edgeRemove(scip, 2, graph, distdata, NULL);
      extreduce_edgeRemove(scip, 4, graph, distdata, NULL);
      extreduce_edgeRemove(scip, 8, graph, distdata, NULL);

      STPTEST_ASSERT(extreduce_distCloseNodesAreValid(scip, graph, distdata));
   }

   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);
   stptest_extreduceTearDown(scip, graph, NULL);

   return SCIP_OKAY;
}


/** test close nodes are computed correctly */
static
SCIP_RETCODE testDistCloseNodesAreValid(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;
   const int nclosenodes = 2;         /**< max. number of close nodes to each node */

   DISTDATA* distdata;

   assert(scip);

   /* 1. build a test graph */

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = root;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.4);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 0.5);
   graph_edge_addBi(scip, graph, 3, 6, 1.6);
   graph_edge_addBi(scip, graph, 3, 7, 1.6);
   graph_edge_addBi(scip, graph, 4, 8, 1.5);
   graph_edge_addBi(scip, graph, 4, 9, 1.5);
   graph_edge_addBi(scip, graph, 5, 10, 1.0);
   graph_edge_addBi(scip, graph, 10, 11, 1.0);
   graph_edge_addBi(scip, graph, 11, 12, 1.0);

   graph_mark(graph);
   SCIP_CALL( graph_initHistory(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );

   /* 2. do the actual test */

   /* test close nodes */
   {
      const RANGE* node_range = distdata->closenodes_range;
      const int* node_idx = distdata->closenodes_indices;
#if 0
      const SCIP_Real* node_dist = distdata->closenodes_distances;

      for( int node = 0; node < 6; node++ )
      {
         printf("node %d: \n", node);

         for( int i = node_range[node].start; i < node_range[node].end; i++ )
         {
            printf("...closenode=%d  dist=%f\n", node_idx[i], node_dist[i]);
         }
      }
#endif

      STPTEST_ASSERT(node_idx[node_range[1].start] == 2);
      STPTEST_ASSERT(node_idx[node_range[1].start + 1] == 5);
      STPTEST_ASSERT(node_idx[node_range[5].start] == 1);
      STPTEST_ASSERT(node_idx[node_range[5].start + 1] == 2);
   }

   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}



/** test root paths are computed correctly */
static
SCIP_RETCODE testDistRootPathsAreValid(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;
   const int nclosenodes = 2;         /**< max. number of close nodes to each node */

   DISTDATA* distdata;

   assert(scip);

   /* 1. build a test graph */

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = root;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.4);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 0.5);
   graph_edge_addBi(scip, graph, 3, 6, 1.6);
   graph_edge_addBi(scip, graph, 3, 7, 1.6);
   graph_edge_addBi(scip, graph, 4, 8, 1.5);
   graph_edge_addBi(scip, graph, 4, 9, 1.5);
   graph_edge_addBi(scip, graph, 5, 10, 1.0);
   graph_edge_addBi(scip, graph, 10, 11, 1.0);
   graph_edge_addBi(scip, graph, 11, 12, 1.0);

   graph_mark(graph);
   SCIP_CALL( graph_initHistory(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );

   /* 2. do the actual test */

   /* test edge root paths */
   {
      const int edge = 2;
      PRSTATE** pathroot_blocks = distdata->pathroot_blocks;
      int* pathroot_blocksizes = distdata->pathroot_blocksizes;


#ifdef SCIP_DEBUG
      printf("edge ");
      graph_edge_printInfo(graph, edge);
#endif

      STPTEST_ASSERT(graph->head[edge] == 2 && graph->tail[edge] == 1);
      STPTEST_ASSERT(pathroot_blocksizes[edge / 2] == 6);

      for( int i = 0; i < pathroot_blocksizes[edge / 2]; i++ )
         SCIPdebugMessage("...root=%d  \n", pathroot_blocks[edge / 2][i].pathroot_id);
   }

   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


/** test that distances are computed correctly */
static
SCIP_RETCODE testDistDistancesAreValid(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;
   const int nclosenodes = 2;         /**< max. number of close nodes to each node */

   DISTDATA* distdata;

   assert(scip);

   /* 1. build a test graph */

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);

   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = root;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.4);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 0.5);
   graph_edge_addBi(scip, graph, 3, 6, 1.6);
   graph_edge_addBi(scip, graph, 3, 7, 1.6);
   graph_edge_addBi(scip, graph, 4, 8, 1.5);
   graph_edge_addBi(scip, graph, 4, 9, 1.5);
   graph_edge_addBi(scip, graph, 5, 10, 1.0);
   graph_edge_addBi(scip, graph, 10, 11, 1.0);
   graph_edge_addBi(scip, graph, 11, 12, 1.0);

   graph_mark(graph);
   SCIP_CALL( graph_initHistory(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );

   /* test distances */
   {
      const int edge = 2;
      const SCIP_Real dist1_2 = extreduce_distDataGetSd(scip, graph, 1, 2, distdata);
      const SCIP_Real dist1_3 = extreduce_distDataGetSd(scip, graph, 1, 3, distdata);
      const SCIP_Real dist1_5 = extreduce_distDataGetSd(scip, graph, 1, 5, distdata);
      const SCIP_Real dist2_1 = extreduce_distDataGetSd(scip, graph, 2, 1, distdata);
      const SCIP_Real dist2_5 = extreduce_distDataGetSd(scip, graph, 2, 5, distdata);
      const SCIP_Real dist2_6 = extreduce_distDataGetSd(scip, graph, 2, 6, distdata);

      STPTEST_ASSERT(dist1_2 == 0.4);
      STPTEST_ASSERT(dist1_3 == -1.0);
      STPTEST_ASSERT(dist1_5 == 0.9);
      STPTEST_ASSERT(dist2_1 == 0.4);
      STPTEST_ASSERT(dist2_5 == 0.5);
      STPTEST_ASSERT(dist2_6 == -1.0);

      assert(graph->head[edge] == 2 && graph->tail[edge] == 1);
      graph_edge_delFull(scip, graph, edge, TRUE);
      extreduce_distDataDeleteEdge(scip, graph, edge, distdata);

      {
         const SCIP_Real dist1_2_b = extreduce_distDataGetSd(scip, graph, 1, 2, distdata);
         const SCIP_Real dist2_6_b = extreduce_distDataGetSd(scip, graph, 2, 6, distdata);
         const SCIP_Real dist2_5_b = extreduce_distDataGetSd(scip, graph, 2, 5, distdata);

         STPTEST_ASSERT(dist1_2_b == -1.0);
         STPTEST_ASSERT(dist2_6_b == -1.0);
         STPTEST_ASSERT(dist2_5_b == 0.5);
      }
   }

   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


/** tests for utilits */
SCIP_RETCODE stptest_extmldists(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testMldistsBuilding(scip) );
   SCIP_CALL( testDistCloseNodesPcAreValid1(scip) );
   SCIP_CALL( testDistCloseNodesPcAreValid2(scip) );
   SCIP_CALL( testDistCloseNodesPcAreValidAfterDeletion(scip) );
   SCIP_CALL( testDistCloseNodesAreValid(scip) );
   SCIP_CALL( testDistRootPathsAreValid(scip) );
   SCIP_CALL( testDistDistancesAreValid(scip) );

   printf("extended utilities test: all ok \n");

   return SCIP_OKAY;
}
