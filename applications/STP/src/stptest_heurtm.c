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

/**@file   stptest_heurtm.c
 * @brief  tests for Steiner tree TM (shortest path based) heuristics
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problem TM (shortest path based) heuristics.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "graph.h"
#include "solstp.h"
#include "heur_tm.h"
#include "portab.h"


/** runs TM heuristic */
static
SCIP_RETCODE runTmPcFull(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   int*                  steinertree_edges   /**< the tree */
)
{
   SCIP_Bool success;

   SCIP_CALL( SCIPStpHeurTMRun(scip, pcmode_fulltree,
      graph, NULL, NULL, steinertree_edges, graph->terms - 1, 2, graph->cost, graph->cost, NULL, NULL, &success) );

   STPTEST_ASSERT(success);

   return SCIP_OKAY;
}


/** tests that PC solution is as expected */
static
SCIP_RETCODE testTmGivesExpectedTreePcFull1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree_edges;
   const int nnodes_org = 5;
   const int nedges_org = 8;
   int nnodes = -1;
   int nedges = -1;
   const SCIP_Real cost_expected = 3.9;
   SCIP_Real cost_tmfull;
   SCIP_Real offset;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0,1
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2,3
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 2.0; /* the pseudo root */
   graph->prize[1] = 1.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 0.9;
   graph->prize[4] = 1.0;

   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, &nnodes, &nedges) );

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_edges, nedges));

   /* actual test */
   SCIP_CALL( runTmPcFull(scip, graph, steinertree_edges) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_tmfull = solstp_getObjBounded(graph, steinertree_edges, offset, nedges);

   STPTEST_ASSERT_MSG_ARGS(EQ(cost_tmfull, cost_expected), "wrong cost: expected=%f, real=%f \n", cost_expected, cost_tmfull);

   stptest_graphTearDown(scip, graph);

   SCIPfreeMemoryArray(scip, &steinertree_edges);

   return SCIP_OKAY;
}


/** tests that PC solution is as expected */
static
SCIP_RETCODE testTmGivesExpectedTreePcFull2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree_edges;
   const int nnodes_org = 5;
   const int nedges_org = 12;
   int nnodes = -1;
   int nedges = -1;
   const SCIP_Real cost_expected = 1.8;
   SCIP_Real cost_tmfull;
   SCIP_Real offset;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM);

   graph->source = 1;

   graph_knot_chg(graph, 0, STP_TERM_NONE);

   graph_edge_addBi(scip, graph, 0, 1, 0.6);
   graph_edge_addBi(scip, graph, 0, 2, 0.6);
   graph_edge_addBi(scip, graph, 0, 3, 0.6);
   graph_edge_addBi(scip, graph, 0, 4, 0.6);
   graph_edge_addBi(scip, graph, 1, 2, 0.3);
   graph_edge_addBi(scip, graph, 3, 4, 0.3);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 0.0;
   graph->prize[1] = 1.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 1.0;
   graph->prize[4] = 1.0;

   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, &nnodes, &nedges) );

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_edges, nedges));

   /* actual test */
   SCIP_CALL( runTmPcFull(scip, graph, steinertree_edges) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_tmfull = solstp_getObjBounded(graph, steinertree_edges, offset, nedges);

   STPTEST_ASSERT_MSG_ARGS(EQ(cost_tmfull, cost_expected), "wrong cost: expected=%f, real=%f \n", cost_expected, cost_tmfull);

   stptest_graphTearDown(scip, graph);

   SCIPfreeMemoryArray(scip, &steinertree_edges);

   return SCIP_OKAY;
}



/** tests that PC solution is as expected */
static
SCIP_RETCODE testTmGivesExpectedTreePcFull3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree_edges;
   const int nnodes_org = 8;
   const int nedges_org = 16;
   int nnodes = -1;
   int nedges = -1;
   const SCIP_Real cost_expected = 6.7;
   SCIP_Real cost_tmfull;
   SCIP_Real offset;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM);

   graph->source = 1;

   graph_knot_chg(graph, 3, STP_TERM_NONE);

   graph_edge_addBi(scip, graph, 0, 1, 1.1);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 6, 7, 1.1);
   graph_edge_addBi(scip, graph, 4, 7, 2.3);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 2.0;
   graph->prize[1] = 2.0;
   graph->prize[2] = 0.1;
   graph->prize[3] = 0.0;
   graph->prize[4] = 1.5;
   graph->prize[5] = 0.1;
   graph->prize[6] = 2.0;
   graph->prize[7] = 2.0;

   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, &nnodes, &nedges) );

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_edges, nedges));

   /* actual test */
   SCIP_CALL( runTmPcFull(scip, graph, steinertree_edges) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_tmfull = solstp_getObjBounded(graph, steinertree_edges, offset, nedges);

   STPTEST_ASSERT_MSG_ARGS(EQ(cost_tmfull, cost_expected), "wrong cost: expected=%f, real=%f \n", cost_expected, cost_tmfull);

   stptest_graphTearDown(scip, graph);

   SCIPfreeMemoryArray(scip, &steinertree_edges);

   return SCIP_OKAY;
}


/** tests that PC solution is improved by pruning */
static
SCIP_RETCODE testPrunedSolIsImprovedPc1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree_edges;
   STP_Bool* steinertree_nodes;
   const int nnodes_org = 5;
   const int nedges_org = 8;
   int nnodes = -1;
   int nedges = -1;
   const SCIP_Real cost_expected = 3.9;
   SCIP_Real cost_pruned;
   SCIP_Real offset;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0,1
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2,3
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 2.0; /* the pseudo root */
   graph->prize[1] = 1.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 0.9;
   graph->prize[4] = 1.0;

   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, &nnodes, &nedges) );

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_edges, nedges));

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = TRUE;

   /* actual test */
   SCIP_CALL( solstp_prune(scip, graph, steinertree_edges, steinertree_nodes) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_pruned = solstp_getObjBounded(graph, steinertree_edges, offset, nedges);

   STPTEST_ASSERT_MSG_ARGS(EQ(cost_pruned, cost_expected), "wrong cost: expected=%f, real=%f \n", cost_expected, cost_pruned);

   stptest_graphTearDown(scip, graph);

   SCIPfreeMemoryArray(scip, &steinertree_edges);
   SCIPfreeMemoryArray(scip, &steinertree_nodes);

   return SCIP_OKAY;
}



/** tests that PC solution is improved by pruning */
static
SCIP_RETCODE testPrunedSolIsImprovedPc2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree_edges;
   STP_Bool* steinertree_nodes;
   const int nnodes_org = 5;
   const int nedges_org = 8;
   int nnodes = -1;
   int nedges = -1;
   const SCIP_Real cost_expected = 4.2;
   SCIP_Real cost_pruned;
   SCIP_Real offset;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 4, 1, 1.0);
   graph_edge_addBi(scip, graph, 4, 2, 1.3);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 0, 1.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 1.1;
   graph->prize[1] = 1.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 1.1;
   graph->prize[4] = 2.0; /* the pseudo-root */

   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, &nnodes, &nedges) );

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_edges, nedges));

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   for( int i = 0; i < nnodes_org; i++ )
      steinertree_nodes[i] = TRUE;

   /* actual test */
   SCIP_CALL( solstp_prune(scip, graph, steinertree_edges, steinertree_nodes) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_pruned = solstp_getObjBounded(graph, steinertree_edges, offset, nedges);

   STPTEST_ASSERT_MSG_ARGS(EQ(cost_pruned, cost_expected), "wrong cost: expected=%f, real=%f \n", cost_expected, cost_pruned);

   stptest_graphTearDown(scip, graph);

   SCIPfreeMemoryArray(scip, &steinertree_edges);
   SCIPfreeMemoryArray(scip, &steinertree_nodes);

   return SCIP_OKAY;
}




/** tests that RMW solution is improved by pruning */
static
SCIP_RETCODE testPrunedSolIsImprovedRmw1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree_edges;
   STP_Bool* steinertree_nodes;
   const int nnodes_org = 5;
   const int nedges_org = 8;
   int nnodes = -1;
   int nedges = -1;
   const SCIP_Real cost_expected = 3.8;
   SCIP_Real cost_pruned;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0,1
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2,3
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = FARAWAY; /* the root */
   graph->prize[1] = 1.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 0.8;
   graph->prize[4] = 1.1;

   SCIP_CALL( stptest_graphSetUpRmwExtended(scip, graph, &nnodes, &nedges) );

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_edges, nedges));

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = TRUE;

   /* actual test */
   SCIP_CALL( solstp_prune(scip, graph, steinertree_edges, steinertree_nodes) );

   cost_pruned = solstp_getObjBounded(graph, steinertree_edges, 0.0, nedges);

   STPTEST_ASSERT_MSG_ARGS(EQ(cost_pruned, cost_expected), "wrong cost: expected=%f, real=%f \n", cost_expected, cost_pruned);

   stptest_graphTearDown(scip, graph);

   SCIPfreeMemoryArray(scip, &steinertree_edges);
   SCIPfreeMemoryArray(scip, &steinertree_nodes);

   return SCIP_OKAY;
}


/** test pruning of solution */
SCIP_RETCODE stptest_testSolPrune(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);
   SCIP_CALL(testPrunedSolIsImprovedRmw1(scip));
   SCIP_CALL(testPrunedSolIsImprovedPc2(scip));
   SCIP_CALL(testPrunedSolIsImprovedPc1(scip));

   return SCIP_OKAY;
}


/** tests TM heuristic */
SCIP_RETCODE stptest_testHeurTm(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL(testTmGivesExpectedTreePcFull3(scip));
   SCIP_CALL(testTmGivesExpectedTreePcFull2(scip) );
   SCIP_CALL(testTmGivesExpectedTreePcFull1(scip) );

   return SCIP_OKAY;
}
