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
#include "heur_tm.h"
#include "portab.h"



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
   SCIP_CALL( graph_solPrune(scip, graph, steinertree_edges, steinertree_nodes) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_pruned = graph_solGetObj(graph, steinertree_edges, offset, nedges);

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
   SCIP_CALL( graph_solPrune(scip, graph, steinertree_edges, steinertree_nodes) );

   offset = graph_pc_getNonLeafTermOffset(scip, graph);
   cost_pruned = graph_solGetObj(graph, steinertree_edges, offset, nedges);

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

   SCIP_CALL(testPrunedSolIsImprovedPc2(scip));
   SCIP_CALL(testPrunedSolIsImprovedPc1(scip));

   return SCIP_OKAY;
}
