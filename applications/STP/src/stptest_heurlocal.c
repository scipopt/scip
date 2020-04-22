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

/**@file   stptest_heurlocal.c
 * @brief  tests for Steiner tree heuristics
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problem local search heuristics.
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
#include "heur_local.h"
#include "heur_tm.h"


/** test key path exchange */
static
SCIP_RETCODE localKeyPathExchange(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   int nnodes = 6;
   int nedges = 6;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);
   graph_edge_add(scip, graph, 4, 5, 1.0, 1.0);
   graph_edge_add(scip, graph, 5, 0, 2.0, 2.0); // 10,11

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   steinertree[0] = CONNECT;
   steinertree[2] = CONNECT;
   steinertree[4] = CONNECT;
   steinertree[6] = CONNECT;

   assert(solstp_isValid(scip, graph, steinertree));

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   assert(steinertree[11] == CONNECT && steinertree[9] == CONNECT);
   assert(steinertree[0] == UNKNOWN && steinertree[2] == UNKNOWN);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test key path exchange */
static
SCIP_RETCODE localKeyPathExchangePc(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 6;
   int nedges = 6;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);
   graph_edge_add(scip, graph, 4, 5, 2.0, 2.0);
   graph_edge_add(scip, graph, 5, 0, 3.0, 3.0); // 10,11

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 5.0;
   graph->prize[2] = 3.0;
   graph->prize[3] = 3.1;
   graph->prize[4] = 3.1;
   graph->prize[5] = 3;

   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 5, 0);

   SCIP_CALL(graph_transPc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   // 4 is the implicit root

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );

   assert(solstp_isValid(scip, graph, steinertree));

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   // 5 is the implicit root

   assert(steinertree[9] == CONNECT && steinertree[10] == CONNECT);
   assert(steinertree[0] == UNKNOWN && steinertree[2] == UNKNOWN);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test key path exchange */
static
SCIP_RETCODE localKeyPathExchangePc2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 7;
   int nedges = 7;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);
   graph_edge_add(scip, graph, 4, 5, 2.0, 2.0);
   graph_edge_add(scip, graph, 5, 6, 1.0, 1.0); // 10,11
   graph_edge_add(scip, graph, 6, 0, 2.0, 2.0); // 12,13

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 5.0;
   graph->prize[2] = 3.0;
   graph->prize[3] = 3.1;
   graph->prize[4] = 3.1;
   graph->prize[5] = 1.0;
   graph->prize[6] = 1.0;


   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 5, 0);
   graph_knot_chg(graph, 6, 0);

   SCIP_CALL(graph_transPc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   // 4 is the implicit root

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );

   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   // 5 is the implicit root

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   if( !SCIPisEQ(scip, cost1 + 1.0, cost0) )
      return SCIP_ERROR;

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test key path exchange */
static
SCIP_RETCODE localKeyPathExchangeMw(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 8;
   int nedges = 8;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 0.0, 0.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 0.0, 0.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 0.0, 0.0);
   graph_edge_add(scip, graph, 2, 4, 0.0, 0.0);
   graph_edge_add(scip, graph, 4, 5, 0.0, 0.0);
   graph_edge_add(scip, graph, 5, 6, 0.0, 0.0); // 10,11
   graph_edge_add(scip, graph, 6, 7, 0.0, 0.0);
   graph_edge_add(scip, graph, 7, 0, 0.0, 0.0); // 14,15

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);


   graph->prize[0] = 5.0;
   graph->prize[1] = -2.0;
   graph->prize[2] = -2.0;
   graph->prize[3] = 3.0;
   graph->prize[4] = 3.0;
   graph->prize[5] = -1.0;
   graph->prize[6] = 0.5;
   graph->prize[7] = -1.0;

   SCIP_CALL(graph_transMw(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   // 4 is the implicit root

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );

   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   // 5 is the implicit root

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   if( !SCIPisEQ(scip, cost1 + 0.5, cost0) )
      return SCIP_ERROR;

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test key vertex elimination */
static
SCIP_RETCODE localKeyVertex(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 7;
   int nedges = 8;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);
   graph_edge_add(scip, graph, 4, 5, 2.0, 2.0);
   graph_edge_add(scip, graph, 5, 0, 3.0, 3.0); // 10,11
   graph_edge_add(scip, graph, 6, 3, 1.0, 1.0); // 12,13
   graph_edge_add(scip, graph, 6, 4, 1.0, 1.0); // 14,15

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   // 4 is the implicit root

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );
   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   assert(steinertree[12] != CONNECT && steinertree[15] != CONNECT);
   assert(steinertree[0] == CONNECT && steinertree[2] == CONNECT);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   // 5 is the implicit root

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);
   if( !SCIPisEQ(scip, cost1 + 1.0, cost0) )
      return SCIP_ERROR;

   assert(steinertree[12] == CONNECT && steinertree[15] == CONNECT);
   assert(steinertree[0] == UNKNOWN && steinertree[2] == UNKNOWN);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test key vertex elimination */
static
SCIP_RETCODE localKeyVertexPc(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 7;
   int nedges = 8;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);
   graph_edge_add(scip, graph, 4, 5, 3.0, 3.0);
   graph_edge_add(scip, graph, 5, 0, 3.0, 3.0); // 10,11
   graph_edge_add(scip, graph, 6, 3, 2.0, 2.0); // 12,13
   graph_edge_add(scip, graph, 6, 4, 2.0, 2.0); // 14,15

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 5.0;
 //  graph->prize[2] = 3.0; // todo
   graph->prize[3] = 3.0;
   graph->prize[4] = 3.0;
   graph->prize[5] = 1.5;
   graph->prize[6] = 1.0;


   graph_knot_chg(graph, 0, 0);
  // graph_knot_chg(graph, 2, 0); todo activate again as non-leaf later
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 5, 0);
   graph_knot_chg(graph, 6, 0);

   SCIP_CALL(graph_transPc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   // 4 is the implicit root

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );
   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   assert(steinertree[12] != CONNECT && steinertree[15] != CONNECT);
   assert(steinertree[1] == CONNECT && steinertree[3] == CONNECT);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );


   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);
   if( !SCIPisEQ(scip, cost1 + 0.5, cost0) )
      return SCIP_ERROR;

   assert(steinertree[0] == UNKNOWN && steinertree[2] == UNKNOWN);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test key vertex elimination */
static
SCIP_RETCODE localKeyVertexPc2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 9;
   int nedges = 10;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   /* the solution */
   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);

   /* the right part */
   graph_edge_add(scip, graph, 4, 5, 2.0, 2.0);
   graph_edge_add(scip, graph, 5, 6, 1.0, 1.0); // 10,11
   graph_edge_add(scip, graph, 6, 7, 1.0, 1.0); // 12,13
   graph_edge_add(scip, graph, 7, 0, 2.0, 2.0); // 14,15

   /* the hat */
   graph_edge_add(scip, graph, 8, 3, 2.0, 2.0); // 16,17
   graph_edge_add(scip, graph, 8, 4, 2.0, 2.0); // 18,19

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 5.0;
 //  graph->prize[2] = 3.0; // todo
   graph->prize[3] = 3.0;
   graph->prize[4] = 3.0;
   graph->prize[6] = 1.5;
   graph->prize[8] = 1.0;


   graph_knot_chg(graph, 0, 0);
  // graph_knot_chg(graph, 2, 0); todo activate again as non-leaf later
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 6, 0);
   graph_knot_chg(graph, 8, 0);


   SCIP_CALL(graph_transPc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   // 4 is the implicit root

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );
   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   assert(steinertree[16] == UNKNOWN && steinertree[18] == UNKNOWN);
   assert(steinertree[1] == CONNECT && steinertree[3] == CONNECT);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   // 6 is the implicit root

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);
   if( !SCIPisEQ(scip, cost1 + 0.5, cost0) )
      return SCIP_ERROR;

   assert(steinertree[0] == UNKNOWN && steinertree[2] == UNKNOWN);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** extension test */
static
SCIP_RETCODE localExtendPc(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 9;
   int nedges = 9;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   /* the solution */
   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);

   /* the right part */
   graph_edge_add(scip, graph, 4, 5, 1.0, 1.0);
   graph_edge_add(scip, graph, 5, 6, 1.0, 1.0); // 10,11
   graph_edge_add(scip, graph, 6, 7, 1.0, 1.0); // 12,13

   /* dummy part */
   graph_edge_add(scip, graph, 7, 8, 5.0, 5.0); // 16,17
   graph_edge_add(scip, graph, 8, 1, 5.0, 5.0); // 18,19

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 5.0;
 //  graph->prize[2] = 3.0; // todo
   graph->prize[3] = 3.0;
   graph->prize[4] = 3.0;
   graph->prize[6] = 1.0;
   graph->prize[7] = 2.5;

   graph->prize[8] = 0.5;


   graph_knot_chg(graph, 0, 0);
  // graph_knot_chg(graph, 2, 0); todo activate again as non-leaf later
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 6, 0);
   graph_knot_chg(graph, 7, 0);
   graph_knot_chg(graph, 8, 0);


   SCIP_CALL(graph_transPc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );
   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalExtendPcMw(scip, graph, graph->cost, steinertree, steinertree_nodes) );

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);
   if( !SCIPisEQ(scip, cost1 + 0.5, cost0) )
      return SCIP_ERROR;

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree_nodes);
   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}

/** test insertion */
static
SCIP_RETCODE localInsertion(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   int nnodes = 6;
   int nedges = 7;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   graph->stp_type = STP_SPG;

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);

   graph_edge_add(scip, graph, 5, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 5, 3, 1.0, 1.0); // 10,11
   graph_edge_add(scip, graph, 5, 4, 1.0, 1.0); // 12,13

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   steinertree[0] = CONNECT;
   steinertree[2] = CONNECT;
   steinertree[4] = CONNECT;
   steinertree[6] = CONNECT;

   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);


   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   if( !SCIPisEQ(scip, cost1 + 1.0, cost0) )
      return SCIP_ERROR;

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test insertion */
static
SCIP_RETCODE localInsertion2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   int nnodes = 6;
   int nedges = 7;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   graph->stp_type = STP_SPG;

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);
   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);

   graph_edge_add(scip, graph, 5, 0, 1.5, 1.5);
   graph_edge_add(scip, graph, 5, 3, 1.5, 1.5); // 10,11
   graph_edge_add(scip, graph, 5, 4, 1.5, 1.5); // 12,13

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   steinertree[0] = CONNECT;
   steinertree[2] = CONNECT;
   steinertree[4] = CONNECT;
   steinertree[6] = CONNECT;

   assert(solstp_isValid(scip, graph, steinertree));

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);


   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   if( !SCIPisEQ(scip, cost1 + 1.5, cost0) )
   {
      printf("localInsertion2 unit test failed \n");
      printf("cost1=%f, cost0=%f \n", cost1, cost0);

      return SCIP_ERROR;
   }

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree);

   return SCIP_OKAY;
}


/** test insertion */
static
SCIP_RETCODE localInsertion2pc(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* steinertree;
   STP_Bool* steinertree_nodes;
   int nnodes = 6;
   int nedges = 7;
   SCIP_Real cost0;
   SCIP_Real cost1;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   graph->stp_type = STP_SPG;

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);


   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0); // 0,1
   graph_edge_add(scip, graph, 1, 2, 2.0, 2.0); // 2,3
   graph_edge_add(scip, graph, 2, 3, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 4, 2.0, 2.0);

   graph_edge_add(scip, graph, 5, 0, 1.5, 1.5);
   graph_edge_add(scip, graph, 5, 3, 1.5, 1.5); // 10,11
   graph_edge_add(scip, graph, 5, 4, 1.5, 1.5); // 12,13

   nnodes = graph->knots;
   nedges = graph->edges;

   graph_pc_initPrizes(scip, graph, nnodes);

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 5.0;
   graph->prize[1] = 0.75;
   graph->prize[2] = 2.1;
   graph->prize[3] = 4.0;
   graph->prize[4] = 4.0;

   graph_knot_chg(graph, 0, 0);
   graph_knot_chg(graph, 1, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 4, 0);

   SCIP_CALL(graph_transPc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree_nodes, nnodes));
   SCIP_CALL(SCIPallocMemoryArray(scip, &steinertree, nedges));

   for( int i = 0; i < nedges; i++ )
      steinertree[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      steinertree_nodes[i] = FALSE;

   steinertree_nodes[0] = TRUE;
   steinertree_nodes[1] = TRUE;
   steinertree_nodes[2] = TRUE;
   steinertree_nodes[3] = TRUE;
   steinertree_nodes[4] = TRUE;

   SCIP_CALL( solstp_prune(scip, graph, steinertree, steinertree_nodes) );

   cost0 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalRun(scip, graph, steinertree) );

   cost1 = solstp_getObjBounded(graph, steinertree, 0.0, nedges);

   if( !SCIPisEQ(scip, cost1 + 0.75, cost0) )
   {
      printf("localInsertion2pc unit test failed \n");
      printf("obj: new=%f, old=%f \n", cost1, cost0);

      return SCIP_ERROR;
   }

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeMemoryArray(scip, &steinertree);
   SCIPfreeMemoryArray(scip, &steinertree_nodes);

   return SCIP_OKAY;
}


/** test local search heuristics */
SCIP_RETCODE stptest_testHeurLocal(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( localInsertion(scip) );
   SCIP_CALL( localInsertion2(scip) );
   SCIP_CALL( localInsertion2pc(scip) );
   SCIP_CALL( localKeyPathExchangeMw(scip) );
   SCIP_CALL( localKeyVertexPc2(scip) );
   SCIP_CALL( localKeyVertex(scip) );
   SCIP_CALL( localKeyVertexPc(scip) );
   SCIP_CALL( localKeyPathExchangePc2(scip) );
   SCIP_CALL( localKeyPathExchangePc(scip) );
   SCIP_CALL( localKeyPathExchange(scip) );

   SCIP_CALL( localExtendPc(scip) );

   printf("stptest_heur_local: all ok \n");

   return SCIP_OKAY;
}
