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

/**@file   stptest_pcreduce.c
 * @brief  tests for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements unit tests for Steiner problem PC/MW reductions.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "graph.h"
#include "reduce.h"
#include "heur_local.h"
#include "heur_tm.h"
#include "portab.h"


static
SCIP_RETCODE checkSdWalk(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             extended,
   GRAPH**               g,
   int*                  nelims
)
{
   DHEAP* dheap;
   int nnodes;
   int nedges;
   int* nodearr_int1;
   int* nodearr_int2;
   int* vbase;
   int* state;
   int* heap;
   SCIP_Real* edgearrreal1;

   SCIP_Real* nodearrreal1;
   SCIP_Bool* isterm;
   STP_Bool* nodearrchar;
   GRAPH* graph = *g;

   /* build PC graph */
   SCIP_CALL( graph_transPc(scip, graph) );

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   graph_mark(graph);

   graph_pc_2org(scip, graph);

   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal1, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal1, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int1, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &isterm, nnodes) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, 4 * nnodes) );

   graph_heap_create(scip, nnodes, NULL, NULL, &dheap);

   graph_get_isTerm(graph, isterm);

   /* actual test */

   if( extended )
   {
      SCIP_CALL( reduce_sdWalkExt(scip, 200, NULL, graph, nodearrreal1, heap, state, vbase, nodearrchar, nelims) );
   }
   else
   {
     SCIP_CALL( reduce_sdWalkTriangle(scip, 5, NULL, graph, nodearr_int1, nodearrreal1, heap, nodearrchar, dheap, nelims));
   }

   /* clean up */

   graph_heap_free(scip, TRUE, TRUE, &dheap);
   graph_path_exit(scip, graph);
   graph_free(scip, g, TRUE);

   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &vbase);

   SCIPfreeBufferArray(scip, &isterm);
   SCIPfreeBufferArray(scip, &nodearr_int2);
   SCIPfreeBufferArray(scip, &nodearr_int1);
   SCIPfreeBufferArray(scip, &nodearrreal1);
   SCIPfreeBufferArray(scip, &edgearrreal1);
   SCIPfreeBufferArray(scip, &heap);

   return SCIP_OKAY;
}


/** tests ANS RMW test */
static
SCIP_RETCODE testRmwAnsDeletesOneNode(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes_org = 5;
   const int nedges_org = 12;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 3, STP_TERM);
   graph->source = 3;

   graph_knot_chg(graph, 4, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 0.0); // 0,1
   graph_edge_addBi(scip, graph, 0, 2, 0.0); // 2,3
   graph_edge_addBi(scip, graph, 0, 3, 0.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.0);
   graph_edge_addBi(scip, graph, 1, 3, 0.0);
   graph_edge_addBi(scip, graph, 2, 4, 0.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = -1.0;
   graph->prize[1] = -1.0;
   graph->prize[2] = -0.5;
   graph->prize[3] = FARAWAY;
   graph->prize[4] = 2.0;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_ans(scip, graph, &nelims) );

   STPTEST_ASSERT(graph->grad[1] == 0);
   STPTEST_ASSERT(graph->grad[0] != 0);
   STPTEST_ASSERT(graph->grad[2] != 0);
   STPTEST_ASSERT(graph->grad[3] != 0);
   STPTEST_ASSERT(graph->grad[4] != 0);


   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** test ANS RMW test */
static
SCIP_RETCODE testRmwAnsDeletesOneEdge(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes_org = 7;
   const int nedges_org = 20;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 4, STP_TERM);
   graph->source = 4;

   graph_knot_chg(graph, 6, STP_TERM);

   graph_edge_addBi(scip, graph, 1, 5, 0.0); // edge to be deleted
   graph_edge_addBi(scip, graph, 0, 1, 0.0);
   graph_edge_addBi(scip, graph, 0, 2, 0.0);
   graph_edge_addBi(scip, graph, 0, 3, 0.0);
   graph_edge_addBi(scip, graph, 0, 6, 0.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.0);
   graph_edge_addBi(scip, graph, 1, 3, 0.0);
   graph_edge_addBi(scip, graph, 2, 4, 0.0);
   graph_edge_addBi(scip, graph, 3, 4, 0.0);
   graph_edge_addBi(scip, graph, 3, 5, 0.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = -1.0;
   graph->prize[1] = -0.6;
   graph->prize[2] = -0.5;
   graph->prize[3] = -2.0;
   graph->prize[4] = FARAWAY;
   graph->prize[5] = -0.4;
   graph->prize[6] = 0.4;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_ans(scip, graph, &nelims) );

   STPTEST_ASSERT(graph->oeat[0] == EAT_FREE);

   for( int i = 0; i < nnodes_org; i++ )
      STPTEST_ASSERT(graph->grad[i] != 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** test ANS RMW test */
static
SCIP_RETCODE testRmwAnsDeletesTwoNodes(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes_org = 7;
   const int nedges_org = 20;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 4, STP_TERM);
   graph->source = 4;

   graph_knot_chg(graph, 6, STP_TERM);

   graph_edge_addBi(scip, graph, 1, 5, 0.0); // edge to be deleted
   graph_edge_addBi(scip, graph, 0, 1, 0.0);
   graph_edge_addBi(scip, graph, 0, 2, 0.0);
   graph_edge_addBi(scip, graph, 0, 3, 0.0);
   graph_edge_addBi(scip, graph, 0, 6, 0.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.0);
   graph_edge_addBi(scip, graph, 1, 3, 0.0);
   graph_edge_addBi(scip, graph, 2, 4, 0.0);
   graph_edge_addBi(scip, graph, 3, 4, 0.0);
   graph_edge_addBi(scip, graph, 3, 5, 0.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = -1.0;
   graph->prize[1] = -1.1;
   graph->prize[2] = -0.5;
   graph->prize[3] = -2.0;
   graph->prize[4] = FARAWAY;
   graph->prize[5] = -1.0;
   graph->prize[6] = 0.4;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_ans(scip, graph, &nelims) );

   STPTEST_ASSERT(graph->grad[1] == 0);
   STPTEST_ASSERT(graph->grad[5] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** test simple RMW test */
static
SCIP_RETCODE testRmwTerminalContraction(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Real offset = 0.0;
   const int nnodes_org = 4;
   const int nedges_org = 10;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 3, STP_TERM);
   graph->source = 0;

   graph_edge_addBi(scip, graph, 3, 1, 0.0); // 0,1
   graph_edge_addBi(scip, graph, 3, 2, 0.0); // 2,3
   graph_edge_addBi(scip, graph, 0, 3, 0.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.0);
   graph_edge_addBi(scip, graph, 2, 3, 0.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = FARAWAY;
   graph->prize[1] = -1.0;
   graph->prize[2] = -5.0;
   graph->prize[3] = 2.0;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_simple_mw(scip, graph, NULL, &offset, &nelims) );

   STPTEST_ASSERT(graph->grad[3] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** test simple RMW test */
static
SCIP_RETCODE testRmwTerminalDeg1Contraction1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Real offset = 0.0;
   const int nnodes_org = 3;
   const int nedges_org = 4;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 0.0);
   graph_edge_addBi(scip, graph, 0, 2, 0.0);
   graph->source = 2;

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = -1.0;
   graph->prize[1] = 1.5;
   graph->prize[2] = FARAWAY;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_simple_mw(scip, graph, NULL, &offset, &nelims) );

   STPTEST_ASSERT(graph->grad[graph->source] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}

/** test simple RMW test */
static
SCIP_RETCODE testRmwTerminalDeg1Contraction2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Real offset = 0.0;
   const int nnodes_org = 3;
   const int nedges_org = 4;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 0.0);
   graph_edge_addBi(scip, graph, 0, 2, 0.0);
   graph->source = 2;

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = -1.0;
   graph->prize[1] = 0.5;
   graph->prize[2] = FARAWAY;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_simple_mw(scip, graph, NULL, &offset, &nelims) );

   STPTEST_ASSERT(EQ(offset, 0.5));
   STPTEST_ASSERT(graph->grad[graph->source] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** test simple RMW test */
static
SCIP_RETCODE testRmwTerminalDeg1Contraction3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Real offset = 0.0;
   const int nnodes_org = 3;
   const int nedges_org = 4;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 0.0);
   graph_edge_addBi(scip, graph, 0, 2, 0.0);
   graph->source = 2;

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = -1.0;
   graph->prize[1] = FARAWAY;
   graph->prize[2] = FARAWAY;

   SCIP_CALL( stptest_graphSetUpRmwOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( reduce_simple_mw(scip, graph, NULL, &offset, &nelims) );

   STPTEST_ASSERT(EQ(offset, 1.0));
   STPTEST_ASSERT(graph->grad[graph->source] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


static
SCIP_RETCODE testSdPcKillsEdge1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nelims;
   const int nnodes = 3;
   const int nedges = 3;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, 2 * nedges, 1) );

   /* build graph */
   graph_knot_add(graph, 0);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 1.0;

   nelims = 0;

   SCIP_CALL( checkSdWalk(scip, FALSE, &graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   assert(graph == NULL);

   return SCIP_OKAY;
}

static
SCIP_RETCODE testSdPcKillsEdge2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nelims;
   const int nnodes = 4;
   const int nedges = 4;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, 2 * nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 3;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 3, 2.0, 2.0);  /* edge to be deleted */
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 2, 3, 1.0, 1.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[3] = 1.0;
   graph->prize[2] = 1.0;

   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 2, 0);

   nelims = 0;

   SCIP_CALL( checkSdWalk(scip, FALSE, &graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);

   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE testSdPcKillsTwoEdges(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nelims;
   const int nnodes = 5;
   const int nedges = 6;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, 2 * nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 3;

   /* square */
   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 3, 3.0, 3.0);   /* edge to be deleted */
   graph_edge_add(scip, graph, 1, 2, 1.1, 1.1);
   graph_edge_add(scip, graph, 2, 3, 1.0, 1.0);

   /* lower hat */
   graph_edge_add(scip, graph, 0, 4, 0.5, 0.5);
   graph_edge_add(scip, graph, 1, 4, 0.5, 0.5);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[3] = 1.0;
   graph->prize[2] = 1.0;
   graph->prize[1] = 0.1;


   graph_knot_chg(graph, 3, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 1, 0);

   nelims = 0;

   SCIP_CALL( checkSdWalk(scip, FALSE, &graph, &nelims) );

   STPTEST_ASSERT(nelims == 2);

   assert(graph == NULL);

   return SCIP_OKAY;
}


/** tests that SD star PC test kills an edge */
static
SCIP_RETCODE testSdStarPcKillsEdge(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DIJK dijkdata;
   GRAPH* graph;
   int* star_base;
   const int nnodes_org = 4;
   const int nedges_org = 8;
   int nnodes = -1;
   int nedges = -1;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org - 1; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_add(graph, STP_TERM); // 3
   graph->source = 3;

   graph_edge_addBi(scip, graph, 0, 1, 2.0); // 0,1
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2,3
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 0.0;
   graph->prize[1] = 0.0;
   graph->prize[2] = 0.0;
   graph->prize[3] = 1.1;  /* the pseudo root */

   SCIP_CALL( stptest_graphSetUpPcOrg(scip, graph, &nnodes, &nedges) );
   SCIP_CALL( graph_dijkLimited_init(scip, graph, &dijkdata) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &star_base, nnodes) );

   /* actual test: edge 0 should have been deleted */
   SCIP_CALL( reduce_sdStarPc2(scip, nedges, NULL, graph, dijkdata.distance, star_base, dijkdata.visitlist, dijkdata.visited, dijkdata.dheap, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->oeat[0] == EAT_FREE);
   for( int e = 2; e < nedges; e++ )
      STPTEST_ASSERT(graph->oeat[e] != EAT_FREE);

   SCIPfreeMemoryArray(scip, &star_base);
   graph_dijkLimited_freeMembers(scip, &dijkdata);
   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests PCMW special distance methods */
SCIP_RETCODE stptest_pcreduce(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testRmwTerminalDeg1Contraction2(scip) );


   SCIP_CALL( testRmwAnsDeletesTwoNodes(scip) );
   SCIP_CALL( testRmwAnsDeletesOneEdge(scip) );
   SCIP_CALL( testRmwAnsDeletesOneNode(scip) );
   SCIP_CALL( testRmwTerminalDeg1Contraction3(scip) );
   SCIP_CALL( testRmwTerminalDeg1Contraction1(scip) );
   SCIP_CALL( testRmwTerminalContraction(scip) );


   SCIP_CALL( testSdStarPcKillsEdge(scip) );
   SCIP_CALL( testSdPcKillsEdge1(scip) );
   SCIP_CALL( testSdPcKillsEdge2(scip) );
   SCIP_CALL( testSdPcKillsTwoEdges(scip) );

   printf("PC/MW reduction test: all ok \n");

   return SCIP_OKAY;
}
