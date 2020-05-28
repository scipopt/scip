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

/**@file   stptest_reducesd.c
 * @brief  tests for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements unit tests for Steiner special distance (Steiner bottleneck distance) reduction methods.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

//#define SCIP_DEBUG

#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "graph.h"
#include "reduce.h"
#include "portab.h"




/** tests SD getter */
static
SCIP_RETCODE testSdGetterReturnsCorrectSds(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SD* sddata;
   GRAPH* graph;
   GRAPH* cliquegraph;
   DIJK* dijkdata;
   int nnodes = 7;
   int nedges = 16;
   int cliqueNodeMap[] = { 3,4,5,6 };

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 3.0); // 0
   graph_edge_addBi(scip, graph, 1, 2, 2.0); // 2

   graph_edge_addBi(scip, graph, 3, 0, 1.0); // 4
   graph_edge_addBi(scip, graph, 4, 0, 1.0); // 6
   graph_edge_addBi(scip, graph, 4, 1, 2.0); // 8
   graph_edge_addBi(scip, graph, 5, 1, 4.0); // 10
   graph_edge_addBi(scip, graph, 6, 1, 4.0); // 12
   graph_edge_addBi(scip, graph, 6, 2, 2.0); // 14

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( graph_buildCompleteGraph(scip, &cliquegraph, 4) );
   SCIP_CALL( graph_path_init(scip, cliquegraph) );
   for( int i = 0; i < 4; i++ )
      cliquegraph->mark[i] = TRUE;

   SCIP_CALL( reduce_sdInit(scip, graph, &sddata) );

   SCIP_CALL( graph_dijkLimited_init(scip, graph, &(dijkdata)) );
   graph_dijkLimited_clean(graph, dijkdata);
   dijkdata->edgelimit = 500;

   SCIP_CALL( reduce_sdGetSdsCliquegraph(scip, graph, cliqueNodeMap, dijkdata, sddata, cliquegraph) );

   for( int e = cliquegraph->outbeg[0]; e != EAT_LAST; e = cliquegraph->oeat[e] )
   {
      const int head = cliquegraph->head[e];
      const SCIP_Real cost = cliquegraph->cost[e];
      if( head == 1 )
      {
         STPTEST_ASSERT(EQ(cost, 1.0));
      }

      if( head == 2 )
      {
         STPTEST_ASSERT(EQ(cost, 3.0));
      }

      if( head == 3 )
      {
         STPTEST_ASSERT(EQ(cost, 3.0));
      }
   }

   for( int e = cliquegraph->outbeg[1]; e != EAT_LAST; e = cliquegraph->oeat[e] )
   {
      const int head = cliquegraph->head[e];
      const SCIP_Real cost = cliquegraph->cost[e];
      if( head == 0 )
      {
         STPTEST_ASSERT(EQ(cost, 1.0));
      }

      if( head == 2 )
      {
         STPTEST_ASSERT(EQ(cost, 3.0));
      }

      if( head == 3 )
      {
         STPTEST_ASSERT(EQ(cost, 2.0));
      }
   }

   graph_dijkLimited_free(scip, &(dijkdata));
   reduce_sdFree(scip, &sddata);
   stptest_graphTearDown(scip, cliquegraph);
   stptest_graphTearDown(scip, graph);


   return SCIP_OKAY;
}



/** tests SD getter */
static
SCIP_RETCODE testSdGraphDistsAreValid(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SD* sddata;
   int nnodes = 5;
   int nedges = 8;
   const SCIP_Real* sddists;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 3, STP_TERM);
   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.9); // 0
   graph_edge_addBi(scip, graph, 0, 2, 2.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 2.0); // 4
   graph_edge_addBi(scip, graph, 1, 4, 1.0); // 4


   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_sdInitBiased(scip, graph, &(sddata)) );
   reduce_sdgraphInitOrderedMstCosts(sddata->sdgraph);

   sddists = reduce_sdgraphGetOrderedMstCosts(sddata->sdgraph);


   STPTEST_ASSERT(EQ(sddists[0], 2.0));
   STPTEST_ASSERT(EQ(sddists[1], 2.0));

   reduce_sdFree(scip, &sddata);
   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests SD getter */
static
SCIP_RETCODE testSdGraphDistsAreValid2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SD* sddata;
   int nnodes = 5;
   int nedges = 10;
   const SCIP_Real* sddists;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 3, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 2.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 2.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 2.0); // 4
   graph_edge_addBi(scip, graph, 2, 4, 3.5);
   graph_edge_addBi(scip, graph, 3, 4, 2.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_sdInitBiased(scip, graph, &sddata) );
   reduce_sdgraphInitOrderedMstCosts(sddata->sdgraph);

   sddists = reduce_sdgraphGetOrderedMstCosts(sddata->sdgraph);

   STPTEST_ASSERT(EQ(sddists[0], 3.5));
   STPTEST_ASSERT(EQ(sddists[1], 2.5));

   reduce_sdFree(scip, &sddata);
   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests that BDk test pseudo-eliminates node of degree 4 */
static
SCIP_RETCODE testBdkTreeDistDeletesNodeDeg4(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 12;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);
   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4

   graph_edge_addBi(scip, graph, 0, 4, 1.1); // 6
   graph_edge_addBi(scip, graph, 1, 5, 1.0); // 8
   graph_edge_addBi(scip, graph, 5, 6, 1.0); // 10

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_bdk(scip, 100, graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->grad[0] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests that BDk test pseudo-eliminates node of degree 4 */
static
SCIP_RETCODE testBdkSdMstDeletesNodeDeg4(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 8;
   int nedges = 20;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);
   graph_knot_chg(graph, 7, STP_TERM);

   graph->source = 5;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 2.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 0, 4, 1.0); // 6

   graph_edge_addBi(scip, graph, 1, 5, 1.0); // 8
   graph_edge_addBi(scip, graph, 3, 5, 1.0); // 10
   graph_edge_addBi(scip, graph, 3, 6, 1.0); // 14
   graph_edge_addBi(scip, graph, 4, 6, 1.0); // 14
   graph_edge_addBi(scip, graph, 3, 7, 1.0); // 14


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_bdk(scip, 100, graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->grad[0] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests that BDk test pseudo-eliminates node of degree 4 */
static
SCIP_RETCODE testBdkSdMstStarDeletesNodeDeg4(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 9;
   int nedges = 22;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);
   graph_knot_chg(graph, 7, STP_TERM);
   graph_knot_chg(graph, 8, STP_TERM);


   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 2.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 0, 4, 1.0); // 6

   graph_edge_addBi(scip, graph, 1, 5, 0.5); // 8
   graph_edge_addBi(scip, graph, 3, 5, 0.5); // 10
   graph_edge_addBi(scip, graph, 3, 6, 1.0); //
   graph_edge_addBi(scip, graph, 4, 6, 1.0); //
   graph_edge_addBi(scip, graph, 3, 7, 1.0); //
   graph_edge_addBi(scip, graph, 7, 8, 1.0); //


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_bdk(scip, 100, graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->grad[0] == 0);

   stptest_graphTearDown(scip, graph);


   return SCIP_OKAY;
}


/** tests that BDk test pseudo-eliminates node of degree 3 */
static
SCIP_RETCODE testBdkSdMstDeletesNodeDeg3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 16;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);

   graph->source = 5;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 4

   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 2.0);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_bdk(scip, 100, graph, &nelims) );

   STPTEST_ASSERT(graph->grad[0] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests that SD star biased test finds edge for deletion */
static
SCIP_RETCODE testSdStarBiasedDeletesEdge(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 10;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_knot_chg(graph, 4, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 2.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 2, 3, 1.0); // 6
   graph_edge_addBi(scip, graph, 3, 4, 1.0); // 6

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_sdStarBiased(scip, 10, NULL, graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->oeat[0] == EAT_FREE);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}

/** tests that SD star biased test finds edge for deletion */
static
SCIP_RETCODE testSdStarBiasedDeletesEdge2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 6;
   int nedges = 12;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_knot_chg(graph, 1, STP_TERM);    /* not necessary ... */
   graph_knot_chg(graph, 4, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 2.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_sdStarBiased(scip, 10, NULL, graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->oeat[0] == EAT_FREE);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests that SD star biased test finds edge for deletion */
static
SCIP_RETCODE testSdStarBiasedDeletesEdge3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 6;
   int nedges = 14;
   int nelims = 0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 4;

   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0); // 2
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 0.9);
   graph_edge_addBi(scip, graph, 3, 5, 0.9);
   graph_edge_addBi(scip, graph, 4, 5, 1.9);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_sdStarBiased(scip, 20, NULL, graph, &nelims) );

   STPTEST_ASSERT(nelims == 1);
   STPTEST_ASSERT(graph->oeat[0] == EAT_FREE);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests clique star correctly identifies adjacency distances for degree 3 node  */
static
SCIP_RETCODE testSdCliqueStarDeg3AdjacencyIsCorrect(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 12;
   const int nsds = 3;
   SCIP_Real sds[3];
   int cliquenodes[] = { 1, 2, 3 };
   SDCLIQUE cliquedata = { .dijkdata = NULL, .cliquenodes = cliquenodes, .ncliquenodes = 3, .sds = sds };

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nsds; i++ )
      sds[i] = FARAWAY - 1.0;

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 4;

   graph_knot_chg(graph, 4, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 0.99); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.8);
   graph_edge_addBi(scip, graph, 1, 3, 1.9);
   graph_edge_addBi(scip, graph, 1, 4, 3.9); // dummy terminal

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( graph_dijkLimited_init(scip, graph, &(cliquedata.dijkdata)) );
   graph_dijkLimited_clean(graph, (cliquedata.dijkdata));
   cliquedata.dijkdata->edgelimit = 50;

   SCIP_CALL( graph_sdComputeCliqueStar(scip, graph, NULL, &cliquedata) );

   STPTEST_ASSERT(EQ(sds[0], 1.8));
   STPTEST_ASSERT(EQ(sds[1], 1.9));
   STPTEST_ASSERT(EQ(sds[2], 1.99));

   graph_dijkLimited_free(scip, &(cliquedata.dijkdata));

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests clique star correctly identifies distances for degree 4 node  */
static
SCIP_RETCODE testSdCliqueStarDeg4IsCorrect(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 16;
   const int nsds = 6;
   SCIP_Real sds[6];
   int cliquenodes[] = { 1, 2, 3, 4 };
   SDCLIQUE cliquedata = { .dijkdata = NULL, .cliquenodes = cliquenodes, .ncliquenodes = 4, .sds = sds };

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nsds; i++ )
      sds[i] = FARAWAY - 1.0;

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 6; // dummy
   graph_knot_chg(graph, 6, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 0.9);
   graph_edge_addBi(scip, graph, 0, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 0.4);
   graph_edge_addBi(scip, graph, 2, 5, 0.4);
   graph_edge_addBi(scip, graph, 2, 3, 1.5);

   graph_edge_addBi(scip, graph, 1, 6, 3.9); // dummy terminal

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( graph_dijkLimited_init(scip, graph, &(cliquedata.dijkdata)) );
   graph_dijkLimited_clean(graph, (cliquedata.dijkdata));
   cliquedata.dijkdata->edgelimit = 50;

   SCIP_CALL( graph_sdComputeCliqueStar(scip, graph, NULL, &cliquedata) );

   STPTEST_ASSERT(EQ(sds[0], 0.8)); // 1-2
   STPTEST_ASSERT(EQ(sds[1], 1.9)); // 1-3
   STPTEST_ASSERT(EQ(sds[2], 2.0)); // 1-4
   STPTEST_ASSERT(EQ(sds[3], 1.5)); // 2-3
   STPTEST_ASSERT(EQ(sds[4], FARAWAY - 1.0)); // 2-4
   STPTEST_ASSERT(EQ(sds[5], 1.9)); // 3-4

   graph_dijkLimited_free(scip, &(cliquedata.dijkdata));

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}

/** tests Bdk methods */
SCIP_RETCODE stptest_reduceBdk(
   SCIP*                 scip                /**< SCIP data structure */
)
{

   SCIP_CALL( testBdkSdMstDeletesNodeDeg3(scip) );
   SCIP_CALL( testBdkTreeDistDeletesNodeDeg4(scip) );
   SCIP_CALL( testBdkSdMstDeletesNodeDeg4(scip) );
   SCIP_CALL( testBdkSdMstStarDeletesNodeDeg4(scip) );




   printf("reduce BDk test: all ok \n");

   return SCIP_OKAY;
}


/** tests SD clique star methods */
SCIP_RETCODE stptest_reduceSdCliqueStar(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testSdCliqueStarDeg3AdjacencyIsCorrect(scip) );
   SCIP_CALL( testSdCliqueStarDeg4IsCorrect(scip) );

   printf("reduce clique star test: all ok \n");

   return SCIP_OKAY;
}


/** tests SD getter methods */
SCIP_RETCODE stptest_reduceSdGetter(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testSdGetterReturnsCorrectSds(scip) );
   SCIP_CALL( testSdGraphDistsAreValid(scip) );
   SCIP_CALL( testSdGraphDistsAreValid2(scip) );

   printf("reduce sd getter test: all ok \n");

   return SCIP_OKAY;
}


/** tests SD biased methods */
SCIP_RETCODE stptest_reduceSdStarBias(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testSdStarBiasedDeletesEdge(scip) );
   SCIP_CALL( testSdStarBiasedDeletesEdge2(scip) );
   SCIP_CALL( testSdStarBiasedDeletesEdge3(scip) );

   printf("reduce SD test: all ok \n");

   return SCIP_OKAY;
}
