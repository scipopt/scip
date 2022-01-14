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

/**@file   stptest_graphpath.c
 * @brief  graph utility tests for Steiner tree problem methods
 * @author Daniel Rehfeldt
 *
 * This file implements graph utility tests for Steiner tree problems.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include "scip/scip.h"
#include "stptest.h"
#include "graph.h"
#include "portab.h"



/** tests that (unbiased) terminal paths are found */
static
SCIP_RETCODE testTerminalPathsTo3NextFound(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_Real dists[4];
   int closeterms[4];
   TPATHS* tpaths;
   GRAPH* graph;
   int nnodes = 6;
   int nedges = 10;
   int ncloseterms = -1;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 2;
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 3, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 3, 2.1); // 2
   graph_edge_addBi(scip, graph, 0, 4, 1.6);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( graph_tpathsInit(scip, graph, &tpaths) );
   graph_tpathsSetAll4(graph, graph->cost, graph->cost, NULL, tpaths);

   graph_tpathsGet4CloseTerms(graph, tpaths, 0, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 3);
   STPTEST_ASSERT(closeterms[0] == 4);
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(closeterms[2] == 3);
   STPTEST_ASSERT(EQ(dists[0], 1.6));
   STPTEST_ASSERT(EQ(dists[1], 2.0));
   STPTEST_ASSERT(EQ(dists[2], 2.1));

   graph_tpathsGet4CloseTerms(graph, tpaths, 1, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 3);
   STPTEST_ASSERT(closeterms[0] == 2);
   STPTEST_ASSERT(closeterms[1] == 4);
   STPTEST_ASSERT(closeterms[2] == 3);
   STPTEST_ASSERT(EQ(dists[0], 1.0));
   STPTEST_ASSERT(EQ(dists[1], 2.6));
   STPTEST_ASSERT(EQ(dists[2], 3.1));

   stptest_graphTearDown(scip, graph);
   graph_tpathsFree(scip, &tpaths);

   return SCIP_OKAY;
}



/** tests that (unbiased) terminal paths are repaired properly (level 1)  */
static
SCIP_RETCODE testTerminalPathsRepair(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_Real dists[4];
   int closeterms[4];
   TPATHS* tpaths;
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 10;
   int ncloseterms = -1;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;
   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.1); // 0
   graph_edge_addBi(scip, graph, 1, 2, 2.1); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0); //
   graph_edge_addBi(scip, graph, 3, 4, 2.2); //
   graph_edge_addBi(scip, graph, 0, 3, 2.15); //

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( graph_tpathsInit(scip, graph, &tpaths) );
   graph_tpathsSetAll4(graph, graph->cost, graph->cost, NULL, tpaths);
   SCIP_CALL( graph_tpathsRepairSetUp(graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 1, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 3);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 1.1));

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 2.1));

   SCIP_CALL( graph_tpathsRepair(scip, 0, FALSE, graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 1, FARAWAY, closeterms, NULL, dists, &ncloseterms);
//   STPTEST_ASSERT(ncloseterms == 1);
   STPTEST_ASSERT(closeterms[0] == 2);
   STPTEST_ASSERT(EQ(dists[0], 2.1));

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 2.15));

   stptest_graphTearDown(scip, graph);
   graph_tpathsFree(scip, &tpaths);

   return SCIP_OKAY;
}


/** tests that (unbiased) terminal paths are repaired properly (level 1 + 2)  */
static
SCIP_RETCODE testTerminalPathsRepair2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_Real dists[4];
   int closeterms[4];
   TPATHS* tpaths;
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 10;
   int ncloseterms = -1;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;
   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);


   graph_edge_addBi(scip, graph, 0, 1, 1.1); // 0
   graph_edge_addBi(scip, graph, 1, 2, 2.0); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0); // 4
   graph_edge_addBi(scip, graph, 3, 4, 2.05); // 6
   graph_edge_addBi(scip, graph, 0, 3, 2.15); // 8

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( graph_tpathsInit(scip, graph, &tpaths) );
   graph_tpathsSetAll4(graph, graph->cost, graph->cost, NULL, tpaths);
   SCIP_CALL( graph_tpathsRepairSetUp(graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 1, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 3);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 1.1));
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(EQ(dists[1], 2.0));

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 4);
   STPTEST_ASSERT(EQ(dists[0], 2.05));
   STPTEST_ASSERT(closeterms[1] == 0);
   STPTEST_ASSERT(EQ(dists[1], 2.1));

   /* delete edge 0-1 */
   SCIP_CALL( graph_tpathsRepair(scip, 0, FALSE, graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 1, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 2);
   STPTEST_ASSERT(EQ(dists[0], 2.0));

   STPTEST_ASSERT(closeterms[1] == 4);
   STPTEST_ASSERT(EQ(dists[1], 3.05));

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 4);
   STPTEST_ASSERT(EQ(dists[0], 2.05));
   STPTEST_ASSERT(closeterms[1] == 0);
   STPTEST_ASSERT(EQ(dists[1], 2.15));

   graph_edge_del(scip, graph, 0, TRUE);

   /* delete edge 0-3 */
   SCIP_CALL( graph_tpathsRepair(scip, 8, FALSE, graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 4);
   STPTEST_ASSERT(EQ(dists[0], 2.05));
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(EQ(dists[1], 3.0));

   stptest_graphTearDown(scip, graph);
   graph_tpathsFree(scip, &tpaths);

   return SCIP_OKAY;
}


/** tests that (unbiased) terminal paths are repaired properly (level 1 + 2 + 3)  */
static
SCIP_RETCODE testTerminalPathsRepair3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_Real dists[4];
   int closeterms[4];
   TPATHS* tpaths;
   GRAPH* graph;
   int nnodes = 6;
   int nedges = 16;
   int ncloseterms = -1;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;
   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 3, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 4, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 5, 1.1); // 4
   graph_edge_addBi(scip, graph, 2, 4, 1.2); // 6
   graph_edge_addBi(scip, graph, 3, 4, 1.0); // 8
   graph_edge_addBi(scip, graph, 3, 5, 1.0); // 10
   graph_edge_addBi(scip, graph, 4, 5, 1.2); // 12
   graph_edge_addBi(scip, graph, 2, 5, 2.1); // 12


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( graph_tpathsInit(scip, graph, &tpaths) );
   graph_tpathsSetAll4(graph, graph->cost, graph->cost, NULL, tpaths);
   SCIP_CALL( graph_tpathsRepairSetUp(graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 3);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 1.0));
   STPTEST_ASSERT(closeterms[1] == 1);
   STPTEST_ASSERT(EQ(dists[1], 2.1));

   graph_tpathsGet4CloseTerms(graph, tpaths, 5, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 1);
   STPTEST_ASSERT(EQ(dists[0], 1.1));
   STPTEST_ASSERT(closeterms[1] == 0);
   STPTEST_ASSERT(EQ(dists[1], 2.0));

   /* delete edge 0-3 */
   SCIP_CALL( graph_tpathsRepair(scip, 0, FALSE, graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 2.0));
   STPTEST_ASSERT(closeterms[1] == 1);
   STPTEST_ASSERT(EQ(dists[1], 2.1));

   graph_tpathsGet4CloseTerms(graph, tpaths, 5, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 1);
   STPTEST_ASSERT(EQ(dists[0], 1.1));
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(EQ(dists[1], 2.1));
   STPTEST_ASSERT(closeterms[2] == 0);
   STPTEST_ASSERT(EQ(dists[2], 2.2));

   graph_tpathsGet4CloseTerms(graph, tpaths, 4, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 0);
   STPTEST_ASSERT(EQ(dists[0], 1.0));
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(EQ(dists[1], 1.2));
   STPTEST_ASSERT(closeterms[2] == 1);
   STPTEST_ASSERT(EQ(dists[2], 2.3));

   graph_edge_del(scip, graph, 0, TRUE);

   /* delete edge 0-4 */
   SCIP_CALL( graph_tpathsRepair(scip, 2, FALSE, graph, tpaths) );

   graph_tpathsGet4CloseTerms(graph, tpaths, 4, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   assert(ncloseterms == 2);
   STPTEST_ASSERT(closeterms[0] == 2);
   STPTEST_ASSERT(EQ(dists[0], 1.2));
   STPTEST_ASSERT(closeterms[1] == 1);
   STPTEST_ASSERT(EQ(dists[1], 2.3));


   graph_tpathsGet4CloseTerms(graph, tpaths, 3, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(closeterms[0] == 1);
   STPTEST_ASSERT(EQ(dists[0], 2.1));
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(EQ(dists[1], 2.2));

   stptest_graphTearDown(scip, graph);
   graph_tpathsFree(scip, &tpaths);

   return SCIP_OKAY;
}


/** tests that SD is repaired properly */
static
SCIP_RETCODE testSdRepair(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SD* sd;
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 22;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 4;
   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 4, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 5, 3.0); // 2
   graph_edge_addBi(scip, graph, 1, 4, 1.1); // 4
   graph_edge_addBi(scip, graph, 2, 5, 2.2); // 6
   graph_edge_addBi(scip, graph, 3, 5, 2.3); // 8
   graph_edge_addBi(scip, graph, 4, 5, 1.0); // 10
   graph_edge_addBi(scip, graph, 1, 5, 1.3); // 12
   graph_edge_addBi(scip, graph, 2, 6, 2.2); // 14
   graph_edge_addBi(scip, graph, 0, 6, 2.4); // 16
   graph_edge_addBi(scip, graph, 5, 6, 2.0); // 18
   graph_edge_addBi(scip, graph, 1, 0, 2.0); // 20

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_sdInit(scip, graph, &sd) );
   SCIP_CALL( reduce_sdRepairSetUp(scip, graph, sd) );


   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 0, 3, FARAWAY, 0.0, sd), 2.3));
   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 1, 2, FARAWAY, 0.0, sd), 2.2));

   SCIP_CALL( reduce_sdRepair(scip, 0, FALSE, graph, sd) );
   graph_edge_del(scip, graph, 0, TRUE);

   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 0, 3, FARAWAY, 0.0, sd), 2.4));
   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 1, 2, FARAWAY, 0.0, sd), 2.2));
   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 0, 3, FARAWAY, 0.0, sd), 2.4));

   SCIP_CALL( reduce_sdRepair(scip, 4, FALSE, graph, sd) );
   graph_edge_del(scip, graph, 4, TRUE);
   SCIP_CALL( reduce_sdRepair(scip, 12, FALSE, graph, sd) );
   graph_edge_del(scip, graph, 12, TRUE);

   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 1, 2, FARAWAY, 0.0, sd), 4.4));
   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 3, 4, FARAWAY, 0.0, sd), 2.3));

   SCIP_CALL( reduce_sdRepair(scip, 8, FALSE, graph, sd) );
   graph_edge_del(scip, graph, 8, TRUE);
   STPTEST_ASSERT(EQ(reduce_sdGetSd(graph, 3, 4, FARAWAY, 0.0, sd), FARAWAY));

   reduce_sdFree(scip, &sd);
   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}

/** tests that (unbiased) terminal paths are found */
static
SCIP_RETCODE testBiasedTerminalPathsTo4NextFound(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_Real dists[4];
   int closeterms[4];
   TPATHS* tpaths;
   SDPROFIT* sdprofit;
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 16;
   int ncloseterms = -1;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 2;
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 3, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 3, 2.1); // 2
   graph_edge_addBi(scip, graph, 0, 4, 1.7);
   graph_edge_addBi(scip, graph, 1, 2, 1.1);
   graph_edge_addBi(scip, graph, 1, 6, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.8);
   graph_edge_addBi(scip, graph, 5, 6, 1.5);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( reduce_sdprofitInit(scip, graph, &sdprofit) );
   SCIP_CALL( graph_tpathsInit(scip, graph, &tpaths) );
   graph_tpathsSetAll4(graph, graph->cost, graph->cost, sdprofit, tpaths);

   graph_tpathsGet4CloseTerms(graph, tpaths, 0, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 4);
   STPTEST_ASSERT(closeterms[0] == 6);
   STPTEST_ASSERT(closeterms[1] == 2);
   STPTEST_ASSERT(closeterms[2] == 4);
   STPTEST_ASSERT(closeterms[3] == 3);
   STPTEST_ASSERT(EQ(dists[0], 1.0));
   STPTEST_ASSERT(EQ(dists[1], 1.6));
   STPTEST_ASSERT(EQ(dists[2], 1.7));
   STPTEST_ASSERT(EQ(dists[3], 2.1));

   graph_tpathsGet4CloseTerms(graph, tpaths, 1, FARAWAY, closeterms, NULL, dists, &ncloseterms);
   STPTEST_ASSERT(ncloseterms == 4);
   STPTEST_ASSERT(closeterms[2] == 4);
   STPTEST_ASSERT(closeterms[3] == 3);
   STPTEST_ASSERT(EQ(dists[0], 1.0));
   STPTEST_ASSERT(EQ(dists[1], 1.1));
   STPTEST_ASSERT(EQ(dists[2], 1.7));
   STPTEST_ASSERT(EQ(dists[3], 3.0));

   stptest_graphTearDown(scip, graph);
   graph_tpathsFree(scip, &tpaths);
   reduce_sdprofitFree(scip, &sdprofit);

   return SCIP_OKAY;
}


/** tests terminal paths */
SCIP_RETCODE stptest_tpaths(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( testSdRepair(scip) );

   SCIP_CALL( testTerminalPathsRepair(scip) );
   SCIP_CALL( testTerminalPathsRepair2(scip) );
   SCIP_CALL( testTerminalPathsRepair3(scip) );

   SCIP_CALL( testTerminalPathsTo3NextFound(scip) );
   SCIP_CALL( testBiasedTerminalPathsTo4NextFound(scip) );

   printf("TPATHS test: all ok \n");

   return SCIP_OKAY;
}
