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

/**@file   stptest_reducesepa.c
 * @brief  tests for Steiner tree node-separator reductions
 * @author Daniel Rehfeldt
 *
 * This file implements unit tests for Steiner tree node-separator based reduction methods.
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
#include "mincut.h"


/** tests terminal separator method */
static
SCIP_RETCODE testTerminalSeparatorsAreFound(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   const int* sepaterms;
   TERMSEPAS* termsepas;
   GRAPH* graph;
   int nnodes = 5;
   int nedges = 10;
   int nsinknodes;
   int sinkterm;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( mincut_termsepasInit(scip, graph, 10, 5, &termsepas) );
   SCIP_CALL( mincut_findTerminalSeparators(scip, NULL, graph, termsepas) );
   sepaterms = mincut_termsepasGetFirst(2, termsepas, &sinkterm, &nsinknodes);

   STPTEST_ASSERT(sepaterms != NULL);
   STPTEST_ASSERT(sinkterm == 4);
   STPTEST_ASSERT(mincut_termsepasGetNall(termsepas) == 1);

   sepaterms = mincut_termsepasGetNext(2, termsepas, &sinkterm, &nsinknodes);
   STPTEST_ASSERT(sepaterms == NULL);

   mincut_termsepasFree(scip, &termsepas);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests terminal separator method */
static
SCIP_RETCODE testTerminalSeparatorsAreFound2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   const int* sepaterms;
   TERMSEPAS* termsepas;
   GRAPH* graph;
   int nnodes = 8;
   int nedges = 22;
   int nsinknodes;
   int sinkterm;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 7, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 6, 1.0);
   graph_edge_addBi(scip, graph, 0, 7, 1.0);

   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);

   graph_edge_addBi(scip, graph, 5, 7, 1.0);
   graph_edge_addBi(scip, graph, 6, 7, 1.0);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( mincut_termsepasInit(scip, graph, 10, 5, &termsepas) );

   SCIP_CALL( mincut_findTerminalSeparators(scip, NULL, graph, termsepas) );
   sepaterms = mincut_termsepasGetFirst(2, termsepas, &sinkterm, &nsinknodes);

   STPTEST_ASSERT(sepaterms != NULL);
   STPTEST_ASSERT(sinkterm == 4 || sinkterm == 5);

   sepaterms = mincut_termsepasGetNext(2, termsepas, &sinkterm, &nsinknodes);

   STPTEST_ASSERT(sepaterms != NULL);
   STPTEST_ASSERT(sinkterm == 4 || sinkterm == 5);

  // printf("%d %d \n", sepaterms[0], sepaterms[1]);

   mincut_termsepasFree(scip, &termsepas);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests terminal separator method */
static
SCIP_RETCODE testTerminalSeparatorsAreFound3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   TERMSEPAS* termsepas;
   GRAPH* graph;
   const int* sepaterms;
   int nsinknodes;
   int sinkterm;
   int nnodes = 9;
   int nedges = 26;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph_knot_chg(graph, 3, STP_TERM);
   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 8, STP_TERM);
   graph_knot_chg(graph, 7, STP_TERM);


   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 0, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 6, 7, 1.0);
   graph_edge_addBi(scip, graph, 2, 8, 1.0);
   graph_edge_addBi(scip, graph, 5, 8, 1.0);
   graph_edge_addBi(scip, graph, 8, 7, 1.0);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( mincut_termsepasInit(scip, graph, 10, 5, &termsepas) );
   SCIP_CALL( mincut_findTerminalSeparators(scip, NULL, graph, termsepas) );
   sepaterms = mincut_termsepasGetNext(3, termsepas, &sinkterm, &nsinknodes);

   STPTEST_ASSERT(sepaterms != NULL);
   STPTEST_ASSERT(sinkterm == 8);
   sepaterms = mincut_termsepasGetNext(3, termsepas, &sinkterm, &nsinknodes);
   STPTEST_ASSERT(sepaterms == NULL);

   mincut_termsepasFree(scip, &termsepas);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests biconnected components method */
static
SCIP_RETCODE testBiconnectedComponentsAreFound(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 8;
   int nedges = 24;
   int nelims = 0;
   SCIP_Real offset = 0.0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 0, STP_TERM);
   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 7, 1.0);
   graph_edge_addBi(scip, graph, 6, 7, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_articulations(scip, graph, &offset, &nelims) );

   STPTEST_ASSERT(graph->grad[5] == 0);
   STPTEST_ASSERT(graph->grad[6] == 0);
   STPTEST_ASSERT(graph->grad[7] == 0);

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests biconnected components method */
static
SCIP_RETCODE testBiconnectedComponentsAreFound2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 11;
   int nedges = 24;
   int nelims = 0;
   SCIP_Real offset = 0.0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 9, STP_TERM);
   graph->source = 9;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 1, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);

   graph_edge_addBi(scip, graph, 0, 7, 1.0);
   graph_edge_addBi(scip, graph, 7, 8, 1.0);
   graph_edge_addBi(scip, graph, 8, 10, 1.0);
   graph_edge_addBi(scip, graph, 8, 9, 1.0);
   graph_edge_addBi(scip, graph, 9, 10, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_articulations(scip, graph, &offset, &nelims) );

   STPTEST_ASSERT(graph->grad[2] == 0);
   STPTEST_ASSERT(graph->grad[3] == 0);
   STPTEST_ASSERT(graph->grad[4] == 0);
   STPTEST_ASSERT(graph->grad[5] == 0);
   STPTEST_ASSERT(graph->grad[6] == 0);
   STPTEST_ASSERT(Is_term(graph->term[0]));
   STPTEST_ASSERT(Is_term(graph->term[7]));
   STPTEST_ASSERT(Is_term(graph->term[8]));

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}





/** tests biconnected components method */
static
SCIP_RETCODE testBiconnectedComponentsAreFound3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 16;
   int nelims = 0;
   SCIP_Real offset = 0.0;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph->source = 2;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 0, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_articulations(scip, graph, &offset, &nelims) );

   STPTEST_ASSERT(graph->grad[1] == 0);
   STPTEST_ASSERT(Is_term(graph->term[0]));
   STPTEST_ASSERT(Is_term(graph->term[4]));

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests biconnected components method */
static
SCIP_RETCODE testBiconnectedDecomposition(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDBASE* redbase;
   GRAPH* graph;
   int nnodes = 8;
   int nedges = 22;
   SCIP_Real obj = 0.0;
   SCIP_Bool wasDecomposed;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 4, STP_TERM);
   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 6, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);
 //  graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 7, 1.0);
   graph_edge_addBi(scip, graph, 6, 7, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );


   SCIP_CALL( reduce_baseInit(scip, graph, &redbase) );

   {
      BIDECPARAMS decparameters = { .depth = 0, .maxdepth = 2, .newLevelStarted = FALSE };
      RPARAMS parameters = { .dualascent = 1, .boundreduce = 1, .nodereplacing = 1, .reductbound_min = 10,
                             .reductbound = 10, .userec = 1, .fullreduce = 1, .usestrongreds = TRUE };
      redbase->redparameters = &parameters;
      redbase->bidecompparams = &decparameters;
      graph->stp_type = STP_SPG;

      SCIP_CALL( reduce_solInit(scip, graph, TRUE, &(redbase->redsol)) );
      SCIP_CALL( reduce_solInitLocal(scip, graph, redbase->redsol, NULL));

      SCIP_CALL( reduce_bidecomposition(scip, graph, redbase, NULL, &wasDecomposed) );

      reduce_solFinalizeLocal(scip, graph, redbase->redsol);
   }

   obj = reduce_solGetUpperBoundWithOffset(redbase->redsol);
   reduce_solFree(scip, &(redbase->redsol));
   reduce_baseFree(scip, &redbase);

   STPTEST_ASSERT(graph->terms == 1);
   STPTEST_ASSERT(EQ(obj, 5.0));
   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}


/** tests biconnected decomposition method */
static
SCIP_RETCODE testBiconnectedDecomposition2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDBASE* redbase;
   GRAPH* graph;
   int nnodes = 11;
   int nedges = 24;
   SCIP_Real obj = 0.0;
   SCIP_Bool wasDecomposed;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 1, STP_TERM);
   graph_knot_chg(graph, 9, STP_TERM);
   graph->source = 9;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 1, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);

   graph_edge_addBi(scip, graph, 0, 7, 1.0);
   graph_edge_addBi(scip, graph, 7, 8, 1.0);
   graph_edge_addBi(scip, graph, 8, 10, 1.0);
   graph_edge_addBi(scip, graph, 8, 9, 1.0);
   graph_edge_addBi(scip, graph, 9, 10, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_baseInit(scip, graph, &redbase) );

   {
      BIDECPARAMS decparameters = { .depth = 0, .maxdepth = 2, .newLevelStarted = FALSE };
      RPARAMS parameters = { .dualascent = 1, .boundreduce = 1, .nodereplacing = 1, .reductbound_min = 10,
                             .reductbound = 10, .userec = 1, .fullreduce = 1, .usestrongreds = TRUE };
      redbase->redparameters = &parameters;
      redbase->bidecompparams = &decparameters;
      graph->stp_type = STP_SPG;

      SCIP_CALL( reduce_solInit(scip, graph, TRUE, &(redbase->redsol)) );
      SCIP_CALL( reduce_solInitLocal(scip, graph, redbase->redsol, NULL));

      SCIP_CALL( reduce_bidecomposition(scip, graph, redbase, NULL, &wasDecomposed) );

      reduce_solFinalizeLocal(scip, graph, redbase->redsol);
   }

   obj = reduce_solGetUpperBoundWithOffset(redbase->redsol);
   reduce_solFree(scip, &(redbase->redsol));
   reduce_baseFree(scip, &redbase);

   STPTEST_ASSERT(graph->terms == 1);
   STPTEST_ASSERT(EQ(obj, 4.0));

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}

/** tests biconnected decomposition method */
static
SCIP_RETCODE testBiconnectedDecomposition3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int nnodes = 7;
   int nedges = 16;
   SCIP_Real obj = 0.0;
   REDBASE* redbase;
   SCIP_Bool wasDecomposed;

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph_knot_chg(graph, 5, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);
   graph->source = 2;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); // 0
   graph_edge_addBi(scip, graph, 0, 2, 1.0); // 2
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 0, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   SCIP_CALL( reduce_baseInit(scip, graph, &redbase) );

   {
      BIDECPARAMS decparameters = { .depth = 0, .maxdepth = 2, .newLevelStarted = FALSE };
      RPARAMS parameters = { .dualascent = 1, .boundreduce = 1, .nodereplacing = 1,
                             .reductbound_min = 10, .reductbound = 10, .userec = 1, .fullreduce = 1, .usestrongreds = TRUE };
      redbase->redparameters = &parameters;
      redbase->bidecompparams = &decparameters;
      graph->stp_type = STP_SPG;

      SCIP_CALL( reduce_solInit(scip, graph, TRUE, &(redbase->redsol)) );
      SCIP_CALL( reduce_solInitLocal(scip, graph, redbase->redsol, NULL));

      SCIP_CALL( reduce_bidecomposition(scip, graph, redbase, NULL, &wasDecomposed) );

      reduce_solFinalizeLocal(scip, graph, redbase->redsol);
   }

   obj = reduce_solGetUpperBoundWithOffset(redbase->redsol);
   reduce_solFree(scip, &(redbase->redsol));
   reduce_baseFree(scip, &redbase);

   STPTEST_ASSERT(graph->terms == 1);
   STPTEST_ASSERT(EQ(obj, 3.0));

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}

/** tests simple cut nodes methods */
SCIP_RETCODE stptest_reduceBiconnected(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( testTerminalSeparatorsAreFound(scip) );
   SCIP_CALL( testTerminalSeparatorsAreFound3(scip) );
   SCIP_CALL( testTerminalSeparatorsAreFound2(scip) );

   SCIP_CALL( testBiconnectedDecomposition(scip) );
   SCIP_CALL( testBiconnectedDecomposition2(scip) );
   SCIP_CALL( testBiconnectedDecomposition3(scip) );

   SCIP_CALL( testBiconnectedComponentsAreFound(scip) );
   SCIP_CALL( testBiconnectedComponentsAreFound2(scip) );
   SCIP_CALL( testBiconnectedComponentsAreFound3(scip) );

   printf("reduce biconnected test: all ok \n");


   return SCIP_OKAY;
}
