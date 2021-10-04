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

/**@file   stptest_da.c
 * @brief  tests for Steiner tree dual-ascent
 * @author Daniel Rehfeldt
 *
 * This file implements unit tests for Steiner tree dual-ascent methods.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

//#define SCIP_DEBUG

#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "graph.h"
#include "dualascent.h"
#include "portab.h"



/** tests that PC DA paths works  */
static
SCIP_RETCODE testDaPathsPcMw3EdgesWorks(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   GRAPH* transgraph;
   SCIP_Real* redcosts;
   const int nnodes_org = 3;
   const int nedges_org = 6;
   SCIP_Real obj;
   SCIP_Real offset = 0.0;

   SCIP_CALL( graph_init(scip, &graph, nnodes_org, nedges_org, 1) );

   for( int i = 0; i < nnodes_org; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_knot_chg(graph, 0, STP_TERM);
   graph_knot_chg(graph, 2, STP_TERM);

   graph_edge_addBi(scip, graph, 0, 1, 4.0); // 0,1
   graph_edge_addBi(scip, graph, 0, 2, 2.0); // 2,3
   graph_edge_addBi(scip, graph, 1, 2, 0.9);

   graph_pc_initPrizes(scip, graph, nnodes_org);
   graph->prize[0] = 2.1; /* the pseudo root */
   graph->prize[1] = 0.0;
   graph->prize[2] = 1.9;

   SCIP_CALL( stptest_graphSetUpPcExtended(scip, graph, NULL, NULL) );
   SCIP_CALL( graph_transPcGetSap(scip, graph, &transgraph, &offset) );
   SCIP_CALL(SCIPallocMemoryArray(scip, &redcosts, transgraph->edges));

   SCIP_CALL( dualascent_pathsPcMw(scip, transgraph, redcosts, &obj, NULL) );

   STPTEST_ASSERT(EQ(redcosts[0], 2.9));
   STPTEST_ASSERT(EQ(redcosts[1], 1.9));
   STPTEST_ASSERT(EQ(redcosts[2], 0.1));
   STPTEST_ASSERT(EQ(redcosts[3], 0.0));

   graph_free(scip, &transgraph, TRUE);
   stptest_graphTearDown(scip, graph);
   SCIPfreeMemoryArray(scip, &redcosts);

   return SCIP_OKAY;
}


/** tests DA paths methods */
SCIP_RETCODE stptest_dapaths(
   SCIP*                 scip                /**< SCIP data structure */
)
{

   SCIP_CALL( testDaPathsPcMw3EdgesWorks(scip) );

   printf("DA paths test: all ok \n");

   return SCIP_OKAY;
}
