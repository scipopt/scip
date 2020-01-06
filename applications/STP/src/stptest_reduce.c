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

/**@file   stptest_reduce.c
 * @brief  tests for Steiner tree reductions
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
#include "graph.h"
#include "reduce.h"
#include "heur_local.h"
#include "heur_tm.h"




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
   SCIP_CALL( graph_pc_2pc(scip, graph) );

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

static
SCIP_RETCODE sdPcMwTest1(
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

   assert(nelims == 1);
   assert(graph == NULL);

   return SCIP_OKAY;
}

static
SCIP_RETCODE sdPcMwTest2(
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

   assert(nelims == 1);

   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE sdPcMwTest4(
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

   assert(nelims == 2);

   assert(graph == NULL);

   return SCIP_OKAY;
}


/** tests PCMW special distance methods */
SCIP_RETCODE stptest_reduce_sdpcmw(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( sdPcMwTest1(scip) );
   SCIP_CALL( sdPcMwTest2(scip) );
   SCIP_CALL( sdPcMwTest4(scip) );

   printf("stptest_reduce_sdpcmw: all ok \n");

   return SCIP_OKAY;
}
