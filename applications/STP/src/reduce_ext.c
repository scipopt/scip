/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_bnd.c
 * @brief  extended reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements extended reduction techniques for several Steiner problems.
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <grph.h>


static
SCIP_RETCODE reduce_extArc(
   SCIP* scip,
   GRAPH* graph,
   SCIP_Real* rootdist,
   SCIP_Real* redcost,
   PATH* termpaths,
   STP_Bool* edgedeleted,
   SCIP_Real cutoff,
   int edge,
   int root,
   SCIP_Bool* deletable
)
{
   const int nnodes = graph->knots;
   int* nodearr_int1;
   int* nodearr_int2;
   int* nodearr_int3;
   int* nodearr_int4;
   int* nodearr_int5;
   SCIP_Bool* isterm;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int1, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int3, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int4, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearr_int5, nnodes) );

   SCIP_CALL( SCIPallocBufferArray(scip, &isterm, nnodes) );

   graph_get_isTerm(graph, isterm);

   /* actual test */
   SCIP_CALL( reduceExtCheckArc(scip, graph, root, redcost, rootdist, termpaths, edgedeleted,
         isterm, cutoff, edge, FALSE, nodearr_int1, nodearr_int2, nodearr_int3, nodearr_int4, nodearr_int5, deletable) );

   /* clean up */
   SCIPfreeBufferArray(scip, &isterm);
   SCIPfreeBufferArray(scip, &nodearr_int5);
   SCIPfreeBufferArray(scip, &nodearr_int4);
   SCIPfreeBufferArray(scip, &nodearr_int3);
   SCIPfreeBufferArray(scip, &nodearr_int2);
   SCIPfreeBufferArray(scip, &nodearr_int1);

   return SCIP_OKAY;
}


SCIP_RETCODE reduce_extTest1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 18;
   const int root = 0;

   SCIP_Real* rootdist;
   SCIP_Real* redcost;
   PATH* termpaths;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff;
   int edge;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   SCIP_CALL( SCIPallocBufferArray(scip, &rootdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termpaths, 3 * nnodes) );


   /* build graph */
   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 4, 1.0, 1.0);

   graph_edge_add(scip, graph, 2, 5, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 6, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 7, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 8, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 9, 1.0, 1.0);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   for( int i = 1; i < nnodes; i++ )
      graph->mark[i] = (graph->grad[i] > 0);


   /* necessary data structures */
   for( int i = 0; i < nnodes; i++ )
   {
      rootdist[i] = 0.0;
      termpaths[i].dist = 0.0;
   }

   for( int i = 0; i < nedges; i++ )
      redcost[i] = 0.0;

   cutoff = 0.0;
   edge = 0;

   graph_edge_printInfo(graph, edge);

   SCIP_CALL(reduce_extArc(scip, graph, rootdist, redcost, termpaths, edgedeleted, cutoff, edge, root, &deletable));

   assert(deletable);

   assert(0);


   /* clean up */

   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}



SCIP_RETCODE reduce_extTest2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;

   SCIP_Real* rootdist;
   SCIP_Real* redcost;
   PATH* termpaths;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff;
   int edge;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   SCIP_CALL( SCIPallocBufferArray(scip, &rootdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termpaths, 3 * nnodes) );


   /* build graph */
   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 4, 1.0, 1.0);

   graph_edge_add(scip, graph, 2, 5, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 6, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 7, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 8, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 9, 1.0, 1.0);

   graph_edge_add(scip, graph, 5, 10, 1.0, 1.0);
   graph_edge_add(scip, graph, 10, 11, 1.0, 1.0);
   graph_edge_add(scip, graph, 11, 12, 1.0, 1.0);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   for( int i = 1; i < nnodes; i++ )
      graph->mark[i] = (graph->grad[i] > 0);


   /* necessary data structures */
   for( int i = 0; i < nnodes; i++ )
   {
      rootdist[i] = 0.0;
      termpaths[i].dist = 0.2;
   }

   for( int i = 0; i < nedges; i++ )
      redcost[i] = 1.0 + (double) i / 20;

   cutoff = 100.0;
   edge = 0;

   graph_edge_printInfo(graph, edge);

   SCIP_CALL(reduce_extArc(scip, graph, rootdist, redcost, termpaths, edgedeleted, cutoff, edge, root, &deletable));

   assert(!deletable);

   graph_knot_del(scip, graph, 12, TRUE);
   graph->mark[12] = FALSE;

   SCIP_CALL(reduce_extArc(scip, graph, rootdist, redcost, termpaths, edgedeleted, cutoff, edge, root, &deletable));
   assert(deletable);


   assert(0);


   /* clean up */

   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}
