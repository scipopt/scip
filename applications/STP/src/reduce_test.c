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

/**@file   reduce_test.c
 * @brief  tests for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problem reductions.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include <stdio.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "heur_local.h"
#include "heur_tm.h"



static
SCIP_RETCODE graphBuildV5E5(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               g,                  /**< the graph */
   SCIP_Bool             pc                  /**< create pc graph? */
)
{
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 5;

   SCIP_CALL(graph_init(scip, g, nnodes, 2 * nedges, 1));
   graph = *g;

   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_add(graph, 0);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 4, 1.0, 1.0);
   graph_edge_add(scip, graph, 2, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 4, 1.0, 1.0);

   if( pc )
   {
      SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

      for( int i = 0; i < nnodes; i++ )
         graph->prize[i] = 0.0;

      graph->prize[0] = 1.0;

      SCIP_CALL( graph_pc_2pc(scip, graph) );
   }

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   return SCIP_OKAY;
}

/** base method for extended edge reduction tests */
static
SCIP_RETCODE extArc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata,        /**< reduced cost data */
   STP_Bool*             edgedeleted,
   int                   edge,
   int                   edgedelete,
   int                   nclosenodes,        /**< max. number of close nodes to each node */
   SCIP_Bool*            deletable,
   SCIP_Bool             equality
)
{
   const int nnodes = graph->knots;
   SCIP_Bool* isterm;
   int* tree_deg;
   SCIP_Real* bottleneckDists;
   SCIP_Real* pcSdToNode = NULL;
   DISTDATA distdata;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( reduce_distDataInit(scip, graph, nclosenodes, FALSE, &distdata) );

   SCIP_CALL( SCIPallocBufferArray(scip, &isterm, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tree_deg, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bottleneckDists, nnodes) );

   if( graph_pc_isPcMw(graph) )
      SCIP_CALL( SCIPallocBufferArray(scip, &pcSdToNode, nnodes) );

   for( int i = 0; i < nnodes; i++ )
   {
      bottleneckDists[i] = -1.0;
      tree_deg[i] = 0;
   }

   graph_get_isTerm(graph, isterm);

   if( edgedelete >= 0 )
   {
      graph_edge_delFull(scip, graph, edgedelete, TRUE);
      reduce_distDataDeleteEdge(scip, graph, edgedelete, &distdata);

      graph_mark(graph);
   }

   /* actual test */
   SCIP_CALL( reduce_extendedCheckArc(scip, graph, redcostdata, edgedeleted,
         isterm, edge, equality, &distdata, bottleneckDists, pcSdToNode, tree_deg, deletable) );

   /* clean up */
   SCIPfreeBufferArrayNull(scip, &pcSdToNode);
   SCIPfreeBufferArray(scip, &bottleneckDists);
   SCIPfreeBufferArray(scip, &tree_deg);
   SCIPfreeBufferArray(scip, &isterm);

   reduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}


static
void ext0InitRedCostArrays(
   const GRAPH*          graph,              /**< the graph */
   SCIP_Real*            redcost,
   SCIP_Real*            rootdist,
   int*                  vbase3,
   PATH*                 termpaths3
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;

   for( int i = 0; i < nnodes; i++ )
      rootdist[i] = 0.0;

   for( int i = 0; i < 3 * nnodes; i++ )
   {
      vbase3[i] = graph->source;
      termpaths3[i].dist = 0.0;
      termpaths3[i].edge = -1;
   }

   for( int i = 0; i < nedges; i++ )
      redcost[i] = 0.0;
}


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
SCIP_RETCODE extTest5_variants(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2 */
)
{
   GRAPH* graph;
   const int nnodes = 85;
   const int nedges = 88;
   const int root = 0;

   int* vbase;
   SCIP_Real* rootdist;
   SCIP_Real* redcost;
   PATH* termpaths;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff;
   int edge;
   SCIP_Bool deletable;

   assert(variant == 1 || variant == 2);

   SCIP_CALL( SCIPallocBufferArray(scip, &rootdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termpaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 4, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 5, 1.0, 1.0);

   graph_edge_add(scip, graph, 2, 6, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 7, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 8, 1.0, 1.0);
   graph_edge_add(scip, graph, 5, 9, 1.0, 1.0);

   graph_edge_add(scip, graph, 6, 10, 1.0, 1.0);
   graph_edge_add(scip, graph, 7, 10, 1.0, 1.0);
   graph_edge_add(scip, graph, 8, 10, 1.0, 1.0);
   graph_edge_add(scip, graph, 9, 10, 1.0, 1.0);

   graph_edge_add(scip, graph, 0, 10, 0.9, 0.9);

   graph_mark(graph);

   ext0InitRedCostArrays(graph, redcost, rootdist, vbase, termpaths);

   for( int i = 0; i < nedges; i++ )
      redcost[i] = 1.0;

   cutoff = 100.0;
   edge = 0;

   graph_knot_chg(graph, 6, 0);
   graph_knot_chg(graph, 7, 0);
   graph_knot_chg(graph, 8, 0);
   graph_knot_chg(graph, 9, 0);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   if( variant == 1 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(deletable);
   }
   else
   {
      int edgedelete = -1;
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};

      assert(variant == 2);

      for( int e = graph->outbeg[0]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( graph->head[e] == 10 )
            edgedelete = e;
      }

      assert(edgedelete >= 0);

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }

   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE extTest4_variants(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2,3 */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;

   SCIP_Real* rootdist;
   int* vbase;
   SCIP_Real* redcost;
   PATH* termpaths;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff;
   int edge;
   SCIP_Bool deletable;

   assert(scip);
   assert(variant == 1 || variant == 2 || variant == 3);

   SCIP_CALL( SCIPallocBufferArray(scip, &rootdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termpaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 2, 4, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 5, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 6, 1.0, 1.0);
   graph_edge_add(scip, graph, 5, 7, 1.0, 1.0);

   graph_edge_add(scip, graph, 0, 6, 0.9, 0.9);
   graph_edge_add(scip, graph, 0, 7, 0.9, 0.9);

   if( variant == 3 )
      graph_edge_add(scip, graph, 7, 8, 0.9, 0.9);

   graph_mark(graph);

   ext0InitRedCostArrays(graph, redcost, rootdist, vbase, termpaths);

   for( int i = 0; i < nedges; i++ )
      redcost[i] = 1.0;

   cutoff = 100.0;
   edge = 0;

   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 5, 0);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   if( variant == 1 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(deletable);
   }
   else if( variant == 2 )
   {
      int edgedelete = -1;
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};

      for( int e = graph->outbeg[0]; e != EAT_LAST; e = graph->oeat[e] )
         if( graph->head[e] == 7 )
            edgedelete = e;

      assert(edgedelete >= 0);

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }
   else
   {
      int edgedelete = -1;
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};

      assert(variant == 3);

      for( int e = graph->outbeg[7]; e != EAT_LAST; e = graph->oeat[e] )
         if( graph->head[e] == 8 )
            edgedelete = e;

      assert(edgedelete >= 0);

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(deletable);
   }

   /* clean up */
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE extTest3_variants(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2,3,4,5 */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;

   int* vbase;
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
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );

   /* build graph */
   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 2, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 2, 4, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 5, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 6, 1.0, 1.0);
   graph_edge_add(scip, graph, 4, 7, 1.0, 1.0);

   graph_edge_add(scip, graph, 7, 9, 0.9, 0.9);
   graph_edge_add(scip, graph, 8, 0, 1.0, 1.0);
   graph_edge_add(scip, graph, 9, 0, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, nnodes - 1, 1.0, 1.0);

   graph_knot_chg(graph, nnodes - 1, 0); /* dummy node */
   graph_mark(graph);

   ext0InitRedCostArrays(graph, redcost, rootdist, vbase, termpaths);

   vbase[5] = nnodes - 1;
   vbase[6] = nnodes - 1;
   termpaths[5].dist = 3.0;
   termpaths[6].dist = 3.0;
   termpaths[5 + nnodes].dist = 3.0;
   termpaths[6 + nnodes].dist = 3.0;

   for( int i = 0; i < nedges; i++ )
      redcost[i] = 1.0;

   cutoff = 7.0;
   edge = 0;

   graph_knot_chg(graph, 7, 0);

   if( variant == 1 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(deletable);
   }
   else if( variant == 2 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));
      assert(deletable);
   }
   else if( variant == 3 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      const int edge2 = 14;
      assert(graph->tail[edge2] == 7 && graph->head[edge2] == 9);

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      graph->cost[edge2] = 1.0;
      graph->cost[flipedge(edge2)] = 1.0;
      graph_knot_chg(graph, 4, 0);

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(!deletable);
   }
   else if( variant == 4 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      const int edgedelete = 14;
      assert(graph->tail[edgedelete] == 7 && graph->head[edgedelete] == 9);

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }
   else
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      const int edgedelete = graph->edges;

      graph_edge_add(scip, graph, 7, 0, 0.1, 0.1);
      assert(graph->tail[edgedelete] == 7 && graph->head[edgedelete] == 0);

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(deletable);
   }

   /* clean up */
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE extTest2_variants(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2,3,4,5 */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;

   SCIP_Real* rootdist;
   SCIP_Real* redcost;
   PATH* termpaths;
   int* vbase;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff;
   int edge;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   SCIP_CALL( SCIPallocBufferArray(scip, &rootdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termpaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );

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

   graph_mark(graph);

   ext0InitRedCostArrays(graph, redcost, rootdist, vbase, termpaths);

   cutoff = 100.0;
   edge = 0;

   if( variant == 1 )
   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 1, graph) );

      for( int e = graph->outbeg[10]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( graph->head[e] == 11  )
         {
            SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, e, 1, graph) );
            break;
         }
      }

      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));
      assert(deletable);
   }

   /* clean up */
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE extTest1(
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
   int* vbase;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff;
   int edge;
   SCIP_Bool deletable;

   SCIP_CALL( SCIPallocBufferArray(scip, &rootdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &redcost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termpaths, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

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

   graph_mark(graph);

   ext0InitRedCostArrays(graph, redcost, rootdist, vbase, termpaths);

   cutoff = 0.0;
   edge = 0;

   {
      REDCOST redcostdata = {redcost, rootdist, termpaths, vbase, cutoff, root};
      SCIP_CALL(extArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));
   }
   assert(deletable);

   /* clean up */
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &termpaths);
   SCIPfreeBufferArray(scip, &redcost);
   SCIPfreeBufferArray(scip, &rootdist);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}


static
SCIP_RETCODE pseudoDelSingle(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 4;
   const int nedges = 3;
   SCIP_Real cutoff[] = {FARAWAY, FARAWAY, FARAWAY};
   SCIP_Bool success;

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_add(graph, 0);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 3, 1.0, 1.0);

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));


   for( int e = 0; e < nedges; e += 2 )
      assert(graph_edge_nPseudoAncestors(graph, e) == 0);

   SCIP_CALL( graph_knot_delPseudo(scip, graph, cutoff, NULL, NULL, 0, &success) );

   assert(success);

   for( int e = 0; e < nedges; e += 2 )
   {
      assert(graph_edge_nPseudoAncestors(graph, e) == 1);
      assert(graph_edge_getPseudoAncestors(graph, e)[0] == 0);
   }


   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}


static
SCIP_RETCODE pseudoDelDouble(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 4;
   SCIP_Real cutoff[] = {FARAWAY, FARAWAY, FARAWAY};
   SCIP_Bool success;

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_add(graph, 0);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 2, 1.0, 1.0);
   graph_edge_add(scip, graph, 0, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 4, 1.0, 1.0);

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   for( int e = 0; e < nedges; e += 2 )
      assert(graph_edge_nPseudoAncestors(graph, e) == 0);

   SCIP_CALL( graph_knot_delPseudo(scip, graph, cutoff, NULL, NULL, 0, &success) );
   assert(success);

   for( int e = graph->outbeg[1]; e != EAT_LAST; e = graph->oeat[e] )
   {
      if( graph->head[e] == 2 )
      {
         graph_edge_del(scip, graph, e, TRUE);

         break;
      }
   }

   SCIP_CALL( graph_knot_delPseudo(scip, graph, cutoff, NULL, NULL, 3, &success) );
   assert(success);

   assert(graph->grad[1] == 1);
   assert(graph->grad[2] == 1);
   assert(graph->grad[4] == 2);


   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}

static
SCIP_RETCODE pseudoAncestorsCreation(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;

   assert(scip);

   SCIP_CALL( graphBuildV5E5(scip, &graph, FALSE) );
   assert(graph->knots == 5 && graph->edges == 10);

   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 2, graph) );
   assert(graph_edge_nPseudoAncestors(graph, 0) == 1);
   assert(graph_edge_nPseudoAncestors(graph, 1) == 1);
   assert(graph_edge_nPseudoAncestors(graph, 2) == 0);

   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 3, graph) );
   assert(graph_edge_nPseudoAncestors(graph, 0) == 2);
   assert(graph_edge_nPseudoAncestors(graph, 1) == 2);

   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 1, 4, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 1, 1, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 1, 0, graph) );
   assert(graph_edge_nPseudoAncestors(graph, 0) == 5);
   assert(graph_edge_nPseudoAncestors(graph, 1) == 5);

   graph_edge_delPseudoAncestors(scip, 1, graph);
   assert(graph_edge_nPseudoAncestors(graph, 0) == 0);
   assert(graph_edge_nPseudoAncestors(graph, 1) == 0);

   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 2, graph) );
   assert(graph_edge_nPseudoAncestors(graph, 0) == 1);
   assert(graph_edge_nPseudoAncestors(graph, 1) == 1);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}


static
SCIP_RETCODE pseudoAncestorsMerge(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Bool conflict;

   assert(scip);

   SCIP_CALL( graphBuildV5E5(scip, &graph, FALSE) );
   assert(graph->knots == 5 && graph->edges == 10);

   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 2, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 4, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 4, 4, graph) );
   SCIP_CALL( graph_pseudoAncestors_appendMoveEdge(scip, 0, 4, FALSE, graph, &conflict)  );
   assert(conflict);
   assert( graph_edge_nPseudoAncestors(graph, 0) == 3 );
   assert( graph_edge_nPseudoAncestors(graph, 4) == 0);


   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 6, 1, graph) );
   SCIP_CALL( graph_pseudoAncestors_appendMoveEdge(scip, 0, 6, FALSE, graph, &conflict)  );
   assert(!conflict);
   assert( graph_edge_nPseudoAncestors(graph, 0) == 4 );
   assert( graph_edge_nPseudoAncestors(graph, 6) == 0);


   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}


static
SCIP_RETCODE pseudoAncestorsMergePc(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Bool conflict;

   assert(scip);

   SCIP_CALL( graphBuildV5E5(scip, &graph, TRUE) );
   assert(graph->knots > 5 && graph->edges > 10);

   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 0, 2, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 0, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 4, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 4, 4, graph) );
   SCIP_CALL( graph_pseudoAncestors_appendMoveNode(scip, 0, 4, FALSE, graph, &conflict)  );
   assert(conflict);
   assert( graph_knot_nPseudoAncestors(graph, 0) == 3 );
   assert( graph_knot_nPseudoAncestors(graph, 4) == 0);


   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 2, 1, graph) );
   SCIP_CALL( graph_pseudoAncestors_appendMoveNode(scip, 0, 2, FALSE, graph, &conflict)  );
   assert(!conflict);
   assert( graph_knot_nPseudoAncestors(graph, 0) == 4 );
   assert( graph_knot_nPseudoAncestors(graph, 2) == 0);


   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}


static
SCIP_RETCODE pseudoAncestorsHash(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Bool conflict;
   int* hasharr;

   assert(scip);

   SCIP_CALL( graphBuildV5E5(scip, &graph, FALSE) );
   assert(graph->knots == 5 && graph->edges == 10);

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, graph->knots) );

   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 2, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 4, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 4, 4, graph) );
   graph_pseudoAncestors_hashEdge(graph->pseudoancestors, 0, hasharr);
   graph_pseudoAncestors_hashEdgeDirty(graph->pseudoancestors, 4, TRUE, &conflict, hasharr);
   assert(conflict);


   SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 6, 1, graph) );
   graph_pseudoAncestors_hashEdgeDirty(graph->pseudoancestors, 6, TRUE, &conflict, hasharr);
   assert(!conflict);

   graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, 6, hasharr);
   graph_pseudoAncestors_unhashEdge(graph->pseudoancestors, 0, hasharr);

   for( int k = 0; k < graph->knots; k++ )
      assert(hasharr[k] == 0);

   SCIPfreeCleanBufferArray(scip, &hasharr);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}

static
SCIP_RETCODE pseudoAncestorsHashPc(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   SCIP_Bool conflict;
   int* hasharr;

   assert(scip);

   SCIP_CALL( graphBuildV5E5(scip, &graph, TRUE) );
   assert(graph->knots > 5 && graph->edges > 10);

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, graph->knots) );

   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 0, 2, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 0, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 4, 3, graph) );
   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 4, 4, graph) );
   graph_pseudoAncestors_hashNode(graph->pseudoancestors, 0, hasharr);
   graph_pseudoAncestors_hashNodeDirty(graph->pseudoancestors, 4, TRUE, &conflict, hasharr);
   assert(conflict);


   SCIP_CALL( graph_pseudoAncestors_addToNode(scip, 2, 1, graph) );
   graph_pseudoAncestors_hashNodeDirty(graph->pseudoancestors, 2, TRUE, &conflict, hasharr);
   assert(!conflict);

   graph_pseudoAncestors_unhashNode(graph->pseudoancestors, 2, hasharr);
   graph_pseudoAncestors_unhashNode(graph->pseudoancestors, 0, hasharr);

   for( int k = 0; k < graph->knots; k++ )
      assert(hasharr[k] == 0);

   SCIPfreeCleanBufferArray(scip, &hasharr);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   return SCIP_OKAY;
}


static
SCIP_RETCODE extDistTest(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;
   const int nclosenodes = 2;         /**< max. number of close nodes to each node */

   DISTDATA distdata;

   assert(scip);

   /* 1. build a test graph */

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = root;

   graph_edge_add(scip, graph, 0, 1, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 2, 0.4, 0.4);
   graph_edge_add(scip, graph, 1, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 1, 4, 1.0, 1.0);

   graph_edge_add(scip, graph, 2, 5, 0.5, 0.5);
   graph_edge_add(scip, graph, 3, 6, 1.6, 1.6);
   graph_edge_add(scip, graph, 3, 7, 1.6, 1.6);
   graph_edge_add(scip, graph, 4, 8, 1.5, 1.5);
   graph_edge_add(scip, graph, 4, 9, 1.5, 1.5);

   graph_edge_add(scip, graph, 5, 10, 1.0, 1.0);
   graph_edge_add(scip, graph, 10, 11, 1.0, 1.0);
   graph_edge_add(scip, graph, 11, 12, 1.0, 1.0);

   graph_mark(graph);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( reduce_distDataInit(scip, graph, nclosenodes, FALSE, &distdata) );

   /* 2. do the actual test */

#ifndef NDEBUG
   /* test close nodes */
   {
      const RANGE* node_range = distdata.closenodes_range;
      const int* node_idx = distdata.closenodes_indices;
#ifdef SCIP_DEBUG
      const SCIP_Real* node_dist = distdata.closenodes_distances;

      for( int node = 0; node < 6; node++ )
      {
         printf("node %d: \n", node);

         for( int i = node_range[node].start; i < node_range[node].end; i++ )
         {
            printf("...closenode=%d  dist=%f\n", node_idx[i], node_dist[i]);
         }
      }
#endif

      assert(node_idx[node_range[1].start] == 2);
      assert(node_idx[node_range[1].start + 1] == 5);
      assert(node_idx[node_range[5].start] == 1);
      assert(node_idx[node_range[5].start + 1] == 2);
   }

   /* test edge root paths */
   {
      const int edge = 2;
      PRSTATE** pathroot_blocks = distdata.pathroot_blocks;
      int* pathroot_blocksizes = distdata.pathroot_blocksizes;


#ifdef SCIP_DEBUG
      printf("edge ");
      graph_edge_printInfo(graph, edge);
#endif

      assert(graph->head[edge] == 2 && graph->tail[edge] == 1);
      assert(pathroot_blocksizes[edge / 2] == 6);

      for( int i = 0; i < pathroot_blocksizes[edge / 2]; i++ )
         SCIPdebugMessage("...root=%d  \n", pathroot_blocks[edge / 2][i].pathroot_id);
   }
#endif


   /* test distances */
#ifndef NDEBUG
   {
      const int edge = 2;
      const SCIP_Real dist1_2 = reduce_distDataGetSD(scip, graph, 1, 2, &distdata);
      const SCIP_Real dist1_3 = reduce_distDataGetSD(scip, graph, 1, 3, &distdata);
      const SCIP_Real dist1_5 = reduce_distDataGetSD(scip, graph, 1, 5, &distdata);
      const SCIP_Real dist2_1 = reduce_distDataGetSD(scip, graph, 2, 1, &distdata);
      const SCIP_Real dist2_5 = reduce_distDataGetSD(scip, graph, 2, 5, &distdata);
      const SCIP_Real dist2_6 = reduce_distDataGetSD(scip, graph, 2, 6, &distdata);

      assert(dist1_2 == 0.4);
      assert(dist1_3 == -1.0);
      assert(dist1_5 == 0.9);
      assert(dist2_1 == 0.4);
      assert(dist2_5 == 0.5);
      assert(dist2_6 == -1.0);

      assert(graph->head[edge] == 2 && graph->tail[edge] == 1);
      graph_edge_delFull(scip, graph, edge, TRUE);
      reduce_distDataDeleteEdge(scip, graph, edge, &distdata);

      {
         const SCIP_Real dist1_2_b = reduce_distDataGetSD(scip, graph, 1, 2, &distdata);
         const SCIP_Real dist2_6_b = reduce_distDataGetSD(scip, graph, 2, 6, &distdata);
         const SCIP_Real dist2_5_b = reduce_distDataGetSD(scip, graph, 2, 5, &distdata);

         assert(dist1_2_b == -1.0);
         assert(dist2_6_b == -1.0);
         assert(dist2_5_b == 0.5);
      }
   }
#endif

   reduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

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


/** test pseudo ancestors */
SCIP_RETCODE pseudoAncestors_test(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( pseudoAncestorsCreation(scip) );
   SCIP_CALL( pseudoAncestorsMerge(scip) );
   SCIP_CALL( pseudoAncestorsHash(scip) );
   SCIP_CALL( pseudoAncestorsMergePc(scip) );
   SCIP_CALL( pseudoAncestorsHashPc(scip) );

   printf("pseudoAncestors_test passed \n");

   return SCIP_OKAY;
}

/** test pseudo deletion */
SCIP_RETCODE pseudoDel_test(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( pseudoDelSingle(scip) );
   SCIP_CALL( pseudoDelDouble(scip) );

   printf("pseudoDeletion test passed \n");

   return SCIP_OKAY;
}


/** tests all */
SCIP_RETCODE testAll(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( pseudoDel_test(scip) );
   SCIP_CALL( reduce_extTest(scip) );
   SCIP_CALL( dheap_Test(scip) );
   SCIP_CALL( reduce_sdPcMwTest(scip) );
   SCIP_CALL( heur_extendPcMwOuterTest(scip) );

   return SCIP_OKAY;
}


/** tests extended reduction techniques */
SCIP_RETCODE reduce_extTest(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);

   SCIP_CALL( extTest5_variants(scip, 1) );
   SCIP_CALL( extTest5_variants(scip, 2) );

   SCIP_CALL( extTest4_variants(scip, 1) );
   SCIP_CALL( extTest4_variants(scip, 2) );
   SCIP_CALL( extTest4_variants(scip, 3) );

   SCIP_CALL( extTest3_variants(scip, 1) );
   SCIP_CALL( extTest3_variants(scip, 2) );
   SCIP_CALL( extTest3_variants(scip, 3) );
   SCIP_CALL( extTest3_variants(scip, 4) );
   SCIP_CALL( extTest3_variants(scip, 5) );

   SCIP_CALL( extTest2_variants(scip, 1) );

   SCIP_CALL( extTest1(scip) );

   SCIP_CALL( extDistTest(scip) );

   printf("reduce_extTest: all ok \n");

   return SCIP_OKAY;
}


/** tests PCMW special distance methods */
SCIP_RETCODE reduce_sdPcMwTest(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( sdPcMwTest1(scip) );
   SCIP_CALL( sdPcMwTest2(scip) );
   SCIP_CALL( sdPcMwTest4(scip) );

   printf("reduce_sdPcMwTest: all ok \n");

   return SCIP_OKAY;
}


SCIP_RETCODE dheap_Test(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DHEAP* heap = NULL;

   int min = -1;
   graph_heap_create(scip, 13, NULL, NULL, &heap);
   assert(heap != NULL);

   graph_heap_correct(1, 2.0, heap);
   graph_heap_correct(2, 1.0, heap);
   graph_heap_correct(0, 1.5, heap);

   assert(heap->size == 3);

   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 2);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 0);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 1);

   assert(heap->size == 0);
   graph_heap_clean(TRUE, heap);


   graph_heap_correct(1, 2.0, heap);
   graph_heap_correct(2, 2.7, heap);
   graph_heap_correct(0, 1.5, heap);
   graph_heap_correct(2, 1.9, heap);
   graph_heap_correct(4, 0.5, heap);

   assert(heap->size == 4);

   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 4);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 0);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 2);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 1);

   assert(heap->size == 0);
   graph_heap_clean(TRUE, heap);

   graph_heap_correct(1, 2.0, heap);
   graph_heap_correct(2, 3.0, heap);
   graph_heap_correct(0, 1.5, heap);
   graph_heap_correct(2, 1.6, heap);
   graph_heap_correct(12, 22.5, heap);
   graph_heap_correct(12, 7.7, heap);
   graph_heap_correct(4, 8.5, heap);


   assert(heap->size == 5);

   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 0);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 2);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 1);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 12);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 4);


   assert(heap->size == 0);


   graph_heap_free(scip, TRUE, TRUE, &heap);
   assert(heap == NULL);

   graph_heap_create(scip, 3, NULL, NULL, &heap);

   graph_heap_correct(1, 2.0, heap);
   graph_heap_correct(2, 3.0, heap);
   graph_heap_correct(0, 1.5, heap);
   graph_heap_correct(2, 2.5, heap);

   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 0);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 1);
   graph_heap_deleteMinGetNode(&min, heap);
   assert(min == 2);

   assert(heap->size == 0);

   graph_heap_free(scip, TRUE, TRUE, &heap);
   assert(heap == NULL);

   return SCIP_OKAY;
}


SCIP_RETCODE heur_extendPcMwOuterTest(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   int* stedge;
   STP_Bool* stvertex;
   int nnodes = 5;
   int nedges = 5;

   assert(scip);

   SCIP_CALL(graph_init(scip, &graph, nnodes, 2 * nedges, 1));

   for( int i = 0; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_add(scip, graph, 0, 1, 2.0, 2.0);
   graph_edge_add(scip, graph, 1, 2, 3.0, 3.0);
   graph_edge_add(scip, graph, 1, 4, 2.0, 2.0);
   graph_edge_add(scip, graph, 2, 3, 1.0, 1.0);
   graph_edge_add(scip, graph, 3, 4, 1.0, 1.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = 0.0;

   graph->prize[0] = 0.1;
   graph->prize[2] = 3.1;
   graph->prize[4] = 3.0;

   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 2, 0);
   graph_knot_chg(graph, 0, 0);

   SCIP_CALL(graph_pc_2pc(scip, graph));

   nnodes = graph->knots;
   nedges = graph->edges;

   SCIP_CALL(graph_init_history(scip, graph));
   SCIP_CALL(graph_path_init(scip, graph));

   graph_mark(graph);

   SCIP_CALL(SCIPallocBufferArray(scip, &stvertex, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &stedge, nedges));

   for( int i = 0; i < nedges; i++ )
      stedge[i] = UNKNOWN;

   for( int i = 0; i < nnodes; i++ )
      stvertex[i] = FALSE;

   stvertex[0] = TRUE;

   SCIP_CALL( SCIPStpHeurTMPrunePc(scip, graph, graph->cost, stedge, stvertex) );

   /* actual test */
   SCIP_CALL( SCIPStpHeurLocalExtendPcMwOut(scip, graph, stedge, stvertex) );

   assert(stvertex[0] && stvertex[1] && stvertex[2] && stvertex[3]);

   /* clean up */

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);

   SCIPfreeBufferArray(scip, &stedge);
   SCIPfreeBufferArray(scip, &stvertex);

   return SCIP_OKAY;
}
