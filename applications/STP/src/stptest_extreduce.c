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
#include "extreduce.h"




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
   DISTDATA distdata;
   EXTPERMA extpermanent;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, graph, edgedeleted, &extpermanent) );

   if( edgedelete >= 0 )
   {
      graph_edge_delFull(scip, graph, edgedelete, TRUE);
      extreduce_distDataDeleteEdge(scip, graph, edgedelete, &distdata);

      graph_mark(graph);
   }

   /* actual test */
   SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, edge, equality, &distdata, &extpermanent, deletable) );

   /* clean up */
   extreduce_extPermaFreeMembers(scip, &extpermanent);
   extreduce_distDataFreeMembers(scip, graph, &distdata);
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
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, &distdata) );

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
      const SCIP_Real dist1_2 = extreduce_distDataGetSD(scip, graph, 1, 2, &distdata);
      const SCIP_Real dist1_3 = extreduce_distDataGetSD(scip, graph, 1, 3, &distdata);
      const SCIP_Real dist1_5 = extreduce_distDataGetSD(scip, graph, 1, 5, &distdata);
      const SCIP_Real dist2_1 = extreduce_distDataGetSD(scip, graph, 2, 1, &distdata);
      const SCIP_Real dist2_5 = extreduce_distDataGetSD(scip, graph, 2, 5, &distdata);
      const SCIP_Real dist2_6 = extreduce_distDataGetSD(scip, graph, 2, 6, &distdata);

      assert(dist1_2 == 0.4);
      assert(dist1_3 == -1.0);
      assert(dist1_5 == 0.9);
      assert(dist2_1 == 0.4);
      assert(dist2_5 == 0.5);
      assert(dist2_6 == -1.0);

      assert(graph->head[edge] == 2 && graph->tail[edge] == 1);
      graph_edge_delFull(scip, graph, edge, TRUE);
      extreduce_distDataDeleteEdge(scip, graph, edge, &distdata);

      {
         const SCIP_Real dist1_2_b = extreduce_distDataGetSD(scip, graph, 1, 2, &distdata);
         const SCIP_Real dist2_6_b = extreduce_distDataGetSD(scip, graph, 2, 6, &distdata);
         const SCIP_Real dist2_5_b = extreduce_distDataGetSD(scip, graph, 2, 5, &distdata);

         assert(dist1_2_b == -1.0);
         assert(dist2_6_b == -1.0);
         assert(dist2_5_b == 0.5);
      }
   }
#endif

   extreduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);

   return SCIP_OKAY;
}




/** tests extended reduction techniques */
SCIP_RETCODE stptest_extreduce(
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

   printf("extreduce test: all ok \n");

   return SCIP_OKAY;
}
