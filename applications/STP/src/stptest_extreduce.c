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

/**@file   stptest_extreduce.c
 * @brief  tests for Steiner tree extended reductions
 * @author Daniel Rehfeldt
 *
 * This file implements tests for Steiner problem extended reductions.
 *
 * A list of all interface methods can be found in stptest.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include <stdio.h>
#include <assert.h>
#include "stptest.h"
#include "portab.h"
#include "graph.h"
#include "reduce.h"
#include "extreduce.h"


#define STPTEST_EXT_MAXNCLOSENODES 64

/** base method for extended edge reduction tests */
static
SCIP_RETCODE extCheckArc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata,        /**< reduced cost data */
   STP_Bool*             edgedeleted,        /**< marks for each edge whether deleted */
   int                   edge,               /**< the edge/arc to check */
   int                   edgedelete,         /**< delete the edge if possible? */
   int                   nclosenodes,        /**< max. number of close nodes to each node */
   SCIP_Bool*            deletable,          /**< is the edge deletable? */
   SCIP_Bool             equality            /**< allow equality for checks */
)
{
   DISTDATA* distdata;
   EXTPERMA* extpermanent;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, nclosenodes, FALSE, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_full, graph, edgedeleted, &extpermanent) );

   if( edgedelete >= 0 )
   {
      extreduce_edgeRemove(scip, edgedelete, graph, distdata, extpermanent);
      assert(graph_isMarked(graph));
   }

   extpermanent->redcostEqualAllow = equality;
   extpermanent->distdata_default = distdata;

   /* actual test */
   SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, edge, extpermanent, deletable) );

   /* clean up */
   extreduce_extPermaFree(scip, &extpermanent);
   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}


/** base method for extended edge reduction tests */
static
SCIP_RETCODE extCheckEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata,        /**< reduced cost data */
   STP_Bool*             edgedeleted,        /**< array to mark */
   int                   edge,               /**< edge to check */
   SCIP_Bool*            deletable,          /**< can edge be deleted?  */
   SCIP_Bool             allowEquality       /**< allow equality for rule-out */
)
{
   DISTDATA* distdata;
   EXTPERMA* extpermanent;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STPTEST_EXT_MAXNCLOSENODES, FALSE, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_full, graph, edgedeleted, &extpermanent) );

   extpermanent->redcostEqualAllow = allowEquality;
   extpermanent->distdata_default = distdata;

   /* actual test */
   SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, edge, extpermanent, deletable) );

   /* clean up */
   extreduce_extPermaFree(scip, &extpermanent);
   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}


/** base method for extended edge reduction tests */
static
SCIP_RETCODE extCheckNode(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata,        /**< reduced cost data */
   STP_Bool*             edgedeleted,        /**< array */
   int                   node,               /**< node to check */
   SCIP_Bool*            deletable,          /**< pointer to mark whether node can be deleted */
   SCIP_Bool             allowEquality       /**< allow equality? */
)
{
   STAR* star;
   DISTDATA* distdata;
   EXTPERMA* extpermanent;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STPTEST_EXT_MAXNCLOSENODES, FALSE, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_full, graph, edgedeleted, &extpermanent) );
   SCIP_CALL( reduce_starInit(scip, graph->grad[node], &star) );

   extpermanent->redcostEqualAllow = allowEquality;
   extpermanent->distdata_default = distdata;

   /* actual test */
   SCIP_CALL( extreduce_checkNode(scip, graph, redcostdata, node, star, extpermanent, deletable) );

   /* clean up */
   reduce_starFree(scip, &star);
   extreduce_extPermaFree(scip, &extpermanent);
   extreduce_distDataFree(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}



/** base method for extended edge reduction tests (with biased SDs) */
static
SCIP_RETCODE extDeleteNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata,        /**< reduced cost data */
   SCIP_Bool             allowEquality       /**< allow equality? */
)
{
   DISTDATA* distdata;
   EXTPERMA* extpermanent;
   int ndeleted;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STPTEST_EXT_MAXNCLOSENODES, FALSE, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, extred_full, graph, NULL, &extpermanent) );

   assert(!extpermanent->distdata_default);

   extpermanent->distdata_default = distdata;
   extpermanent->redcostdata = redcostdata;
   extpermanent->redcostEqualAllow = allowEquality;

   SCIP_CALL( extreduce_pseudoDeleteNodes(scip, NULL, extpermanent, graph, NULL, &ndeleted) );

   /* clean up */
   extreduce_distDataFree(scip, graph, &distdata);

   if( extpermanent->distdata_biased )
      extreduce_distDataFree(scip, graph, &(extpermanent->distdata_biased));

   extreduce_extPermaFree(scip, &extpermanent);

   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}

/** initializes to default */
static
void extInitRedCostArrays(
   const GRAPH*          graph,              /**< the graph */
   REDCOST*              redcostdata         /**< reduced costs data */
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   SCIP_Real* const rootdist = redcosts_getRootToNodeDistTop(redcostdata);
   SCIP_Real* const redcost = redcosts_getEdgeCostsTop(redcostdata);
   PATH* const termpaths3 = redcosts_getNodeToTermsPathsTop(redcostdata);
   int* const vbase3 = redcosts_getNodeToTermsBasesTop(redcostdata);

   assert(redcosts_getNnodes(redcostdata) == nnodes);
   assert(redcosts_getNedges(redcostdata) == nedges);

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


/** initializes to default for PC */
static
void extInitRedCostArraysPc(
   const GRAPH*          graph,              /**< the graph */
   REDCOST*              redcostdata         /**< reduced costs data */
)
{
   const int nnodes = graph->knots;
   int* const vbase3 = redcosts_getNodeToTermsBasesTop(redcostdata);

   assert(graph_pc_isPc(graph));

   extInitRedCostArrays(graph, redcostdata);

   for( int i = 0; i < 3 * nnodes; i++ )
   {
      vbase3[i] = 0;
   }
}

#ifdef SCIP_DISABLED
/** initializes to default for PC */
static
void extInitRedCostArraysPcWithBase(
   const GRAPH*          graph,              /**< the graph */
   int                   base,               /**< the base */
   REDCOST*              redcostdata         /**< reduced costs data */
)
{
   const int nnodes = graph->knots;
   int* const vbase3 = redcosts_getNodeToTermsBasesTop(redcostdata);

   assert(graph_pc_isPc(graph));

   extInitRedCostArrays(graph, redcostdata);

   for( int i = 0; i < 3 * nnodes; i++ )
   {
      vbase3[i] = base;
   }
}
#endif


static
SCIP_RETCODE testEdgeDeletion5_deprecated(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2 */
)
{
   GRAPH* graph;
   const int nnodes = 85;
   const int nedges = 88;
   const int root = 0;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 100.0;
   int testedge;
   SCIP_Bool deletable;
   REDCOST* redcostdata;

   assert(variant == 1 || variant == 2);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, STP_TERM);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph_knot_chg(graph, 6, STP_TERM);
   graph_knot_chg(graph, 7, STP_TERM);
   graph_knot_chg(graph, 8, STP_TERM);
   graph_knot_chg(graph, 9, STP_TERM);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);

   graph_edge_addBi(scip, graph, 2, 6, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 8, 1.0);
   graph_edge_addBi(scip, graph, 5, 9, 1.0);

   graph_edge_addBi(scip, graph, 6, 10, 1.0);
   graph_edge_addBi(scip, graph, 7, 10, 1.0);
   graph_edge_addBi(scip, graph, 8, 10, 1.0);
   graph_edge_addBi(scip, graph, 9, 10, 1.0);

   graph_edge_addBi(scip, graph, 0, 10, 0.9);

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   for( int i = 0; i < graph->edges; i++ )
      redcosts_getEdgeCostsTop(redcostdata)[i] = 1.0;

   testedge = 0;

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   if( variant == 1 )
   {
      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, testedge, -1, nnodes, &deletable, TRUE));
      assert(deletable);
   }
   else
   {
      int edgedelete = -1;

      assert(variant == 2);

      for( int e = graph->outbeg[0]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( graph->head[e] == 10 )
            edgedelete = e;
      }

      assert(edgedelete >= 0);

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, testedge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


static
SCIP_RETCODE testEdgeDeletion4_deprecated(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2,3 */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 100.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST* redcostdata;

   assert(scip);
   assert(variant == 1 || variant == 2 || variant == 3);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 7, 1.0);

   graph_edge_addBi(scip, graph, 0, 6, 0.9);
   graph_edge_addBi(scip, graph, 0, 7, 0.9);

   if( variant == 3 )
      graph_edge_addBi(scip, graph, 7, 8, 0.9);

   graph_mark(graph);

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   for( int i = 0; i < graph->edges; i++ )
      redcosts_getEdgeCostsTop(redcostdata)[i] = 1.0;

   edge = 0;

   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 5, 0);

   SCIP_CALL( graph_initHistory(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   if( variant == 1 )
   {
      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));

      assert(deletable);
   }
   else if( variant == 2 )
   {
      int edgedelete = -1;

      for( int e = graph->outbeg[0]; e != EAT_LAST; e = graph->oeat[e] )
         if( graph->head[e] == 7 )
            edgedelete = e;

      assert(edgedelete >= 0);

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }
   else
   {
      int edgedelete = -1;

      assert(variant == 3);

      for( int e = graph->outbeg[7]; e != EAT_LAST; e = graph->oeat[e] )
         if( graph->head[e] == 8 )
            edgedelete = e;

      assert(edgedelete >= 0);

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(deletable);
   }

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


static
SCIP_RETCODE testEdgeDeletion3_deprecated(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2,3,4 */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 7.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST* redcostdata;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build graph */
   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 6, 1.0);
   graph_edge_addBi(scip, graph, 4, 7, 1.0);

   graph_edge_addBi(scip, graph, 7, 9, 0.9);
   graph_edge_addBi(scip, graph, 8, 0, 1.0);
   graph_edge_addBi(scip, graph, 9, 0, 1.0);
   graph_edge_addBi(scip, graph, 0, nnodes - 1, 1.0);

   graph_knot_chg(graph, nnodes - 1, 0); /* dummy node */
   graph_mark(graph);

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   redcosts_getNodeToTermsBasesTop(redcostdata)[5] = nnodes - 1;
   redcosts_getNodeToTermsBasesTop(redcostdata)[6] = nnodes - 1;
   redcosts_getNodeToTermsPathsTop(redcostdata)[5].dist = 3.0;
   redcosts_getNodeToTermsPathsTop(redcostdata)[6].dist = 3.0;
   redcosts_getNodeToTermsPathsTop(redcostdata)[5 + nnodes].dist = 3.0;
   redcosts_getNodeToTermsPathsTop(redcostdata)[6 + nnodes].dist = 3.0;

   for( int i = 0; i < graph->edges; i++ )
      redcosts_getEdgeCostsTop(redcostdata)[i] = 1.0;

   edge = 0;

   graph_knot_chg(graph, 7, 0);

   if( variant == 1 )
   {
      SCIP_CALL( graph_initHistory(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(deletable);
   }
   else if( variant == 2 )
   {
      SCIP_CALL( graph_initHistory(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));
      assert(deletable);
   }
   else if( variant == 3 )
   {
      const int edge2 = 14;
      assert(graph->tail[edge2] == 7 && graph->head[edge2] == 9);

      SCIP_CALL( graph_initHistory(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      graph->cost[edge2] = 0.99;
      graph->cost[flipedge(edge2)] = 0.99;
      graph_knot_chg(graph, 4, 0);

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
   }
   else if( variant == 4 )
   {
      const int edgedelete = 14;
      assert(graph->tail[edgedelete] == 7 && graph->head[edgedelete] == 9);

      SCIP_CALL( graph_initHistory(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


static
SCIP_RETCODE testEdgeDeletion2_deprecated(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   variant             /**< 1,2,3,4,5 */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 28;
   const int root = 0;

   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 100.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST* redcostdata;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);

   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 6, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 8, 1.0);
   graph_edge_addBi(scip, graph, 4, 9, 1.0);

   graph_edge_addBi(scip, graph, 5, 10, 1.0);
   graph_edge_addBi(scip, graph, 10, 11, 1.0);
   graph_edge_addBi(scip, graph, 11, 12, 1.0);

   graph_mark(graph);

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   edge = 0;

   if( variant == 1 )
   {
      int pseudoancestor;
      SCIP_CALL( graph_initHistory(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      graph_addPseudoAncestor(graph, &pseudoancestor);
      graph_addPseudoAncestor(graph, &pseudoancestor);
      SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, 0, 1, graph) );

      for( int e = graph->outbeg[10]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( graph->head[e] == 11  )
         {
            SCIP_CALL( graph_pseudoAncestors_addToEdge(scip, e, 1, graph) );
            break;
         }
      }

      SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));

      STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");
   }

   stptest_extreduceTearDown(scip, graph, &redcostdata);


   return SCIP_OKAY;
}


static
SCIP_RETCODE testEdgeDeletion1_deprecated(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 55;
   const int nedges = 18;
   const int root = 0;
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 0.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST* redcostdata;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   graph_knot_add(graph, 0);

   /* also add dummy nodes to avoid full stack */
   for( int i = 1; i < nnodes; i++ )
      graph_knot_add(graph, -1);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);

   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 6, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 8, 1.0);
   graph_edge_addBi(scip, graph, 4, 9, 1.0);

   SCIP_CALL( graph_initHistory(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   graph_mark(graph);

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   edge = 0;

   SCIP_CALL(extCheckArc(scip, graph, redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that edge can be deleted by exploiting that two vertices have the
 *  same reduced costs closest terminal */
static
SCIP_RETCODE testEdgeDeletedByCommonRedCostsTargets(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 15;
   const int nedges = 64; /* just some upper bound */
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;
   const int firstdummynode = 6;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */

   /* add dummy nodes 6-14 */
   for( int i = firstdummynode; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.9);
   graph_edge_addBi(scip, graph, 0, 4, 1.9);

   /* also add dummy edges to avoid extension */
   for( int i = firstdummynode; i < nnodes; i++ )
   {
      graph_edge_addBi(scip, graph, 3, i, 1.0);
      graph_edge_addBi(scip, graph, 4, i, 1.0);
   }

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   /* put 5 as first target for both nodes 3 and 4 and set the distance to the second target too high */
   redcosts_getNodeToTermsBasesTop(redcostdata)[3] = 5;
   redcosts_getNodeToTermsBasesTop(redcostdata)[4] = 5;
   redcosts_getNodeToTermsPathsTop(redcostdata)[3 + nnodes].dist = cutoff + 0.1;
   redcosts_getNodeToTermsPathsTop(redcostdata)[4 + nnodes].dist = cutoff + 0.1;

   SCIP_CALL(extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}




/** tests that edge can be deleted by exploiting multi-level reduced costs */
static
SCIP_RETCODE testEdgeDeletedByMultiRedCosts(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 15;
   const int nedges = 80; /* just some upper bound */
   SCIP_Real* redcosts = NULL;
   int testedge = 0;
   SCIP_Bool deletable;
   const int firstdummynode = 7;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */

   for( int i = firstdummynode; i < nnodes; i++ )
      graph_knot_add(graph, STP_TERM_NONE);

   graph->source = 6;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); /* 0 */
   graph_edge_addBi(scip, graph, 1, 2, 1.0); /* 2 */
   graph_edge_addBi(scip, graph, 1, 3, 1.0); /* 4 */
   graph_edge_addBi(scip, graph, 2, 4, 2.0); /* 6 */
   graph_edge_addBi(scip, graph, 3, 5, 2.0); /* 8 */

   graph_edge_addBi(scip, graph, 0, 4, 1.0); /* 10 */
   graph_edge_addBi(scip, graph, 0, 6, 1.0); /* 12 */


   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 2, 1.9);
   graph_edge_addBi(scip, graph, 0, 3, 1.9);

   /* also add dummy edges to avoid extension */
   for( int i = firstdummynode; i < nnodes; i++ )
   {
      graph_edge_addBi(scip, graph, 0, i, 3.0);
      graph_edge_addBi(scip, graph, 2, i, 3.0);
      graph_edge_addBi(scip, graph, 3, i, 3.0);
   }

   SCIP_CALL( stptest_graphSetUp(scip, graph) );

   {
       RCPARAMS rcparams = { .cutoff = 7.0, .nLevels = 2, .nCloseTerms = 3,
                             .nnodes = graph->knots, .nedges = graph->edges, .redCostRoot = graph->source };
       SCIP_CALL( redcosts_initFromParams(scip, &rcparams, &redcostdata) );
   }

   /* 1st level */
   redcosts = redcosts_getEdgeCostsTop(redcostdata);
   for( int i = 0; i < graph->edges; i++ )
      redcosts[i] = graph->cost[i];

   redcosts[0] = 0;

   SCIP_CALL( redcosts_initializeDistancesTop(scip, graph, redcostdata) );

   /* 2nd level */
   redcosts_addLevel(redcostdata);
   redcosts = redcosts_getEdgeCostsTop(redcostdata);
   for( int i = 0; i < graph->edges; i++ )
      redcosts[i] = graph->cost[i];

   redcosts_setCutoffTop(6.5, redcostdata);
   redcosts_setRootTop(4, redcostdata);
   SCIP_CALL( redcosts_initializeDistancesTop(scip, graph, redcostdata) );

   /* call actual reduction method */
   SCIP_CALL(extCheckEdge(scip, graph, redcostdata, NULL, testedge, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}

/** tests that edge can be deleted by using SD MST argument */
static
SCIP_RETCODE testEdgeDeletedByMst1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 16;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 2, 1.4);
   graph_edge_addBi(scip, graph, 0, 3, 1.4);
   graph_edge_addBi(scip, graph, 0, 4, 1.4);
   graph_edge_addBi(scip, graph, 2, 3, 1.1);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL(extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge can be deleted by using SD MST argument */
static
SCIP_RETCODE testEdgeDeletedByMst2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 8;
   const int nedges = 22;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */
   graph_knot_add(graph, STP_TERM);       /* node 7 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.2);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 6, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 3, 1.5);
   graph_edge_addBi(scip, graph, 0, 4, 2.5);
   graph_edge_addBi(scip, graph, 0, 5, 2.5);
   graph_edge_addBi(scip, graph, 4, 5, 1.1);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "edge was deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that edge can be deleted by using equality bottlneck test*/
static
SCIP_RETCODE testEdgeDeletedByEqBottleneck(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 4;
   const int nedges = 12;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 2, 1.9);
   graph_edge_addBi(scip, graph, 0, 3, 1.9);
   graph_edge_addBi(scip, graph, 2, 3, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL(extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge cannot be deleted */
static
SCIP_RETCODE testEdgeNotDeleted1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 8;
   const int nedges = 22;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( redcosts_init(scip, nnodes, nedges, cutoff, root, &redcostdata) );
   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */
   graph_knot_add(graph, STP_TERM);       /* node 7 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 6, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 3, 1.5);
   graph_edge_addBi(scip, graph, 0, 4, 2.5);
   graph_edge_addBi(scip, graph, 0, 5, 2.5);
   graph_edge_addBi(scip, graph, 4, 5, 1.1);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(!deletable, "edge was deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge can be deleted (with equality argument) */
static
SCIP_RETCODE testEdgeDeletedByEqBottleneck2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 12;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);

   graph_edge_addBi(scip, graph, 0, 4, 3.0);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that edge can be deleted by extended SPG rule-out */
static
SCIP_RETCODE testEdgeDeletedBy3LeafSpg(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 12;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM);  /* node 2 */
   graph_knot_add(graph, STP_TERM);  /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.1);
   graph_edge_addBi(scip, graph, 1, 2, 1.1);
   graph_edge_addBi(scip, graph, 1, 3, 1.1);

   graph_edge_addBi(scip, graph, 0, 4, 1.1);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge can be deleted by reduced costs argument */
static
SCIP_RETCODE testNode3PseudoDeletedByRedCosts1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 6;
   const int nedges = 14;
   const int root = 0;
   SCIP_Real cutoff = 5.9;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;
   SCIP_Real* redEdgeCosts;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0); /* 0 */
   graph_edge_addBi(scip, graph, 0, 2, 1.0); /* 2 */
   graph_edge_addBi(scip, graph, 0, 3, 1.0); /* 4 */
   graph_edge_addBi(scip, graph, 1, 4, 1.0); /* 6 */
   graph_edge_addBi(scip, graph, 2, 4, 1.0); /* 8 */
   graph_edge_addBi(scip, graph, 2, 5, 1.0); /* 10 */
   graph_edge_addBi(scip, graph, 3, 5, 2.2); /* 12 */

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   redEdgeCosts = redcosts_getEdgeCostsTop(redcostdata);

   redEdgeCosts[2] = 3.0;
   redEdgeCosts[3] = 3.0;
   redEdgeCosts[10] = 3.0;
   redEdgeCosts[11] = 3.0;

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that node can be deleted */
static
SCIP_RETCODE testNode3PseudoDeletedByContraction(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 9;
   const int nedges = 22;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 1 */
   graph_knot_add(graph, STP_TERM);        /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 5 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 6 */
   graph_knot_add(graph, STP_TERM);       /* node 7 */
   graph_knot_add(graph, STP_TERM);       /* node 8 */


   graph->source = 7;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 0.5);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 5, 7, 1.0);
   graph_edge_addBi(scip, graph, 6, 8, 1.0);

   graph_edge_addBi(scip, graph, 1, 2, 1.05);
   graph_edge_addBi(scip, graph, 3, 5, 1.9);
   graph_edge_addBi(scip, graph, 3, 6, 1.9);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that node can be pseudo-deleted */
static
SCIP_RETCODE testNode3PseudoDeletedBySdBiasedSimple(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   DISTDATA* distdata_default;
   DISTDATA* distdata_biased;
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 7;
   const int nedges = 16;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   int testnode = 0;
   SCIP_Bool deletable = FALSE;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */


   graph->source = 5;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 2.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);

   /* dummy */
   graph_edge_addBi(scip, graph, 5, 6, 12.0);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STPTEST_EXT_MAXNCLOSENODES, FALSE, FALSE, &distdata_default) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STPTEST_EXT_MAXNCLOSENODES, TRUE, TRUE, &distdata_biased) );

   SCIP_CALL(  extreduce_mstbiasedCheck3NodeSimple(scip, graph, testnode, distdata_default, distdata_biased, &deletable) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   extreduce_distDataFree(scip, graph, &distdata_biased);
   extreduce_distDataFree(scip, graph, &distdata_default);
   graph_free_dcsr(scip, graph);
   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that node can be pseudo-deleted */
static
SCIP_RETCODE testNode3PseudoDeletedBySdBiased(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 10;
   const int nedges = 22;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   int testnode = 0;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 7 */
   graph_knot_add(graph, STP_TERM);       /* node 8 */
   graph_knot_add(graph, STP_TERM);       /* node 9 */

   graph->source = 5;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 7, 1.0);
   graph_edge_addBi(scip, graph, 7, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.5);
   graph_edge_addBi(scip, graph, 3, 4, 1.5);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);

   /* dummy */
   graph_edge_addBi(scip, graph, 5, 6, 12.0);
   graph_edge_addBi(scip, graph, 1, 8, 1.0);
   graph_edge_addBi(scip, graph, 3, 9, 1.0);
   graph_edge_addBi(scip, graph, 2, 6, 2.1);


   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extDeleteNodes(scip, graph, redcostdata, TRUE) );

   STPTEST_ASSERT_MSG(graph->grad[testnode] == 0, "node was not deleted! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that node can be deleted */
static
SCIP_RETCODE testNode3PseudoDeletedBySd1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 6;
   const int nedges = 14;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that node can be deleted */
static
SCIP_RETCODE testNode3PseudoDeletedBySd2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 6;
   const int nedges = 14;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */

   graph->source = 3;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.5);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}

/** tests that node can be deleted */
static
SCIP_RETCODE testNode3PseudoDeletedBySd3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 8;
   const int nedges = 20;
   const int root = 7;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */
   graph_knot_add(graph, STP_TERM);       /* node 7 */

   graph->source = 7;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.5);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 6, 2.0);
   graph_edge_addBi(scip, graph, 2, 5, 2.0);
   graph_edge_addBi(scip, graph, 2, 6, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 5, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that node can be deleted */
static
SCIP_RETCODE testNode4PseudoDeletedBySd1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 7;
   const int nedges = 18;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.5);
   graph_edge_addBi(scip, graph, 0, 3, 2.0);
   graph_edge_addBi(scip, graph, 0, 4, 2.0);
   graph_edge_addBi(scip, graph, 2, 5, 2.0);
   graph_edge_addBi(scip, graph, 2, 6, 2.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 6, 3, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}




/** tests that node cannot be deleted */
static
SCIP_RETCODE testNode4PseudoNotDeletedBySd1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 7;
   const int nedges = 16;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.5);
   graph_edge_addBi(scip, graph, 0, 3, 2.0);
   graph_edge_addBi(scip, graph, 0, 4, 2.0);
   graph_edge_addBi(scip, graph, 2, 5, 2.0);
   graph_edge_addBi(scip, graph, 2, 6, 2.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 6, 3, 1.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(!deletable, "node was marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge can be deleted by general star test */
static
SCIP_RETCODE testGeneralStarDeletedEdge1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   int nelims = 0;
   const int nnodes = 8;
   const int nedges = 20;
   const int root = 0;
   SCIP_Real cutoff = 100.0;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 1 */
   graph_knot_add(graph, STP_TERM);            /* node 2 */
   graph_knot_add(graph, STP_TERM);            /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */
   graph_knot_add(graph, STP_TERM);            /* node 5 */
   graph_knot_add(graph, STP_TERM);            /* node 6 */
   graph_knot_add(graph, STP_TERM_NONE);            /* node 7 */

   graph->source = 2;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);

   graph_edge_addBi(scip, graph, 4, 6, 1.1);
   graph_edge_addBi(scip, graph, 4, 7, 1.0);

   /* short cuts */
   graph_edge_addBi(scip, graph, 2, 3, 2.0);
   graph_edge_addBi(scip, graph, 3, 6, 2.0);
   graph_edge_addBi(scip, graph, 2, 5, 2.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   /* actual test */
   SCIP_CALL( extreduce_deleteGeneralStars(scip, redcostdata, NULL, graph, NULL, &nelims) );

   STPTEST_ASSERT_MSG(graph_edge_isDeleted(graph, 0), "edge 0 was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge can be deleted by general star test */
static
SCIP_RETCODE testGeneralStarDeletedEdge2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   int nelims = 0;
   const int nnodes = 6;
   const int nedges = 20;
   const int root = 0;
   SCIP_Real cutoff = 100.0;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 1 */
   graph_knot_add(graph, STP_TERM);            /* node 2 */
   graph_knot_add(graph, STP_TERM);            /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 5 */

   graph->source = 2;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);

   /* redundant */
   graph_edge_addBi(scip, graph, 0, 5, 1.0);

   /* short cuts */
   graph_edge_addBi(scip, graph, 2, 3, 3.0);
   graph_edge_addBi(scip, graph, 3, 4, 3.0);
   graph_edge_addBi(scip, graph, 4, 5, 3.0);
   graph_edge_addBi(scip, graph, 5, 2, 3.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArrays(graph, redcostdata);

   /* actual test */
   SCIP_CALL( extreduce_deleteGeneralStars(scip, redcostdata, NULL, graph, NULL, &nelims) );

   STPTEST_ASSERT_MSG(graph_edge_isDeleted(graph, 0), "edge 0 was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that edge can be deleted by general star test */
static
SCIP_RETCODE testGeneralStarDeletedEdge3(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   int nelims = 0;
   const int nnodes = 10;
   const int nedges = 26;
   const int root = 0;
   SCIP_Real cutoff = 100.0;

   assert(scip);

   SCIP_CALL( redcosts_init(scip, nnodes, nedges, cutoff, root, &redcostdata) );
   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 0 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 1 */
   graph_knot_add(graph, STP_TERM);            /* node 2 */
   graph_knot_add(graph, STP_TERM);            /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 4 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 5 */
   graph_knot_add(graph, STP_TERM);            /* node 6 */
   graph_knot_add(graph, STP_TERM);            /* node 7 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 8 */
   graph_knot_add(graph, STP_TERM);       /* node 9 */

   graph->source = 2;

   graph_edge_addBi(scip, graph, 0, 1, 2.0);
   graph_edge_addBi(scip, graph, 0, 2, 2.0);
   graph_edge_addBi(scip, graph, 0, 3, 2.0);
   graph_edge_addBi(scip, graph, 1, 4, 2.0);
   graph_edge_addBi(scip, graph, 1, 5, 3.0);

   graph_edge_addBi(scip, graph, 4, 7, 2.0);
   graph_edge_addBi(scip, graph, 4, 8, 2.0);
   graph_edge_addBi(scip, graph, 8, 9, 3.0);
   graph_edge_addBi(scip, graph, 5, 6, 10.0);

   /* short cuts */
   graph_edge_addBi(scip, graph, 2, 3, 2.0);
   graph_edge_addBi(scip, graph, 2, 5, 7.0);
   graph_edge_addBi(scip, graph, 3, 7, 2.0);
   graph_edge_addBi(scip, graph, 5, 8, 3.0);

   SCIP_CALL( stptest_graphSetUp(scip, graph) );
   extInitRedCostArrays(graph, redcostdata);

   /* actual test */
   SCIP_CALL( extreduce_deleteGeneralStars(scip, redcostdata, NULL, graph, NULL, &nelims) );

   STPTEST_ASSERT_MSG(graph_edge_isDeleted(graph, 0), "edge 0 was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that edge can be deleted by using SD MST argument */
static
SCIP_RETCODE testPcEdgeDeletedByMst1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 16;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 2, 1.4);
   graph_edge_addBi(scip, graph, 0, 3, 1.4);
   graph_edge_addBi(scip, graph, 0, 4, 1.4);
   graph_edge_addBi(scip, graph, 2, 3, 1.1);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = FARAWAY;

   graph->prize[1] = 0.09;

   SCIP_CALL( stptest_graphSetUpRpcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArraysPc(graph, redcostdata);

   SCIP_CALL(extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}



/** tests that edge cannot be deleted (by using SD MST argument) */
static
SCIP_RETCODE testPcEdgeNotDeleted(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 16;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */

   graph->source = 0;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 1.0);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);

   /* add shortcut edges */
   graph_edge_addBi(scip, graph, 0, 2, 1.4);
   graph_edge_addBi(scip, graph, 0, 3, 1.4);
   graph_edge_addBi(scip, graph, 0, 4, 1.4);
   graph_edge_addBi(scip, graph, 2, 3, 1.1);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = FARAWAY;

   graph->prize[1] = 0.11;

   SCIP_CALL( stptest_graphSetUpRpcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArraysPc(graph, redcostdata);

   SCIP_CALL(extCheckEdge(scip, graph, redcostdata, edgedeleted, testedge, &deletable, FALSE));

   STPTEST_ASSERT_MSG(!deletable, "edge was deleted \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that node can be deleted */
static
SCIP_RETCODE testPcNode3PseudoDeletedBySd1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 6;
   const int nedges = 14;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = FARAWAY;

   graph->prize[0] = 1.0;
   graph->prize[2] = 0.9;

   SCIP_CALL( stptest_graphSetUpRpcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArraysPc(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}


/** tests that node can be deleted */
static
SCIP_RETCODE testPcNode4PseudoDeletedBySd1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 7;
   const int nedges = 18;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   STP_Bool* edgedeleted = NULL;
   int testnode = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */
   graph_knot_add(graph, STP_TERM);       /* node 6 */

   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.5);
   graph_edge_addBi(scip, graph, 0, 3, 2.0);
   graph_edge_addBi(scip, graph, 0, 4, 2.0);
   graph_edge_addBi(scip, graph, 2, 5, 2.0);
   graph_edge_addBi(scip, graph, 2, 6, 2.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 6, 3, 2.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = FARAWAY;

   graph->prize[0] = 0.0;
   graph->prize[2] = 1.0;
   graph->prize[3] = 0.6;

   SCIP_CALL( stptest_graphSetUpRpcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArraysPc(graph, redcostdata);

   SCIP_CALL( extCheckNode(scip, graph, redcostdata, edgedeleted, testnode, &deletable, FALSE) );

   STPTEST_ASSERT_MSG(deletable, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}

/* requires interface change! */
#ifdef SCIP_DISABLED
/** tests correctness of pseudo deletion of all nodes */
static
SCIP_RETCODE testPcNodesPseudoDeletedBySd1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   REDCOST* redcostdata;
   GRAPH* graph;
   const int nnodes = 6;
   const int nedges = 14;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   int nelims = 0;
   SCIP_Real offset = 0.0;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM);       /* node 0 */
   graph_knot_add(graph, STP_TERM);       /* node 1 */
   graph_knot_add(graph, STP_TERM);       /* node 2 */
   graph_knot_add(graph, STP_TERM);       /* node 3 */
   graph_knot_add(graph, STP_TERM);       /* node 4 */
   graph_knot_add(graph, STP_TERM);       /* node 5 */

   graph->source = 4;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 4, 1.0);
   graph_edge_addBi(scip, graph, 2, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);

   SCIP_CALL( graph_pc_initPrizes(scip, graph, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      graph->prize[i] = FARAWAY;

   graph->prize[0] = 1.0;
   graph->prize[2] = 0.9;

   SCIP_CALL( stptest_graphSetUpRpcOrg(scip, graph, NULL, NULL) );

   SCIP_CALL( redcosts_init(scip, graph->knots, graph->edges, cutoff, root, &redcostdata) );
   extInitRedCostArraysPcWithBase(graph, nnodes - 2, redcostdata);

   SCIP_CALL( extreduce_pseudoDeleteNodes(scip, NULL, redcostdata, graph, &offset, &nelims) );

   STPTEST_ASSERT_MSG(nelims == 1, "not enough eliminations \n");
   STPTEST_ASSERT_MSG(graph->grad[0] == 0, "node was not marked as deleteable! \n");

   stptest_extreduceTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}
#endif




/** tests that path replacement deletes an edge */
static
SCIP_RETCODE testPathReplaceDeletesEdge(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 9;
   const int nedges = 28;
   int nelims = 0;
   int testedge = 4;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);  /* node 1 */
   graph_knot_add(graph, STP_TERM);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 4 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 5 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 6 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 7 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 8 */


   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 7, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 4, 8, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 7, 8, 0.1);

   stptest_graphSetUp(scip, graph);

   SCIP_CALL( reduce_pathreplace(scip, graph, &nelims) );
   STPTEST_ASSERT_MSG(graph_edge_isDeleted(graph, testedge), "edge was not deleted! \n");

   stptest_graphTearDown(scip, graph);

   return SCIP_OKAY;
}



/** tests that path replacement deletes an edge */
static
SCIP_RETCODE testPathReplaceDeletesEdge2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 10;
   const int nedges = 28;
   int nelims = 0;
   int testedge = 4;

   assert(scip);

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 0 */
   graph_knot_add(graph, STP_TERM);  /* node 1 */
   graph_knot_add(graph, STP_TERM);  /* node 2 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 3 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 4 */
   graph_knot_add(graph, STP_TERM_NONE);       /* node 5 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 6 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 7 */
   graph_knot_add(graph, STP_TERM_NONE);  /* node 8 */
   graph_knot_add(graph, STP_TERM);  /* node 9 */


   graph->source = 1;

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 0, 2, 1.0);
   graph_edge_addBi(scip, graph, 0, 3, 2.0);
   graph_edge_addBi(scip, graph, 1, 5, 1.0);
   graph_edge_addBi(scip, graph, 2, 7, 1.0);
   graph_edge_addBi(scip, graph, 3, 4, 1.0);
   graph_edge_addBi(scip, graph, 3, 5, 1.0);
   graph_edge_addBi(scip, graph, 3, 7, 1.0);
   graph_edge_addBi(scip, graph, 4, 6, 1.0);
   graph_edge_addBi(scip, graph, 4, 8, 1.0);
   graph_edge_addBi(scip, graph, 5, 6, 1.0);
   graph_edge_addBi(scip, graph, 7, 8, 0.1);

   graph_edge_addBi(scip, graph, 9, 2, 0.5);
   graph_edge_addBi(scip, graph, 0, 9, 0.5);


   stptest_graphSetUp(scip, graph);

   SCIP_CALL( reduce_pathreplace(scip, graph, &nelims) );
   STPTEST_ASSERT_MSG(graph_edge_isDeleted(graph, testedge), "edge was not deleted! \n");

   stptest_graphTearDown(scip, graph);


   return SCIP_OKAY;
}


/** frees, etc. */
void stptest_extreduceTearDown(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST**             redcostdata         /**< reduced cost data */
   )
{
   assert(scip && graph);

   if( redcostdata )
      redcosts_free(scip, redcostdata);

   stptest_graphTearDown(scip, graph);
}


/** tests extended reduction techniques */
SCIP_RETCODE stptest_extreduce(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   assert(scip);
   SCIP_CALL( testPathReplaceDeletesEdge2(scip) );

   SCIP_CALL( testPathReplaceDeletesEdge(scip) );


   SCIP_CALL( testNode3PseudoDeletedBySdBiased(scip) );


   SCIP_CALL( testNode3PseudoDeletedBySdBiasedSimple(scip) );


   SCIP_CALL( testNode3PseudoDeletedByContraction(scip) );

   SCIP_CALL( testEdgeDeletedBy3LeafSpg(scip) );

   SCIP_CALL( testEdgeDeletedByMultiRedCosts(scip) );

   SCIP_CALL( testGeneralStarDeletedEdge1(scip) );
   SCIP_CALL( testGeneralStarDeletedEdge2(scip) );
   SCIP_CALL( testGeneralStarDeletedEdge3(scip) );

   SCIP_CALL( testEdgeDeletedByEqBottleneck2(scip) );
   SCIP_CALL( testEdgeDeletedByEqBottleneck(scip) );


   SCIP_CALL( testNode4PseudoDeletedBySd1(scip) );
   SCIP_CALL( testNode4PseudoNotDeletedBySd1(scip) );

   SCIP_CALL( testNode3PseudoDeletedByRedCosts1(scip) );

   SCIP_CALL( testPcNode4PseudoDeletedBySd1(scip) );
  // SCIP_CALL( testPcNodesPseudoDeletedBySd1(scip) );
   SCIP_CALL( testPcNode3PseudoDeletedBySd1(scip) );
   SCIP_CALL( testPcEdgeDeletedByMst1(scip) );
   SCIP_CALL( testPcEdgeNotDeleted(scip) );

   //SCIP_CALL( testPcNode3NotPseudoDeletedBySd1(scip) );
   SCIP_CALL( testNode3PseudoDeletedBySd1(scip) );
   SCIP_CALL( testNode3PseudoDeletedBySd2(scip) );
   SCIP_CALL( testNode3PseudoDeletedBySd3(scip) );


   SCIP_CALL( testEdgeDeletedByCommonRedCostsTargets(scip) );
   SCIP_CALL( testEdgeDeletedByMst2(scip) );
   SCIP_CALL( testEdgeDeletedByMst1(scip) );
   SCIP_CALL( testEdgeNotDeleted1(scip) );
   SCIP_CALL( testEdgeDeletion5_deprecated(scip, 1) );
   SCIP_CALL( testEdgeDeletion5_deprecated(scip, 2) );
   SCIP_CALL( testEdgeDeletion4_deprecated(scip, 1) );
   SCIP_CALL( testEdgeDeletion4_deprecated(scip, 2) );
   SCIP_CALL( testEdgeDeletion4_deprecated(scip, 3) );
   SCIP_CALL( testEdgeDeletion3_deprecated(scip, 1) );
   SCIP_CALL( testEdgeDeletion3_deprecated(scip, 2) );
   SCIP_CALL( testEdgeDeletion3_deprecated(scip, 3) );
   SCIP_CALL( testEdgeDeletion3_deprecated(scip, 4) );
   SCIP_CALL( testEdgeDeletion2_deprecated(scip, 1) );
   SCIP_CALL( testEdgeDeletion1_deprecated(scip) );



   printf("extreduce test: all ok \n");

   return SCIP_OKAY;
}
