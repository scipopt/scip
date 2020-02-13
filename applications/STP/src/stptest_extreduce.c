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
      extreduce_removeEdge(scip, edgedelete, graph, &distdata);
      assert(graph_isMarked(graph));
   }

   /* actual test */
   SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, edge, equality, &distdata, &extpermanent, deletable) );

   /* clean up */
   extreduce_extPermaFreeMembers(scip, &extpermanent);
   extreduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}


/** base method for extended edge reduction tests */
static
SCIP_RETCODE extCheckEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata,        /**< reduced cost data */
   STP_Bool*             edgedeleted,
   int                   edge,
   SCIP_Bool*            deletable,
   SCIP_Bool             allowEquality
)
{
   DISTDATA distdata;
   EXTPERMA extpermanent;

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STPTEST_EXT_MAXNCLOSENODES, FALSE, &distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, graph, edgedeleted, &extpermanent) );

   /* actual test */
   SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, edge, allowEquality, &distdata, &extpermanent, deletable) );

   /* clean up */
   extreduce_extPermaFreeMembers(scip, &extpermanent);
   extreduce_distDataFreeMembers(scip, graph, &distdata);
   graph_free_dcsr(scip, graph);

   return SCIP_OKAY;
}


/** frees, etc. */
static
void extTearDown(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   REDCOST*              redcostdata         /**< reduced cost data */
   )
{
   reduce_redcostdataFreeMembers(scip, redcostdata);

   graph_path_exit(scip, graph);
   graph_free(scip, &graph, TRUE);
   assert(graph == NULL);
}


/** helper */
static
void mldistsAddLevel(
   int                   nslots,
   int                   ntargets,
   const int             bases[nslots],
   const SCIP_Real       dists[nslots][ntargets],
   MLDISTS*              mldists
 )
{
   assert(nslots > 0 && ntargets > 0);

   extreduce_mldistsLevelAddTop(nslots, ntargets, mldists);

   for( int i = 0; i < nslots; i++ )
   {
      SCIP_Real* emptydists;

      extreduce_mldistsEmptySlotSetBase(bases[i], mldists);
      emptydists = extreduce_mldistsEmptySlotTargetDists(mldists);

      BMScopyMemoryArray(emptydists, dists[i], ntargets);

      extreduce_mldistsEmptySlotSetFilled(mldists);
   }
}

/** helper */
static
SCIP_Bool mldistsEqualDists(
   const MLDISTS*        mldists,
   const SCIP_Real       dists[],
   int                   level,
   int                   baseid,
   int                   ntargets
   )
{
   const SCIP_Real* dists_ml = extreduce_mldistsTargetDists(mldists, level, baseid);

   assert(ntargets == extreduce_mldistsLevelNTargets(mldists, level));

   for( int i = 0; i < ntargets; i++ )
   {
      if( !EQ(dists_ml[i], dists[i]) )
      {
         SCIPdebugMessage("wrong distance at %d: %f!=%f \n", i, dists_ml[i], dists[i]);
         return FALSE;
      }
   }

   return TRUE;
}


/** helper */
static
SCIP_Bool mldistsContainsBases(
   const MLDISTS*        mldists,
   const int             bases[],
   int                   level,
   int                   nslots
   )
{
   assert(nslots == extreduce_mldistsLevelNSlots(mldists, level));

   for( int i = 0; i < nslots; i++ )
   {
      if( !extreduce_mldistsLevelContainsBase(mldists, level, bases[i]) )
      {
         SCIPdebugMessage("bases %d not contained! \n", bases[i]);
         return FALSE;
      }
   }

   return TRUE;
}


/** initializes to default */
static
void initRedCostArrays(
   const GRAPH*          graph,              /**< the graph */
   REDCOST*              redcostdata         /**< reduced costs data */
)
{
   const int nnodes = graph->knots;
   const int nedges = graph->edges;
   SCIP_Real* const rootdist = redcostdata->rootToNodeDist;
   SCIP_Real* const redcost = redcostdata->redEdgeCost;
   PATH* const termpaths3 = redcostdata->nodeTo3TermsPaths;
   int* const vbase3 = redcostdata->nodeTo3TermsBases;

   assert(redcostdata->nnodes >= nnodes);
   assert(redcostdata->nedges >= nedges);

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


/** tests building and un-building */
static
SCIP_RETCODE mldistsTest1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   const int maxnlevels = 3;
   const int maxnslots = 12;
   const int maxntargets = 7;
#ifndef NDEBUG
   const SCIP_Real* dists;
#endif

   MLDISTS* mldists;

   SCIP_CALL( extreduce_mldistsInit(scip, maxnlevels, maxnslots, maxntargets, 0, TRUE, &mldists) );

   /* add */

   extreduce_mldistsLevelAddTop(2, 7, mldists);

#ifndef NDEBUG
   dists = extreduce_mldistsEmptySlotTargetDists(mldists);
   assert(dists);
#endif

   for( int i = 0; i < 2; i++ )
   {
      extreduce_mldistsEmptySlotSetBase(i + 10, mldists);
      extreduce_mldistsEmptySlotSetFilled(mldists);
   }

   extreduce_mldistsLevelAddTop(12, 5, mldists);

   for( int i = 0; i < 12; i++ )
   {
      extreduce_mldistsEmptySlotSetBase(i, mldists);
      extreduce_mldistsEmptySlotSetFilled(mldists);
   }

   extreduce_mldistsLevelAddTop(2, 7, mldists);

   for( int i = 0; i < 2; i++ )
   {
      extreduce_mldistsEmptySlotSetBase(i, mldists);
      extreduce_mldistsEmptySlotSetFilled(mldists);
   }

   /* remove */

   extreduce_mldistsLevelRemoveTop(mldists);
   extreduce_mldistsLevelRemoveTop(mldists);
   extreduce_mldistsLevelRemoveTop(mldists);

   extreduce_mldistsFree(scip, &mldists);

   return SCIP_OKAY;
}



/** tests correct storing */
static
SCIP_RETCODE mldistsTest2(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   const int maxnlevels = 4;
   const int maxnslots = 3;
   const int maxntargets = 4;

   const int nslots1 = 2;
   const int ntargets1 = 3;
   const int bases1[] = { 1, 4 };
   const SCIP_Real dists1[2][3] = { {1.0, 2.5,   4.5},
                                  {1.3, 0.5, 225.7} };

   const int nslots2 = 2;
   const int ntargets2 = 4;
   const int bases2[] = { 2, 22 };
   const SCIP_Real dists2[2][4] = { {15.0, 22.5, 455.2, 5.3},
                                  {65.0,  7.5,   5.2, 5.7} };

   const int nslots3 = 3;
   const int ntargets3 = 2;
   const int bases3[] = { 2, 22, 6 };
   const SCIP_Real dists3[3][2] = { { 0.1, 0.05},
                                  { 4.1, 2.11},
                                  {44.1, 9.14} };

   MLDISTS* mldists;

   SCIP_CALL( extreduce_mldistsInit(scip, maxnlevels, maxnslots, maxntargets, 0, TRUE, &mldists) );

   /* add first */
   mldistsAddLevel(nslots1, ntargets1, bases1, dists1, mldists);

   /* add second */
   mldistsAddLevel(nslots2, ntargets2, bases2, dists2, mldists);

   if( !mldistsEqualDists(mldists, dists1[1], 0, 4, ntargets1) )
   {
      return SCIP_ERROR;
   }

   if( !mldistsContainsBases(mldists, bases1, 0, nslots1) )
   {
      SCIPdebugMessage("bases fail \n");
      return SCIP_ERROR;
   }

   /* add third */
   mldistsAddLevel(nslots3, ntargets3, bases3, dists3, mldists);

   if( !mldistsEqualDists(mldists, dists3[2], 2, 6, ntargets3) )
   {
      return SCIP_ERROR;
   }

   /* remove third */

   extreduce_mldistsLevelRemoveTop(mldists);


   if( !mldistsEqualDists(mldists, dists2[0], 1, 2, ntargets2) )
   {
      return SCIP_ERROR;
   }

   if( !mldistsContainsBases(mldists, bases2, 1, nslots2) )
   {
      SCIPdebugMessage("bases fail \n");
      return SCIP_ERROR;
   }

   /* remove second */
   extreduce_mldistsLevelRemoveTop(mldists);

   if( !mldistsEqualDists(mldists, dists1[0], 0, 1, ntargets1) )
   {
      return SCIP_ERROR;
   }

   /* remove first */
   extreduce_mldistsLevelRemoveTop(mldists);

   extreduce_mldistsFree(scip, &mldists);

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
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 100.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST redcostdata;

   assert(variant == 1 || variant == 2);

   SCIP_CALL( reduce_redcostdataInit(scip, nnodes, nedges, cutoff, root, &redcostdata) );

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

   graph_mark(graph);

   initRedCostArrays(graph, &redcostdata);

   for( int i = 0; i < nedges; i++ )
      redcostdata.redEdgeCost[i] = 1.0;

   edge = 0;

   graph_knot_chg(graph, 6, 0);
   graph_knot_chg(graph, 7, 0);
   graph_knot_chg(graph, 8, 0);
   graph_knot_chg(graph, 9, 0);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   if( variant == 1 )
   {
      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
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

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }

   extTearDown(scip, graph, &redcostdata);

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
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 100.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST redcostdata;

   assert(scip);
   assert(variant == 1 || variant == 2 || variant == 3);

   SCIP_CALL( reduce_redcostdataInit(scip, nnodes, nedges, cutoff, root, &redcostdata) );

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

   initRedCostArrays(graph, &redcostdata);

   for( int i = 0; i < nedges; i++ )
      redcostdata.redEdgeCost[i] = 1.0;

   edge = 0;

   graph_knot_chg(graph, 4, 0);
   graph_knot_chg(graph, 5, 0);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   if( variant == 1 )
   {
      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));

      assert(deletable);
   }
   else if( variant == 2 )
   {
      int edgedelete = -1;

      for( int e = graph->outbeg[0]; e != EAT_LAST; e = graph->oeat[e] )
         if( graph->head[e] == 7 )
            edgedelete = e;

      assert(edgedelete >= 0);

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
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

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(deletable);
   }

   extTearDown(scip, graph, &redcostdata);

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
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 7.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST redcostdata;

   assert(scip);

   SCIP_CALL( reduce_redcostdataInit(scip, nnodes, nedges, cutoff, root, &redcostdata) );

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

   initRedCostArrays(graph, &redcostdata);

   redcostdata.nodeTo3TermsBases[5] = nnodes - 1;
   redcostdata.nodeTo3TermsBases[6] = nnodes - 1;
   redcostdata.nodeTo3TermsPaths[5].dist = 3.0;
   redcostdata.nodeTo3TermsPaths[6].dist = 3.0;
   redcostdata.nodeTo3TermsPaths[5 + nnodes].dist = 3.0;
   redcostdata.nodeTo3TermsPaths[6 + nnodes].dist = 3.0;

   for( int i = 0; i < nedges; i++ )
      redcostdata.redEdgeCost[i] = 1.0;

   edge = 0;

   graph_knot_chg(graph, 7, 0);

   if( variant == 1 )
   {
      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(deletable);
   }
   else if( variant == 2 )
   {
      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));
      assert(deletable);
   }
   else if( variant == 3 )
   {
      const int edge2 = 14;
      assert(graph->tail[edge2] == 7 && graph->head[edge2] == 9);

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      graph->cost[edge2] = 1.0;
      graph->cost[flipedge(edge2)] = 1.0;
      graph_knot_chg(graph, 4, 0);

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, TRUE));
      assert(!deletable);
   }
   else if( variant == 4 )
   {
      const int edgedelete = 14;
      assert(graph->tail[edgedelete] == 7 && graph->head[edgedelete] == 9);

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(!deletable);
   }
   else
   {
      const int edgedelete = graph->edges;

      graph_edge_addBi(scip, graph, 7, 0, 0.1);
      assert(graph->tail[edgedelete] == 7 && graph->head[edgedelete] == 0);

      SCIP_CALL( graph_init_history(scip, graph) );
      SCIP_CALL( graph_path_init(scip, graph) );

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, edgedelete, nnodes, &deletable, TRUE));
      assert(deletable);
   }

   extTearDown(scip, graph, &redcostdata);

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

   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 100.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST redcostdata;

   assert(scip);

   SCIP_CALL( reduce_redcostdataInit(scip, nnodes, nedges, cutoff, root, &redcostdata) );

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

   initRedCostArrays(graph, &redcostdata);

   edge = 0;

   if( variant == 1 )
   {
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

      SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));

      STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");
   }

   extTearDown(scip, graph, &redcostdata);


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
   STP_Bool* edgedeleted = NULL;
   SCIP_Real cutoff = 0.0;
   int edge;
   SCIP_Bool deletable;
   REDCOST redcostdata;

   assert(scip);

   SCIP_CALL( reduce_redcostdataInit(scip, nnodes, nedges, cutoff, root, &redcostdata) );

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

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   graph_mark(graph);

   initRedCostArrays(graph, &redcostdata);

   edge = 0;

   SCIP_CALL(extCheckArc(scip, graph, &redcostdata, edgedeleted, edge, -1, nnodes, &deletable, FALSE));

   STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");

   extTearDown(scip, graph, &redcostdata);

   return SCIP_OKAY;
}

static
SCIP_RETCODE testEdgeDeletedByMst1(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   GRAPH* graph;
   const int nnodes = 5;
   const int nedges = 16;
   const int root = 0;
   SCIP_Real cutoff = 100.0;
   REDCOST redcostdata;

   STP_Bool* edgedeleted = NULL;
   int testedge = 0;
   SCIP_Bool deletable;

   assert(scip);

   SCIP_CALL( reduce_redcostdataInit(scip, nnodes, nedges, cutoff, root, &redcostdata) );

   SCIP_CALL( graph_init(scip, &graph, nnodes, nedges, 1) );

   /* build tree*/
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

   graph_mark(graph);

   initRedCostArrays(graph, &redcostdata);

   SCIP_CALL( graph_init_history(scip, graph) );
   SCIP_CALL( graph_path_init(scip, graph) );

   SCIP_CALL(extCheckEdge(scip, graph, &redcostdata, edgedeleted, testedge, &deletable, FALSE));

 //  STPTEST_ASSERT_MSG(deletable, "edge was not deleted \n");


   extTearDown(scip, graph, &redcostdata);

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

   graph_edge_addBi(scip, graph, 0, 1, 1.0);
   graph_edge_addBi(scip, graph, 1, 2, 0.4);
   graph_edge_addBi(scip, graph, 1, 3, 1.0);
   graph_edge_addBi(scip, graph, 1, 4, 1.0);

   graph_edge_addBi(scip, graph, 2, 5, 0.5);
   graph_edge_addBi(scip, graph, 3, 6, 1.6);
   graph_edge_addBi(scip, graph, 3, 7, 1.6);
   graph_edge_addBi(scip, graph, 4, 8, 1.5);
   graph_edge_addBi(scip, graph, 4, 9, 1.5);

   graph_edge_addBi(scip, graph, 5, 10, 1.0);
   graph_edge_addBi(scip, graph, 10, 11, 1.0);
   graph_edge_addBi(scip, graph, 11, 12, 1.0);

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
      const SCIP_Real dist1_2 = extreduce_distDataGetSd(scip, graph, 1, 2, &distdata);
      const SCIP_Real dist1_3 = extreduce_distDataGetSd(scip, graph, 1, 3, &distdata);
      const SCIP_Real dist1_5 = extreduce_distDataGetSd(scip, graph, 1, 5, &distdata);
      const SCIP_Real dist2_1 = extreduce_distDataGetSd(scip, graph, 2, 1, &distdata);
      const SCIP_Real dist2_5 = extreduce_distDataGetSd(scip, graph, 2, 5, &distdata);
      const SCIP_Real dist2_6 = extreduce_distDataGetSd(scip, graph, 2, 6, &distdata);

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
         const SCIP_Real dist1_2_b = extreduce_distDataGetSd(scip, graph, 1, 2, &distdata);
         const SCIP_Real dist2_6_b = extreduce_distDataGetSd(scip, graph, 2, 6, &distdata);
         const SCIP_Real dist2_5_b = extreduce_distDataGetSd(scip, graph, 2, 5, &distdata);

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

/** tests for mldists */
SCIP_RETCODE stptest_extmldists(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( mldistsTest1(scip) );
   SCIP_CALL( mldistsTest2(scip) );


   printf("mldists test: all ok \n");

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

   SCIP_CALL( testEdgeDeletedByMst1(scip) );

   SCIP_CALL( extTest2_variants(scip, 1) );

   SCIP_CALL( extTest1(scip) );

   SCIP_CALL( extDistTest(scip) );

   printf("extreduce test: all ok \n");

   return SCIP_OKAY;
}
