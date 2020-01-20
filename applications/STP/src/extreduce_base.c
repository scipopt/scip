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

/**@file   extreduce_base.c
 * @brief  extended reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements interface methods for extended reduction techniques for several Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "extreduce.h"

#ifndef NDEBUG
/** all good with the graph->mark array? */
static
SCIP_Bool graphmarkIsClean(
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const GRAPH*          graph              /**< graph data structure */
)
{
   const int nnodes = graph_get_nNodes(graph);

   for( int k = 0; k < nnodes; k++ )
   {
      if( graph->grad[k] == 0 && k != redcostdata->redCostRoot && !Is_term(graph->term[k]) )
      {
         if( graph->mark[k] )
            return FALSE;
      }
   }

   return TRUE;
}
#endif


/** deletes an edge and makes corresponding adaptations */
static
void removeEdge(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to delete */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   DISTDATA*             distdata            /**< distance data (in/out) */
)
{
   const int tail = graph->tail[edge];
   const int head = graph->head[edge];

#ifdef SCIP_DEBUG
   SCIPdebugMessage("removing edge ");
   graph_edge_printInfo(graph, edge);
#endif

   assert(extreduce_distCloseNodesAreValid(scip, graph, distdata));

   graph_edge_delFull(scip, graph, edge, TRUE);
   extreduce_distDataDeleteEdge(scip, graph, edge, distdata);

   if( graph->grad[tail] == 0 )
   {
      if( Is_term(graph->term[tail])  )
      {
         assert(graph_pc_isPcMw(graph) || tail == graph->source);
      }
      else
      {
         graph->mark[tail] = FALSE;
      }
   }

   if( graph->grad[head] == 0 )
   {
      if( Is_term(graph->term[head]) || head == graph->source )
      {
         assert(graph_pc_isPcMw(graph));
      }
      else
      {
         graph->mark[head] = FALSE;
      }
   }

   assert(extreduce_distCloseNodesAreValid(scip, graph, distdata));
}

/** initialize */
static
SCIP_RETCODE extInit(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   DISTDATA*             distdata,           /**< distance data (out) */
   EXTPERMA*             extpermanent        /**< permanent extension data (out) */
)
{
   assert(!graph_pc_isPcMw(graph) || !graph->extended);

   graph_mark(graph);

   SCIP_CALL( graph_init_dcsr(scip, graph) );
   SCIP_CALL( extreduce_distDataInit(scip, graph, STP_EXT_CLOSENODES_MAXN, FALSE, distdata) );
   SCIP_CALL( extreduce_extPermaInit(scip, graph, edgedeletable, extpermanent) );

   return SCIP_OKAY;
}


/** free */
static
void extFree(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   DISTDATA*             distdata,           /**< distance data (in/out) */
   EXTPERMA*             extpermanent        /**< permanent extension data (in/out) */
)
{
   extreduce_extPermaFreeMembers(scip, extpermanent);
   extreduce_distDataFreeMembers(scip, graph, distdata);
   graph_free_dcsr(scip, graph);
}


/** Extended reduction test for arcs.
 * This method will also set edgedeletable[a] to TRUE if arc 'a' can be deleted, but its anti-parallel arc not. */
SCIP_RETCODE extreduce_deleteArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed */
   int*                  nelims              /**< number of eliminations */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const redcost = redcostdata->redEdgeCost;
   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, graph, edgedeletable, &distdata, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, e) )
      {
         const int erev = e + 1;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) && SCIPisZero(scip, redcost[erev]) )
            continue;

         if( !edgedeletable[e] )
         {
            SCIP_Bool deletable;
            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, e, allowequality, &distdata, &extpermanent,
                  &deletable) );

            if( deletable )
               edgedeletable[e] = TRUE;
         }

         if( !edgedeletable[erev] )
         {
            SCIP_Bool erevdeletable = FALSE;

            SCIP_CALL( extreduce_checkArc(scip, graph, redcostdata, erev, allowequality, &distdata, &extpermanent,
                  &erevdeletable) );

            if( erevdeletable )
               edgedeletable[erev] = TRUE;
         }

         if( edgedeletable[e] && edgedeletable[erev] )
         {
            assert(edgedeletable[e] && edgedeletable[erev]);

            removeEdge(scip, e, graph, &distdata);

            (*nelims)++;
         }
      }
   }

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}


/** extended reduction test for edges */
SCIP_RETCODE extreduce_deleteEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const int nedges = graph_get_nEdges(graph);
   const SCIP_Real* const redcost = redcostdata->redEdgeCost;
   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);
   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, graph, edgedeletable, &distdata, &extpermanent) );

   /* main loop */
   for( int e = 0; e < nedges; e += 2 )
   {
      if( extreduce_edgeIsValid(graph, e) )
      {
         const int erev = e + 1;
         SCIP_Bool deletable = TRUE;
         const SCIP_Bool allowequality = (result != NULL && result[e] != CONNECT && result[erev] != CONNECT);

         assert(flipedge(e) == erev && SCIPisEQ(scip, graph->cost[e], graph->cost[erev]));

         if( SCIPisZero(scip, redcost[e]) && SCIPisZero(scip, redcost[erev]) )
            continue;

         SCIP_CALL( extreduce_checkEdge(scip, graph, redcostdata, e, allowequality, &distdata, &extpermanent, &deletable) );

         if( deletable )
         {
            removeEdge(scip, e, graph, &distdata);

            (*nelims)++;
         }
      }
   }

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}


/** extended reduction test for pseudo-eliminating nodes */
SCIP_RETCODE extreduce_pseudodeleteNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const REDCOST*        redcostdata,        /**< reduced cost data */
   const int*            result,             /**< solution array or NULL */
   GRAPH*                graph,              /**< graph data structure (in/out) */
   STP_Bool*             edgedeletable,      /**< edge array to mark which (directed) edge can be removed (in/out) */
   int*                  nelims              /**< number of eliminations (out) */
)
{
   const int nnodes = graph_get_nNodes(graph);
   SCIP_Bool* pseudoDeletable;
   DISTDATA distdata;
   EXTPERMA extpermanent;

   assert(scip && redcostdata && edgedeletable);
   assert(redcostdata->redCostRoot >= 0 && redcostdata->redCostRoot < graph->knots);

   *nelims = 0;

   if( SCIPisZero(scip, redcostdata->cutoff) )
      return SCIP_OKAY;

   SCIP_CALL( extInit(scip, graph, edgedeletable, &distdata, &extpermanent) );

   SCIP_CALL( SCIPallocBufferArray(scip, &pseudoDeletable, nnodes) );

   /* main loop */
   for( int i = 0; i < nnodes; ++i )
   {
      SCIP_Bool nodeisDeletable = FALSE;

      SCIP_CALL( extreduce_checkNode(scip, graph, redcostdata, i, FALSE, &distdata, &extpermanent, &nodeisDeletable) );

      pseudoDeletable[i] = nodeisDeletable;
   }

   // todo: perform the actual pseudo-elimination based on pseudoDeletable



   SCIPfreeBufferArray(scip, &pseudoDeletable);

   extFree(scip, graph, &distdata, &extpermanent);

   assert(graphmarkIsClean(redcostdata, graph));

   return SCIP_OKAY;
}


/** get maximum allowed stack size */
int extreduce_getMaxStackSize(void)
{
   return STP_EXT_MAXSTACKSIZE;
}


/** get maximum allowed depth for extended tree in given graph */
int extreduce_getMaxTreeDepth(
   const GRAPH*          graph               /**< graph data structure */
)
{
   const int maxdepth = (graph->edges > STP_EXT_EDGELIMIT) ? STP_EXT_MINDFSDEPTH : STP_EXT_MAXDFSDEPTH;

   assert(maxdepth > 0);

   return maxdepth;
}
