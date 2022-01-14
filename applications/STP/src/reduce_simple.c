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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_simple.c
 * @brief  several basic reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements simple reduction techniques for several Steiner problems.
 * Mosts tests are described in "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "portab.h"
#include "scip/scip.h"


#ifdef SCIP_DISABLED
/** deletes entire graph */
static
void deleteAllEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int nedges = graph_get_nEdges(g);

   for( int e = 0; e < nedges; e += 2 )
   {
      if( g->oeat[e] != EAT_FREE )
         graph_edge_del(scip, g, e, TRUE);
   }
}
#endif


/** prune */
static
SCIP_RETCODE cutEdgePrune(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   start,              /**< the node to start from */
   int                   end,                /**< note to ignore */
   int                   cutedge,            /**< the edge */
   GRAPH*                g                   /**< graph data structure */
)
{
   if( Is_term(g->term[start]) )
   {
      int e = g->outbeg[start];

      while( e >= 0 )
      {
         const int enext = g->oeat[e];
         const int head = g->head[e];
         if( head != end )
         {
            assert(!Is_term(g->term[head]));
            graph_edge_del(scip, g, e, TRUE);
         }

         e = enext;
      }
   }
   else
   {
      graph_edge_del(scip, g, cutedge, TRUE);
   }

   SCIP_CALL( reduce_unconnected(scip, g) );

   return SCIP_OKAY;
}


/** probe */
static
SCIP_RETCODE cutEdgeProbe(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< the node to start from */
   int                   end,                /**< note to ignore */
   SCIP_Bool* RESTRICT   nodes_visited,      /**< marks for each node whether visited or not */
   SCIP_Bool*            terminalFound       /**< found? */
)
{
   int* RESTRICT stackarr;
   int stacksize = 0;
   const int nnodes = graph_get_nNodes(g);
   *terminalFound = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &stackarr, nnodes) );

   nodes_visited[start] = TRUE;
   nodes_visited[end] = TRUE;
   stackarr[stacksize++] = start;

   while( stacksize > 0 )
   {
      const int node = stackarr[--stacksize];

      for( int e = g->outbeg[node]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];

         if( !nodes_visited[head] )
         {
            if( Is_term(g->term[head]) )
            {
               *terminalFound = TRUE;
               stacksize = 0;
               break;
            }

            nodes_visited[head] = TRUE;
            stackarr[stacksize++] = head;
         }
      }
   }

   SCIPfreeBuffer(scip, &stackarr)

   return SCIP_OKAY;
}


/** removes nodes of degree one and unmarks */
void reduce_nodesDeg1(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                graph               /**< graph data structure (in/out) */
)
{
   const int nnodes = graph_get_nNodes(graph);
   SCIP_Bool rerun = TRUE;

   while( rerun )
   {
      rerun = FALSE;
      for( int i = 0; i < nnodes; ++i )
      {
         if( graph->grad[i] == 1 && !Is_term(graph->term[i]) )
         {
            const int sibling = graph->head[graph->outbeg[i]];
            assert(graph_knot_isInRange(graph, sibling));

            graph_knot_del(scip, graph, i, TRUE);
            graph->mark[i] = FALSE;

            if( Is_term(graph->term[sibling]) )
               continue;

            if( graph->grad[sibling] == 0 )
               graph->mark[sibling] = FALSE;
            else if( graph->grad[sibling] == 1 )
               rerun = TRUE;
         }
      }
   }
}


/** basic reduction tests for the STP */
SCIP_RETCODE reduce_simple(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                  or NULL */
   int*                  nelims,             /**< pointer to number of reductions */
   int*                  edgestate           /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   )
{
   SCIP_Bool rerun = TRUE;
   int elimscount = 0;
   const int nnodes = g->knots;
   const SCIP_Bool checkstate = (edgestate != NULL);
   SCIP_Bool replaceConflict = FALSE;

   assert(scip && g && fixed && nelims);

   SCIPdebugMessage("Degree Test: ");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( int i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);

         if( g->grad[i] == 1 && g->terms > 1 )
         {
            const int e1  = g->outbeg[i];
            const int i1  = g->head[e1];

            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->inpbeg[i]));
            assert(g->oeat[e1] == EAT_LAST);
            assert(g->ieat[g->inpbeg[i]] == EAT_LAST);

            if( checkstate )
            {
               if( edgestate[e1] == EDGE_BLOCKED || Is_term(g->term[i]) )
                  continue;
               else
               {
                  SCIPdebugMessage("delete degree 1 node: %d \n", i);
                  graph_edge_del(scip, g, e1, TRUE);
               }
            }
            else
            {
               if( Is_term(g->term[i]) )
               {
                  *fixed += g->cost[e1];
                  SCIP_CALL( graph_fixed_addEdge(scip, e1, g) );
               }

               SCIPdebugMessage("contract degree 1 terminal: %d \n", i);
               SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
            }
            elimscount++;

            assert(g->grad[i] == 0);

            /* the last node in the graph? */
            if( g->grad[i1] == 0 )
            {
               rerun = FALSE;
               break;
            }
            if( (i1 < i) && (g->grad[i1] < 3) )
               rerun = TRUE;

            continue;
         }

         if( g->grad[i] == 2 && !checkstate  )
         {
            const int e1 = g->outbeg[i];
            const int e2 = g->oeat[e1];
            const int i1 = g->head[e1];
            const int i2 = g->head[e2];
            SCIP_Bool eliminationDone = TRUE;

            assert(e1 >= 0);
            assert(e2 >= 0);

            do
            {
               if( !Is_term(g->term[i]) )
               {
                  SCIP_Bool conflict;
                  SCIPdebugMessage("replace degree 2 node: %d \n", i);
                  SCIP_CALL( graph_knot_replaceDeg2(scip, i, g->cost[e1] + g->cost[e2], -1, g, &conflict) );

                  if( conflict)
                     replaceConflict = TRUE;

                  elimscount++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if( Is_term(g->term[i1]) && Is_term(g->term[i2]) )
               {
                  SCIPdebugMessage("contract degree 2 terminal (with terminal neighbors): %d \n", i);

                  if( SCIPisLT(scip, g->cost[e1], g->cost[e2]) )
                  {
                     *fixed += g->cost[e1];

                     SCIP_CALL( graph_fixed_addEdge(scip, e1, g) );
                     SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
                  }
                  else
                  {
                     *fixed += g->cost[e2];

                     SCIP_CALL( graph_fixed_addEdge(scip, e2, g) );
                     SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
                  }
                  elimscount++;

                  break;
               }
               if( Is_term(g->term[i1]) && !Is_term(g->term[i2]) && SCIPisLE(scip, g->cost[e1], g->cost[e2]) )
               {
                  SCIPdebugMessage("contract degree 2 terminal (with one terminal neighbor): %d \n", i);

                  *fixed += g->cost[e1];

                  SCIP_CALL( graph_fixed_addEdge(scip, e1, g) );
                  SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
                  elimscount++;
                  break;
               }
               if( Is_term(g->term[i2]) && !Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e2], g->cost[e1]) )
               {
                  SCIPdebugMessage("contract degree 2 terminal (with one terminal neighbor): %d \n", i);

                  *fixed += g->cost[e2];

                  SCIP_CALL( graph_fixed_addEdge(scip, e2, g) );
                  SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
                  elimscount++;
                  break;
               }
               eliminationDone = FALSE;
            }
            while( FALSE );

            if( eliminationDone && (((i1 < i) && (g->grad[i1] < 3)) || ((i2 < i) && (g->grad[i2] < 3))) )
               rerun = TRUE;
         }

         if( Is_term(g->term[i]) && g->grad[i] > 2 && !checkstate )
         {
            SCIP_Real mincost = FARAWAY;
            int ett = UNKNOWN;
            for( int e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            {
               const int i1 = g->head[e1];

               if( SCIPisLT(scip, g->cost[e1], mincost) )
               {
                  mincost = g->cost[e1];
                  if( Is_term(g->term[i1]) )
                     ett = e1;
               }
               else if( Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e1], mincost) )
               {
                  ett = e1;
               }
            }
            if( ett != UNKNOWN && SCIPisLE(scip, g->cost[ett], mincost) )
            {
               SCIPdebugMessage("contract terminal into terminal: %d \n", i);

               *fixed += g->cost[ett];

               SCIP_CALL( graph_fixed_addEdge(scip, ett, g) );
               SCIP_CALL( graph_knot_contractLowdeg2High(scip, g, solnode, i, g->head[ett]) );

               rerun = TRUE;
            }
         }
      }
   }

   /* NOTE: in this case we might have made the graph unconnected... */
   if( replaceConflict )
   {
      SCIP_CALL( reduce_unconnected(scip, g) );
   }

   /* todo: seems to hurt performance in ascend-prune, not sure why.... */
#ifdef SCIP_DISABLED
   if( g->terms == 1 )
   {
      deleteAllEdges(scip, g);
   }
#endif

   SCIPdebugMessage(" %d Knots deleted\n", elimscount);
   assert(graph_valid(scip, g));

   *nelims += elimscount;
   return SCIP_OKAY;
}


/** basic reduction tests for the SAP */
SCIP_RETCODE reduce_simple_sap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_QUEUE* queue;
   int e;
   int i1;
   int i2;
   int e1;
   int e2;
   const int nnodes = graph_get_nNodes(g);
   char rerun;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);

   rerun = TRUE;
   *count = 0;

   SCIPdebugMessage("Degree Test: \n");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( int i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);

         if( g->grad[i] == 1 )
         {
            e1  = g->inpbeg[i];
            i1  = g->tail[e1];

            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->outbeg[i]));
            assert(g->ieat[e1] == EAT_LAST);
            assert(g->oeat[g->outbeg[i]] == EAT_LAST);

            if( Is_term(g->term[i]) )
            {
               if( g->terms == 1 )
                  continue;

               if( i == g->source )
               {
                  e2 = flipedge(e1);
                  *fixed += g->cost[e2];

                  SCIPdebugMessage("contract grad-1 terminal (source) \n");
                  SCIP_CALL( graph_knot_contractFixed(scip, g, NULL, e2, i1, i) );
               }
               else
               {
                  *fixed += g->cost[e1];
#ifdef SCIP_DEBUG
                  SCIPdebugMessage("contract grad-1 terminal ");
                  graph_knot_printInfo(g, i);
                  SCIPdebugMessage("with edge ");
                  graph_edge_printInfo(g, e1);
                  SCIPdebugMessage("into ");
                  graph_knot_printInfo(g, i1);
#endif

                  SCIP_CALL( graph_knot_contractFixed(scip, g, NULL, e1, i1, i) );
               }
            }
            else
            {
               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete grad-1 node \n");
            }

            assert(g->grad[i] == 0);

            if ((i1 < i) && (g->grad[i1] < 3))
               rerun = TRUE;

            (*count)++;

            continue;
         }

         if( g->grad[i] == 2 )
         {
            e1 = g->outbeg[i];
            e2 = g->oeat[e1];
            i1 = g->head[e1];
            i2 = g->head[e2];

            assert(e1 >= 0);
            assert(e2 >= 0);

            if( !Is_term(g->term[i]) )
            {
               if( (!Is_term(g->term[i2]) && !Is_term(g->term[i1])) )
               {
                  g->cost[e1] += g->cost[Edge_anti(e2)];
                  g->cost[Edge_anti(e1)] += g->cost[e2];

                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e1]), g->ancestors[Edge_anti(e2)], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(e1)]), g->ancestors[e2], NULL) );

                  if( SCIPisGT(scip, g->cost[e1], FARAWAY) )
                     g->cost[e1] = FARAWAY;
                  if( SCIPisGT(scip, g->cost[Edge_anti(e1)], FARAWAY) )
                     g->cost[Edge_anti(e1)] = FARAWAY;

                  SCIP_CALL( graph_knot_contract(scip, g, NULL, i2, i) );

                  SCIPdebugMessage("replace grad-2 node %d \n", i);

                  (*count)++;
                  if( ((i1 < i) && (g->grad[i1] < 3))
                     || ((i2 < i) && (g->grad[i2] < 3)) )
                     rerun = TRUE;
               }
            }
            /* CONSTCOND */
            /*lint -save -e717 */
            /*lint -restore */
         }
      }
   }

   /* delete all arcs in \delta^-(root) */
   for( e = g->inpbeg[g->source]; e != EAT_LAST; e = g->ieat[e] )
      g->cost[e] = FARAWAY;

   /* delete all nodes not reachable from the root by forward arcs */

   /* BFS  */
   SCIP_CALL(SCIPqueueCreate(&queue, nnodes, 1.1));

   for (int i = 0; i < nnodes; i++)
      g->mark[i] = FALSE;

   g->mark[g->source] = TRUE;
   SCIP_CALL(SCIPqueueInsert(queue, &(g->source)));

   while (!SCIPqueueIsEmpty(queue))
   {
      int* pnode = (SCIPqueueRemove(queue));
      for (e = g->outbeg[*pnode]; e != EAT_LAST; e = g->oeat[e])
      {
         if( !g->mark[g->head[e]] && LT(g->cost[e], FARAWAY) )
         {
            g->mark[g->head[e]] = TRUE;
            SCIP_CALL(SCIPqueueInsert(queue, &(g->head[e])));
         }
      }
   }

   for( int i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] && g->grad[i] > 0 )
      {
         assert(!Is_term(g->term[i]));
         SCIPdebugMessage("deleting unreachable (forward) node %d \n", i);
         graph_knot_del(scip, g, i, TRUE);
      }
   }

   /* delete all nodes that cannot reach a terminal other than the root by forward arcs (not using the root) */

   /* BFS */

   assert(SCIPqueueIsEmpty(queue));

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != g->source )
      {
         g->mark[i] = TRUE;
         SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[g->outbeg[i]])) );
      }
      else
      {
         g->mark[i] = FALSE;
      }
   }

   g->mark[g->source] = TRUE;

   while( !SCIPqueueIsEmpty(queue) )
   {
      int* pnode = (SCIPqueueRemove(queue));
      for( e = g->inpbeg[*pnode]; e != EAT_LAST; e = g->ieat[e] )
      {
         if( !g->mark[g->tail[e]] && LT(g->cost[e], FARAWAY) )
         {
            g->mark[g->tail[e]] = TRUE;
            SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[e])) );
         }
      }
   }

   SCIPqueueFree(&queue);

   for( int i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] && g->grad[i] > 0 )
      {
         assert(!Is_term(g->term[i]));
         SCIPdebugMessage("deleting unreachable (backward) node %d \n", i);
         graph_knot_del(scip, g, i, TRUE);
      }
   }

   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);
   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/** basic reduction tests for the DCSTP...call only once! */
void reduce_simple_dc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
   )
{
   const int nnodes = graph_get_nNodes(g);
   const int* const maxdeg = g->maxdeg;

   assert(g->stp_type == STP_DCSTP);
   assert(scip);
   assert(g->maxdeg);

   if( g->terms <= 2 )
      return;

   for( int i = 0; i < nnodes; i++ )
   {
      int enext;

      if( maxdeg[i] != 1 )
         continue;

      for( int e = g->outbeg[i]; e != EAT_LAST; e = enext )
      {
         const int head = g->head[e];
         enext = g->oeat[e];
         if( maxdeg[head] == 1 )
         {
            graph_edge_del(scip, g, e, TRUE);
         }
      }
   }
}


/** basic reduction tests for the HCDSTP */
SCIP_RETCODE reduce_simple_hc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int i;
   int e;
   int e2;
   int nnodes;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(g->stp_type == STP_DHCSTP);

   nnodes = g->knots;

   SCIPdebugMessage("basic HC test: \n");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      /* delete incoming arcs of the root */
      e = g->inpbeg[g->source];
      while( e != EAT_LAST )
      {
         e2 = g->ieat[e];

         if( SCIPisGE(scip, g->cost[flipedge(e)], FARAWAY) )
         {
            SCIPdebugMessage("delete incoming root arc \n");
            (*count)++;
            graph_edge_del(scip, g, e, TRUE);
         }
         else if( SCIPisLT(scip, g->cost[e], FARAWAY) )
         {
            SCIPdebugMessage("delete anti-parallel root arcs \n");
            g->cost[e] = FARAWAY;
         }

         e = e2;
      }

      /* delete outgoing arcs of the terminals (other than the root) */
      for( i = 0; i < nnodes; i++ )
      {
         if( Is_term(g->term[i]) && i != g->source )
         {
            e = g->outbeg[i];
            while( e != EAT_LAST )
            {
               e2 = g->oeat[e];

               if( SCIPisGE(scip, g->cost[flipedge(e)], FARAWAY) )
               {
                  SCIPdebugMessage("delete anti-parallel terminal arcs \n");
                  (*count)++;
                  graph_edge_del(scip, g, e, TRUE);
               }

               e = e2;
            }
         }
      }
   }

   SCIPdebugMessage("HC basic reduction package has deleted %d edges\n", *count);

   return SCIP_OKAY;
}


/** contract edges of weight zero */
SCIP_RETCODE reduce_contract0Edges(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< solution node mark or NULL */
   SCIP_Bool             savehistory         /**< save the history? */
   )
{
   int count;
   int nedges;

   assert(g != NULL);
   assert(scip != NULL);

   nedges = g->edges;

   do
   {
      count = 0;
      for( int e = 0; e < nedges; e += 2 )
      {
         if( g->oeat[e] != EAT_FREE && SCIPisZero(scip, g->cost[e]) && SCIPisZero(scip, g->cost[flipedge(e)]) )
         {
            if( savehistory )
            {
               if( graph_pc_isPcMw(g) )
               {
                  IDX* ans = graph_edge_getAncestors(g, e);
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[g->head[e]]), ans, NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[g->tail[e]]), ans, NULL) );
                  assert(0 && "currently not implemented");
               }
               else
               {
                  SCIP_CALL( graph_fixed_addEdge(scip, e, g) );
               }

            }

            SCIP_CALL( graph_knot_contract(scip, g, solnode, g->tail[e], g->head[e]) );
            count++;
         }
      }
   } while( count > 0 );

   return SCIP_OKAY;
}


/* removes parallel edges */
SCIP_RETCODE reduce_deleteMultiedges(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int nnodes = g->knots;
   int* count;

   assert(scip != NULL);
   assert(g != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &count, nnodes) );

   for( int k = 0; k < nnodes; k++ )
      count[k] = 0;

   for( int k = 0; k < nnodes; k++ )
   {
      int enext;
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         count[head]++;
      }

      for( int e = g->outbeg[k]; e != EAT_LAST; e = enext )
      {
         const int head = g->head[e];
         enext = g->oeat[e];

         if( count[head] > 1 )
         {
            graph_edge_del(scip, g, e, TRUE);
            return SCIP_ERROR;
         }
         count[head]--;

      }

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         assert(count[g->head[e]] == 0);
   }

   SCIPfreeBufferArray(scip, &count);

   return SCIP_OKAY;
}


/** basic reductions */
SCIP_RETCODE reduce_fixedConflicts(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   int*                  countnew            /**< pointer to number of new reductions (will be initially set to 0) */
   )
{
   int* hasharr;
   const int* fixednodes = graph_getFixpseudonodes(scip, g);
   const int nedges = graph_get_nEdges(g);
   const int arrsize = graph_pseudoAncestorsGetHashArraySize(g->pseudoancestors);
   const int nfixednodes = graph_getNfixpseudonodes(g);

   *countnew = 0;

   assert(scip && g && g->pseudoancestors);

   if( nfixednodes == 0 )
      return SCIP_OKAY;

   if( g->is_packed )
      return SCIP_OKAY;

   assert(fixednodes);

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, arrsize) );

   /* hash fixed nodes */
   for( int k = 0; k < nfixednodes; k++  )
   {
      const int node = fixednodes[k];
      assert(node >= 0 && node < arrsize);
      assert(hasharr[node] == 0 );

      hasharr[node] = 1;
   }

   for( int e = 0; e < nedges; e += 2 )
   {
      if( edgestate && edgestate[e] == EDGE_BLOCKED )
         continue;

      if( g->oeat[e] != EAT_FREE )
      {
         const int* pseudoancestors = graph_edge_getPseudoAncestors(g, e);
         const int nPseudoancestors = graph_edge_nPseudoAncestors(g, e);

         assert(g->oeat[e + 1] != EAT_FREE);
         assert(nPseudoancestors == 0 || pseudoancestors != NULL);

         for( int k = 0; k < nPseudoancestors; k++  )
         {
            const int ancestor = pseudoancestors[k];
            assert(ancestor >= 0 && ancestor < arrsize);

            if( hasharr[ancestor] != 0 )
            {
               graph_edge_del(scip, g, e, TRUE);
               (*countnew)++;

               SCIPdebugMessage("conflict deleted edge %d \n", e);
               break;
            }
         }
      }
   }

   /* unhash fixed nodes */
   for( int k = 0; k < nfixednodes; k++  )
   {
      const int node = fixednodes[k];
      assert(hasharr[node] == 1 );

      hasharr[node] = 0;
   }

   SCIPfreeCleanBufferArray(scip, &hasharr);


   return SCIP_OKAY;
}


/** root proximity terminal test (SAP) */
SCIP_RETCODE reduce_rpt(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   STP_Bool* nodes_isForbidden;
   SCIP_Real* dijkdist;
   const int nnodes = graph_get_nNodes(g);
   int* dijkedge;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   assert(count != NULL);
   assert(!graph_typeIsUndirected(g));
   assert(g->stp_type != STP_NWPTSPG || !graph_knotIsNWLeaf(g, g->source) || g->terms == 1);

   *count = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_isForbidden, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dijkdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dijkedge, nnodes) );

   graph_path_execX(scip, g, g->source, g->cost, dijkdist, dijkedge);

   for( int i = 0; i < nnodes; i++ )
      nodes_isForbidden[i] = FALSE;

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != g->source && g->grad[i] > 0 && !nodes_isForbidden[i] )
      {
         int e;
         const int e1 = dijkedge[i];
         const SCIP_Real pathcost = dijkdist[i];

         for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         {
            if( e == e1 )
               continue;

            if( GT(pathcost, g->cost[e]) )
               break;
         }

         if( e == EAT_LAST )
         {
            const int i1 = g->tail[e1];
            const int old = g->grad[i] + g->grad[i1] - 1;

            *fixed += g->cost[e1];

            assert(g->grad[i1] > 0);

            for( int e2 = g->outbeg[i]; e2 != EAT_LAST; e2 = g->oeat[e2] )
               nodes_isForbidden[g->head[e2]] = TRUE;

#ifdef SCIP_DEBUG
            SCIPdebugMessage("contracting (a) -> (b) with (a) root-dist=%f \n", pathcost);
            graph_edge_printInfo(g, e1);
            graph_knot_printInfo(g, i);
            graph_knot_printInfo(g, i1);
#endif

            SCIP_CALL( graph_knot_contractFixed(scip, g, NULL, e1, i1, i) );

            *count += old - g->grad[i1];
            SCIPdebugMessage("contract degree=%d\n", old - g->grad[i] - g->grad[i1]);
            assert(old - g->grad[i1] > 0);
         }

      }
   }

   SCIPfreeBufferArray(scip, &dijkedge);
   SCIPfreeBufferArray(scip, &dijkdist);
   SCIPfreeBufferArray(scip, &nodes_isForbidden);

   return SCIP_OKAY;
}


/** try to remove cute edge and prune one side of the graph */
SCIP_RETCODE reduce_cutEdgeTryPrune(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   cutedge,            /**< the edge */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Bool*            success             /**< could we prune the edge? */
)
{
   SCIP_Bool* RESTRICT nodes_visited;
   const int nnodes = graph_get_nNodes(g);
   SCIP_Bool terminalFound;
   const int tail = g->tail[cutedge];
   const int head = g->head[cutedge];

   assert(scip && success);
   assert(graph_edge_isInRange(g, cutedge));

   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_visited, nnodes) );

   for( int i = 0; i < nnodes; i++ )
      nodes_visited[i] = FALSE;

   *success = FALSE;

   /* first side */
   SCIP_CALL( cutEdgeProbe(scip, g, tail, head, nodes_visited, &terminalFound) );
   if( !terminalFound )
   {
      SCIP_CALL( cutEdgePrune(scip, tail, head, cutedge, g) );
      assert(g->grad[tail] == 0 || Is_term(g->term[tail]));
      *success = TRUE;
   }

   /* try second side? */
   if( !(*success) )
   {
      SCIP_CALL( cutEdgeProbe(scip, g, head, tail, nodes_visited, &terminalFound) );
      if( !terminalFound )
      {
         SCIP_CALL( cutEdgePrune(scip, head, tail, cutedge, g) );
         assert(g->grad[head] == 0 || Is_term(g->term[head]));
         *success = TRUE;
      }
   }

   SCIPfreeBufferArray(scip, &nodes_visited);

   return SCIP_OKAY;
}


/** identify non-leaf terminals and remove extensions */
void reduce_identifyNonLeafTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int nnodes = graph_get_nNodes(g);

   assert(scip && g);
   assert(graph_pc_isPc(g));
   assert(!g->extended);
   assert(g->prize && g->term2edge);
   assert(graph_valid(scip, g));
   assert(graph_pc_term2edgeIsConsistent(scip, g));

   /* make sure to not delete the last terminal connected to the root */
   if( g->stp_type == STP_PCSPG && g->grad[g->source] <= 2 )
   {
      return;
   }

   for( int k = 0; k < nnodes; ++k )
   {
      if( Is_term(g->term[k]) && !graph_pc_knotIsFixedTerm(g, k) && !graph_pc_termIsNonLeafTerm(g, k) )
      {
         assert(SCIPisGE(scip, g->prize[k], 0.0));
         assert(k != g->source);

         if( graph_pc_evalTermIsNonLeaf(scip, g, k) )
         {
            SCIPdebugMessage("transform term %d to non-leaf term \n", k);

            graph_pc_termToNonLeafTerm(scip, g, k, FALSE);

            assert(Is_term(g->term[k]));
            assert(graph_pc_knotIsNonLeafTerm(g, k));
         }
      }
   }

   assert(graph_pc_term2edgeIsConsistent(scip, g));
}



/* remove unconnected vertices */
SCIP_RETCODE reduce_unconnected(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int nnodes = g->knots;
   SCIP_Bool* nodevisited;

   assert(scip && g);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, nnodes) );

   SCIP_CALL( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = nnodes - 1; k >= 0 ; k-- )
   {
      if( !nodevisited[k] && (g->grad[k] > 0) )
      {
         if( Is_term(g->term[k]) )
         {
            assert(graph_pc_isPc(g));
            assert(graph_pc_termIsNonLeafTerm(g, k));

            continue;
         }

         graph_knot_del(scip, g, k, TRUE);
      }
   }

   SCIPfreeBufferArray(scip, &nodevisited);

   return SCIP_OKAY;
}


/* remove unconnected vertices for directed problems */
SCIP_RETCODE reduce_unconnectedForDirected(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int nnodes = graph_get_nNodes(g);
   SCIP_Bool* nodevisited;

   assert(scip && g);
   assert(!graph_typeIsUndirected(g));

   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, nnodes) );
   SCIP_CALL( graph_trail_costAware(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
   {
      if( !nodevisited[k] && (g->grad[k] > 0) )
      {
         assert(!Is_term(g->term[k]));
         graph_knot_del(scip, g, k, TRUE);
      }
   }

   SCIPfreeBufferArray(scip, &nodevisited);

   return SCIP_OKAY;
}


/** remove unconnected vertices and checks whether problem is infeasible  */
SCIP_RETCODE reduce_unconnectedInfeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             beVerbose,          /**< be verbose? */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Bool*            infeas              /**< is problem infeasible? */
)
{
   const int nnodes = graph_get_nNodes(g);
   SCIP_Bool* nodevisited;

   assert(scip && g && infeas);

   SCIP_CALL( SCIPallocBufferArray(scip, &nodevisited, nnodes) );

   *infeas = FALSE;

   SCIP_CALL( graph_trail_arr(scip, g, g->source, nodevisited) );

   for( int k = 0; k < nnodes; k++ )
      if( !nodevisited[k] && Is_term(g->term[k]) )
      {
         assert(k != g->source);
         *infeas = TRUE;
         if( beVerbose )
         {
            printf("terminal node %d (original index %d) lies in separate connected component \n", k, k + 1);
         }

         break;
      }

   for( int k = 0; k < nnodes; k++ )
   {
      if( !nodevisited[k] && (g->grad[k] > 0) )
      {
         while( g->inpbeg[k] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[k], TRUE);
      }
   }

   SCIPfreeBufferArray(scip, &nodevisited);

   return SCIP_OKAY;
}
