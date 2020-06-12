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


static
void getArticulationPoints(
   const GRAPH*          g,                  /**< graph data structure */
   int                   k,                  /**< vertex */
   int                   d,                  /**< depth */
   int*                  artpoints,          /**< artpoints list */
   int*                  nartpointsp,        /**< artpoints list counter */
   int*                  depth,              /**< depth */
   int*                  lowpoint,           /**< lowpoint */
   int*                  parent,             /**< parent */
   SCIP_Bool*            visited             /**< visited */
   )
{
   int nchildren = 0;
   SCIP_Bool isart = FALSE;

   visited[k] = TRUE;
   depth[k] = d;
   lowpoint[k] = d;

   for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int head = g->head[e];

      if( !g->mark[head] )
         continue;

      if( !visited[head] )
      {
         parent[head] = k;
         getArticulationPoints(g, head, d + 1, artpoints, nartpointsp, depth, lowpoint, parent, visited);
         nchildren++;
         if( lowpoint[head] >= depth[k] )
            isart = TRUE;
         lowpoint[k] = MIN(lowpoint[k], lowpoint[head]);
      }
      else if( head != parent[k] )
      {
         assert(depth[head] >= 0);
         lowpoint[k] = MIN(lowpoint[k], depth[head]);
      }
   }

   if( (parent[k] != -1 && isart) || (parent[k] == -1 && nchildren > 1) )
      artpoints[(*nartpointsp)++] = k;
}


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

   SCIP_CALL( reduceLevel0(scip, g) );

   return SCIP_OKAY;
}


/** probe */
static
SCIP_RETCODE cutEdgeProbe(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   start,              /**< the node to start from */
   int                   end,                /**< note to ignore */
   SCIP_Bool* RESTRICT   nodes_visited,
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
                  SCIP_CALL( graph_knot_replaceDeg2(scip, i, g, solnode, &conflict) );

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
      SCIP_CALL( reduceLevel0(scip, g) );
   }

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
   int i;
   int e;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;
   int* pnode;
   char rerun;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);

   rerun = TRUE;
   nnodes = g->knots;

   *count = 0;
   SCIPdebugMessage("Degree Test: ");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( i = 0; i < nnodes; i++ )
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
               if( i == g->source )
               {
                  e2 = flipedge(e1);

                  *fixed += g->cost[e2];

                  SCIP_CALL( graph_knot_contractFixed(scip, g, NULL, e2, i1, i) );
               }
               else
               {
                  *fixed += g->cost[e1];

                  SCIP_CALL( graph_knot_contractFixed(scip, g, NULL, e1, i1, i) );
               }
            }
            else
            {
               graph_edge_del(scip, g, e1, TRUE);
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

   for (i = 0; i < nnodes; i++)
      g->mark[i] = FALSE;

   g->mark[g->source] = TRUE;
   SCIP_CALL(SCIPqueueInsert(queue, &(g->source)));

   while (!SCIPqueueIsEmpty(queue))
   {
      pnode = (SCIPqueueRemove(queue));
      for (e = g->outbeg[*pnode]; e != EAT_LAST; e = g->oeat[e])
      {
         if( !g->mark[g->head[e]] && SCIPisLT(scip, g->cost[e], FARAWAY) )
         {
            g->mark[g->head[e]] = TRUE;
            SCIP_CALL(SCIPqueueInsert(queue, &(g->head[e])));
         }
      }
   }

   for (i = 0; i < nnodes; i++)
      if( !g->mark[i] )
         while (g->inpbeg[i] != EAT_LAST)
            graph_edge_del(scip, g, g->inpbeg[i], TRUE);

   /* delete all nodes that cannot reach a terminal other than the root by forward arcs (not using the root) */

   /* BFS */

   assert(SCIPqueueIsEmpty(queue));

   for( i = 0; i < nnodes; i++ )
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
      pnode = (SCIPqueueRemove(queue));
      for( e = g->inpbeg[*pnode]; e != EAT_LAST; e = g->ieat[e] )
      {
         if( !g->mark[g->tail[e]] && SCIPisLT(scip, g->cost[e], FARAWAY) )
         {
            g->mark[g->tail[e]] = TRUE;
            SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[e])) );
         }
      }
   }

   SCIPqueueFree(&queue);

   for( i = 0; i < nnodes; i++ )
      if( !g->mark[i] )
         while( g->inpbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[i], TRUE);

   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);
   assert(graph_valid(scip, g));

   return SCIP_OKAY;
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
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[g->head[e]]), g->ancestors[e], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->pcancestors[g->tail[e]]), g->ancestors[e], NULL) );
                  assert(0 && "currently not implemented");
               }
               else
               {
                  SCIP_CALL( graph_fixed_addEdge(scip, e, g) );
               }

            }

            SCIP_CALL( graph_knot_contract(scip, g, NULL, g->tail[e], g->head[e]) );
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


/** articulation points based reduction */
SCIP_RETCODE reduce_aritculations(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixedp,             /**< pointer to offset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   const int root = g->source;
   const int nnodes = g->knots;
   int* depth;
   int* lowpoint;
   int* parent;
   int* artpoints;
   int nartpoints;

   SCIP_Bool* visited;

   SCIP_CALL( SCIPallocBufferArray(scip, &artpoints, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &depth, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lowpoint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &parent, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &visited, nnodes) );

   for( int k = 0; k < nnodes; k++ )
   {
      visited[k] = FALSE;
      depth[k] = -1;
      parent[k] = -1;
      lowpoint[k] = -1;
   }
   nartpoints = 0;

   getArticulationPoints(g, root, 0, artpoints, &nartpoints, depth, lowpoint, parent, visited);

      printf("aritculation points found %d \n", nartpoints);

   SCIPfreeBufferArray(scip, &visited);
   SCIPfreeBufferArray(scip, &parent);
   SCIPfreeBufferArray(scip, &lowpoint);
   SCIPfreeBufferArray(scip, &depth);
   SCIPfreeBufferArray(scip, &artpoints);

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

   const int nedges = g->edges;
   const int nnodes = g->knots;
   const int nfixednodes = graph_get_nFixpseudonodes(scip, g);
   const int* fixednodes = graph_get_fixpseudonodes(scip, g);

   *countnew = 0;

   assert(scip && g && g->pseudoancestors);

   if( nfixednodes == 0 )
      return SCIP_OKAY;

   if( g->is_packed )
      return SCIP_OKAY;

   assert(fixednodes);

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &hasharr, nnodes) );

   /* hash fixed nodes */
   for( int k = 0; k < nfixednodes; k++  )
   {
      const int node = fixednodes[k];
      assert(node >= 0 && node < nnodes);
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
            assert(ancestor >= 0 && ancestor < nnodes);

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
   SCIP_Real pathcost;
   SCIP_Real* dijkdist;
   int e;
   int i1;
   int old;
   int nnodes;
   int* dijkedge;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   assert(count != NULL);
   assert(g->stp_type == STP_SAP);

   nnodes = g->knots;
   *count = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &dijkdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dijkedge, nnodes) );

   graph_path_execX(scip, g, g->source, g->cost, dijkdist, dijkedge);

   for( int i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != g->source && g->grad[i] > 0 )
      {
         const int e1 = dijkedge[i];
         pathcost = dijkdist[i];

         for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         {
            if( e == e1 )
               continue;

            if( SCIPisGT(scip, pathcost, g->cost[e]) )
               break;
         }
         if( e == EAT_LAST )
         {
            i1 = g->tail[e1];
            old = g->grad[i] + g->grad[i1] - 1;

            *fixed += g->cost[e1];

            SCIP_CALL( graph_knot_contractFixed(scip, g, NULL, e1, i1, i) );

            assert(old - g->grad[i1] > 0);
            *count += old - g->grad[i1];
            SCIPdebugMessage("contract %d\n", old - g->grad[i] - g->grad[i1]);
         }

      }
   }

   SCIPfreeBufferArray(scip, &dijkedge);
   SCIPfreeBufferArray(scip, &dijkdist);

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

  // int c = 0;
  // SCIP_CALL( reduce_aritculations(scip, g, NULL, &c ) );

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
