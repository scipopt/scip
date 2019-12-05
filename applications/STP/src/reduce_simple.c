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
 * This file implements basic reduction techniques for several Steiner problems.
 * All tests are described in "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "portab.h"
#include "scip/scip.h"

#ifndef NDEBUG
/** check whether problem has adjacent terminals */
static
SCIP_Bool hasAdjacentTerminals(
   const GRAPH*          g                   /**< graph data structure */
)
{
   for( int e = 0; e < g->edges; e++ )
   {
      if( g->oeat[e] != EAT_FREE )
      {
         const int tail = g->tail[e];
         const int head = g->head[e];
         if( Is_term(g->term[tail]) && Is_term(g->term[head]) && g->mark[head] && g->mark[tail] )
            return TRUE;
      }
   }

   return FALSE;
}
#endif

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

/** count numbers of chains */
static
int nChains(
   const GRAPH*          g                   /**< graph data structure */
   )
{
   int ccount = 0;

   assert(graph_pc_isMw(g));

   for( int e = 0; e < g->edges; e++ )
   {
      if( g->oeat[e] != EAT_FREE )
      {
         const int tail = g->tail[e];
         const int head = g->head[e];

         if( !Is_term(g->term[tail]) && !Is_term(g->term[head]) && g->grad[head] == 2 && g->grad[tail] == 2 )
            ccount++;
      }
   }

   return ccount;
}


/** traverse one side of a chain (MWCSP) */
static
SCIP_RETCODE traverseChain(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  length,             /**< pointer to store length of chain */
   int*                  final,              /**< pointer to store final vertex */
   int                   i,                  /**< start vertex */
   int                   i1,                 /**< first vertex */
   int                   i2,                 /**< last vertex */
   int                   e1                  /**< first edge */
   )
{
   IDX* ancestors = NULL;
   IDX* revancestors = NULL;
   SCIP_Real sum;
   int k;
   int e;

   assert(g != NULL);
   assert(scip != NULL);
   assert(length != NULL);
   assert(graph_pc_isMw(g));

   k = i1;
   e = e1;
   sum = 0.0;

   while( g->grad[k] == 2 && !Is_term(g->term[k]) && k != i2 )
   {
      assert(g->mark[k]);

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e], NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[flipedge(e)], NULL) );

      if( e != e1 )
         graph_edge_del(scip, g, e, TRUE);

      e = g->outbeg[k];
      sum += g->prize[k];
      (*length)++;

      if( e == flipedge(e1) )
         e = g->oeat[e];

      assert(e != EAT_LAST);
      assert(SCIPisLE(scip, g->prize[k], 0.0));

      k = g->head[e];
   }
   if( k != i1 )
   {
      int ne;

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e], NULL) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[flipedge(e)], NULL) );

      graph_edge_del(scip, g, e, TRUE);

      g->prize[i] += sum;
      ne = graph_edge_redirect(scip, g, e1, i, k, 1.0, TRUE, TRUE);

      if( ne != -1 )
      {
         e1 = ne;

         graph_edge_delHistory(scip, g, e1);

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e1]), ancestors, NULL) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(e1)]), revancestors, NULL) );

      }
      else
      {
         for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            if( g->head[e1] == k )
               break;
         assert(e1 != EAT_LAST);
      }

      SCIPintListNodeFree(scip, &(ancestors));
      SCIPintListNodeFree(scip, &(revancestors));

      if( SCIPisGE(scip, g->prize[k], 0.0) )
         g->cost[e1] = 0.0;
      else
         g->cost[e1] = -g->prize[k];

      assert(SCIPisLE(scip, g->prize[i], 0.0));
   }

   *final = k;

   return SCIP_OKAY;
}


/** is there no vertex of higher prize? */
static
SCIP_Bool isMaxprizeTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   i,                  /**< the terminal to be checked */
   SCIP_Real*            maxprize            /**< stores incumbent prize (can be updated) */
   )
{
   int t = -1;
   SCIP_Real max;
   const int nnodes = graph_get_nNodes(g);
   const int root = g->source;

   assert(i >= 0 && i < nnodes);
   assert(Is_term(g->term[i]) && g->prize[i] > 0.0);

   if( g->stp_type == STP_RPCSPG )
   {
      return (i == root);
   }

   max = *maxprize;

   if( max > g->prize[i] )
      return FALSE;

   max = -1.0;

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && k != root )
      {
         assert(g->mark[k]);

         if( g->prize[k] > max )
         {
            max = g->prize[k];
            t = k;
         }
         else if( t == i && g->prize[k] >= max )
         {
            t = k;
         }
      }
   }

   *maxprize = max;

   assert(t >= 0);

   SCIPdebugMessage("maxprize: %f (from %d) \n", g->prize[t], t );

   return (t == i);
}


/** reduces non-terminal of degree 1 for a (rooted) PC problem */
static
void pcReduceKnotDeg1(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i,                  /**< index of the terminal */
   SCIP_Bool*            rerun               /**< further eliminations possible? */
   )
{
   const int e1 = g->inpbeg[i];
   const int i1 = g->tail[e1];

   assert(!Is_term(g->term[i]));
   assert(e1 >= 0);
   assert(e1 == Edge_anti(g->outbeg[i]));
   assert(g->ieat[e1] == EAT_LAST);
   assert(g->oeat[g->outbeg[i]] == EAT_LAST);

   graph_edge_del(scip, g, e1, TRUE);
   SCIPdebugMessage("delete non-terminal of degree 1 %d\n ",  i);

   assert(g->grad[i] == 0);

   /* no more reductions possible from i1? */
   if( g->grad[i1] == 0 )
      *rerun = FALSE;
   else if( (i1 < i) && (g->grad[i1] < 3 || Is_term(g->term[i1])) )
      *rerun = TRUE;
}


/** reduces non-terminal of degree 2 for a (rooted) PC problem (by replacement) */
static
SCIP_RETCODE pcReduceKnotDeg2(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i,                  /**< index of the terminal */
   int*                  solnode,            /**< solution nodes or NULL */
   SCIP_Bool*            rerun               /**< further eliminations possible? */
   )
{
   const int e1 = g->outbeg[i];
   const int e2 = g->oeat[e1];
   const int i1 = g->head[e1];
   const int i2 = g->head[e2];

   assert(!Is_term(g->term[i]));

   SCIPdebugMessage("replace degree 2 non-terminal %d \n ", i);

   SCIP_CALL( graph_knot_replaceDeg2(scip, i, g, solnode) );

   if( (Is_term(g->term[i2]) && (i2 < i)) || (Is_term(g->term[i1]) && (i1 < i)) )
      *rerun = TRUE;

   return SCIP_OKAY;
}

/** adjust for a (rooted) PC or MW problem */
static
void pcmwReduceTerm0Prize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   assert(!g->extended);
   assert(Is_term(g->term[i]));
   assert(g->source != i);
   assert(SCIPisZero(scip, g->prize[i]));
   assert(!graph_pc_knotIsFixedTerm(g, i));

   if( graph_pc_termIsNonLeafTerm(g, i) )
   {
      assert(graph_pc_isPcMw(g));

      graph_pc_knotToNonTermProperty(g, i);
   }
   else
   {
      int t = UNKNOWN;
      int e2 = UNKNOWN;

      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int i1 = g->head[e];
         if( Is_pseudoTerm(g->term[i1]) && g->source != i1 )
            t = i1;
         else if( g->source == i1 )
            e2 = e;
      }

      assert(t != UNKNOWN);
      assert(g->head[g->term2edge[i]] == t);

      /* i is not a terminal anymore */
      graph_pc_knotToNonTermProperty(g, i);

      if( g->stp_type != STP_RPCSPG )
      {
         assert(e2 != UNKNOWN);
         graph_edge_del(scip, g, e2, TRUE);
      }

      /* delete artificial terminal */
      graph_pc_knotToNonTermProperty(g, t);
      graph_knot_del(scip, g, t, TRUE);
   }

   assert(!Is_term(g->term[i]));
}


/** try to eliminate a terminal of degree one */
static
SCIP_RETCODE pcmwReduceTermDeg1(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  count,              /**< pointer storing number of eliminated edges */
   int                   i,                  /**< the terminal to be checked */
   int                   iout,               /**< outgoing arc */
   SCIP_Bool*            rerun,              /**< further eliminations possible? */
   SCIP_Real*            maxprize            /**< stores incumbent prize (can be updated) */
   )
{
   int i1;
   int degsum;

   assert(scip && g && count);
   assert(Is_term(g->term[i]));

   if( isMaxprizeTerm(scip, g, i, maxprize) )
      return SCIP_OKAY;

   i1 = g->head[iout];

   if( SCIPisLE(scip, g->prize[i], g->cost[iout]) && g->stp_type != STP_MWCSP )
   {
      /* delete terminal */

      assert(!graph_pc_knotIsFixedTerm(g, i));

      if( (i1 < i) && (Is_term(g->term[i1]) || g->grad[i1] == 2 || g->grad[i1] == 3) )
         (*rerun) = TRUE;

      SCIPdebugMessage("Delete (degree 1) terminal %d \n", i);

      (*count) += graph_pc_deleteTerm(scip, g, i, offset);
   }
   else if( edgestate == NULL || edgestate[iout] != EDGE_BLOCKED )
   {
      /* contract terminal */

      (*rerun) = TRUE;
      assert(SCIPisGT(scip, g->prize[i], 0.0 ));

      if( g->stp_type == STP_MWCSP )
      {
         if( SCIPisLE(scip, g->prize[i], -g->prize[i1]) )
            *offset += g->prize[i];
         else
            *offset -= g->prize[i1];
      }
      else
      {
         *offset += g->cost[iout];
      }

      degsum = g->grad[i] + g->grad[i1];

      if( Is_term(g->term[i1]) && !graph_pc_termIsNonLeafTerm(g, i1) )
      {
         SCIP_CALL(graph_pc_contractEdge(scip, g, solnode, i1, i, i1));
         degsum -= g->grad[i1];
      }
      else
      {
         SCIP_CALL(graph_pc_contractEdge(scip, g, solnode, i, i1, i));
         degsum -= g->grad[i];
      }

      assert(degsum >= 1);

      if( g->stp_type == STP_MWCSP )
      {
         int e;
         int t = UNKNOWN;
         int e2 = UNKNOWN;
         if( SCIPisLE(scip, g->prize[i], 0.0) )
         {
            for (e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
            {
               i1 = g->head[e];
               if( Is_pseudoTerm(g->term[i1]) && g->source != i1 )
                  t = i1;
               else if( g->source == i1 )
                  e2 = e;
            }

            assert(t != UNKNOWN);
            assert(e2 != UNKNOWN);

            /* delete artificial terminal */
            graph_pc_knotToNonTermProperty(g, t);
            while (g->outbeg[t] != EAT_LAST)
            {
               e = g->outbeg[t];
               g->cost[e] = 0.0;
               g->cost[flipedge(e)] = 0.0;
               graph_edge_del(scip, g, e, TRUE);
               (*count)++;
            }

            assert(g->grad[t] == 0);

            /* i is not a terminal anymore */
            graph_pc_knotToNonTermProperty(g, i);
            graph_edge_del(scip, g, e2, TRUE);

            for (e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e])
               if( g->mark[g->tail[e]] )
                  g->cost[e] = -g->prize[i];

            for (e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
            {
               i1 = g->head[e];
               if( g->mark[i1] )
               {
                  if( !Is_term(g->term[i1]) )
                  {
                     g->cost[e] = -g->prize[i1];
                  }
                  else
                  {
                     g->cost[e] = 0.0;
                  }
               }
            }
         }
         else
         {
            for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
               if( g->mark[g->tail[e]] )
                  g->cost[e] = 0.0;

            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
               i1 = g->head[e];
               if( g->mark[i1] )
               {
                  if( !Is_term(g->term[i1]) )
                  {
                     assert(SCIPisLE(scip, g->prize[i1], 0.0));
                     g->cost[e] = -g->prize[i1];
                  }
                  else
                  {
                     assert(SCIPisGE(scip, g->prize[i1], 0.0));
                     g->cost[e] = 0.0;
                  }
               }
               else if( Is_pseudoTerm(g->term[i1]) && g->source != i1 )
               {
                  t = i1;
               }

            }
            assert(t != UNKNOWN);

            for (e = g->inpbeg[t]; e != EAT_LAST; e = g->ieat[e])
               if( g->tail[e] == g->source )
                  break;
            assert(e != EAT_LAST);
            g->cost[e] = g->prize[i];
         }
      }
      (*count) += degsum;
   }
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

   assert(scip && g && fixed && nelims);

   SCIPdebugMessage("Degree Test: ");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( int i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);

         if( g->grad[i] == 1 )
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
                  graph_edge_del(scip, g, e1, TRUE);
            }
            else
            {
               if( Is_term(g->term[i]) )
               {
                  *fixed += g->cost[e1];
                  SCIP_CALL( graph_fixed_addEdge(scip, e1, g) );
               }

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
                  SCIP_CALL( graph_knot_replaceDeg2(scip, i, g, solnode) );

                  elimscount++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if( Is_term(g->term[i1]) && Is_term(g->term[i2]) )
               {

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
                  *fixed += g->cost[e1];

                  SCIP_CALL( graph_fixed_addEdge(scip, e1, g) );
                  SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
                  elimscount++;
                  break;
               }
               if( Is_term(g->term[i2]) && !Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e2], g->cost[e1]) )
               {
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
               *fixed += g->cost[ett];

               SCIP_CALL( graph_fixed_addEdge(scip, ett, g) );
               SCIP_CALL( graph_knot_contractLowdeg2High(scip, g, solnode, i, g->head[ett]) );

               rerun = TRUE;
            }
         }
      }
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


/** basic reduction tests for the MWCS problem */
SCIP_RETCODE reduce_simple_mw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real maxprize = -1.0;
   const int nnodes = g->knots;
   int localcount = 0;

   SCIP_Bool rerun;
   SCIP_Bool contracted;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);
   assert(g->stp_type == STP_MWCSP);

   SCIPdebugMessage("MW degree test: \n");

   contracted = TRUE;
   while( contracted )
   {
      contracted = FALSE;

      /* contract adjacent positive vertices */
      for( int i = 0; i < nnodes; i++ )
      {
         int i1 = -1;
         int grad;
         SCIP_Bool hit;

         if( !Is_term(g->term[i]) || !(g->mark[i]) )
            continue;

         grad = g->grad[i];
         hit = FALSE;

         for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];

            if( Is_term(g->term[head]) && head != g->source )
            {
               assert(g->mark[head]);

               if( (g->grad[head] <= grad) )
               {
                  grad = g->grad[head];
                  i1 = head;
               }
               else if( head < i )
               {
                  hit = TRUE;
               }
            }
         }

         while( i1 >= 0 )
         {
            int i2 = -1;

            assert(g->mark[i1]);
            assert(g->grad[i1] > 0);
            assert(Is_term(g->term[i1]));

            grad = g->grad[i];
            hit = FALSE;
            for( int e = g->outbeg[i1]; e != EAT_LAST; e = g->oeat[e] )
            {
               const int head = g->head[e];
               if( Is_term(g->term[head]) && head != i && head != g->source )
               {
                  assert(g->mark[head]);

                  if( (g->grad[head] <= grad) )
                  {
                     i2 = head;
                     grad = g->grad[head];
                  }
                  else if( head < i )
                  {
                     hit = TRUE;
                  }
               }
            }

            SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i, i1, i));

            localcount++;

            i1 = i2;
         }
         if( hit )
            contracted = TRUE;
      }
   }

   /* contract adjacent 0 vertices */
   for( int i = 0; i < nnodes; i++ )
   {
      if( !(g->mark[i]) || !SCIPisZero(scip, g->prize[i]) )
         continue;

      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int i2 = g->head[e];

         if( g->mark[i2] && SCIPisGE(scip, g->prize[i2], 0.0) )
         {
            if( Is_term(g->term[i2]) )
            {
               SCIP_CALL(graph_pc_contractEdge(scip, g, solnode, i2, i, i2));
            }
            else
            {
               SCIP_CALL( graph_pc_contractNodeAncestors(scip, g, i2, i, flipedge_Uint(e)) );
               SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
            }

            localcount++;
            break;
         }
      }
   }

   SCIPdebugMessage("chains before: %d \n", nChains(g));

   rerun = TRUE;

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      /* main loop for remaining tests */
      for( int i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);
         if( !g->mark[i] || g->grad[i] == 0 )
            continue;

         assert(!Is_pseudoTerm(g->term[i]));

         /* non-positive vertex? */
         if( !Is_term(g->term[i]) )
         {
            if( g->grad[i] == 1 )
            {
               const int e1 = g->inpbeg[i];
               const int i1 = g->tail[e1];

               assert(e1 >= 0);
               assert(e1 == Edge_anti(g->outbeg[i]));
               assert(g->ieat[e1] == EAT_LAST);
               assert(g->oeat[g->outbeg[i]] == EAT_LAST);
               assert(SCIPisLE(scip, g->prize[i], 0.0));

               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete negative vertex of degree 1 (%d)\n ",  i);
               assert(g->grad[i] == 0);

               if( (i1 < i) && (g->grad[i1] < 3 || (g->grad[i1] == 3 && Is_term(g->term[i1]))) )
                  rerun = TRUE;

               localcount++;
               continue;
            }

            /* contract non-positive chains */
            if( g->grad[i] == 2 )
            {
               int f1 = -1;
               int f2 = -1;
               int length = 0;

               const int e1 = g->outbeg[i];
               const int e2 = g->oeat[e1];
               const int i1 = g->head[e1];
               const int i2 = g->head[e2];

               assert(e1 >= 0);
               assert(e2 >= 0);
               assert(i1 != i2);
               assert(g->mark[i1]);
               assert(g->mark[i2]);

               SCIP_CALL( traverseChain(scip, g, &length, &f1, i, i1, i2, e1) );
               SCIP_CALL( traverseChain(scip, g, &length, &f2, i, i2, i1, e2) );

               if( f1 == f2 )
               {
                  while( g->outbeg[i] != EAT_LAST )
                     graph_edge_del(scip, g, g->outbeg[i], TRUE);
               }
               else if( length > 0 )
               {
                  assert(g->grad[i] <= 2);

                  for( int e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
                     g->cost[e] = -g->prize[i];

                  localcount += length;
               }
            }
            continue;
         }

         /* node i is of positive weight (terminal): */

         /* terminal of 0-prize? */
         if( SCIPisLE(scip, g->prize[i], 0.0) )
         {
            pcmwReduceTerm0Prize(scip, g, i);
            localcount += 2;
            continue;
         }

         /* terminal of (real) degree 0? */
         if( g->grad[i] == 2 )
         {
            /* if terminal node i is not the one with the highest prize, delete */
            if( !isMaxprizeTerm(scip, g, i, &maxprize) )
            {
               SCIPdebugMessage("delete degree 0 term %d prize: %f count:%d\n ", i, g->prize[i], localcount);

               localcount += graph_pc_deleteTerm(scip, g, i, fixed);
            }
         }
         /* terminal of (real) degree 1? */
         else if( g->grad[i] == 3 )
         {
            int e;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] )
                  break;

            assert(e != EAT_LAST);
            assert(g->head[e] != g->source);

            if( !Is_term(g->term[g->head[e]]) )
            {
               SCIP_CALL( pcmwReduceTermDeg1(scip, NULL, g, fixed, NULL, count, i, e, &rerun, &maxprize) );
               continue;
            }
         }
      } /* i = 1 ... nnodes */
   } /* main loop */

   /* contract adjacent positive vertices */
   for( int i = 0; i < nnodes; i++ )
   {
      int i1;

      if( !(g->mark[i]) || !Is_term(g->term[i]) )
         continue;

      i1 = i;

      do
      {
         assert(g->mark[i1]);
         assert(g->grad[i1] > 0);
         assert(Is_term(g->term[i1]));

         contracted = FALSE;

         for( int e = g->outbeg[i1]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int i2 = g->head[e];
            if( g->mark[i2] && Is_term(g->term[i2]) )
            {
               SCIPdebugMessage("contract tt after (local) main loop %d->%d\n ", i1, i2);
               SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i1, i2, i1) );
               localcount++;
               contracted = TRUE;
               break;
            }
         }
      }
      while( contracted );
   }

   (*count) += localcount;
   SCIPdebugMessage("MW basic reduction package has deleted %d edges\n", *count);

   SCIPdebugMessage("chains after: %d \n", nChains(g));
   assert(!hasAdjacentTerminals(g));

   assert(graph_valid(scip, g));

   SCIP_CALL( level0(scip, g) );

   return SCIP_OKAY;
}


/** basic reductions for RPCSTP and PCSPG */
SCIP_RETCODE reduce_simple_pc(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  countnew,           /**< pointer to number of new reductions (will be initially set to 0) */
   int*                  countall,           /**< pointer to number of all reductions or NULL */
   int*                  solnode             /**< solution nodes */
   )
{
   int edges2[2];
   int nodes2[2];
   SCIP_Real maxprize = -1.0;
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool pc = (g->stp_type == STP_PCSPG);
   const SCIP_Bool checkstate = (edgestate != NULL);
   SCIP_Bool rerun = TRUE;

   assert(scip && fixed && countnew);
   assert(graph_pc_isPc(g));
   assert(!g->extended);

   *countnew = 0;

   SCIPdebugMessage("Degree Test: ");

   if( !pc )
      g->mark[g->source] = FALSE;

   /* main loop */
   while( rerun )
   {
      SCIP_Bool fixedterm = FALSE;
      rerun = FALSE;

      for( int i = 0; i < nnodes; i++ )
      {
         assert(!(g->mark[i] && Is_pseudoTerm(g->term[i])));

         if( (!g->mark[i] || g->grad[i] == 0) && !graph_pc_knotIsNonLeafTerm(g, i) )
         {
            assert(!Is_term(g->term[i]) || i == g->source);
            continue;
         }

         assert(!Is_pseudoTerm(g->term[i]) && i != g->source);

         if( !Is_term(g->term[i]) )
         {
            assert(0.0 == g->prize[i]);

            if( g->grad[i] == 1 )
            {
               pcReduceKnotDeg1(scip, g, i, &rerun);
               (*countnew)++;

               continue;
            }

            if( g->grad[i] == 2 )
            {
               SCIP_CALL( pcReduceKnotDeg2(scip, g, i, solnode, &rerun) );
               (*countnew)++;
            }

            continue;
         }

         /*
          * node i is a terminal:
          */

         assert(Is_term(g->term[i]));
         fixedterm = (!pc && graph_pc_knotIsFixedTerm(g, i));

         /* terminal of 0-prize? */
         if( SCIPisLE(scip, g->prize[i], 0.0) && i != g->source )
         {
            pcmwReduceTerm0Prize(scip, g, i);
            (*countnew) += 2;

            continue;
         }

         /* terminal of (real) degree 0? */
         if( graph_pc_realDegree(g, i, fixedterm) == 0 )
         {
            assert(!fixedterm);

            /* if terminal node i is not the one with the highest prize, delete */
            if( !isMaxprizeTerm(scip, g, i, &maxprize) )
            {
               SCIPdebugMessage("delete 0 term %d prize: %f countnew:%d\n ", i, g->prize[i], *countnew);

               (*countnew) += graph_pc_deleteTerm(scip, g, i, fixed);
            }
         }
         /* terminal of (real) degree 1? */
         else if( graph_pc_realDegree(g, i, fixedterm) == 1 )
         {
            int e;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] || (!pc && g->head[e] == g->source) )
                  break;

            assert(e != EAT_LAST);
            assert(g->head[e] != g->source || !pc);

            SCIP_CALL( pcmwReduceTermDeg1(scip, edgestate, g, fixed, solnode, countnew, i, e, &rerun, &maxprize) );
         }
         /* terminal of (real) degree 2? */
         else if( graph_pc_realDegree(g, i, fixedterm) == 2 )
         {
            if( !isMaxprizeTerm(scip, g, i, &maxprize) )
            {
               int edgecount = 0;
               for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               {
                  const int i1 = g->head[e];
                  if( g->mark[i1] || (!pc && i1 == g->source) )
                  {
                     assert(edgecount < 2);

                     edges2[edgecount] = e;
                     nodes2[edgecount++] = i1;
                  }
               }

               assert(edgecount >= 2);
               if( SCIPisLE(scip, g->prize[i], g->cost[edges2[0]]) && SCIPisLE(scip, g->prize[i], g->cost[edges2[1]]) )
               {
                  SINGLETONANS ancestors0;
                  SINGLETONANS ancestors1;
                  int newedge;
                  const int e0 = edges2[0];
                  const int e1 = edges2[1];
                  SCIP_Bool conflict;

                  SCIP_CALL( graph_singletonAncestors_init(scip, g, e0, &(ancestors0)) );
                  SCIP_CALL( graph_singletonAncestors_init(scip, g, e1, &(ancestors1)) );

                  assert(!fixedterm);

                  SCIPdebugMessage("delete - term - %d\n ", i);

                  SCIP_CALL( graph_edge_reinsert(scip, g, e0, nodes2[1], nodes2[0], g->cost[e0] + g->cost[e1] - g->prize[i],
                        i, &ancestors1, &ancestors0, &newedge, &conflict) );

                  (*countnew) += graph_pc_deleteTerm(scip, g, i, fixed);

                  graph_singletonAncestors_freeMembers(scip, &(ancestors0));
                  graph_singletonAncestors_freeMembers(scip, &(ancestors1));

                  if( conflict )
                  {
                     assert(newedge >= 0);
                     graph_edge_del(scip, g, newedge, TRUE);
                     (*countnew)++;
                  }
               }
            }
         }

         /* try to contract adjacent terminals */
         if( g->grad[i] > 0 )
         {
            SCIP_Real mincost = FARAWAY;
            int ett = UNKNOWN;

            for( int e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            {
               const int i1 = g->head[e1];

               if( !g->mark[i1] && (pc || i1 != g->source) )
                  continue;

               if( SCIPisLT(scip, g->cost[e1], mincost) )
               {
                  mincost = g->cost[e1];
                  if( Is_term(g->term[i1]) )
                     ett = e1;
               }
               else if( Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e1], mincost) )
               {
                  assert(SCIPisLT(scip, g->cost[e1], FARAWAY));
                  assert(SCIPisEQ(scip, g->cost[e1], mincost));
                  ett = e1;
               }
            }

            if( ett != UNKNOWN && SCIPisLE(scip, g->cost[ett], mincost) && SCIPisLE(scip, g->cost[ett], g->prize[i])
               && SCIPisLE(scip, g->cost[ett], g->prize[g->head[ett]]) )
            {
               const int i1 = g->head[ett];
               if( checkstate && edgestate[ett] == EDGE_BLOCKED )
                  continue;

               SCIPdebugMessage("contract tt %d->%d\n ", i, i1);
               assert(SCIPisLT(scip, mincost, FARAWAY));
               *fixed += g->cost[ett];
               (*countnew)++;

               SCIP_CALL( graph_pc_contractEdgeUnordered(scip, g, solnode, i, i1) );

               rerun = TRUE;
            }
         } /* contract adjacent terminals */
      } /* for i = 1, ..., nnodes */
   } /* main loops */

   if( !pc )
      g->mark[g->source] = TRUE;

   SCIPdebugMessage("degree test pc: %d nodes deleted\n", *countnew);

   reduce_identifyNonLeafTerms(scip, g);

   if( countall != NULL )
      (*countall) += (*countnew);

   assert(graph_valid(scip, g));

   SCIP_CALL( level0(scip, g) );

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

/** reduction test for PCSPG */
void
reduce_removeDeg0NonLeafTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offsetp             /**< pointer to offset value */
   )
{
   SCIP_Real maxprize = -1.0;
   const int nnodes = graph_get_nNodes(g);

   assert(scip && offsetp);
   assert(graph_pc_isPc(g));
   assert(!g->extended);

   for( int k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 && graph_pc_knotIsNonLeafTerm(g, k) && !isMaxprizeTerm(scip, g, k, &maxprize) )
      {
         graph_pc_deleteTerm(scip, g, k, offsetp);
      }
   }
}
