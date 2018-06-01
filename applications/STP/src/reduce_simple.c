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

/**@file   reduce_simple.c
 * @brief  several basic reductions for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements basic reduction techniques for several Steiner problems.
 * All tests are described in "A Generic Approach to Solving the Steiner Tree Problem and Variants" by Daniel Rehfeldt.
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"

#ifndef NDEBUG
/** check whether problem has adjacent terminals */
static
SCIP_Bool adjterms(
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

/** count numbers of chains */
static
unsigned nchains(
   const GRAPH*          g                   /**< graph data structure */
   )
{
   unsigned ccount = 0;
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

/** is there no vertex of higher prize? */
static
SCIP_Bool is_maxprize(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   i,                  /**< the terminal to be checked */
   SCIP_Real*            maxprize            /**< stores incumbent prize (can be updated) */
   )
{
   int t = -1;
   SCIP_Real max;

   assert(i >= 0 && Is_term(g->term[i]) && g->prize[i] > 0.0);

   if( g->stp_type == STP_RPCSPG && i != g->source )
      return FALSE;
   else if( g->stp_type == STP_RPCSPG && i == g->source )
      return TRUE;

   max = *maxprize;

   if( max > g->prize[i] )
      return FALSE;

   max = -1.0;

   for( int k = 0; k < g->knots; k++ )
   {
      if( Is_term(g->term[k]) && g->mark[k] && g->grad[k] > 0 )
      {
         assert(k != g->source);
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

   SCIPdebugMessage("maxprize: %f (from %d) \n", g->prize[t], t );
   return (t == i);
}

/** try to eliminate a terminal of degree one */
static
SCIP_RETCODE trydg1edgepc(
   SCIP*                 scip,               /**< SCIP data structure */
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

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(count  != NULL);
   assert(Is_term(g->term[i]));

   if( is_maxprize(scip, g, i, maxprize) )
      return SCIP_OKAY;

   i1 = g->head[iout];

   if( SCIPisLE(scip, g->prize[i], g->cost[iout]) && g->stp_type != STP_MWCSP )
   {
      /* delete terminal */

      if( (i1 < i) && (Is_term(g->term[i1]) || g->grad[i1] == 2 || g->grad[i1] == 3) )
         (*rerun) = TRUE;
      SCIPdebugMessage("Delete (degree 1) terminal %d \n", i);
      (*offset) += g->prize[i];
      (*count) += graph_pc_deleteTerm(scip, g, i);
   }
   else
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

      if( g->source == i1 )
      {
         assert(g->stp_type == STP_RPCSPG );

         if( g->pcancestors[i] != NULL )
         {
            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->pcancestors[i], NULL) );
            SCIPintListNodeFree(scip, &(g->pcancestors[i]));
         }
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[iout], NULL) );
         (*count) += graph_pc_deleteTerm(scip, g, i);
         return SCIP_OKAY;
      }

      degsum = g->grad[i] + g->grad[i1];

      SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i, i1, i) );

      degsum -= g->grad[i];

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
               if( Is_pterm(g->term[i1]) && g->source != i1 )
                  t = i1;
               else if( g->source == i1 )
                  e2 = e;
            }

            assert(t != UNKNOWN);
            assert(e2 != UNKNOWN);

            /* delete artificial terminal */
            graph_pc_knot2nonTerm(g, t);
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
            graph_pc_knot2nonTerm(g, i);
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
               else if( Is_pterm(g->term[i1]) && g->source != i1 )
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
      ne = graph_edge_redirect(scip, g, e1, i, k, 1.0, TRUE);

      if( ne != -1 )
      {
         e1 = ne;

         SCIPintListNodeFree(scip, &(g->ancestors[e1]));
         SCIPintListNodeFree(scip, &(g->ancestors[flipedge(e1)]));

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
      assert(SCIPisLE(scip, g->prize[i], 0.0) );
   }

   *final = k;

   return SCIP_OKAY;
}


/** adjust for a (rooted) PC or MW problem */
static
void adjust0term(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int t;
   int e2;

   assert(Is_term(g->term[i]));
   assert(g->source != i);
   assert(SCIPisZero(scip, g->prize[i]));

   t = UNKNOWN;
   e2 = UNKNOWN;

   if( !Is_term(g->term[i]) || SCIPisGT(scip, g->prize[i], 0.0) )
      return;

   for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int i1 = g->head[e];
      if( Is_pterm(g->term[i1]) && g->source != i1 )
         t = i1;
      else if( g->source == i1 )
         e2 = e;
   }

   assert(t != UNKNOWN);
   assert(g->head[g->term2edge[i]] == t);

   /* i is not a terminal anymore */
   graph_pc_knot2nonTerm(g, i);

   if( g->stp_type != STP_RPCSPG )
   {
      assert(e2 != UNKNOWN);
      graph_edge_del(scip, g, e2, TRUE);
   }

   /* delete artificial terminal */
   graph_pc_knot2nonTerm(g, t);

   while( g->outbeg[t] != EAT_LAST )
      graph_edge_del(scip, g, g->outbeg[t], TRUE);

   assert(g->grad[t] == 0);
}


/** contract edges of weight zero */
SCIP_RETCODE reduce_contractZeroEdges(
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
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e], NULL) );
            }
            SCIP_CALL( graph_knot_contract(scip, g, NULL, g->tail[e], g->head[e]) );
            count++;
         }
      }
   } while( count > 0 );

   return SCIP_OKAY;
}


/** basic reduction tests for the STP */
SCIP_RETCODE reduce_simple(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  solnode,            /**< node array to mark whether an node is part of a given solution (CONNECT),
                                                  or NULL */
   int*                  nelims,              /**< pointer to number of reductions */
   int*                  edgestate           /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   )
{
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;
   int rerun = TRUE;
   int done  = TRUE;
   int count = 0;
   SCIP_Bool checkstate = (edgestate != NULL);

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(nelims != NULL);

   nnodes = g->knots;

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
            e1  = g->outbeg[i];
            i1  = g->head[e1];

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
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1], NULL) );
               }

               SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
            }
            count++;

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
            e1 = g->outbeg[i];
            e2 = g->oeat[e1];
            i1 = g->head[e1];
            i2 = g->head[e2];

            assert(e1 >= 0);
            assert(e2 >= 0);

            do
            {
               done = TRUE;

               if( !Is_term(g->term[i]) )
               {
                  assert(EQ(g->cost[e2], g->cost[Edge_anti(e2)]));

                  g->cost[e1]            += g->cost[e2];
                  g->cost[Edge_anti(e1)] += g->cost[e2];

                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e1]), g->ancestors[flipedge(e2)], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(e1)]), g->ancestors[e2], NULL) );

                  SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
                  count++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if( Is_term(g->term[i1]) && Is_term(g->term[i2]) )
               {

                  if( SCIPisLT(scip, g->cost[e1], g->cost[e2]) )
                  {
                     *fixed += g->cost[e1];
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1], NULL) );
                     SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
                  }
                  else
                  {
                     *fixed += g->cost[e2];
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e2], NULL) );
                     SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
                  }
                  count++;

                  break;
               }
               if( Is_term(g->term[i1]) && !Is_term(g->term[i2]) && SCIPisLE(scip, g->cost[e1], g->cost[e2]) )
               {
                  *fixed += g->cost[e1];
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1], NULL) );
                  SCIP_CALL( graph_knot_contract(scip, g, solnode, i1, i) );
                  count++;
                  break;
               }
               if( Is_term(g->term[i2]) && !Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e2], g->cost[e1]) )
               {
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e2], NULL) );
                  *fixed += g->cost[e2];
                  SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
                  count++;
                  break;
               }
               done = FALSE;
            }
            while( FALSE );

            if (done
               && (((i1 < i) && (g->grad[i1] < 3))
                  || ((i2 < i) && (g->grad[i2] < 3))))
               rerun = TRUE;
         }
         if( Is_term(g->term[i]) && g->grad[i] > 2 && !checkstate )
         {
            SCIP_Real mincost = FARAWAY;
            int ett = UNKNOWN;
            for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            {
               i1 = g->head[e1];

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
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[ett], NULL) );
               SCIP_CALL( graph_knot_contractLowdeg2High(scip, g, solnode, i, g->head[ett]) );

               rerun = TRUE;
            }
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", count);
   assert(graph_valid(g));

   *nelims += count;
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
   int root;
   int nnodes;
   int* pnode;
   char rerun;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);

   root = g->source;
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
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1], NULL) );
               *fixed += g->cost[e1];
               SCIP_CALL( graph_knot_contract(scip, g, NULL, i1, i) );
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
   for( e = g->inpbeg[root]; e != EAT_LAST; e = g->ieat[e] )
      g->cost[e] = FARAWAY;

   /* delete all nodes not reachable from the root by forward arcs */

   /* BFS  */
   SCIP_CALL(SCIPqueueCreate(&queue, nnodes, 1.1));

   for (i = 0; i < nnodes; i++)
      g->mark[i] = FALSE;

   g->mark[root] = TRUE;
   SCIP_CALL(SCIPqueueInsert(queue, &(root)));

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
      if( Is_term(g->term[i]) && i != root )
      {
         g->mark[i] = TRUE;
         SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[g->outbeg[i]])) );
      }
      else
      {
         g->mark[i] = FALSE;
      }
   }

   g->mark[root] = TRUE;

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
   assert(graph_valid(g));

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
   int i;
   int e;
   int i1;
   int e1;
   int old;
   int root;
   int nnodes;
   int* dijkedge;

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);
   assert(count != NULL);

   root = g->source;
   nnodes = g->knots;
   *count = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &dijkdist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dijkedge, nnodes) );

   graph_path_execX(scip, g, root, g->cost, dijkdist, dijkedge);

   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(g->term[i]) && i != root && g->grad[i] > 0 )
      {
         e1 = dijkedge[i];
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

            SCIP_CALL(SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1], NULL));
            *fixed += g->cost[e1];
            SCIP_CALL(graph_knot_contract(scip, g, NULL, i1, i));

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
   const int root = g->source;
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

            if( Is_term(g->term[head]) && head != root )
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
               if( Is_term(g->term[head]) && head != i && head != root )
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
               SCIP_CALL( graph_pc_contractEdgeAncestors(scip, g, i2, i, flipedge_Uint(e)) );
               SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
            }

            localcount++;
            break;
         }
      }
   }

   SCIPdebugMessage("chains before: %d \n", nchains(g));

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

         assert(!Is_pterm(g->term[i]));

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
            adjust0term(scip, g, i);
            localcount += 2;
            continue;
         }

         /* terminal of (real) degree 0? */
         if( g->grad[i] == 2 )
         {
            /* if terminal node i is not the one with the highest prize, delete */
            if( !is_maxprize(scip, g, i, &maxprize) )
            {
               SCIPdebugMessage("delete degree 0 term %d prize: %f count:%d\n ", i, g->prize[i], localcount);
               (*fixed) += g->prize[i];
               localcount += graph_pc_deleteTerm(scip, g, i);
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
               SCIP_CALL( trydg1edgepc(scip, g, fixed, NULL, count, i, e, &rerun, &maxprize) );
               continue;
            }
         }
      }
   }

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

   SCIPdebugMessage("chains after: %d \n", nchains(g));
   assert(!adjterms(g));

   assert(graph_valid(g));

   SCIP_CALL( level0save(scip, g) );

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
   int root;
   int nnodes;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(g->stp_type == STP_DHCSTP);

   nnodes = g->knots;
   root = g->source;

   SCIPdebugMessage("basic HC test: \n");

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      /* delete incoming arcs of the root */
      e = g->inpbeg[root];
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
         if( Is_term(g->term[i]) && i != root )
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


/** basic reductions for RPCSTP and PCSPG */
SCIP_RETCODE reduce_simple_pc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  count,              /**< pointer to number of reductions */
   int*                  solnode,            /**< solution nodes */
   SCIP_Bool             contractroot        /**< contract vertices into root (for rooted prize-collecting) */
   )
{
   int edges2[2];
   int nodes2[2];
   SCIP_Real maxprize = -1.0;
   const int root = g->source;
   const int nnodes = g->knots;
   const SCIP_Bool pc = (g->stp_type == STP_PCSPG);
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG);
   assert(!g->extended);

   *count = 0;

   SCIPdebugMessage("Degree Test: ");

   if( !pc )
      g->mark[root] = FALSE;

   /* main loop */
   while( rerun )
   {
      rerun = FALSE;

      for( int i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);
         assert(!(g->mark[i] && Is_pterm(g->term[i])));

         /* last condition should never be true, but just in case ... */
         if( !g->mark[i] || g->grad[i] == 0 || Is_pterm(g->term[i]) )
            continue;

         if( !Is_term(g->term[i]) )
         {
            assert(SCIPisZero(scip, g->prize[i]));

            if( g->grad[i] == 1 )
            {
               const int e1 = g->inpbeg[i];
               const int i1 = g->tail[e1];

               assert(e1 >= 0);
               assert(e1 == Edge_anti(g->outbeg[i]));
               assert(g->ieat[e1] == EAT_LAST);
               assert(g->oeat[g->outbeg[i]] == EAT_LAST);

               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete non-terminal of degree 1 %d\n ",  i);
               (*count)++;

               assert(g->grad[i] == 0);

               /* the last node? */
               if( g->grad[i1] == 0 )
               {
                  rerun = FALSE;
                  break;
               }
               if( (i1 < i) && (g->grad[i1] < 3 || Is_term(g->term[i1])) )
                  rerun = TRUE;

               continue;
            }

            /* contract non terminals of degree 2 */
            if( g->grad[i] == 2 )
            {
               const int e1 = g->outbeg[i];
               const int e2 = g->oeat[e1];
               const int i1 = g->head[e1];
               const int i2 = g->head[e2];

               assert(e1 >= 0);
               assert(e2 >= 0);

               assert(g->mark[i1] || i1 == g->source);
               assert(g->mark[i2] || i2 == g->source);
               assert(SCIPisEQ(scip, g->cost[e2], g->cost[flipedge(e2)]));

               g->cost[e1]            += g->cost[e2];
               g->cost[flipedge(e1)]  += g->cost[e2];

               SCIPdebugMessage("contract non-terminals %d %d \n ", i2, i);
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e1]), g->ancestors[flipedge(e2)], NULL) );
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(e1)]), g->ancestors[e2], NULL) );

               SCIP_CALL( graph_knot_contract(scip, g, solnode, i2, i) );
               (*count)++;

               if( (Is_term(g->term[i2]) && (i2 < i)) || (Is_term(g->term[i1]) && (i1 < i)) )
                  rerun = TRUE;
            }
            continue;
         }

         /*
          * node i is a terminal:
          */

         /* terminal of 0-prize? */
         if( SCIPisLE(scip, g->prize[i], 0.0) && i != root )
         {
            adjust0term(scip, g, i);
            (*count) += 2;
            continue;
         }

         assert(Is_term(g->term[i]));

         /* terminal of (real) degree 0? */
         if( ( (g->grad[i] == 2 && pc) || (g->grad[i] == 1 && !pc) ) )
         {
            /* if terminal node i is node the one with the highest prize, delete */
            if( !is_maxprize(scip, g, i, &maxprize) )
            {
               SCIPdebugMessage("delete 0 term %d prize: %f count:%d\n ", i, g->prize[i], *count);
               (*fixed) += g->prize[i];
               (*count) += graph_pc_deleteTerm(scip, g, i);
            }
         }
         /* terminal of (real) degree 1? */
         else if( ( (g->grad[i] == 3 && pc) || (g->grad[i] == 2 && !pc) ) )
         {
            int e;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] || (!pc && g->head[e] == root) )
                  break;

            assert(e != EAT_LAST);
            assert(g->head[e] != root || !pc);

            SCIP_CALL( trydg1edgepc(scip, g, fixed, solnode, count, i, e, &rerun, &maxprize) );
         }
         /* terminal of (real) degree 2? */
         else if( ( (g->grad[i] == 4 && pc) || (g->grad[i] == 3 && !pc)  )  )
         {
            if( !is_maxprize(scip, g, i, &maxprize) )
            {
               int i2 = 0;
               for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               {
                  const int i1 = g->head[e];
                  if( g->mark[i1] || (!pc && i1 == root) )
                  {
                     assert(i2 < 2);

                     edges2[i2] = e;
                     nodes2[i2++] = i1;
                  }
               }

               assert(i2 >= 2);
               if( SCIPisLE(scip, g->prize[i], g->cost[edges2[0]]) && SCIPisLE(scip, g->prize[i], g->cost[edges2[1]]) )
               {
                  IDX* ancestors = NULL;
                  IDX* revancestors = NULL;
                  int n1;
                  const int e = edges2[0];
                  const int e1 = edges2[1];

                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[Edge_anti(e1)], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[Edge_anti(e)], NULL) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[e1], NULL) );
                  SCIPdebugMessage("delete - term - %d\n ", i);

                  /* contract edge */
                  n1 = graph_edge_redirect(scip, g, e, nodes2[1], nodes2[0], g->cost[e] + g->cost[e1] - g->prize[i], TRUE);

                  /* new edge inserted? */
                  if( n1 >= 0)
                  {
                     /* add ancestors */
                     SCIPintListNodeFree(scip, &(g->ancestors[n1]));
                     SCIPintListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors, NULL) );
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors, NULL) );
                  }
                  (*count) += graph_pc_deleteTerm(scip, g, i);
                  (*fixed) += g->prize[i];
                  SCIPintListNodeFree(scip, &(ancestors));
                  SCIPintListNodeFree(scip, &(revancestors));
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

               if( !g->mark[i1] && (pc || i1 != root || !contractroot) )
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
               SCIPdebugMessage("contract tt %d->%d\n ", i, i1);
               assert(SCIPisLT(scip, mincost, FARAWAY));
               *fixed += g->cost[ett];
               (*count)++;

               if( i1 == root )
               {
                  int j;
                  int e;

                  /* get edge from i to its artificial terminal */
                  for (e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e])
                     if( Is_pterm(g->term[g->head[e]]) && g->head[e] != root )
                        break;

                  assert(e != EAT_LAST);
                  SCIPdebugMessage("contract rt %d->%d \n", i, i1);

                  if( g->pcancestors[i] != NULL )
                  {
                     SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->pcancestors[i], NULL));
                     SCIPintListNodeFree(scip, &(g->pcancestors[i]));
                  }
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[ett], NULL));

                  /* artificial terminal to i */
                  j = g->head[e];
                  assert(!g->mark[j]);

                  /* delete edge and unmark artificial terminal */
                  graph_pc_knot2nonTerm(g, j);
                  graph_edge_del(scip, g, e, TRUE);

                  /* delete remaining incident edge of artificial terminal */
                  e = g->inpbeg[j];

                  assert(e != EAT_LAST);
                  assert(g->source == g->tail[e] || g->source == j);
                  assert(SCIPisEQ(scip, g->prize[i], g->cost[e]));

                  graph_edge_del(scip, g, e, TRUE);

                  assert(g->inpbeg[j] == EAT_LAST);

                  SCIP_CALL(graph_knot_contract(scip, g, solnode, i1, i));
                  graph_pc_knot2nonTerm(g, i);
               } /* i1 == root */
               else
               {
                  if( g->grad[i] >= g->grad[i1] )
                     SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i, i1, i) );
                  else
                     SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i1, i, i1) );
               }

               rerun = TRUE;
            }
         } /* contract adjacent terminals */
      } /* for i = 1, ..., nnodes */
   } /* main loops */

   if( !pc )
      g->mark[root] = TRUE;
   SCIPdebugMessage("degree test pc: %d nodes deleted\n", *count);

   assert(graph_valid(g));

   SCIP_CALL( level0save(scip, g) );

   return SCIP_OKAY;
}
