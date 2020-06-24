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

/**@file   reduce_pcsimple.c
 * @brief  several basic reductions for Steiner tree PC/MW problems
 * @author Daniel Rehfeldt
 *
 * This file implements basic reduction techniques for prize-collecting Steiner tree and maximum-weight connected subgraph.
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
#include "enumeration.h"
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

   if( graph_pc_isRootedPcMw(g) )
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


/** count numbers of chains */
static inline
int mwGetNchains(
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
SCIP_RETCODE mwTraverseChain(
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


/** contracts non-positive chains (path of degree 2 vertices) for (R)MWCS */
static
SCIP_RETCODE mwContractNonPositiveChain(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   i,                  /**< the start node */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< pointer to number of reductions */
   )
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
   assert(graph_pc_isMw(g));

   SCIP_CALL( mwTraverseChain(scip, g, &length, &f1, i, i1, i2, e1) );
   SCIP_CALL( mwTraverseChain(scip, g, &length, &f2, i, i2, i1, e2) );

   if( f1 == f2 )
   {
      while( g->outbeg[i] != EAT_LAST )
         graph_edge_del(scip, g, g->outbeg[i], TRUE);

      SCIPdebugMessage("deleting chain from vertex %d \n", i);

      *nelims += 1;
   }
   else if( length > 0 )
   {
      assert(g->grad[i] <= 2);

      for( int e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
         g->cost[e] = -g->prize[i];

      SCIPdebugMessage("deleting chain from vertex %d \n", i);

      *nelims += length;
   }

   return SCIP_OKAY;
}

/** contract 0-weight vertices for the MWCS problem */
static
SCIP_RETCODE mwContract0WeightVertices(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   const int nnodes = graph_get_nNodes(g);
   int nfixed_local = 0;

   assert(graph_pc_isMw(g));

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

            SCIPdebugMessage("contracted 0->positive %d->%d \n", i, i2);

            assert(g->grad[i] == 0);
            g->mark[i] = FALSE;

            nfixed_local++;
            break;
         }
      }
   }

   *nelims += nfixed_local;

   return SCIP_OKAY;
}


/** contract positive vertices for the MWCS problem */
static
SCIP_RETCODE mwContractTerminalsChainWise(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   const int* const gTerm = g->term;
   const int* const gOutbeg = g->outbeg;
   const int* const gOeat = g->oeat;
   const int* const gHead = g->head;
   const int nnodes = graph_get_nNodes(g);
   int nfixed_local = 0;
   SCIP_Bool contracted = TRUE;
   const SCIP_Bool isRooted = graph_pc_isRootedPcMw(g);

   assert(graph_pc_isMw(g));

   while( contracted )
   {
      contracted = FALSE;

      /* contract adjacent positive vertices */
      for( int i = 0; i < nnodes; i++ )
      {
         int i1 = -1;
         int grad;
         SCIP_Bool hit;

         if( !Is_term(gTerm[i]) || !(g->mark[i]) )
         {
            assert(LE(g->prize[i], 0.0) || graph_pc_knotIsDummyTerm(g, i));
            continue;
         }

         if( isRooted && graph_pc_knotIsFixedTerm(g, i) )
            continue;

         grad = g->grad[i];
         hit = FALSE;

         for( int e = gOutbeg[i]; e >= 0; e = gOeat[e] )
         {
            const int head = gHead[e];

            if( Is_term(gTerm[head]) && !graph_pc_knotIsFixedTerm(g, head) )
            {
               assert(g->mark[head]);
               assert(head != g->source);

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
            assert(Is_term(gTerm[i1]));

            grad = g->grad[i];
            hit = FALSE;
            for( int e = gOutbeg[i1]; e >= 0; e = gOeat[e] )
            {
               const int head = gHead[e];
               if( Is_term(gTerm[head]) && head != i && !graph_pc_knotIsFixedTerm(g, head) )
               {
                  assert(g->mark[head]);
                  assert(head != g->source);

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

            SCIPdebugMessage("contracted positive->positive %d->%d \n", i1, i);

            assert(g->grad[i1] == 0);
            g->mark[i1] = FALSE;
            i1 = i2;

            nfixed_local++;
         }
         if( hit )
            contracted = TRUE;
      }
   }

   *nelims += nfixed_local;

   return SCIP_OKAY;
}


/** contract positive vertices for the MWCS problem */
static
SCIP_RETCODE mwContractTerminalsSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   const int nnodes = graph_get_nNodes(g);
   int nfixed_local = 0;
   SCIP_Bool contracted = TRUE;
   const SCIP_Bool isRooted = graph_pc_isRootedPcMw(g);

   assert(graph_pc_isMw(g));


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
         assert(g->grad[i1] > 0 || (i1 == g->source) );
         assert(Is_term(g->term[i1]));

         contracted = FALSE;

         for( int e = g->outbeg[i1]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int i2 = g->head[e];
            if( g->mark[i2] && Is_term(g->term[i2]) && i2 != g->source )
            {
               /* NOTE necessary, because this would currently not work with the edge contraction */
               if( isRooted && graph_pc_knotIsFixedTerm(g, i2) && !graph_pc_knotIsFixedTerm(g, i1) )
                  continue;

               SCIP_CALL( graph_pc_contractEdge(scip, g, solnode, i1, i2, i1) );

               SCIPdebugMessage("contracted simple positive->positive %d->%d \n", i2, i1);

               assert(g->grad[i2] == 0);
               g->mark[i2] = FALSE;
               nfixed_local++;
               contracted = TRUE;
               break;
            }
         }
      }
      while( contracted );
   }

   *nelims += nfixed_local;

   return SCIP_OKAY;
}


/** try to eliminate a terminal of degree one */
static
SCIP_RETCODE mwReduceTermDeg1(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  nelims,             /**< pointer storing number of eliminated edges */
   int                   i,                  /**< the terminal to be checked */
   int                   iout,               /**< outgoing arc */
   SCIP_Bool*            rerun,              /**< further eliminations possible? */
   SCIP_Real*            maxprize            /**< stores incumbent prize (can be updated) */
   )
{
   assert(scip && g && nelims);
   assert(Is_term(g->term[i]));
   assert(g->tail[iout] == i);

   if( isMaxprizeTerm(scip, g, i, maxprize) )
      return SCIP_OKAY;

   /* can we contract? */
   if( edgestate == NULL || edgestate[iout] != EDGE_BLOCKED )
   {
      const int i1 = g->head[iout];
#ifndef NDEBUG
      const SCIP_Real newprize = g->prize[i] + g->prize[i1];
#endif

      (*rerun) = TRUE;
      assert(SCIPisGT(scip, g->prize[i], 0.0 ));
      assert(!Is_term(g->term[i1]));

      if( graph_pc_knotIsFixedTerm(g, i) )
      {
         *offset -= g->prize[i1];
      }
      else
      {
         if( SCIPisLE(scip, g->prize[i], -g->prize[i1]) )
            *offset += g->prize[i];
         else
            *offset -= g->prize[i1];
      }

      SCIP_CALL(graph_pc_contractEdge(scip, g, solnode, i, i1, i));
      assert(g->grad[i1] == 0);
      assert(EQ(g->prize[i1], 0.0));

      SCIPdebugMessage("degree-1 contraction: %d->%d new weight of %d: %f  \n", i1, i, i, g->prize[i]);

#ifndef NDEBUG
      assert(EQ(g->prize[i], newprize) || graph_pc_knotIsFixedTerm(g, i));

      for( int e = g->inpbeg[i]; e >= 0; e = g->ieat[e] )
      {
         const int tail = g->tail[e];
         if( g->mark[tail] )
         {
            assert(EQ(g->cost[e], MAX(-newprize, 0.0)));
         }
      }

      /* we also need to adapt the outgoing arcs because the contraction might have destroyed something  */
      for( int e = g->outbeg[i]; e >= 0; e = g->oeat[e] )
      {
         const int head = g->head[e];
         assert((!graph_pc_knotIsDummyTerm(g, head)) == g->mark[head]);

         if( g->mark[head] )
         {
            if( !Is_term(g->term[head]) )
            {
               assert(SCIPisLE(scip, g->prize[head], 0.0));
               assert(EQ(g->cost[e], -g->prize[head]));
            }
            else
            {
               assert(SCIPisGE(scip, g->prize[head], 0.0));
               assert(EQ(g->cost[e], 0.0));
            }
         }
      }
#endif

      (*nelims) += 1;
   }
   return SCIP_OKAY;
}


#if 1
/** tries to reduce full graph for a rooted PC problem */
static
void rpcTryFullReduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g                   /**< graph data structure */
)
{
   const int root = g->source;

   assert(g->stp_type == STP_RPCSPG);

   if( g->grad[root] == 0 )
      return;

   if( graph_pc_nProperPotentialTerms(g) == 0 && graph_pc_nFixedTerms(g) == 1 )
   {
      const int nnodes = g->knots;

      for( int i = 0; i < nnodes; i++ )
      {
         if( i != root )
         {
            graph_knot_del(scip, g, i, TRUE);
         }
      }
   }
}
#endif


/** reduces non-terminal of degree 1 for a (rooted) PC/MW problem */
static
void pcmwReduceKnotDeg1(
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
   assert(SCIPisLE(scip, g->prize[i], 0.0));

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
   SCIP_Bool conflict;

   assert(!Is_term(g->term[i]));

   SCIPdebugMessage("replace degree 2 non-terminal %d \n ", i);

   SCIP_CALL( graph_knot_replaceDeg2(scip, i, g, solnode, &conflict) );

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


/** check for possible enumeration */
static
void pcmwDeleteNonSolEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            result,             /**< the result */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< number of eliminations */
)
{
   const int nnodes = graph_get_nNodes(g);

   assert(scip && nelims);
   assert(*nelims >= 0);
   graph_mark(g);

   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      if( !g->mark[i] )
         continue;

      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         const int enext = g->oeat[e];
         const int head = g->head[e];

         if( !g->mark[head] )
         {
            e = enext;
            continue;
         }

         if( result[e] == CONNECT || result[flipedge(e)] == CONNECT )
         {
            e = enext;
            continue;
         }

         assert(LT(g->cost[e], FARAWAY));
         assert(LT(g->cost[flipedge(e)], FARAWAY));

         graph_edge_del(scip, g, e, TRUE);
         (*nelims)++;

         e = enext;
      }
   }
}


/** check for possible enumeration */
static
SCIP_RETCODE pcmwEnumerationTry(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims,             /**< number of eliminations */
   SCIP_Bool*            rerun               /**< further eliminations possible? */
   )
{
   if( g->grad[g->source] == 0 )
   {
      return SCIP_OKAY;
   }

   if( enumeration_isPossible(g) )
   {
      int nelims_new = 0;
      const int nedges = graph_get_nEdges(g);
      int* result = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );

      SCIP_CALL( enumeration_findSolPcMw(scip, g, result) );

      pcmwDeleteNonSolEdges(scip, result, g, &nelims_new);

      if( nelims_new > 0 )
      {
         *rerun = TRUE;
         *nelims += nelims_new;
      }

#ifdef SCIP_DEBUG
      graph_printInfo(g);
      printf("enumeriation elimination: %d \n", nelims_new);
#endif

      SCIPfreeBufferArray(scip, &result);
   }

   if( graph_pc_isRootedPcMw(g) )
      g->mark[g->source] = FALSE;

   return SCIP_OKAY;
}


/** try to eliminate a terminal of degree one */
static inline
SCIP_RETCODE pcReduceTermDeg1(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  nelims,             /**< pointer storing number of eliminated edges */
   int                   i,                  /**< the terminal to be checked */
   SCIP_Bool*            rerun,              /**< further eliminations possible? */
   SCIP_Real*            maxprize            /**< stores incumbent prize (can be updated) */
   )
{
   const SCIP_Bool isUnrooted = (g->stp_type == STP_PCSPG);
   int iout;

   assert(scip && nelims && offset);
   assert(Is_term(g->term[i]));

   if( isMaxprizeTerm(scip, g, i, maxprize) )
      return SCIP_OKAY;

   for( iout = g->outbeg[i]; iout != EAT_LAST; iout = g->oeat[iout] )
      if( g->mark[g->head[iout]] || (!isUnrooted && g->head[iout] == g->source) )
         break;

   assert(iout != EAT_LAST);
   assert(g->head[iout] != g->source || !isUnrooted);
   assert(g->tail[iout] == i);

   /* can we just delete the terminal? */
   if( SCIPisLE(scip, g->prize[i], g->cost[iout]) )
   {
      const int i1 = g->head[iout];

      assert(!graph_pc_knotIsFixedTerm(g, i));

      if( (i1 < i) && (Is_term(g->term[i1]) || g->grad[i1] == 2 || g->grad[i1] == 3) )
         (*rerun) = TRUE;

      SCIPdebugMessage("Delete (degree 1) terminal %d \n", i);

      (*nelims) += graph_pc_deleteTerm(scip, g, i, offset);
   }
   /* cannot be deleted, so contract terminal if not blocked */
   else if( edgestate == NULL || edgestate[iout] != EDGE_BLOCKED )
   {
      const int i1 = g->head[iout];
      int degsum = g->grad[i] + g->grad[i1];

      (*rerun) = TRUE;
      assert(SCIPisGT(scip, g->prize[i], 0.0 ));
      *offset += g->cost[iout];

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

      (*nelims) += degsum;
   }
   return SCIP_OKAY;
}


/* try to eliminate a terminal of degree 2 */
static inline
SCIP_RETCODE pcReduceTermDeg2(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  nelims,             /**< pointer storing number of eliminated edges */
   int                   i,                  /**< the terminal to be checked */
   SCIP_Bool*            rerun,              /**< further eliminations possible? */
   SCIP_Real*            maxprize            /**< stores incumbent prize (can be updated) */
)
{
   assert(SCIPisGT(scip, g->prize[i], 0.0));

   if( !isMaxprizeTerm(scip, g, i, maxprize) )
   {
      int edges2[2];
      int edgecount = 0;
      const SCIP_Bool isUnrooted = (g->stp_type == STP_PCSPG);

      for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int i1 = g->head[e];
         if( g->mark[i1] || (!isUnrooted && i1 == g->source) )
         {
            assert(edgecount < 2);
            edges2[edgecount++] = e;
         }
      }

      assert(edgecount == 2);

      /* is i a non-proper potential terminal? */
      if( SCIPisLE(scip, g->prize[i], g->cost[edges2[0]]) && SCIPisLE(scip, g->prize[i], g->cost[edges2[1]]) )
      {
         SINGLETONANS ancestors0;
         SINGLETONANS ancestors1;
         int newedge;
         const int e0 = edges2[0];
         const int e1 = edges2[1];
         const int i0 = g->head[e0];
         const int i1 = g->head[e1];
         SCIP_Bool conflict;

         SCIP_CALL( graph_singletonAncestors_init(scip, g, e0, &(ancestors0)) );
         SCIP_CALL( graph_singletonAncestors_init(scip, g, e1, &(ancestors1)) );

         assert(!graph_pc_knotIsFixedTerm(g, i));

         SCIPdebugMessage("delete degree 2 terminal %d\n ", i);

         SCIP_CALL( graph_edge_reinsert(scip, g, e0, i1, i0, g->cost[e0] + g->cost[e1] - g->prize[i],
               i, &ancestors1, &ancestors0, &newedge, &conflict) );

         (*nelims) += graph_pc_deleteTerm(scip, g, i, offset);

         graph_singletonAncestors_freeMembers(scip, &(ancestors0));
         graph_singletonAncestors_freeMembers(scip, &(ancestors1));

         if( conflict )
         {
            assert(newedge >= 0);
            graph_edge_del(scip, g, newedge, TRUE);
            (*nelims)++;
         }
      }
      else if( SCIPisLE(scip, g->prize[i], g->cost[edges2[0]]) || SCIPisLE(scip, g->prize[i], g->cost[edges2[1]]) )
      {
         /* i is semi-proper! */

         const int smalledge = (g->cost[edges2[0]] < g->cost[edges2[1]]) ? edges2[0] : edges2[1];

         assert(SCIPisGE(scip, g->prize[i], g->cost[smalledge]));

         if( edgestate == NULL || edgestate[smalledge] != EDGE_BLOCKED )
         {
            const int i1 = g->head[smalledge];
            SCIPdebugMessage("contract semi-proper degree 2 terminal %d with node %d \n", i, i1);

            (*rerun) = TRUE;
            *offset += g->cost[smalledge];

            if( Is_term(g->term[i1]) && !graph_pc_termIsNonLeafTerm(g, i1) )
            {
               SCIP_CALL(graph_pc_contractEdge(scip, g, solnode, i1, i, i1));
            }
            else
            {
               SCIP_CALL(graph_pc_contractEdge(scip, g, solnode, i, i1, i));
            }

            (*nelims) += 1;
         }
      }
   }

   return SCIP_OKAY;
}


/** try to contract terminal with adjacent one */
static inline
SCIP_RETCODE pcContractWithAdjacentTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< for propagation or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            offset,             /**< pointer to store the offset */
   int*                  solnode,            /**< solution nodes or NULL */
   int*                  nelims,             /**< pointer storing number of eliminated edges */
   int                   i,                  /**< the terminal to be checked */
   SCIP_Bool*            rerun               /**< further eliminations possible? */
)
{
   SCIP_Real mincost = FARAWAY;
   int edge_i2t = UNKNOWN;
   const SCIP_Bool isUnrooted = (g->stp_type == STP_PCSPG);

   for( int e1 = g->outbeg[i]; e1 >= 0; e1 = g->oeat[e1] )
   {
      const int i1 = g->head[e1];

      if( !g->mark[i1] && (isUnrooted || i1 != g->source) )
         continue;

      /* do we have to update edge_i2t? */
      if( SCIPisLT(scip, g->cost[e1], mincost) )
      {
         mincost = g->cost[e1];
         if( Is_term(g->term[i1]) )
            edge_i2t = e1;
      }
      else if( SCIPisLE(scip, g->cost[e1], mincost) )
      {
         /* we don't have to, but would it make sense to update edge_i2t? */
         if( Is_term(g->term[i1]) && SCIPisLE(scip, g->cost[e1], g->prize[i1])
                                  && SCIPisLE(scip, g->cost[e1], g->prize[i]) )
         {
            assert(SCIPisLT(scip, g->cost[e1], FARAWAY));
            assert(SCIPisEQ(scip, g->cost[e1], mincost));
            edge_i2t = e1;
         }
      }
   }

   if( edge_i2t != UNKNOWN && SCIPisLE(scip, g->cost[edge_i2t], mincost)
                           && SCIPisLE(scip, g->cost[edge_i2t], g->prize[i])
                           && SCIPisLE(scip, g->cost[edge_i2t], g->prize[g->head[edge_i2t]]) )
   {
      const int i1 = g->head[edge_i2t];
      const SCIP_Bool checkstate = (edgestate != NULL);

      if( checkstate && edgestate[edge_i2t] == EDGE_BLOCKED )
         return SCIP_OKAY;

      SCIPdebugMessage("contract tt %d->%d\n ", i, i1);
      assert(Is_term(g->term[i1]));
      assert(SCIPisLT(scip, mincost, FARAWAY));

      *offset += g->cost[edge_i2t];
      (*nelims)++;

      SCIP_CALL( graph_pc_contractEdgeUnordered(scip, g, solnode, i, i1) );

      *rerun = TRUE;
   }

   return SCIP_OKAY;
}


/** basic reduction tests for the MWCS problem */
SCIP_RETCODE reduce_simple_mw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int*                  solnode,            /**< array to indicate whether a node is part of the current solution (==CONNECT) */
   SCIP_Real*            fixed,              /**< pointer to offset value */
   int*                  nelims              /**< pointer to number of reductions */
   )
{
   SCIP_Real maxprize = -1.0;
   const int nnodes = graph_get_nNodes(g);
   int nelims_local = 0;
   SCIP_Bool rerun;
   const SCIP_Bool isRooted = graph_pc_isRootedPcMw(g);

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(nelims != NULL);
   assert(g->stp_type == STP_MWCSP || g->stp_type == STP_RMWCSP);

   SCIPdebugMessage("MW degree test: \n");

   graph_mark(g);

   SCIP_CALL( mwContractTerminalsChainWise(scip, g, fixed, solnode, &nelims_local) );
   SCIP_CALL( mwContract0WeightVertices(scip, g, solnode, &nelims_local) );

   SCIPdebugMessage("chains before: %d \n", mwGetNchains(g));

   rerun = TRUE;

   /* main loop */
   while( rerun )
   {
      SCIP_Bool fixedterm;
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
            assert(LE(g->prize[i], 0.0));

            if( g->grad[i] == 1 )
            {
               pcmwReduceKnotDeg1(scip, g, i, &rerun);
               nelims_local++;
            }
            else if( g->grad[i] == 2 )
            {
               SCIP_CALL( mwContractNonPositiveChain(scip, i, g, &nelims_local) );
            }
            continue;
         }

         /* node i is of positive weight (terminal): */

         assert(GE(g->prize[i], 0.0));

         /* terminal of 0-prize? */
         if( SCIPisLE(scip, g->prize[i], 0.0) )
         {
            pcmwReduceTerm0Prize(scip, g, i);
            nelims_local += 2;
            continue;
         }

         fixedterm = (isRooted && graph_pc_knotIsFixedTerm(g, i));

         /* terminal of (real) degree 0? */
         if( graph_pc_realDegree(g, i, fixedterm) == 0 )
         {
            /* if terminal node i is not the one with the highest prize, delete */
            if( !isMaxprizeTerm(scip, g, i, &maxprize) )
            {
               SCIPdebugMessage("delete degree 0 term %d prize: %f count:%d\n ", i, g->prize[i], nelims_local);

               nelims_local += graph_pc_deleteTerm(scip, g, i, fixed);
            }
         }
         /* terminal of (real) degree 1? */
         else if( graph_pc_realDegree(g, i, fixedterm) == 1 )
         {
            int e;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] )
                  break;

            if( e == EAT_LAST )
            {
               assert(isRooted);
               assert(i == g->source);
               continue;
            }

            assert(!graph_pc_knotIsDummyTerm(g, g->head[e]));

            if( !Is_term(g->term[g->head[e]]) )
            {
               SCIP_CALL( mwReduceTermDeg1(scip, NULL, g, fixed, NULL, nelims, i, e, &rerun, &maxprize) );
               continue;
            }
         }
      } /* i = 1 ... nnodes */

      SCIP_CALL( pcmwEnumerationTry(scip, g, &nelims_local, &rerun) );
   } /* main loop */

   SCIP_CALL( mwContractTerminalsSimple(scip, g, fixed, solnode, &nelims_local) );

   (*nelims) += nelims_local;
   SCIPdebugMessage("MW basic reduction package has done %d eliminations \n", *nelims);
   SCIPdebugMessage("chains after: %d \n", mwGetNchains(g));

   assert(!hasAdjacentTerminals(g));
   assert(graph_valid(scip, g));

   SCIP_CALL( reduceLevel0(scip, g) );

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
   SCIP_Real maxprize = -1.0;
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool isRooted = (g->stp_type == STP_RPCSPG);
   SCIP_Bool rerun = TRUE;

   assert(scip && fixed && countnew);
   assert(graph_pc_isPc(g));
   assert(!g->extended);

   *countnew = 0;

   SCIPdebugMessage("Degree Test: ");

   if( isRooted )
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
               pcmwReduceKnotDeg1(scip, g, i, &rerun);
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
         fixedterm = (isRooted && graph_pc_knotIsFixedTerm(g, i));

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
            SCIP_CALL( pcReduceTermDeg1(scip, edgestate, g, fixed, solnode, countnew, i, &rerun, &maxprize) );
         }
         /* terminal of (real) degree 2? */
         else if( graph_pc_realDegree(g, i, fixedterm) == 2 )
         {
            SCIP_CALL( pcReduceTermDeg2(scip, edgestate, g, fixed, solnode, countnew, i, &rerun, &maxprize) );
         }

         /* try to contract adjacent terminals */
         if( g->grad[i] > 0 )
         {
            SCIP_CALL( pcContractWithAdjacentTerm(scip, edgestate, g, fixed, solnode, countnew, i, &rerun) );
         }

      } /* for i = 1, ..., nnodes */

      SCIP_CALL( pcmwEnumerationTry(scip, g, countnew, &rerun) );
   } /* main loop */

   if( isRooted )
      g->mark[g->source] = TRUE;

   SCIPdebugMessage("degree test pc: %d nodes deleted\n", *countnew);

#if 1
   if( isRooted )
      rpcTryFullReduce(scip, g);
#endif

   reduce_identifyNonLeafTerms(scip, g);

   if( countall != NULL )
      (*countall) += (*countnew);

   assert(graph_valid(scip, g));

   SCIP_CALL( reduceLevel0(scip, g) );

   return SCIP_OKAY;
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
