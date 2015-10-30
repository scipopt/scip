/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dirreduce.c
 * @brief  several simple reductions for Steiner tree problems
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 *
 *
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

/** delete a terminal for a (rooted) prize collecting problem */
int deleteterm(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< index of the terminal */
   )
{
   int e;
   int t;
   int i1;
   int count;

   assert(g != NULL);
   assert(scip != NULL);
   assert(Is_term(g->term[i]));

   t = UNKNOWN;
   count = g->grad[i] + 2;

   /* delete terminal */

   graph_knot_chg(g, i, -1);
   g->mark[i] = FALSE;

   while( (e = g->outbeg[i]) != EAT_LAST )
   {
      i1 = g->head[e];

      if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
         t = g->head[e];
      graph_edge_del(scip, g, e, TRUE);
   }

   assert(t != UNKNOWN);

   /* delete artifical terminal */

   graph_knot_chg(g, t, -1);

   while( g->outbeg[t] != EAT_LAST )
      graph_edge_del(scip, g, g->outbeg[t], TRUE);


   return count;
}

static
SCIP_Bool maxprize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   int                   i                   /**< the terminal to be checked */
   )
{
   int k;
   int t = -1;
   SCIP_Real max = -1.0;

   for( k = 0; k < g->knots; k++ )
   {
      if( Is_term(g->term[k]) && g->mark[k] && g->grad[k] > 0 )
      {
	 assert(k != g->source[0]);
	 if( SCIPisGT(scip, g->prize[k], max) )
	 {
            max = g->prize[k];
            t = k;
	 }
	 else if( t == i && SCIPisGE(scip, g->prize[k], max) )
	 {
            t = k;
	 }
      }
   }

   SCIPdebugMessage("maxprize: %f (from %d) \n", g->prize[t], t );
   return (t == i);
}

static
SCIP_RETCODE trydg1edgepc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real* offset,
   int* count,
   int i,
   int iout,
   SCIP_Bool* rerun
   )
{
   int i1;
   int degsum;
   assert(scip  != NULL);
   assert(g      != NULL);
   assert(count != NULL);
   assert(Is_term(g->term[i]));

   if( maxprize(scip, g, i) )
      return SCIP_OKAY;

   i1 = g->head[iout];

   if( SCIPisLE(scip, g->prize[i], g->cost[iout]) && g->stp_type != STP_MAX_NODE_WEIGHT )
   {
      if( (i1 < i) && (Is_term(g->term[i1]) || g->grad[i1] == 2 || g->grad[i1] == 3) )
         (*rerun) = TRUE;
      SCIPdebugMessage("Delete (degree 1) terminal %d \n", i);
      (*offset) += g->prize[i];
      *count += deleteterm(scip, g, i);
   }
   else
   {
      (*rerun) = TRUE;
      assert(SCIPisGT(scip, g->prize[i], 0.0 ));
#if 0
      printf("prize before p1: %f, p2: %f off %f \n", g->prize[i], g->prize[i1], *offset);

      int et;
      for( et = g->outbeg[i1]; et != EAT_LAST; et = g->oeat[et] )
         printf("i1: outedge: %d %d %f, %f \n", g->tail[et], g->head[et], g->cost[et], g->cost[flipedge(et)]);


      for( et = g->outbeg[i]; et != EAT_LAST; et = g->oeat[et] )
         printf("i: outedge: %d %d %f, %f \n", g->tail[et], g->head[et], g->cost[et], g->cost[flipedge(et)]);
#endif
      if( g->stp_type == STP_MAX_NODE_WEIGHT )
      {
	 if( SCIPisLT(scip, g->prize[i], -g->prize[i1]) )
            *offset += g->prize[i];
	 else
            *offset -= g->prize[i1];
      }
      else
      {
         *offset += g->cost[iout];
      }

      degsum = g->grad[i] + g->grad[i1];
      SCIP_CALL( graph_knot_contractpc(scip, g, i, i1, i) );
      degsum = degsum - g->grad[i];
      assert(degsum >= 1);
      if( g->stp_type == STP_MAX_NODE_WEIGHT )
      {
	 int e;
	 int e2 = UNKNOWN;
	 int t = UNKNOWN;

	 if( SCIPisLT(scip, g->prize[i], 0.0 ) )
	 {
            i1 = UNKNOWN;
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
               i1 = g->head[e];
               if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
                  t = i1;
	       else if( g->source[0] == i1 )
                  e2 = e;
            }
            assert(e2 != UNKNOWN);
            assert(t != UNKNOWN);

            /* delete artifical terminal */
            graph_knot_chg(g, t, -1);

            while( g->outbeg[t] != EAT_LAST )
            {
	       e = g->outbeg[t];
	       g->cost[e] = 0.0;
	       g->cost[flipedge(e)] = 0.0;
               graph_edge_del(scip, g, e, TRUE);
               count++;
            }
            assert(g->grad[t] == 0);
	    /* i is not a terminal anymore */
	    graph_knot_chg(g, i, -1);
	    graph_edge_del(scip, g, e2, TRUE);

            for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
               if( g->mark[g->tail[e]] )
                  g->cost[e] = -g->prize[i];

            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
	    {
               i1 = g->head[e];
               if( g->mark[i1]  )
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
               if( g->mark[i1]  )
               {
                  if( !Is_term(g->term[i1]) )
                  {
                     assert(SCIPisLE(scip, g->prize[i1], 0.0 ));
                     g->cost[e] = -g->prize[i1];
                  }
                  else
                  {
                     assert(SCIPisGE(scip, g->prize[i1], 0.0 ));
                     g->cost[e] = 0.0;
                  }
               }
               else if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
	       {
		  t = i1;
	       }

	    }
	    assert(t != UNKNOWN);

	    for( e = g->inpbeg[t]; e != EAT_LAST; e = g->ieat[e] )
               if( g->tail[e] == g->source[0] )
                  break;
            assert(e != EAT_LAST);
	    g->cost[e] = g->prize[i];
	 }
      }
      *count += degsum;
   }
   return SCIP_OKAY;
}





#if 0

/** contract a chain */
static
SCIP_RETCODE traverseChain2(
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
   int k;
   int e;
   int sum;

   assert(g != NULL);
   assert(scip != NULL);
   assert(length != NULL);

   k = i1;
   e = e1;
   return SCIP_OKAY;
}

#endif

/** contract a chain */
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
   sum = 0;
   while( g->grad[k] == 2 && !Is_term(g->term[k]) && k != i2 )
   {
      assert(g->mark[k]);

      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e]) );
      SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[flipedge(e)]) );

      if( e != e1 )
         graph_edge_del(scip, g, e, TRUE);
      sum += g->prize[k];
      e = g->outbeg[k];

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

      graph_edge_del(scip, g, e, TRUE);
      g->prize[i] += sum;
      ne = graph_edge_redirect(scip, g, e1, i, k, 1.0);
      if( ne != -1 )
      {
         e1 = ne;

         SCIPintListNodeFree(scip, &(g->ancestors[e1]));
         SCIPintListNodeFree(scip, &(g->ancestors[flipedge(e1)]));

         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[e1]), ancestors) );
         SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->ancestors[flipedge(e1)]), revancestors) );

      }
      else
      {
	 for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            if(  g->head[e1] == k )
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

SCIP_RETCODE degree_test_dir(
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

   root = g->source[0];
   rerun = TRUE;
   nnodes = g->knots;

   *count = 0;
   SCIPdebugMessage("Degree Test: ");

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
#if 0
            if( !g->mark[i1] )
               continue;
#endif
            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->outbeg[i]));
            assert(g->ieat[e1] == EAT_LAST);
            assert(g->oeat[g->outbeg[i]] == EAT_LAST);

            if( Is_term(g->term[i]) )
            {
               SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
               *fixed += g->cost[e1];
               SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
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

         /* Note on costs in the directed graph
          * g->outbeg[i] is the outgoing directed edge
          * g->inpbeg[i] is the incoming directed edge */
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
               /* both the edges are outgoing from node i
                * need to ensure that the flow of the edge costs is correct
                * Edge_anti(e2) -> e1 and Edge_anti(e1) -> e2  */
               if( (!Is_term(g->term[i2]) && !Is_term(g->term[i1])) )
               {
                  g->cost[e1] += g->cost[Edge_anti(e2)];
                  g->cost[Edge_anti(e1)] += g->cost[e2];
		  if( SCIPisGT(scip, g->cost[e1], FARAWAY) )
                     g->cost[e1] = FARAWAY;
		  if( SCIPisGT(scip, g->cost[Edge_anti(e1)], FARAWAY) )
                     g->cost[Edge_anti(e1)] = FARAWAY;
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
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

   /* delete all arcs in not connected to a terminal other than the root by forward arcs */

   /* BFS until all terminals are reached */
   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );

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
         if( !g->mark[g->tail[e]] )
         {
            g->mark[g->tail[e]] = TRUE;
            SCIP_CALL( SCIPqueueInsert(queue, &(g->tail[e])) );
         }
      }
   }

   SCIPqueueFree(&queue);

   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
      {
	 while( g->inpbeg[i] != EAT_LAST )
	 {
	    printf("remove edge to node %d \n", i);
            graph_edge_del(scip, g, g->inpbeg[i], TRUE);
	 }
      }
      for( e = g->outbeg[i]; e != EAT_LAST && 0; e = g->oeat[e] )
      {
	 if( SCIPisGE(scip, g->cost[e], FARAWAY) &&  SCIPisGE(scip, g->cost[flipedge(e)], FARAWAY) )
	 {
	    printf("remove high cost edge to node %d \n", i);
            graph_edge_del(scip, g, e, TRUE);
	 }
      }
   }

   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);
   assert(graph_valid(g));

   return SCIP_OKAY;
}


/** Root proximity terminal */
SCIP_RETCODE rptReduction(
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

   root = g->source[0];
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

            SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
            *fixed += g->cost[e1];
            SCIP_CALL( graph_knot_contract(scip, g, i1, i) );

            assert(old - g->grad[i1] > 0);
            *count += old - g->grad[i1];
            printf("contract %d\n", old - g->grad[i] - g->grad[i1]);
	 }

      }
   }

   SCIPfreeBufferArray(scip, &dijkedge);
   SCIPfreeBufferArray(scip, &dijkdist);

   return SCIP_OKAY;
}

/** adjacent neighbourhood reduction */
SCIP_RETCODE ansReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   int k;
   int j;
   int e;
   int e2;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighbourhood of all nodes */
   for( k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
	 continue;

      /* mark adjacent vertices and k*/
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
	 marked[g->head[e]] = TRUE;
      marked[k] = TRUE;

      if( SCIPisLT(scip, g->prize[k], 0.0) )
         min = g->prize[k];
      else
         min = 0.0;

      maxgrad = g->grad[k];

      /* check all neighbours of k */
      e = g->outbeg[k];
      while( e != EAT_LAST )
      {
	 j = g->head[e];
	 e = g->oeat[e];
	 /* valid candidate? */
	 if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) )
	 {
	    for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
	       if( !marked[g->head[e2]] )
                  break;

	    /* neighbours of j subset of those of k? */
	    if( e2 == EAT_LAST )
	    {
	       while( g->outbeg[j] != EAT_LAST )
               {
                  e2 = g->outbeg[j];
                  (*count)++;
                  graph_edge_del(scip, g, e2, TRUE);
               }
               g->mark[j] = FALSE;
	       marked[j] = FALSE;
	    }
	 }
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
	 marked[g->head[e]] = FALSE;
      marked[k] = FALSE;

      for( k = 0; k < nnodes; k++ )
         assert(marked[k] == FALSE);

      for( k = 0; k < nnodes; k++ )
         if( marked[k] )
            printf("FAIL %d \n \n", k );
   }

   return SCIP_OKAY;
}


/** advanced adjacent neighbourhood reduction */
SCIP_RETCODE ansadvReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   int* neighbarr;
   int k;
   int j;
   int e;
   int k2;
   int e2;
   int l;
   int nn;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &neighbarr, 15) );

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   /* check neighbourhood of all nodes */
   for( k = 0; k < nnodes; k++ )
   {

      if( (!(g->mark[k])) || (g->grad[k] < 2) || Is_term(g->term[k]) )
         continue;

      maxgrad = 0;
      nn = 0;

      /* mark adjacent vertices and k*/
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
	 j = g->head[e];
	 marked[j] = TRUE;
	 if( Is_term(g->term[j]) && nn < 15 )
            neighbarr[nn++] = j;
      }

      marked[k] = TRUE;

      maxgrad = g->grad[k];
      for( l = 0; l < nn; l++ )
      {
         for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = TRUE;
         maxgrad += g->grad[neighbarr[l]];
      }

      assert(SCIPisLE(scip, g->prize[k], 0.0));

      min = g->prize[k];

      /* check all neighbours of k */
      e = g->outbeg[k];
      while( e != EAT_LAST )
      {
	 j = g->head[e];
	 e = g->oeat[e];
	 /* valid candidate? */
	 if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) )
	 {
	    for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
	       if( !marked[g->head[e2]] )
                  break;

	    /* neighbours of j subset of those of k? */
	    if( e2 == EAT_LAST )
	    {
	       while( g->outbeg[j] != EAT_LAST )
               {
                  e2 = g->outbeg[j];
                  (*count)++;
                  graph_edge_del(scip, g, e2, TRUE);
               }
               g->mark[j] = FALSE;
	       marked[j] = FALSE;
	    }
	 }
      }

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         marked[g->head[e]] = FALSE;

      for( l = 0; l < nn; l++ )
         for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = FALSE;

      marked[k] = FALSE;

      for( k2 = 0; k2 < nnodes; k2++ )
         assert(marked[k2] == FALSE);
   }

   SCIPfreeBufferArray(scip, &neighbarr);
   return SCIP_OKAY;
}


/** advanced adjacent neighbourhood reduction */
SCIP_RETCODE ansadv2Reduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  marked,             /**< nodes array */
   int*                  count               /**< pointer to number of reductions */
   )
{
   SCIP_Real min;
   SCIP_Real maxprize;
   int* neighbarr;
   int k;
   int j;
   int e;
   int k2;
   int e2;
   int l;
   int run;
   int nn;
   int nnodes;
   int maxgrad;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &neighbarr, 15) );

   /* unmark all nodes */
   for( k = 0; k < nnodes; k++ )
      marked[k] = FALSE;

   for( run = 0; run < 2; run++ )
   {
      /* check neighbourhood of all nodes */
      for( k = 0; k < nnodes; k++ )
      {
         if( (!(g->mark[k])) || (g->grad[k] < 2) || Is_term(g->term[k]) )
            continue;

         maxgrad = 0;
         nn = 0;

         maxprize = g->prize[k];

         k2 = UNKNOWN;

         /* mark adjacent vertices and k*/
         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
         {
            j = g->head[e];
            marked[j] = TRUE;
            if( Is_term(g->term[j]) && nn < 14 )
            {
               neighbarr[nn++] = j;
            }
            else if( SCIPisGT(scip, g->prize[j], maxprize) && g->mark[j] )
            {
               maxprize = g->prize[j];
               k2 = j;
            }
         }

         marked[k] = TRUE;

         maxgrad = g->grad[k];

         if( run == 0 && k2 != UNKNOWN )
            neighbarr[nn++] = k2;

         for( l = 0; l < nn; l++ )
         {
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
            {
               j = g->head[e];
               if( run == 1 && g->mark[j] && !Is_term(g->term[j]) && SCIPisGT(scip, g->prize[j], maxprize) )
               {
                  maxprize = g->prize[j];
                  k2 = j;
               }
               marked[j] = TRUE;
            }
            maxgrad += g->grad[neighbarr[l]];
         }

         if( run == 1 && k2 != UNKNOWN )
         {
            neighbarr[nn++] = k2;
            for( e = g->outbeg[k2]; e != EAT_LAST; e = g->oeat[e] )
               marked[g->head[e]] = TRUE;
            maxgrad += g->grad[k2];
         }

         assert(SCIPisLE(scip, g->prize[k], 0.0));

         min = g->prize[k];
         if( k2 != UNKNOWN )
            min += g->prize[k2];

         /* check all neighbours of k */
         e = g->outbeg[k];
         while( e != EAT_LAST )
         {
            j = g->head[e];
            e = g->oeat[e];
            /* valid candidate? */
            if( g->grad[j] <= maxgrad && g->mark[j] && SCIPisLE(scip, g->prize[j], min) )
            {
               for( e2 = g->outbeg[j]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                  if( !marked[g->head[e2]] )
                     break;

               /* neighbours of j subset of those of k? */
               if( e2 == EAT_LAST )
               {
                  while( g->outbeg[j] != EAT_LAST )
                  {
                     e2 = g->outbeg[j];
                     (*count)++;
                     graph_edge_del(scip, g, e2, TRUE);
                  }
                  g->mark[j] = FALSE;
                  marked[j] = FALSE;
                  // printf("eliminate edge \n" );
               }
            }
         }

         for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
            marked[g->head[e]] = FALSE;

         for( l = 0; l < nn; l++ )
            for( e = g->outbeg[neighbarr[l]]; e != EAT_LAST; e = g->oeat[e] )
               marked[g->head[e]] = FALSE;

         marked[k] = FALSE;

         for( k2 = 0; k2 < nnodes; k2++ )
            assert(marked[k2] == FALSE);

         for( k2 = 0; k2 < nnodes; k2++ )
            if( marked[k2] )
               printf("FAIIL \n \n" );
      }
   }

   SCIPfreeBufferArray(scip, &neighbarr);
   return SCIP_OKAY;
}


/** NPV reduction */
SCIP_RETCODE npvReduction(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit
   )
{
   GRAPH* auxg;
   PATH* mst;
   SCIP_Real prize;
   SCIP_Real sdist0;
   SCIP_Real sdist1;
   SCIP_Real sdist2;
   int* adjverts;
   int i;
   int k;
   int k2;
   int s;
   int l;
   int e;
   int nnodes;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(nelims != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &adjverts, 5) );

   *nelims = 0;
   nnodes = g->knots;

   /* initialize arrays */
   for( i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }


   /* --- NPV3 test --- */

   /* try to eliminate non-positive vertices of degree 3 */
   for( i = 0; i < nnodes; i++ )
   {
      assert(g->grad[i] >= 0);
      /* only non-positive vertices of degree 3 */
      if( !g->mark[i] || g->grad[i] != 3 || Is_term(g->term[i]) )
         continue;

      k = 0;
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
	 adjverts[k++] = g->head[e];

      assert(k == 3);

      g->mark[i] = FALSE;
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist0, heap, statetail, statehead, memlbltail, memlblhead, adjverts[0], adjverts[1], limit, FALSE, TRUE) );
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist1, heap, statetail, statehead, memlbltail, memlblhead, adjverts[1], adjverts[2], limit, FALSE, TRUE) );
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist2, heap, statetail, statehead, memlbltail, memlblhead, adjverts[2], adjverts[0], limit, FALSE, TRUE) );
      prize = g->prize[i];

      /* can vertex be deleted? */
      if( (SCIPisGE(scip, -sdist0 - sdist1, prize) && SCIPisGE(scip, -sdist2, prize))
         || (SCIPisGE(scip, -sdist1 - sdist2, prize) && SCIPisGE(scip, -sdist0, prize))
         || (SCIPisGE(scip, -sdist2 - sdist0, prize) && SCIPisGE(scip, -sdist1, prize))
         )
      {
         SCIPdebugMessage("npv3Reduction delete: %d (prize: %f) \n", i,  g->prize[i]);
	 (*nelims) +=  g->grad[i];
         while( g->outbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->outbeg[i], TRUE);
      }
      else
      {
         g->mark[i] = TRUE;
      }
   }


   /* --- NPV4 test --- */

   /* initialize mst struct and new graph for further tests */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, 5) );
   SCIP_CALL( graph_init(scip, &auxg, 5, 40, 1, 0) );

   for( k = 0; k < 4; k++ )
      graph_knot_add(auxg, -1);
   for( k = 0; k < 4; k++ )
      for( k2 = k + 1; k2 < 4; k2++ )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   SCIP_CALL( graph_path_init(scip, auxg) );

   /* try to eliminate non-positive vertices of degree 4 */
   for( i = 0; i < nnodes; i++ )
   {
      /* only non-positive vertices of degree 4 */
      if( !g->mark[i] || g->grad[i] != 4 || Is_term(g->term[i]) )
         continue;
      k = 0;

      /* store neighbours */
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         adjverts[k++] = g->head[e];

      assert(k == 4);

      g->mark[i] = FALSE;
      prize = g->prize[i];

      /* compute mw bottleneck distance to each pair of neighbours */
      for( k = 0; k < 4; k++ )
      {
         auxg->mark[k] = TRUE;
         for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
         {
            k2 = auxg->head[e];
            if( k2 > k )
            {
               SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sdist0), heap, statetail, statehead, memlbltail, memlblhead, adjverts[k], adjverts[k2], limit, FALSE, TRUE) );
               auxg->cost[e] = sdist0;
	       //printf("%d %d cost %f ", k, k2, auxg->cost[e]);
               if( SCIPisGT(scip, prize, -auxg->cost[e]) )
                  break;

               auxg->cost[flipedge(e)] = auxg->cost[e];
            }
         }
         if( e != EAT_LAST )
            break;
      }

      k = UNKNOWN;
      if( e == EAT_LAST )
      {
	 /* compute mst on all neighbours */
         graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);

	 /* calculate mst cost */
         sdist0 = 0.0;
         for( l = 1; l < 4; l++ )
            sdist0 += mst[l].dist;

	 //printf("prize: %f  mst4cost: %f \n", prize, -sdist0);
         if( SCIPisLE(scip, prize, -sdist0) )
         {
            /* compute subset msts on all neighbours */
            for( k = 0; k < 4; k++ )
            {
               auxg->mark[k] = FALSE;
               sdist0 = 0.0;
               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( l = 1; l < 4; l++ )
                     if( auxg->mark[l] )
                        sdist0 += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( l = 2; l < 4; l++ )
                     sdist0 += mst[l].dist;
               }
               //printf("  mst3cost: %f \n", -sdist0);
               auxg->mark[k] = TRUE;
               if( SCIPisGT(scip, prize, -sdist0) )
                  break;
            }
         }
      }

      /* can node be eliminated? */
      if( k == 4 )
      {
         SCIPdebugMessage("npv4Reduction delete: %d (prize: %f) \n", i,  g->prize[i]);
	 (*nelims) += g->grad[i];
         while( g->outbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->outbeg[i], TRUE);
      }
      else
      {
         g->mark[i] = TRUE;
      }
   }

   /* --- NPV5 test --- */

   /* enlarge graph for NPV5 test*/
   graph_knot_add(auxg, -1);
   for( k = 0; k < 4; k++ )
      graph_edge_add(scip, auxg, k, 4, 1.0, 1.0);
   graph_path_exit(scip, auxg);
   SCIP_CALL( graph_path_init(scip, auxg) );

   /* try to eliminate non-positive vertices of degree 5 */
   for( i = 0; i < nnodes; i++ )
   {
      /* only non-positive vertices of degree 5 */
      if( !g->mark[i] || g->grad[i] != 5 || Is_term(g->term[i]) )
         continue;
      k = 0;

      /* store neighbours */
      for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         adjverts[k++] = g->head[e];

      assert(k == 5);

      g->mark[i] = FALSE;
      prize = g->prize[i];

      for( k = 0; k < 5; k++ )
      {
         auxg->mark[k] = TRUE;
         for( e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
         {
            k2 = auxg->head[e];
            if( k2 > k )
            {
               SCIP_CALL( getSD(scip, g, pathtail, pathhead, &(sdist0), heap, statetail, statehead, memlbltail, memlblhead, adjverts[k], adjverts[k2], limit, FALSE, TRUE) );
               auxg->cost[e] = sdist0;
	       //printf("%d %d cost %f ", k, k2, auxg->cost[e]);
               if( SCIPisGT(scip, prize, -auxg->cost[e]) )
                  break;

               auxg->cost[flipedge(e)] = auxg->cost[e];
            }
         }
         if( e != EAT_LAST )
            break;
      }
      //printf("\n");
      k = UNKNOWN;
      if( e == EAT_LAST )
      {
         //printf("\n  possible 5 delete: %d \n", i);

         graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
         sdist0 = 0.0;
         for( l = 1; l < 5; l++ )
            sdist0 += mst[l].dist;

	 //printf("prize: %f  mst5cost: %f \n", prize, -sdist0);
         if( SCIPisLE(scip, prize, -sdist0) )
         {
            // printf("possible 5 delete: %d \n", i);
            for( k = 0; k < 5; k++ )
            {
               auxg->mark[k] = FALSE;
               sdist0 = 0.0;
               if( k != 0 )
               {
                  graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                  for( l = 1; l < 5; l++ )
                     if( auxg->mark[l] )
                        sdist0 += mst[l].dist;
               }
               else
               {
                  graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                  for( l = 2; l < 5; l++ )
                     sdist0 += mst[l].dist;
               }
               // printf("  mst4cost: %f \n", -sdist0);
               k2 = UNKNOWN;
	       if( SCIPisLE(scip, prize, -sdist0) )
	       {
                  for( k2 = k + 1; k2 < 5; k2++ )
                  {
                     if( k2 == k )
                        continue;
                     auxg->mark[k2] = FALSE;
                     sdist0 = 0.0;
                     if( k2 != 0 && k != 0)
                     {
                        graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                        for( l = 1; l < 5; l++ )
                           if( auxg->mark[l] )
                              sdist0 += mst[l].dist;
                     }
                     else
                     {
		        if( k != 1 && k2 != 1 )
                           s = 1;
			else
                           s = 2;
                        graph_path_exec(scip, auxg, MST_MODE, s, auxg->cost, mst);
                        for( l = 0; l < 5; l++ )
                           if( auxg->mark[l] && l != s  )
                              sdist0 += mst[l].dist;
                     }
                     //      printf("  mst3cost: %f \n", -sdist0);
                     auxg->mark[k2] = TRUE;
                     if( SCIPisGT(scip, prize, -sdist0) )
                        break;
                  }
	       }
               auxg->mark[k] = TRUE;
               if( k2 != 5 )
		  break;
               //if( SCIPisGT(scip, prize, -sdist0) )
               // break;
            }
         }
      }

      if( k == 5 )
      {
         printf(" \n npv5Reduction delete: %d (prize: %f) \n", i,  g->prize[i]);
         (*nelims) += g->grad[i];
         while( g->outbeg[i] != EAT_LAST )
            graph_edge_del(scip, g, g->outbeg[i], TRUE);

      }
      else
      {
         g->mark[i] = TRUE;
      }
   }

   /* free memory*/
   graph_path_exit(scip, auxg);
   graph_free(scip, auxg, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &adjverts);

   return SCIP_OKAY;
}


/** chain reduction */
SCIP_RETCODE chain2Reduction(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 pathtail,
   PATH*                 pathhead,
   int*                  heap,
   int*                  statetail,
   int*                  statehead,
   int*                  memlbltail,
   int*                  memlblhead,
   int*                  nelims,
   int                   limit
   )
{
   SCIP_Real sdist;
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(nelims != NULL);

   *nelims = 0;
   nnodes = g->knots;

   /* initialize arrays */
   for( i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

   for( i = 0; i < nnodes; i++ )
   {
      assert(g->grad[i] >= 0);
      if( !g->mark[i] || g->grad[i] == 0 || Is_term(g->term[i]) || g->grad[i] != 2 )
         continue;

      /* non-positive chains */
      e1 = g->outbeg[i];
      e2 = g->oeat[e1];
      i1 = g->head[e1];
      i2 = g->head[e2];

      assert(e1 >= 0);
      assert(e2 >= 0);
      assert(g->mark[i1]);
      assert(g->mark[i2]);
      g->mark[i] = FALSE;
      SCIP_CALL( getSD(scip, g, pathtail, pathhead, &sdist, heap, statetail, statehead, memlbltail, memlblhead, i1, i2, limit, FALSE, TRUE) );
      if( SCIPisGE(scip, -sdist, g->prize[i]) )
      {
         SCIPdebugMessage("delete : %d prize: %f sd: %f \n", i,  g->prize[i], -sdist );
         graph_edge_del(scip, g, e1, TRUE);
         graph_edge_del(scip, g, e2, TRUE);
         (*nelims) += 2;
      }
      else
      {
         g->mark[i] = TRUE;
      }
   }
   return SCIP_OKAY;
}
/** NNP reduction */
SCIP_RETCODE nnpReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  mem,                /**< nodes array */
   int*                  marked,             /**< nodes array */
   int*                  visited,            /**< nodes array */
   int*                  count,              /**< pointer to number of reductions */
   int                   maxniter,           /**< max number of edges to check */
   char*                 positive            /**< nodes array */
   )
{
   SCIP_QUEUE* queue;
   int* pnode;
   int i;
   int j;
   int k;
   int e;
   int e2;
   int erev;
   int enext;
   int nnodes;
   int nvisited1;
   int nvisited2;
   int iterations;
   SCIP_Bool success;

   assert(scip   != NULL);
   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count  != NULL);
   assert(mem != NULL);
   assert(positive != NULL);
   assert(visited != NULL);
   assert(marked != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   *count = 0;
   nnodes = g->knots;

   /* unmark all nodes */
   for( i = 0; i < nnodes; i++ )
   {
      if( g->mark[i] && SCIPisGE(scip, g->prize[i], 0.0) )
         positive[i] = TRUE;
      else
         positive[i] = FALSE;
      marked[i] = FALSE;
      visited[i] = FALSE;
   }

   SCIP_CALL( SCIPqueueCreate(&queue, nnodes, 2.0) );

   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
	 continue;

      /* mark adjacent vertices and i*/
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
	 j = g->head[e];
	 enext = g->oeat[e];
	 if( g->mark[j] )
	 {
	    nvisited1 = 0;
            iterations = 0;
	    success = FALSE;
            assert(SCIPqueueIsEmpty(queue));

            SCIP_CALL( SCIPqueueInsert(queue, &i) );
	    marked[i] = TRUE;
	    mem[nvisited1++] = i;

            /* BFS */
            while( !SCIPqueueIsEmpty(queue) )
            {
               pnode = SCIPqueueRemove(queue);
               for( e2 = g->outbeg[*pnode]; e2 != EAT_LAST; e2 = g->oeat[e2] )
               {
                  if( e2 == e )
                     continue;
                  k = g->head[e2];
		  if( iterations++ >= maxniter || k == j )
		  {
		     if( k == j )
		        success = TRUE;
		     SCIPqueueClear(queue);
                     break;
		  }

                  if( positive[k] && !marked[k] )
                  {
		     marked[k] = TRUE;
		     assert(nvisited1 < nnodes);
                     mem[nvisited1++] = k;
                     SCIP_CALL( SCIPqueueInsert(queue, &g->head[e2]) );
                  }
               }
            }

            nvisited2 = 0;
            /* vertex j not reached yet? */
            if( !success )
            {
	       assert(SCIPqueueIsEmpty(queue));
               SCIP_CALL( SCIPqueueInsert(queue, &j) );
               visited[j] = TRUE;
	       assert(nvisited2 + nvisited1 < nnodes);
               mem[nvisited1 + nvisited2++] = j;
               erev = flipedge(e);

               iterations = 0;
               while( !SCIPqueueIsEmpty(queue) )
               {
                  pnode = (SCIPqueueRemove(queue));
                  for( e2 = g->outbeg[*pnode]; e2 != EAT_LAST; e2 = g->oeat[e2] )
                  {
                     if( e2 == erev )
                        continue;
                     k = g->head[e2];

                     if( iterations++ >= maxniter || marked[k] )
                     {
                        if( marked[k] )
                           success = TRUE;
                        SCIPqueueClear(queue);
                        break;
                     }

                     if( positive[k] && !visited[k] )
                     {
                        visited[k] = TRUE;
                        assert(nvisited2 + nvisited1 < nnodes);
                        mem[nvisited1 + nvisited2++] = k;
                        SCIP_CALL( SCIPqueueInsert(queue, &g->head[e2]) );
                     }
                  }
               }
            }

            if( success )
            {
               (*count)++;
               graph_edge_del(scip, g, e, TRUE);
            }
            for( j = 0; j < nvisited1; j++ )
               marked[mem[j]] = FALSE;

            for( j = nvisited1; j < nvisited1 + nvisited2; j++ )
               visited[mem[j]] = FALSE;

            for( j = 0; j < nnodes; j++ )
	    {
               assert(marked[j] == FALSE);
               assert(visited[j] == FALSE);
            }
         }
         e = enext;
      }
   }
   SCIPqueueFree(&queue);

   return SCIP_OKAY;
}

/** simple reduction tests for the MWCS problem */
SCIP_RETCODE degree_test_mw(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count,              /**< pointer to number of reductions */
   int                   maxnvisits          /**< maximal number of neighbours to visit per node */
   )
{
   int i;
   int e;
   int i1;
   int i2;
   int e1;
   int e2;
   int nnodes;
   int nedges;
   int edgecount;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);
   assert(g->stp_type == STP_MAX_NODE_WEIGHT);

   nnodes = g->knots;
   nedges = g->edges;

   SCIPdebugMessage("MW degree test: ");

   while( rerun )
   {
      rerun = FALSE;
      for( e = 0; e < nedges; e += 2 )
      {
	 i1 = g->tail[e];
	 i2 = g->head[e];
	 if( g->mark[i1] && g->mark[i2] && Is_term(g->term[i1]) && Is_term(g->term[i2]) )
	 {
            SCIPdebugMessage("contract tt %d->%d\n ", i1, i2);
            (*count)++;
            SCIP_CALL( graph_knot_contractpc(scip, g, i1, i2, i1) );
	 }
      }
      for( i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);
         if( !g->mark[i] || g->grad[i] == 0 )
            continue;

         assert( !SCIPisEQ(scip, g->prize[i], 0.0) );

         if( !Is_term(g->term[i]) )
         {
            if( g->grad[i] == 1 )
            {
               e1 = g->inpbeg[i];
               i1 = g->tail[e1];
               assert(e1 >= 0);
               assert(e1 == Edge_anti(g->outbeg[i]));
               assert(g->ieat[e1] == EAT_LAST);
               assert(g->oeat[g->outbeg[i]] == EAT_LAST);
	       assert(SCIPisLE(scip, g->prize[i], 0.0));

               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete NT %d\n ",  i);
               assert(g->grad[i] == 0);

               if( (i1 < i) && (g->grad[i1] < 3 || (g->grad[i1] == 3 && Is_term(g->term[i1]))) )
                  rerun = TRUE;

               (*count)++;
               continue;
            }

            /* contract non-positive chains */
            if( g->grad[i] == 2 )
            {
#if 1
	       int length = 0;
	       int f1 = -1;
	       int f2 = -1;
               e1 = g->outbeg[i];
               e2 = g->oeat[e1];
               i1 = g->head[e1];
               i2 = g->head[e2];

               assert(e1 >= 0);
               assert(e2 >= 0);
               assert(g->mark[i1]);
               assert(g->mark[i2]);
               assert(i1 != i2);

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
                  for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
                     g->cost[e] = -g->prize[i];
                  e1 = g->outbeg[i];
                  e2 = g->oeat[e1];


                  //  printf("contract NT %d %d %d length: %d \n ", i1, i, i2, length);
                  (*count) += length;
	       }
# if 0
               assert(g->grad[i] <= 2);
	       for( e = g->inpbeg[i]; e != EAT_LAST; e = g->ieat[e] )
                  g->cost[e] = -g->prize[i];
               if( length > 1 )
		  SCIPdebugMessage("contract NT %d %d %d length: %d \n ", i1, i, i2, length);
               (*count) += length;
#endif
#endif
	    }
            continue;
         }

         /* node i is a terminal: */

         /* terminal of (real) degree 0? */
         if( g->grad[i] == 2 )
         {
            /* if terminal node i is not the one with the highest prize, delete */
            if( !maxprize(scip, g, i) )
            {
               SCIPdebugMessage("MW delete 0 term %d prize: %f count:%d\n ", i, g->prize[i], *count);
               (*fixed) += g->prize[i];
               (*count) += deleteterm(scip, g, i);
            }
         }
         /* terminal of (real) degree 1? */
         else if( g->grad[i] == 3 )
         {
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] )
                  break;
            assert(e != EAT_LAST);
            assert(g->head[e] != g->source[0]);
	    if( !Is_term(g->term[g->head[e]]) )
	    {
               SCIP_CALL( trydg1edgepc(scip, g, fixed, count, i, e, &rerun) );
               continue;
	    }
         }
#if 0
         /* contract adjacent terminals */
         if( g->grad[i] >= 3 )
         {
	    edgecount = 0;
            for( e1 = g->outbeg[i]; e1 != EAT_LAST && edgecount++ <= maxnvisits; e1 = g->oeat[e1] )
            {
               i1 = g->head[e1];
               if( !g->mark[i1] || !Is_term(g->term[i1]) )
                  continue;

               SCIPdebugMessage("contract tt %d->%d\n ", i, i1);

               (*count)++;
               SCIP_CALL( graph_knot_contractpc(scip, g, i, i1, i) );
	       break;
            }
         }
#endif
      }
   }

   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);

   return SCIP_OKAY;
}


SCIP_RETCODE degree_test_pc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int* edges2;
   int* nodes2;
   int i;
   int i1;
   int i2;
   int e;
   int e1;
   int e2;
   int nnodes;
   SCIP_Bool pc;
   SCIP_Bool rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);
   assert(g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING);

   pc = (g->stp_type == STP_PRIZE_COLLECTING);

   nnodes = g->knots;
   *count = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edges2, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes2, 2) );

   SCIPdebugMessage("Degree Test: ");

   if( !pc )
      g->mark[g->source[0]] = FALSE;

   while( rerun )
   {
      rerun = FALSE;

      for( i = 0; i < nnodes; i++ )
      {
         assert(g->grad[i] >= 0);
         if( !g->mark[i] || g->grad[i] == 0 )
            continue;

         if( !Is_term(g->term[i]) )
         {
            if( g->grad[i] == 1 )
            {
               e1 = g->inpbeg[i];
               i1 = g->tail[e1];
               assert(e1 >= 0);
               assert(e1 == Edge_anti(g->outbeg[i]));
               assert(g->ieat[e1] == EAT_LAST);
               assert(g->oeat[g->outbeg[i]] == EAT_LAST);

               graph_edge_del(scip, g, e1, TRUE);
               SCIPdebugMessage("delete NT %d\n ",  i);
               assert(g->grad[i] == 0);

               /* the last node? */
               if( g->grad[i1] == 0 )
               {
                  rerun = FALSE;
                  break;
               }
               if( (i1 < i) && (g->grad[i1] < 3 || Is_term(g->term[i1])) )
                  rerun = TRUE;

               (*count)++;
               continue;
            }

            /* contract non terminals of degree 2 */
            if( g->grad[i] == 2 )
            {
               e1 = g->outbeg[i];
               e2 = g->oeat[e1];
               i1 = g->head[e1];
               i2 = g->head[e2];

               assert(e1 >= 0);
               assert(e2 >= 0);
               assert(g->mark[i1] || i1 == g->source[0]);
               assert(g->mark[i2] || i2 == g->source[0]);
               assert(EQ(g->cost[e2], g->cost[Edge_anti(e2)]));

               g->cost[e1]            += g->cost[e2];
               g->cost[Edge_anti(e1)] += g->cost[e2];

               SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
               (*count)++;
               SCIPdebugMessage("contract NT %d %d \n ", i2, i);
               if( (Is_term(g->term[i2]) && (i2 < i)) || (Is_term(g->term[i1]) && (i1 < i)) )
                  rerun = TRUE;
            }
            continue;
         }

         /* node i is a terminal: */

         /* terminal of (real) degree 0? */
         if( (g->grad[i] == 2 && pc) || (g->grad[i] == 1 && !pc) )
         {
            /* if terminal node i is node the one with the highest prize, delete*/
            if( !maxprize(scip, g, i) )
            {
               SCIPdebugMessage("delete 0 term %d prize: %f count:%d\n ", i, g->prize[i], *count);
               (*fixed) += g->prize[i];
               (*count) += deleteterm(scip, g, i);
            }
         }
         /* terminal of (real) degree 1? */
         else if( (g->grad[i] == 3 && pc) || (g->grad[i] == 2 && !pc)  )
         {
            for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               if( g->mark[g->head[e]] || (!pc && g->head[e] == g->source[0]) )
                  break;
            assert(e != EAT_LAST);
            assert(g->head[e] != g->source[0] || !pc);
            SCIP_CALL( trydg1edgepc(scip, g, fixed, count, i, e, &rerun) );
            SCIPdebugMessage("delete 1 term %d\n ", i);
         }
         /* terminal of (real) degree 2? */
         else if( (g->grad[i] == 4 && pc) || (g->grad[i] == 3 && !pc) )
         {
            if( !maxprize(scip, g, i) )
            {
               i2 = 0;
               for( e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
               {
                  i1 = g->head[e];
                  if( g->mark[i1] )
                  {
                     if( i2 >= 2 )
                        assert(i2 < 2);
                     edges2[i2] = e;
                     nodes2[i2++] = i1;
                  }
               }
               if( SCIPisLE(scip, g->prize[i], g->cost[edges2[0]]) && SCIPisLE(scip, g->prize[i], g->cost[edges2[1]]) )
               {
                  int n1;
                  IDX* ancestors = NULL;
                  IDX* revancestors = NULL;

                  e = edges2[0];
                  e1 = edges2[1];
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[e]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(ancestors), g->ancestors[Edge_anti(e1)]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[Edge_anti(e)]) );
                  SCIP_CALL( SCIPintListNodeAppendCopy(scip, &(revancestors), g->ancestors[e1]) );
                  SCIPdebugMessage("delete - term - %d\n ", i);
                  /* contract edge */
                  n1 = graph_edge_redirect(scip, g, e, nodes2[1], nodes2[0], g->cost[e] + g->cost[e1] - g->prize[i]);
                  /* new edge inserted? */
                  if( n1 >= 0)
                  {
                     /* add ancestors */
                     SCIPintListNodeFree(scip, &(g->ancestors[n1]));
                     SCIPintListNodeFree(scip, &(g->ancestors[Edge_anti(n1)]));
                     SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[n1]), ancestors) );
                     SCIP_CALL(  SCIPintListNodeAppendCopy(scip, &(g->ancestors[Edge_anti(n1)]), revancestors) );
                  }
                  (*count) += deleteterm(scip, g, i);
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

            for( e1 = g->outbeg[i]; e1 != EAT_LAST; e1 = g->oeat[e1] )
            {
               i1 = g->head[e1];
               if( !g->mark[i1] )
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
               i1 = g->head[ett];
               SCIPdebugMessage("contract tt %d->%d\n ", i, i1);
               assert(SCIPisLT(scip, mincost, FARAWAY));
               *fixed += g->cost[ett];
               (*count)++;
               SCIP_CALL( graph_knot_contractpc(scip, g, i, i1, i) );
               rerun = TRUE;
            }
         }
      }
   }

   if( !pc )
      g->mark[g->source[0]] = TRUE;
   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);
   /* free memory */
   SCIPfreeBufferArray(scip, &nodes2);
   SCIPfreeBufferArray(scip, &edges2);

   return SCIP_OKAY;
}
