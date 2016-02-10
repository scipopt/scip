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


static
int deleteterm(
   SCIP* scip,
   GRAPH* g,
   int i
   )
{
   int e;
   int t;
   int i1;
   int etemp;
   int count;

   assert(g      != NULL);
   assert(scip      != NULL);
   assert(Is_term(g->term[i]));
   t = UNKNOWN;
   count = 0;
   /* delete terminal */
   graph_knot_chg(g, i, -1);
   g->mark[i] = FALSE;
   e = g->outbeg[i];
   while( e != EAT_LAST )
   {
      i1 = g->head[e];
      if( Is_pterm(g->term[i1]) && g->source[0] != i1 )
         t = g->head[e];
      assert(e >= 0);
      etemp = g->oeat[e];
      count++;
      graph_edge_del(scip, g, e, TRUE);
      e = etemp;
   }
   assert(t != UNKNOWN);

   /* delete artifical terminal */
   graph_knot_chg(g, t, -1);
   e = g->outbeg[t];
   while( e != EAT_LAST )
   {
      assert(e >= 0);
      etemp = g->oeat[e];
      count++;
      graph_edge_del(scip, g, e, TRUE);
      e = etemp;
   }
   return count;
}

static
SCIP_Bool maxprize(
   SCIP* scip,
   GRAPH* g,
   int i
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
   SCIP* scip,
   GRAPH* g,
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

   if( SCIPisLE(scip, g->prize[i], g->cost[iout]) )
   {
      if( (i1 < i) && (Is_term(g->term[i1]) || g->grad[i1] == 2 || g->grad[i1] == 3) )
         (*rerun) = TRUE;
      SCIPdebugMessage("DEL (1 edge) terminal %d \n", i);
      (*offset) += g->prize[i];
      *count += deleteterm(scip, g, i);
   }
   else
   {
      (*rerun) = TRUE;
      printf("contract (1 edge) terminal %d \n", i);
      *offset += g->cost[iout];
      degsum = g->grad[i] + g->grad[i1];
      SCIP_CALL( graph_knot_contractpc(scip, g, i, i1, i) );
      degsum = degsum - g->grad[i];
      assert(degsum >= 1);
      *count += degsum;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE degree_test_dir(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to offfset value */
   int*                  count               /**< pointer to number of reductions */
   )
{
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int rerun = TRUE;
   int nnodes;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);

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
            if( !g->mark[i1] )
               continue;
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

            /* the last node? */
            if( g->grad[i1] == 0 )
            {
               rerun = FALSE;
               break;
            }
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
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  (*count)++;

                  if (!Is_term(g->term[i2]) && (((i1 < i) && (g->grad[i1] < 3))
                        || ((i2 < i) && (g->grad[i2] < 3))))
                     rerun = TRUE;

               }
            }
            /* CONSTCOND */
            /*lint -save -e717 */
            /*lint -restore */
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", *count);
   SCIPdebugMessage("dirdeg %d Knots deleted\n", *count);
   assert(graph_valid(g));

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

   pc = g->stp_type == STP_PRIZE_COLLECTING;

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
	       /*
                 printf("knot: %d \n ",  i);
                 printf("eq %d %d (%d)\n ", g->tail[e2], g->head[e2], e2);
                 printf("term? %d %d \n ", Is_term(g->term[g->tail[e2]]),  Is_term(g->term[g->head[e2]]));
	       */
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
