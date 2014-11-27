/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: reduce.c                                                      */
/*   Name....: Steiner Tree Reduction                                        */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h)           */

#define REDUCE_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "heur_tm.h"
#include "portab.h"
#include "misc_stp.h"
#include "scip/scip.h"

/* Moeglichkeiten:
 *    a 1 b 2 c
 *    *---*       : Kontraction entlang 1, a -> b
 *    *---o---*   : Kontraktion entlang 2, b -> c
 *    t---t---t   : Kontraktion entlang min(1, 2), b -> min(d(a),d(c))
 *    o---t---o   : Nichts
 *    t---t---o   : Kontraktion entlang 1, b -> a, wenn c(1) <= c(2)
 *    o---t---t   : Kontraktion entlang 2, b -> c, wenn c(2) <= c(1)
 */
static int degree_test(
   GRAPH*  g,
   double* fixed)
{
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int rerun = TRUE;
   int done  = TRUE;
   int count = 0;

   //int e;
   //IDX* curr;
   assert(g      != NULL);
   assert(fixed  != NULL);

   SCIPdebugMessage("Degree Test: ");
   fflush(stdout);

   while(rerun)
   {
      rerun = FALSE;

      SCIPdebug(fputc('.', stdout));
      SCIPdebug(fflush(stdout));

      for(i = 0; i < g->knots; i++)
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

            if( Is_term(g->term[i]) )
	    {
               *fixed += g->cost[e1];
	       //printf("fix1 %d, %d\n", i1, i);
	       SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e1]);
	    }
            graph_knot_contract(g, i1, i);
	    /*
              printf("con1 %d, %d\n", i1, i);
              for( e = g->outbeg[i1]; e != EAT_LAST; e = g->oeat[e] )
              {
              curr = g->ancestors[e];
              while( curr != NULL )
              {
              printf("%d->%d ancestor: %d->%d \n", g->tail[e], g->head[e], g->tail[curr->index], g->head[curr->index] );
              curr = curr->parent;
              }
              }
            */
            assert(g->grad[i] == 0);

            /* Ist es etwa der Letzte gewesen ?
             */
            if (g->grad[i1] == 0)
            {
               rerun = FALSE;
               break;
            }
            if ((i1 < i) && (g->grad[i1] < 3))
               rerun = TRUE;

            count++;
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

            do
            {
               done = TRUE;

               if( !Is_term(g->term[i]) )
               {
                  assert(EQ(g->cost[e2], g->cost[Edge_anti(e2)]));

                  g->cost[e1]            += g->cost[e2];
                  g->cost[Edge_anti(e1)] += g->cost[e2];
                  //printf("con0 %d, %d\n", i2, i);
                  graph_knot_contract(g, i2, i);
                  /*
                    for( e = g->outbeg[i2]; e != EAT_LAST; e = g->oeat[e] )
                    {
                    curr = g->ancestors[e];
                    while( curr != NULL )
                    {
                    printf("%d->%d ancestor: %d->%d \n", g->tail[e], g->head[e], g->tail[curr->index], g->head[curr->index] );
                    curr = curr->parent;
                    }
                    }*/
                  count++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if (Is_term(g->term[i1]) && Is_term(g->term[i2]))
               {
                  if (LT(g->cost[e1], g->cost[e2]))
                  {
                     //printf("con2 %d, %d\n", i1, i);
                     *fixed += g->cost[e1];
		     SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e1]);
                     graph_knot_contract(g, i1, i);
                     /*
                       for( e = g->outbeg[i1]; e != EAT_LAST; e = g->oeat[e] )
                       {
                       curr = g->ancestors[e];
                       while( curr != NULL )
                       {
                       printf("%d->%d ancestor: %d->%d \n", g->tail[e], g->head[e], g->tail[curr->index], g->head[curr->index] );
                       curr = curr->parent;
                       }
                       }
                     */
                  }
                  else
                  {
                     //printf("con3 \n");
                     *fixed += g->cost[e2];
		     SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e2]);
                     graph_knot_contract(g, i2, i);
                  }
                  count++;

                  break;
               }
               if (Is_term(g->term[i1]) && !Is_term(g->term[i2]) && LE(g->cost[e1], g->cost[e2]))
               {
		  //printf("con4 \n");
                  *fixed += g->cost[e1];
		  SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e1]);
                  graph_knot_contract(g, i1, i);
                  count++;

                  break;
               }
               if (Is_term(g->term[i2]) && !Is_term(g->term[i1]) && LE(g->cost[e2], g->cost[e1]))
               {
                  //printf("con5 %d, %d\n", i2, i);
		  SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e2]);
                  *fixed += g->cost[e2];
		  /*
                    curr = g->fixedges;
                    while( curr != NULL )
                    {
                    printf("fixed (con5): %d->%d \n",  g->orgtail[curr->index], g->orghead[curr->index] );
                    curr = curr->parent;
                    }
                  */
                  graph_knot_contract(g, i2, i);
                  /*
                    for( e = g->outbeg[i2]; e != EAT_LAST; e = g->oeat[e] )
                    {
                    curr = g->ancestors[e];
                    while( curr != NULL )
                    {
                    printf("%d->%d ancestor(con5): %d->%d \n", g->tail[e], g->head[e], g->orgtail[curr->index], g->orghead[curr->index] );
                    curr = curr->parent;
                    }
                    } */
                  count++;

                  break;
               }
               done = FALSE;
            }
            /* CONSTCOND */
            /*lint -save -e717 */
            while(FALSE);
            /*lint -restore */

            if (done
               && (((i1 < i) && (g->grad[i1] < 3))
                  || ((i2 < i) && (g->grad[i2] < 3))))
               rerun = TRUE;
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", count);
   printf(" %d Knots deleted\n", count);
   assert(graph_valid(g));

   return count;
}

/* A. Balakrishnan and N. R. Patel
 *
 * "Problem Reduction Methods and a Tree Generation Algorithm
 *             for the Steiner Network Problem"
 *
 * NETWORKS, Vol 17 (1985) 65-85
 *
 * Demand Node Aggregation, Page 68
 */
static int tt_aggregation(
   GRAPH*  g,
   double* fixed)
{
   PATH*  mst;
   int    i;
   int    head;
   int    tail;
   int    count = 0;
   int    retry;

   assert(g      != NULL);
   assert(fixed  != NULL);

   SCIPdebugMessage("T-T Aggregation: ");
   fflush(stdout);

   mst = malloc((size_t)g->knots * sizeof(PATH));

   assert(mst != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = TRUE;

   do
   {
      retry = FALSE;

      graph_path_exec(g, MST_MODE, g->source[0], g->cost, mst);

      for(i = 0; !retry && (i < g->knots); i++)
      {
         assert((mst[i].edge >= 0) || (g->grad[i] == 0) || (i == g->source[0]));

         if (mst[i].edge >= 0)
         {
            head = g->head[mst[i].edge];
            tail = g->tail[mst[i].edge];

            /* Ist es eine T-T Verbindung ?
             */
            if (Is_term(g->term[head]) && Is_term(g->term[tail]))
            {
               graph_knot_contract(g, head, tail);

               *fixed += mst[i].dist;

               count++;

               /* Neuen MST berechnen und nochmal.
                */
               retry = TRUE;
            }
         }
      }
   } while(retry);

   free(mst);

   SCIPdebugMessage("%d Knots deleted\n", count);

   assert(graph_valid(g));

   return(count);
}

/* A. Balakrishnan and N. R. Patel
 *
 * "Problem Reduction Methods and a Tree Generation Algorithm
 *             for the Steiner Network Problem"
 *
 * NETWORKS, Vol 17 (1985) 65-85
 *
 * R-R Edge Deletion, Page 69
 */
static int tt_deletion(
   GRAPH* g)
{
   PATH*   mst;
   double* cost;
   int     count = 0;
   int     i;
   int     k;
   int     e;
   int     f;
   double  max = 0.0;

   assert(g != NULL);

   SCIPdebugMessage("T-T Edge deletion: ");
   fflush(stdout);

   for(i = 0; i < g->edges; i++)
      if (GT(g->cost[i], max))
         max = g->cost[i];

   assert(GT(max, 0.0));

   /* Groesser als alle anderen.
    */
   max = max * 2.0 + 1.0;

   cost = malloc((size_t)g->edges * sizeof(double));

   assert(cost != NULL);

   for(i = 0; i < g->edges; i++)
      cost[i] = (Is_term(g->term[g->head[i]]) && Is_term(g->term[g->tail[i]])) ? g->cost[i] : max;

   mst = malloc((size_t)g->knots * sizeof(PATH));

   assert(mst != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = (g->grad[i] > 0);

   graph_path_exec(g, MST_MODE, g->source[0], cost, mst);

   for(i = 0; i < g->knots; i++)
   {
      if (!g->mark[i])
         continue;

      assert(g->grad[i] > 0);

      if (!Is_term(g->term[i]))
         continue;

      e = g->outbeg[i];

      while(e != EAT_LAST)
      {
         k = g->head[e];
         f = e;
         e = g->oeat[e];

         if (!Is_term(g->term[k]))
            continue;

         assert(Is_term(g->term[i]));
         assert(Is_term(g->term[k]));
         assert(g->tail[f] == i);
         assert(g->head[f] == k);

         if ((mst[k].edge == f) || (mst[k].edge == Edge_anti(f)))
            continue;

         if ((mst[i].edge == f) || (mst[i].edge == Edge_anti(f)))
            continue;

         graph_edge_del(g, f);

         /* Weil ja zwei Kanten geloescht wurden
          */
         count += 2;
      }
   }
#ifndef NDEBUG
   /* Es darf keine Kante geloescht worden sein, die Teil des MST war.
    */
   for(i = 0; i < g->knots; i++)
   {
      if ((e = mst[i].edge) < 0)
         continue;

      assert(g->oeat[e] != EAT_FREE);
      assert(g->ieat[e] != EAT_FREE);
   }
#endif

   free(mst);
   free(cost);

   /* Vielleicht sogar noch mehr Kanten, bei der Knoten Kontraktion
    */
   SCIPdebugMessage("%d Edges deleted\n", count);

   assert(graph_valid(g));

   return(count);
}

static
int bound_test(
   SCIP*  scip,
   GRAPH* graph
   )
{
   //SCIP_Real** pathdist;
   SCIP_Real dist;
   SCIP_Real obj;
   SCIP_Real* rad;
   PATH* vnoi;
   PATH** path;

   int k;
   int r;
   int i;
   int e;
   int nnodes;
   //int** pathedge;

   int* vbase;
   int* result;
   char* base;

   assert(scip != NULL);
   assert(graph != NULL);
   nnodes = graph->knots;

   SCIP_CALL( SCIPallocBufferArray(scip, &result, graph->edges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rad, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &base, nnodes) );
   /*SCIP_CALL( SCIPallocBufferArray(scip, &pathdist, nnodes) );
     SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
     BMSclearMemoryArray(pathdist, nnodes);
     BMSclearMemoryArray(pathedge, nnodes);*/
   for( k = nnodes - 1; k >= 0; k-- )
   {
      if( Is_term(graph->term[k]) )
	 base[k] = TRUE;
      else
	 base[k] = FALSE;
   }


   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
   BMSclearMemoryArray(path, nnodes);

   SCIP_CALL( SCIPtmHeur(scip, graph, path, graph->cost, graph->cost, result) );

   obj = 0.0;

   for( e = 0; e < graph->edges; e++ )
      if( result[e] == CONNECT )
         obj += graph->cost[e];

   voronoi_radius(scip, graph, rad, graph->cost, graph->cost, base, vbase, vnoi);
   SCIPsortReal(rad, nnodes);
   /* test each node k */
   for( k = 0; k < nnodes; k++ )
   {
      dist = FARAWAY;
      r = UNKNOWN;
      for( i = 0; i < nnodes; i++ )
      {
	 if( Is_term(graph->term[i]) && SCIPisLT(scip, path[i][k].dist, dist ) )
	 {
	    r = i;
	    dist = path[i][k].dist;
	 }
      }
      i = k;
      while(i != r)
      {
	 if( i == vbase[k] )
            break;
         e = path[r][i].edge;
         i = graph->tail[e];
      }
      if( i != r )
         continue;
      for( i = 0; i < graph->terms - 2; i++ )
	 dist += rad[i];
      if( SCIPisGT(scip, dist + path[vbase[k]][k].dist, obj) )
	 printf("node to eliminate!! \n");


   }



   for( k = graph->edges - 1; k >= 0; k-- )
      if( result[k] == 0 )
         printf("%d->%d \n", graph->tail[k] + 1, graph->head[k] + 1);
   for( k = 0; k < nnodes; k++ )
   {
      printf("vbase[%d] = %d \n", k+1, vbase[k] + 1);
      printf("RAD: %f \n", rad[k] );

   }
   for( k = nnodes - 1; k >= 0; k-- )
   {

      assert(path[k] == NULL || graph->term[k] == 0);
      SCIPfreeBufferArrayNull(scip, &(path[k]));

   }
   SCIPfreeBufferArray(scip, &path);

   SCIPfreeBufferArray(scip, &base);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &rad);
   assert(0);
   return 1;
}


/* P. Winter and J. MacGregor Smith
 *
 * "Path-Distance Heuristics for the Steiner Problem in Undirected Networks"
 *
 * Algorithmica, 1992, 309-327
 *
 * Closest Z-Vertices Reductions, Page 316
 */
static int czv_reduction(
   GRAPH*  g,
   double* fixed)
{
   int     i;
   int     k;
   int     e1;
   int     e2;
   int     count = 0;
   double  cost;
   double  min1;
   double  min2;
   double  min3;

   assert(g     != NULL);
   assert(fixed != NULL);

   SCIPdebugMessage("Lazy closest T-Edge Reduction: ");
   fflush(stdout);

   for(i = 0; i < g->knots; i++)
   {
      if (!Is_term(g->term[i]))
         continue;

      if (g->grad[i] < 2)
         continue;

      min1 = FARAWAY;
      min2 = FARAWAY;
      e1   = -1;
      e2   = -1;

      /* Die beiden kuerzesten Kanten finden
       */
      for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
      {
         cost = g->cost[k];

         if (LT(cost, min1))
         {
            min2 = min1;
            e2   = e1;
            min1 = cost;
            e1   = k;
         }
         else
         {
            if (LT(cost, min2))
            {
               min2 = cost;
               e2   = k;
            }
         }
      }
      /* Da ja Grad > 1 war muesste es die Kanten geben
       */
      if( e2 < 0 )
         printf("Error in czv_reduction (will be ignored) \n");
      assert(e1 >= 0);
      assert(e2 >= 0);
      assert(e1 != e2);

      /* Wenn der naechste ein Terminal war, brauchen wir nicht mehr suchen.
       */
      if (Is_term(g->term[g->head[e1]]))
         min3 = 0.0;
      else
      {
         min3 = FARAWAY;

         for(k = g->outbeg[g->head[e1]]; k != EAT_LAST; k = g->oeat[k])
         {
            if (!Is_term(g->term[g->head[k]]))
               continue;

            if (g->head[k] == i)
               continue;

            cost = g->cost[k];

            if (LT(cost, min3))
               min3 = cost;
         }
      }
      /* Gab es denn ein Terminal ?
       */
      if (EQ(min3, FARAWAY))
         continue;

      if (GE(min1 + min3, min2))
         continue;

      assert(g->tail[e1] == i);

      graph_knot_contract(g, i, g->head[e1]);

      *fixed += min1;

      count++;
   }
   SCIPdebugMessage("%d Knots deleted\n", count);

   assert(graph_valid(g));

   return(count);
}

/* P. Winter and J. MacGregor Smith
 *
 * "Path-Distance Heuristics for the Steiner Problem in Undirected Networks"
 *
 * Algorithmica, 1992, 309-327
 *
 * Longest Edge Reductions, Page 316
 */
static int lle_reduction(
   GRAPH* g)
{
   int    i;
   int    k;
   int    j;
   int    l;
   int    m;
   int    i1;
   int    i2;
   int    count = 0;

   assert(g      != NULL);

   SCIPdebugMessage("Lazy Longest Edge Reduction: ");
   fflush(stdout);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = FALSE;

   for(i = 0; i < g->knots; i++)
   {
      if (!Is_term(g->term[i]))
         continue;

      if (g->grad[i] < 2)
         continue;

      assert(Is_term(g->term[i]));

      /* Erst mal markieren wo wir ueberall hinkommen
       */
      for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
         g->mark[g->head[k]] = TRUE;

      for(k = g->outbeg[i]; k != EAT_LAST; k = g->oeat[k])
      {
         i1 = g->head[k];
         j  = g->outbeg[i1];

         while(j != EAT_LAST)
         {
            m  = j;
            j  = g->oeat[j];
            i2 = g->head[m];

            if (g->mark[i2])
            {
               for(l = g->outbeg[i]; l != EAT_LAST; l = g->oeat[l])
                  if (g->head[l] == i2)
                     break;

               assert(l != EAT_LAST);

               if (LE(Max(g->cost[l], g->cost[k]), g->cost[m]))
               {
                  graph_edge_del(g, m);
                  count++;
               }
            }
         }
         g->mark[i1] = FALSE;
      }
   }
   SCIPdebugMessage("%d Edges deleted\n", count * 2);

   assert(graph_valid(g));

   return(count);
}

/* P. Winter and J. MacGregor Smith
 *
 * "Path-Distance Heuristics for the Steiner Problem in Undirected Networks"
 *
 * Algorithmica, 1992, 309-327
 *
 * Longest Edge Reductions, Page 316
 */
static int le_reduction(
   GRAPH* g)
{
   PATH** path;
   char*  used;
   int*   term;
   int    terms;
   int    i;
   int    k;
   int    j;
   int    m;
   int    i1;
   int    i2;
   int    count = 0;

   assert(g      != NULL);

   SCIPdebugMessage("Longest Edge Reduction: ");
   fflush(stdout);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = (g->grad[i] > 0);

   terms = 0;
   term  = malloc((size_t)g->knots * sizeof(int));

   assert(term != NULL);

   used = malloc((size_t)(g->edges / 2) * sizeof(char));

   assert(used != NULL);

   for(i = 0; i < g->edges / 2; i++)
      used[i] = FALSE;

   path = malloc((size_t)g->knots * sizeof(PATH*));

   assert(path != NULL);

   for(i = 0; i < g->knots; i++)
   {
      if ((g->grad[i] < 1) || (!Is_term(g->term[i])))
         path[i] = NULL;
      else
      {
         term[terms++] = i;

         path[i] = malloc((size_t)g->knots * sizeof(PATH));

         assert(path[i] != NULL);

         graph_path_exec(g, FSP_MODE, i, g->cost, path[i]);

         /* Alle benutzen Kanten markieren
          */
         for(k = 0; k < g->knots; k++)
         {
            assert((k == i) || !g->mark[k] || (path[i][k].edge >= 0));

            if (path[i][k].edge >= 0)
               used[path[i][k].edge / 2] = TRUE;
         }
#ifndef NDEBUG
         for(k = 0; k < g->knots; k++)
            assert(!g->mark[k] || (path[i][k].dist < FARAWAY));
#endif
      }
   }
   for(i = 0; i < g->knots; i++)
   {
      if (Is_term(g->term[i]))
         continue;

      /* Nicht Terminals mit Grad 2 macht die Degree-Reduktion weg
       */
      if (g->grad[i] < 3)
         continue;

      assert(!Is_term(g->term[i]));

      k = g->outbeg[i];

      while(k != EAT_LAST)
      {
         m  = k;
         k  = g->oeat[k];

         if (used[m / 2])
            continue;

         i1 = g->head[m];

         for(j = 0; j < terms; j++)
         {
            i2 = term[j];

            assert(Is_term(g->term[i2]));
            assert(g->grad[i2] > 0);
            assert(path[i2] != NULL);
            assert(LT(path[i2][i ].dist, FARAWAY));
            assert(LT(path[i2][i1].dist, FARAWAY));

            if (LE(Max(path[i2][i].dist, path[i2][i1].dist), g->cost[k]))
            {
               graph_edge_del(g, m);
               count++;
               break;
            }
         }
      }
   }
   for(i = 0; i < g->knots; i++)
      if (path[i] != NULL)
         free(path[i]);

   free(path);
   free(used);
   free(term);

   SCIPdebugMessage("%d Edges deleted\n", count * 2);

   assert(graph_valid(g));

   return(count);
}
#if 0
static int hops_test(
   GRAPH* g)
{
   const int max_hops = param_get("MAX_HOPS")->i;

   PATH*   root;
   PATH*   path;
   double* cost;
   int*    mark;
   int     min_hops  = 0;
   double  term_dist;
   int     count     = 0;
   int     retry     = TRUE;
   int     i;
   int     k;
   int     e;

   assert(g != NULL);

   if (!param_get("HOPS_SEP")->i)
      return(0);

   assert(g->layers == 1);

   SCIPdebugMessage("Max-Hops reachability Test: ");
   fflush(stdout);

   cost  = malloc((size_t)g->edges * sizeof(*cost));
   path  = malloc((size_t)g->knots * sizeof(*path));
   root  = malloc((size_t)g->knots * sizeof(*root));
   mark  = malloc((size_t)g->knots * sizeof(*mark));

   assert(cost != NULL);
   assert(path != NULL);
   assert(root != NULL);
   assert(mark != NULL);

   for(i = 0; i < g->knots; i++)
      g->mark[i] = TRUE;

   for(i = 0; i < g->edges; i++)
      cost[i] = 1.0;

   while(retry)
   {
      retry = FALSE;

      for(i = 0; i < g->knots; i++)
      {
         if (!Is_term(g->term[i]) && g->mark[i] && (g->grad[i] == 1))
         {
            graph_edge_del(g, g->outbeg[i]);
            g->mark[i] = FALSE;
            retry      = TRUE;
            count++;
         }
      }
   }

   graph_path_exec(g, FSP_MODE, g->source[0], cost, root);

   for(i = 0; i < g->knots; i++)
   {
      if (!g->mark[i])
         continue;

      if (Is_term(g->term[i]))
      {
         if (root[i].dist > min_hops)
            min_hops = root[i].dist;
      }
      else
      {
         if (root[i].dist >= max_hops)
         {
            for(e = g->outbeg[i]; e != EAT_LAST; e = g->outbeg[i])
               graph_edge_del(g, e);

            g->mark[i] = FALSE;
            count++;
         }
      }
   }
   for(i = 0; i < g->knots; i++)
      mark[i] = g->mark[i];
   while(retry)
   {
      retry = FALSE;

      for(i = 0; i < g->knots; i++)
      {
         if (!Is_term(g->term[i]) && g->mark[i] && (g->grad[i] == 1))
         {
            graph_edge_del(g, g->outbeg[i]);
            g->mark[i] = FALSE;
            retry      = TRUE;
            count++;
         }
      }
   }


   free(mark);
   free(root);
   free(path);
   free(cost);

   SCIPdebugMessage("%d Knots removed (%d/%d)\n", count, min_hops, max_hops);

   assert(graph_valid(g));

   return(count);
}
#endif

static double level1(
   SCIP* scip,
   GRAPH* g)
{
   double fixed   = 0.0;

   assert(g != NULL);

#if 0
   hops_test(g);
#endif

   degree_test(g, &fixed);
   if( 0 )
   {
      bound_test(scip, g);

      degree_test(g, &fixed);

      tt_deletion(g);

      czv_reduction(g, &fixed);

      tt_aggregation(g, &fixed);

      lle_reduction(g);

      degree_test(g, &fixed);

      assert(graph_valid(g));
   }
   SCIPdebugMessage("Reduction Level 1: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}

static double level2(
   SCIP*  scip,
   GRAPH* g)
{
   double fixed   = 0.0;
   int    rerun   = TRUE;

   assert(g      != NULL);

   while(rerun)
   {
      rerun = FALSE;

      if (tt_aggregation(g, &fixed))
         rerun = TRUE;

      while(czv_reduction(g, &fixed))
         rerun = TRUE;
      /*
        while(sd_reduction(g))
        ;
      */
      if (le_reduction(g) > 0)
         rerun = TRUE;

      if (bd3_reduction(g))
         rerun = TRUE;

      if (degree_test(g, &fixed) > 0)
         rerun = TRUE;

      if (nsv_reduction(scip, g, &fixed))
         rerun = TRUE;

   }
   SCIPdebugMessage("Reduction Level 2: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}

static double level3(
   GRAPH* g)
{
   double fixed   = 0.0;

   assert(g != NULL);

   degree_test(g, &fixed);

   // sd_reduction(g);

   degree_test(g, &fixed);
#if 0
   nsv_reduction(g, &fixed);
#else
   bd3_reduction(g);
#endif
   // sd_reduction(g);

   degree_test(g, &fixed);

   //sd_reduction(g);

   degree_test(g, &fixed);

   SCIPdebugMessage("Reduction Level 3: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}

static double level4(
   SCIP* scip,
   GRAPH* g
   )
{
   SCIP_Real timelimit;
   double fixed   = 0.0;
   char    rerun   = TRUE;
   int    i;
   //int    edgebound;
   int    nodebound;
   double*  sddist;
   double*  sdtrans;
   double*  sdrand;
   double* cost;
   double* random;
   int*    heap;
   int*    state;
   char    sd = TRUE;
   char    bd3 = FALSE;
   char    nsv = TRUE;
   char    timebreak = FALSE;
   assert(g != NULL);
   //bound_test(scip, g);

   /* define the miimial number of edge/node eleminations for a reduction test to be continued */
   //edgebound = MAX(g->edges / 100, 5 );
   nodebound = MAX(g->knots / 500, 10);
   //printf("edgebound: %d \n", edgebound );
   printf("nodebound: %d \n", nodebound );

   degree_test(g, &fixed);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   heap  = malloc((size_t)g->knots * sizeof(int));
   state = malloc((size_t)g->knots * sizeof(int));
   sddist = malloc((size_t)g->knots * sizeof(double));
   sdtrans = malloc((size_t)g->knots * sizeof(double));
   sdrand = malloc((size_t)g->knots * sizeof(double));
   cost  = malloc((size_t)g->edges * sizeof(double));
   random  = malloc((size_t)g->edges * sizeof(double));

   while(rerun && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      printf("new presolving run \n");
      rerun = FALSE;

      if( sd )
      {
         for( i = 0; i < 4; i++ ) //TODO 6
         {
            if( sd_reduction(scip, g, sddist, sdtrans, sdrand, cost, random, heap, state) > nodebound )
               rerun = TRUE;

            if( SCIPgetTotalTime(scip) > timelimit )
            {
               timebreak = TRUE;
               break;
            }
         }
         sd = rerun;
      }

      if( timebreak )
         break;

      if( degree_test(g, &fixed) > 0.5 * nodebound )
         rerun = TRUE;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( nsv )
      {
         if( !(nsv_reduction(scip, g, &fixed) > nodebound) )
            nsv = FALSE;
         else
            rerun = TRUE;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( bd3 )
      {
         if( !(bd3_reduction(g) > nodebound) )
            bd3 = FALSE;
         else
            rerun = TRUE;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( degree_test(g, &fixed) > 0.5 * nodebound )
         rerun = TRUE;
   }
   /*
     IDX* curr;
     int es;
     for(es = 0; es < g->edges; es++)
     {
     curr = g->ancestors[es];
     while( curr != NULL )
     {
     printf("afterred %d->%d ancestor: %d->%d \n", g->tail[es], g->head[es], g->orgtail[curr->index], g->orghead[curr->index] );
     curr = curr->parent;
     }
     }
   */
   SCIPdebugMessage("Reduction Level 4: Fixed Cost = %.12e\n",
      fixed);
   free(sddist);
   free(sdtrans);
   free(sdrand);
   free(heap);
   free(state);
   free(cost);
   free(random);

   return(fixed);
}

static double levelm4(
   SCIP* scip,
   GRAPH* g
   )
{
   SCIP_Real timelimit;
   double fixed   = 0.0;
   int    rerun   = TRUE;
   int    i;
   int    numelim;
   int    redbound;
   double*  sddist;
   double*  sdtrans;
   double*  sdrand;
#if 0
   /* These are only used for the sd_reduction_dir.
    * Since the HOP constrained problems are not being solved, then this function is currently redundant.
    */
   double** sd_indist;
   double** sd_intran;
   double** sd_outdist;
   double** sd_outtran;
#endif
   double* cost;
   double* random;
   int*    heap;
   int*    state;
   int*     outterms;
   char    sd = TRUE;
   char    nsv = TRUE;
   char    timebreak = FALSE;

   assert(g != NULL);
   redbound = MAX(g->knots / 500, 8);
   printf("redbound: %d \n", redbound );
   heap        = malloc((size_t)g->knots * sizeof(int));
   state       = malloc((size_t)g->knots * sizeof(int));
   sddist      = malloc((size_t)g->knots * sizeof(double));
   sdtrans     = malloc((size_t)g->knots * sizeof(double));
   sdrand      = malloc((size_t)g->knots * sizeof(double));
#if 0
   sd_indist   = malloc((size_t)g->knots * sizeof(double*));
   sd_intran   = malloc((size_t)g->knots * sizeof(double*));
   sd_outdist  = malloc((size_t)g->knots * sizeof(double*));
   sd_outtran  = malloc((size_t)g->knots * sizeof(double*));
#endif
   cost        = malloc((size_t)g->edges * sizeof(double));
   random        = malloc((size_t)g->edges * sizeof(double));
   outterms    = malloc((size_t)g->knots * sizeof(int));

#if 0
   assert(sd_indist  != NULL);
   assert(sd_intran  != NULL);
   assert(sd_outdist != NULL);
   assert(sd_outtran != NULL);

   for( i = 0; i < g->knots; i++ )
   {
      sd_indist[i]   = malloc((size_t)g->knots * sizeof(double));
      sd_intran[i]   = malloc((size_t)g->knots * sizeof(double));
      sd_outdist[i]  = malloc((size_t)g->knots * sizeof(double));
      sd_outtran[i]  = malloc((size_t)g->knots * sizeof(double));
   }
#endif

   //voronoi_inout(g);

   degree_test_dir(g, &fixed);

   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   while(rerun && !SCIPisStopped(scip))
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      rerun = FALSE;
      if( sd )
      {
         sd = FALSE;
         for(i = 0; i < 2; i++)
         {
#if 0
            if( g->stp_type == STP_HOP_CONS )
               numelim = sd_reduction_dir(g, sd_indist, sd_intran, sd_outdist, sd_outtran, cost, heap, state, outterms);
            else
#endif
               numelim = sd_reduction(scip, g, sddist, sdtrans, sdrand, cost, random, heap, state);
            printf("SD Reduction %d: %d\n", i, numelim);
            if( numelim > redbound )
            {
               rerun = TRUE;
               sd = TRUE;
            }

            if( SCIPgetTotalTime(scip) > timelimit )
            {
               timebreak = TRUE;
               break;
            }
         }
      }

      if( timebreak )
         break;

      if( degree_test_dir(g, &fixed) > redbound / 2 )
         rerun = TRUE;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( nsv )
      {
         nsv = FALSE;
         if( 1 || g->stp_type != STP_HOP_CONS )
         {
            for (i = 0; i < 4; i++)
            {
               numelim = nv_reduction_optimal(g, &fixed);
               printf("NV Reduction %d: %d\n", i, numelim);

               if( SCIPgetTotalTime(scip) > timelimit )
               {
                  timebreak = TRUE;
                  break;
               }

               if( numelim > redbound )
               {
                  rerun = TRUE;
                  nsv = TRUE;
               }
               else
                  break;
            }
         }
      }

      if( timebreak )
         break;

      //if (bd3_reduction(g))
      //rerun = TRUE;

      if (degree_test_dir(g, &fixed) > redbound / 2 )
         rerun = TRUE;
   }
   SCIPdebugMessage("Reduction Level 4: Fixed Cost = %.12e\n", fixed);

#if 0
   for( i = 0; i < g->knots; i++ )
   {
      free(sd_indist[i]);
      free(sd_intran[i]);
      free(sd_outdist[i]);
      free(sd_outtran[i]);
   }
#endif


   free(sddist);
   free(sdtrans);
   free(sdrand);
#if 0
   free(sd_indist);
   free(sd_intran);
   free(sd_outdist);
   free(sd_outtran);
#endif
   free(heap);
   free(state);
   free(cost);
   free(random);
   free(outterms);

   return(fixed);
}

static double level5(
   GRAPH* g)
{
   double fixed   = 0.0;
   int    rerun   = TRUE;
   int    i;
   int    maxruns = 20;

   assert(g != NULL);

   degree_test(g, &fixed);
   tt_aggregation(g, &fixed);

   while(rerun && (--maxruns > 0))
   {
      rerun = FALSE;

      for(i = 0; i < 4; i++)
         /* if (sd_reduction(g))
            rerun = TRUE;
         */
         if (degree_test(g, &fixed) > 0)
            rerun = TRUE;

      if (bd3_reduction(g))
         rerun = TRUE;
   }
   tt_aggregation(g, &fixed);
   degree_test(g, &fixed);

   SCIPdebugMessage("Reduction Level 5: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}


double reduce(
   GRAPH* g,
   int    level,
   SCIP*  scip
   )
{
   double fixed = 0.0;
   printf("Level: %d\n", level);

   assert(g      != NULL);
   assert(level  >= 0 || level == -4);

   if( g->layers != 1 )
      return(0);
   assert(g->fixedges == NULL);
   graph_init_history(g, &(g->orgtail), &(g->orghead), &(g->ancestors));
#if 0
   for( i = 0; i < g->edges; i++ )
   {
      printf("%d->%d ancestor: %d->%d \n", g->tail[i], g->head[i], g->tail[(g->ancestors[i])->index], g->head[(g->ancestors[i])->index] );
      assert((g->ancestors[i])->parent == NULL );
   }
#endif

   /* only use reduction for undirected STP's in graphs */
   printf("type: %d\n", g->stp_type);
   if( 0 && g->stp_type != STP_UNDIRECTED )
      return fixed;

   if( g->stp_type == STP_GRID )
      return fixed;

   if( g->stp_type == STP_HOP_CONS )
      return fixed;

   if( g->stp_type == STP_DEG_CONS )
      return fixed;

   if( g->stp_type != STP_UNDIRECTED )
      level = -4;

   assert(g->layers == 1);

   if (level == 1)
      fixed = level1(scip, g);

   if (level == 2)
      fixed = level1(scip, g) + level2(scip, g);

   if (level == 3)
      fixed = level3(g);

   if (level == 4)
      fixed = level4(scip, g);

   if (level == 5)
      fixed = level5(g);

   if (level == -4)
      fixed = levelm4(scip, g);

   return(fixed);
}
