/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce.c
 * @brief  Reduction tests for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 *
 * This file includes several easy reduction techniques ('degree_test'), bound-based reductions (e.g. 'bound_reduce')
 * and several packages of reduction techniques for different Steiner problem variants.
 *
 * A list of all interface methods can be found in grph.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h)           */

#define REDUCE_C
#define SDSP_BOUND    400          /**< visited edges bound for SDSP test  */
#define BD3_BOUND     400          /**< visited edges bound for BD3 test  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "heur_tm.h"
#include "misc_stp.h"
#include "scip/scip.h"
#include "probdata_stp.h"
#include "prop_stp.h"

/** set entries of (char) array to FALSE */
static
void setTrue(
   char**                arr,
   int                   arrsize
   )
{
   int i;

   assert(arr != NULL);
   assert(arrsize >= 0);

   for( i = 0; i < arrsize; i++ )
   {
      assert(arr[i] != NULL);
      *(arr[i]) = TRUE;
   }
}

/** iterate NV and SL test while at least minelims many contractions are being performed */
static
SCIP_RETCODE nvsl_reduction(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 vnoi,
   SCIP_Real*            nodearrreal,
   SCIP_Real*            fixed,
   int*                  edgearrint,
   int*                  heap,
   int*                  state,
   int*                  vbase,
   int*                  neighb,
   int*                  distnode,
   char*                 visited,
   int*                  nelims,
   int                   minelims
   )
{
   int elims;
   int nvelims;
   int slelims;
   int degelims;
   int totalelims;

   assert(g != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(vnoi != NULL);
   assert(nodearrreal != NULL);
   assert(visited != NULL);
   assert(minelims >= 0);

   *nelims = 0;
   totalelims = 0;

   do
   {
      elims = 0;
      degelims = 0;

      /* NV-reduction */
      SCIP_CALL( nv_reductionAdv(scip, g, vnoi, nodearrreal, fixed, edgearrint, heap, state, vbase, neighb, distnode, &nvelims) );
      elims += nvelims;

      SCIPdebugMessage("NV-reduction (in NVSL): %d \n", nvelims);

      /* SL-reduction */
      SCIP_CALL( sl_reduction(scip, g, vnoi, fixed, heap, state, vbase, neighb, visited, &slelims) );
      elims += slelims;

      SCIPdebugMessage("SL-reduction (in NVSL): %d \n", slelims);

      /* trivial reductions */
      if( elims > 0 )
      {
         if( g->stp_type == STP_PRIZE_COLLECTING || g->stp_type == STP_ROOTED_PRIZE_COLLECTING )
            SCIP_CALL( degree_test_pc(scip, g, fixed, &degelims) );
         else
            SCIP_CALL( degree_test(scip, g, fixed, &degelims) );
      }
      else
      {
         degelims = 0;
      }

      elims += degelims;

      SCIPdebugMessage("Degree Test-reduction (in NVSL): %d \n", degelims);

      totalelims += elims;
   }while( elims > minelims );

   *nelims = totalelims;
   return SCIP_OKAY;
}

#if 0

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
   SCIP* scip,
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
               SCIP_CALL( graph_knot_contract(scip, g, head, tail) );

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

         graph_edge_del(NULL, g, f, FALSE);

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
#endif





#if 0
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

      SCIP_CALL( graph_knot_contract(scip, g, i, g->head[e1]) );

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
                  graph_edge_del(NULL, g, m, FALSE);
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
               graph_edge_del(scip, g, m);
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
#endif

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
            graph_edge_del(scip, g, g->outbeg[i]);
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
               graph_edge_del(scip, g, e);

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
            graph_edge_del(scip, g, g->outbeg[i]);
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

void level0(
   SCIP* scip,
   GRAPH* g)
{
   int k;
   assert(scip != NULL);
   assert(g != NULL);

   for( k = 0; k < g->knots; k++ )
      g->mark[k] = FALSE;

   graph_trail(g, g->source[0]);

   for( k = 0; k < g->knots; k++ )
   {
      if( !g->mark[k] && (g->grad[k] > 0) )
      {
         assert(!Is_term(g->term[k]));

         while( g->inpbeg[k] != EAT_LAST )
            graph_edge_del(scip, g, g->inpbeg[k], TRUE);
      }
   }
}

/** basic reduction package for the STP */
static
SCIP_RETCODE reduceStp(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   SCIP_Bool             dualascent          /**< perform dualascent reductions? */
   )
{
   PATH* vnoi;
   PATH* path;
   SCIP_Real timelimit;
   GRAPH* g = *graph;
   GNODE** gnodearr;
   SCIP_Real*  nodearrreal;
#if 0
   SCIP_Real*  nodearrreal2;
   SCIP_Real*  nodearrreal3;
#endif
   SCIP_Real*  edgearrreal;
   SCIP_Real*  edgearrreal2;
   SCIP_Real   upperbound;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     i;
   int     nnodes;
   int     nedges;
   int     nterms;
   int     danelims;
   int     sdnelims;
   int     lenelims;
   int     sdcnelims;
   int     bd3nelims;
   int     nvslnelims;
   int     brednelims;
   int     degtnelims;
   int     reductbound;
   char* nodearrchar;

   SCIP_Bool    le = TRUE;
   SCIP_Bool    sd = TRUE;
   SCIP_Bool    da = dualascent;
   SCIP_Bool    sdc = TRUE;
   SCIP_Bool    bd3 = TRUE;
   SCIP_Bool    nvsl = TRUE;
   SCIP_Bool    bred = FALSE;
   SCIP_Bool    rerun = TRUE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nterms = g->terms;
   nnodes = g->knots;
   nedges = g->edges;
   degtnelims = 0;

   if( SCIPisLE(scip, (double) g->terms / (double) nnodes, 0.03 ) )
      bred = TRUE;

#if 0
da = FALSE
bd3 = FALSE;
sdc = FALSE;

le = FALSE;
nvsl = FALSE;
bred = FALSE;
#endif

   /* get timelimit parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   if( da )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[i]) ); /*lint !e866*/
      }
   }
   else
   {
      gnodearr = NULL;
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes) );
#if 0
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal3, nnodes) );
#endif
   SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   if( bred || da )
   {
        SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal2, nedges) );
   }
   else
   {
        edgearrreal2 = NULL;
   }

   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = MAX(nnodes / 1000, minelims);

   SCIP_CALL( degree_test(scip, g, fixed, &degtnelims) );

   while( rerun && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      danelims = 0;
      lenelims = 0;
      sdnelims = 0;
      sdcnelims = 0;
      bd3nelims = 0;
      nvslnelims = 0;
      degtnelims = 0;
      brednelims = 0;
      upperbound = -1.0;

      if( le )
      {
         SCIP_CALL( ledge_reduction(scip, g, vnoi, heap, state, vbase, &lenelims) );

         if( lenelims <= reductbound )
            le = FALSE;
         else
            SCIP_CALL( degree_test(scip, g, fixed, &degtnelims) );

         SCIPdebugMessage("le: %d \n", lenelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd )
      {
         SCIP_CALL( sd_red(scip, g, vnoi, edgearrreal, nodearrreal, heap, state, vbase, nodearrint, nodearrint2, edgearrint, &sdnelims) );

         if( sdnelims <= reductbound )
            sd = FALSE;

         SCIPdebugMessage("sd: %d, \n", sdnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sdc )
      {
         SCIP_CALL( sdsp_reduction(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims, SDSP_BOUND) );

         if( sdcnelims <= reductbound )
            sdc = FALSE;

	 SCIPdebugMessage("sdsp: %d \n", sdcnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd || sdc )
         SCIP_CALL( degree_test(scip, g, fixed, &degtnelims) );

      if( bd3 )
      {
         SCIP_CALL( bd3_reduction(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &bd3nelims, BD3_BOUND) );
         if( bd3nelims <= reductbound )
            bd3 = FALSE;
         else
            SCIP_CALL( degree_test(scip, g, fixed, &degtnelims) );

         SCIPdebugMessage("bd3: %d \n", bd3nelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( nvsl )
      {
         SCIP_CALL( nvsl_reduction(scip, g, vnoi, nodearrreal, fixed, edgearrint, heap, state, vbase, nodearrint, NULL, nodearrchar, &nvslnelims, reductbound) );

         if( nvslnelims <= reductbound )
            nvsl = FALSE;

         SCIPdebugMessage("nvsl: %d \n", nvslnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( bred )
      {
         SCIP_CALL( bound_reduce(scip, g, vnoi, edgearrreal, NULL, nodearrreal, edgearrreal2, fixed, &upperbound, heap, state, vbase, &brednelims) );

	 if( brednelims <= 2 * reductbound )
            bred = FALSE;

         SCIPdebugMessage("bnd: %d \n", brednelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( da )
      {
         SCIP_CALL( da_reduce(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, edgearrint, vbase, heap, state, nodearrint2, nodearrchar, &danelims) );

	 if( danelims <= 5 * reductbound )
            da = FALSE;

         SCIPdebugMessage("da: %d \n", danelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      SCIP_CALL( degree_test(scip, g, fixed, &degtnelims) );

      if( (danelims + sdnelims + bd3nelims + nvslnelims + lenelims + brednelims + sdcnelims) <= 2 * reductbound  )
      {
         rerun = FALSE;
#if 0
         if( advanced )
         {
            int     nelims;
                     int     runnum;
   unsigned int seed;
      seed = 0;
   runnum = 0;
            advanced = FALSE;
            sdnelims = 0;
            for( i = 0; i < nnodes; i++ )
               nodearrint[i] = -1;
            for( i = 0; i < 5; i++ )
            {
               SCIP_CALL( sd_reduction(scip, g, nodearrreal, nodearrreal2, nodearrreal3, edgearrreal, edgearrreal2, heap, state, nodearrint, &nelims, runnum, &seed) );
               runnum++;
               sdnelims += nelims;
               if( nelims <= 5 * reductbound )
                  break;
               SCIP_CALL( degree_test(scip, g, fixed, &degtnelims) );
            }

            printf("final sd: %d, \n", sdnelims);
            if( sdnelims > 5 * reductbound )
            {
               rerun = TRUE;
               le = TRUE;
               sd = TRUE;
               sdc = TRUE;
               bd3 = TRUE;
               nvsl = TRUE;
	       da = TRUE;
            }
         }
#endif
      }
   }

   SCIPdebugMessage("Reduction Level 1: Fixed Cost = %.12e\n", *fixed);

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &edgearrreal2);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &edgearrreal);
#if 0
   SCIPfreeBufferArray(scip, &nodearrreal3);
   SCIPfreeBufferArray(scip, &nodearrreal2);
#endif
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &edgearrint);

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBuffer(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }

   return SCIP_OKAY;
}

/** basic reduction package for the (R)PCSTP */
static
SCIP_RETCODE reducePc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   char                  dualascent          /**< perform dual ascent reductions? */
   )
{
   PATH* vnoi;
   PATH* path;
   GRAPH* g = *graph;
   GNODE** gnodearr;
   SCIP_Real* exedgearrreal;
   SCIP_Real* nodearrreal;
   SCIP_Real* exedgearrreal2;
   SCIP_Real timelimit;
   SCIP_Real upperbound;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    edgearrint;
   int*    nodearrint2;
   int     i;
   int     nelims;
   int     nnodes;
   int     nterms;
   int     nedges;
   int     danelims;
   int     sdnelims;
   int     extnedges;
   int     sdcnelims;
   int     bd3nelims;
   int     degnelims;
   int     nvslnelims;
   int     brednelims;
   int     reductbound;
   char*   nodearrchar;
   char    sd = TRUE;
   char    da = dualascent;
   char    sdc = TRUE;
   char    bd3 = TRUE;
   char    nvsl = TRUE;
   char    bred = FALSE;
   char    rerun = TRUE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nterms = g->terms;
   nnodes = g->knots;
   nedges = g->edges;

   /* for PCSPG more memory is necessary */
   if( g->stp_type == STP_ROOTED_PRIZE_COLLECTING || !da )
      extnedges = nedges;
   else
      extnedges = nedges + 2 * (g->terms - 1);

   /* get timelimit parameter*/
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &exedgearrreal, extnedges ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 4 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes + 1) );

   if( SCIPisLE(scip, (double) g->terms / (double) nnodes, 0.03) )
      bred = TRUE;

   if( bred || da )
   {
        SCIP_CALL( SCIPallocBufferArray(scip, &exedgearrreal2, extnedges) );
   }
   else
   {
        exedgearrreal2 = NULL;
   }

   if( da )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, extnedges) );
      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[i]) ); /*lint !e866*/
      }
   }
   else
   {
      gnodearr = NULL;
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, nedges) );
   }


   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = MAX(nnodes / 1000, minelims);

   SCIP_CALL( pcgraphorg(scip, g) );

   SCIP_CALL( degree_test_pc(scip, g, fixed, &degnelims) );

   while( rerun && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      danelims = 0;
      sdnelims = 0;
      sdcnelims = 0;
      bd3nelims = 0;
      nvslnelims = 0;
      degnelims = 0;
      brednelims = 0;
      upperbound = -1.0;

      if( sdc )
      {
         SCIP_CALL( sdsp_reduction(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdcnelims, SDSP_BOUND) );

         if( sdcnelims <= reductbound )
            sdc = FALSE;

         SCIPdebugMessage("SDsp: %d \n", sdcnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd )
      {
         SCIP_CALL( sdpc_reduction(scip, g, vnoi, exedgearrreal, heap, state, vbase, nodearrint, nodearrint2, &sdnelims) );
         if( sdnelims <= reductbound )
            sd = FALSE;

         SCIPdebugMessage("SDpc: %d \n", sdnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      SCIP_CALL( degree_test_pc(scip, g, fixed, &nelims) );
      degnelims += nelims;

      if( bd3 )
      {
         SCIP_CALL( bd3_reduction(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &bd3nelims, BD3_BOUND) );
         if( bd3nelims <= reductbound )
            bd3 = FALSE;

         SCIPdebugMessage("bd3: %d, ", bd3nelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( nvsl )
      {
         SCIP_CALL( nvsl_reduction(scip, g, vnoi, nodearrreal, fixed, edgearrint, heap, state, vbase, nodearrint, nodearrint2, nodearrchar, &nvslnelims, reductbound) );

         if( nvslnelims <= 0.5 * reductbound )
            nvsl = FALSE;

         SCIPdebugMessage("nvsl: %d \n", nvslnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( bred )
      {
         SCIP_CALL( bound_reduce(scip, g, vnoi, exedgearrreal, g->prize, nodearrreal, exedgearrreal2, fixed, &upperbound, heap, state, vbase, &brednelims) );
         if( brednelims <= reductbound )
	    bred = FALSE;

         SCIPdebugMessage("bnd %d \n", brednelims);
	 if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( da )
      {
	 if( g->stp_type == STP_ROOTED_PRIZE_COLLECTING )
	 {
            SCIP_CALL( da_reduce(scip, g, vnoi, gnodearr, exedgearrreal, exedgearrreal2, nodearrreal, edgearrint, vbase, heap, state, nodearrint2, nodearrchar, &danelims) );
	 }
	 else
	 {
	    SCIP_CALL( daPc_reduce(scip, g, vnoi, gnodearr, exedgearrreal, exedgearrreal2, nodearrreal, vbase, heap, edgearrint, state, nodearrchar, &danelims) );
	 }

	 if( danelims <= 2 * reductbound )
            da = FALSE;

         SCIPdebugMessage("da: %d \n", danelims);
      }

      if( degnelims + sdnelims + sdcnelims + bd3nelims + danelims + brednelims + nvslnelims <= reductbound )
         rerun = FALSE;
   }

   SCIP_CALL( pcgraphtrans(scip, g) );
   SCIPdebugMessage("Reduction Level PC 1: Fixed Cost = %.12e\n", *fixed);

   /* free memory */

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBuffer(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }
   SCIPfreeBufferArray(scip, &edgearrint);
   SCIPfreeBufferArrayNull(scip, &exedgearrreal2);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &exedgearrreal);
   SCIPfreeBufferArray(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);

   return SCIP_OKAY;
}

/** reduction package for the MWCSP */
static
SCIP_RETCODE reduceMwcs(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims,           /**< minimal number of edges to be eliminated in order to reiterate reductions */
   char                  dualascent          /**< perform dual ascent reductions? */
   )
{
   GRAPH* g = *graph;
   PATH* vnoi;
   PATH* path;
   GNODE** gnodearr;
   SCIP_Real*  nodearrreal;
   SCIP_Real* edgearrreal;
   SCIP_Real* edgearrreal2;
   SCIP_Real timelimit;
   SCIP_Real upperbound;
   int* state;
   int* vbase;
   int* edgearrint;
   int* nodearrint;
   int* nodearrint2;
   int* nodearrint3;
   int i;
   int nterms;
   int nnodes;
   int nedges;
   int daelims;
   int anselims;
   int nnpelims;
   int degelims;
   int npvelims;
   int redbound;
   int extnedges;
   int bredelims;
   int ansadelims;
   int ansad2elims;
   int chain2elims;
   char* nodearrchar;
   char da = dualascent;
   char ans = TRUE;
   char nnp = TRUE;
   char npv = TRUE;
   char bred = TRUE;
   char rerun = TRUE;
   char ansad = TRUE;
   char ansad2 = TRUE;
   char chain2 = TRUE;
   char* boolarr[16];

   assert(scip != NULL);
   assert(g != NULL);
   assert(fixed != NULL);

   boolarr[0] = &ans;
   boolarr[1] = &nnp;
   boolarr[2] = &npv;
   boolarr[3] = &bred;
   boolarr[4] = &ansad;
   boolarr[5] = &ansad2;
   boolarr[6] = &chain2;
   boolarr[7] = &da;

   nnodes = g->knots;
   nedges = g->edges;
   nterms = g->terms;
   redbound = MAX(nnodes / 1000, minelims);

   if( da )
   {
      extnedges = nedges + 2 * (g->terms - 1);


      SCIP_CALL( SCIPallocBufferArray(scip, &gnodearr, nterms - 1) );
      for( i = 0; i < nterms - 1; i++ )
      {
         SCIP_CALL( SCIPallocBuffer(scip, &gnodearr[i]) ); /*lint !e866*/
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &edgearrint, extnedges) );
   }
   else
   {
      extnedges = nedges;
      edgearrint = NULL;
      gnodearr = NULL;
   }
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint3, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );

   if( bred || da )
   {
           SCIP_CALL( SCIPallocBufferArray(scip, &nodearrreal, nnodes + 1) );
       SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal, extnedges) );
        SCIP_CALL( SCIPallocBufferArray(scip, &edgearrreal2, extnedges) );
   }
   else
   {
        nodearrreal = NULL;
        edgearrreal = NULL;
        edgearrreal2 = NULL;
   }


   SCIP_CALL( pcgraphorg(scip, g) );

   degelims = 0;

   SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );

   while( rerun && !SCIPisStopped(scip) )
   {
      daelims = 0;
      anselims = 0;
      nnpelims = 0;
      degelims = 0;
      npvelims = 0;
      bredelims = 0;
      ansadelims = 0;
      ansad2elims = 0;
      chain2elims = 0;

      upperbound = -1.0;

      if( SCIPgetTotalTime(scip) > timelimit )
         break;
#if 0
      if( chain2 )
      {
         SCIP_CALL( chain2Reduction(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &chain2elims, 300) );

         if( chain2elims <= redbound )
            chain2 = FALSE;

         printf("chain2 delete: %d \n", chain2elims);
      }

      if( npv )
      {
         SCIP_CALL( npvReduction(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &npvelims, 400) );

         if( npvelims <= redbound )
            npv = FALSE;

         printf("npv delete: %d \n", npvelims);
      }
#endif

      if( ans )
      {
         SCIP_CALL( ansReduction(scip, g, fixed, nodearrint2, &anselims) );

         if( anselims <= redbound )
            ans = FALSE;

         SCIPdebugMessage("ans deleted: %d \n", anselims);
      }

      if( ansad )
      {
         SCIP_CALL( ansadvReduction(scip, g, fixed, nodearrint2, &ansadelims) );

	 if( ansadelims <= redbound )
            ansad = FALSE;

         SCIPdebugMessage("ans advanced deleted: %d \n", ansadelims);
      }
#if 0
      if( ansad2 )
      {
         SCIP_CALL( ansadv2Reduction(scip, g, fixed, nodearrint2, &ansad2elims) );

         if( ansad2elims <= redbound )
            ansad2 = FALSE;
	 else
	    SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );

         printf("ans advanced 2 deleted: %d \n", degelims);
      }
#endif
      if( ans || ansad )
	 SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );

#if 1
      if( chain2 )
      {
         SCIP_CALL( chain2Reduction(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &chain2elims, 300) );

         if( chain2elims <= redbound )
            chain2 = FALSE;

         SCIPdebugMessage("chain2 delete: %d \n", chain2elims);
      }

      if( npv )
      {
         SCIP_CALL( npvReduction(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &npvelims, 400) );

         if( npvelims <= redbound )
            npv = FALSE;

         SCIPdebugMessage("npv delete: %d \n", npvelims);
      }
#endif
#if 1 /* org */
      if( nnp )
      {
         SCIP_CALL( nnpReduction(scip, g, fixed, nodearrint, nodearrint2, nodearrint3, &nnpelims, 300, nodearrchar) );

         if( nnpelims <= redbound )
            nnp = FALSE;

         SCIPdebugMessage("nnp deleted: %d \n", nnpelims);
      }
#endif
      if( nnp )
      {
         SCIP_CALL( chain2Reduction(scip, g, vnoi, path, state, vbase, nodearrint, nodearrint2, nodearrint3, &chain2elims, 500) );

         if( chain2elims <= redbound )
            chain2 = FALSE;

         SCIPdebugMessage("chain2 delete: %d \n", chain2elims);

	 if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );
#if 1 /* org */
      if( ansad2 )
      {
         SCIP_CALL( ansadv2Reduction(scip, g, fixed, nodearrint2, &ansad2elims) );

         if( ansad2elims <= redbound )
            ansad2 = FALSE;
	 else
	    SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );

         SCIPdebugMessage("ans advanced 2 deleted: %d \n", degelims);
      }
#endif

      if( bred )
      {
         SCIP_CALL( bound_reduce(scip, g, vnoi, edgearrreal, g->prize, nodearrreal, edgearrreal2, fixed, &upperbound, nodearrint, state, vbase, &bredelims) );

	 if( bredelims <= redbound )
            bred = FALSE;
         else
            SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );

         SCIPdebugMessage("bound_reduce: %d \n", bredelims);
      }

      if( da )
      {
         SCIP_CALL( daPc_reduce(scip, g, vnoi, gnodearr, edgearrreal, edgearrreal2, nodearrreal, vbase, nodearrint, edgearrint, nodearrint2, nodearrchar, &daelims) );

         if( daelims <= 2 * redbound )
            da = FALSE;
         else
            SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );
      }

      if( anselims + nnpelims + chain2elims + bredelims + npvelims + ansadelims + ansad2elims + daelims <= redbound )
	 rerun = FALSE;

      if( !rerun && dualascent )
      {
	 int cnsadvelims = 0;
	 SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );
	 SCIP_CALL( cnsAdvReduction(scip, g, nodearrint2, &cnsadvelims) );

         if( cnsadvelims > 2 * redbound )
	 {
            setTrue(boolarr, 8);
            rerun = TRUE;
            dualascent = FALSE;
            SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );
            SCIPdebugMessage("RELAOD! %d\n\n ", cnsadvelims);
	 }
      }
   }

   SCIP_CALL( degree_test_mw(scip, g, fixed, &degelims) );

   /* go back to the extended graph */
   SCIP_CALL( pcgraphtrans(scip, g) );

   level0(scip, g);

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &edgearrreal2);
   SCIPfreeBufferArrayNull(scip, &edgearrreal);
   SCIPfreeBufferArrayNull(scip, &nodearrreal);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &nodearrchar);
   SCIPfreeBufferArray(scip, &nodearrint3);
   SCIPfreeBufferArray(scip, &nodearrint2);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArrayNull(scip, &edgearrint);

   if( gnodearr != NULL )
   {
      for( i = nterms - 2; i >= 0; i-- )
         SCIPfreeBuffer(scip, &gnodearr[i]);
      SCIPfreeBufferArray(scip, &gnodearr);
   }

   return SCIP_OKAY;
}

/** basic reduction package for the HCDSTP */
static
SCIP_RETCODE reduceHc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   GRAPH* g = *graph;
   PATH* vnoi;
   SCIP_Real*  cost;
   SCIP_Real*  radius;
   SCIP_Real*  costrev;
   SCIP_Real timelimit;
   SCIP_Real upperbound;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    pathedge;
   int     nnodes;
   int     nedges;
   int     redbound;
#if 0
   int     danelims;
#endif
   int     degnelims;
   int     brednelims;
   int     hbrednelims;
   int     hcrnelims;
   int     hcrcnelims;
   char*   nodearrchar;
#if 0
   DOES NOT WORK for HC!
      char    da = !TRUE;
#endif
   char    bred = TRUE;
   char    hbred = TRUE;
   char    rbred = TRUE;
   char    rcbred = TRUE;

   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nnodes = g->knots;
   nedges = g->edges;
   degnelims = 0;
   upperbound = -1.0;
   redbound = MAX(g->knots / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrchar, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &radius, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, 3 * nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pathedge, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, 3 * nnodes) );

   SCIP_CALL( degree_test_hc(scip, g, fixed, &degnelims) );

   while( (bred || hbred || rbred || rcbred) && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( rbred )
      {
         SCIP_CALL( hcrbound_reduce(scip, g, vnoi, cost, costrev, radius, heap, state, vbase, &hcrnelims, pathedge) );
         if( hcrnelims <= redbound )
            rbred = FALSE;
      }

      if( rcbred )
      {
         SCIP_CALL( hcrcbound_reduce(scip, g, vnoi, cost, costrev, radius, -1.0, heap, state, vbase, &hcrcnelims, pathedge, FALSE) );
         if( hcrcnelims <= redbound )
            rcbred = FALSE;
      }

      if( bred )
      {
         SCIP_CALL( bound_reduce(scip, g, vnoi, cost, NULL, radius, costrev, fixed, &upperbound, heap, state, vbase, &brednelims) );
         if( brednelims <= redbound )
            bred = FALSE;
      }

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( hbred )
      {
         SCIP_CALL( hopbound_reduce(scip, g, vnoi, cost, radius, costrev, heap, state, vbase, &hbrednelims) );
         if( hbrednelims <= redbound )
            hbred = FALSE;
	 if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }
#if 0
      if( da )
      {
	 SCIP_CALL( da_reduce(scip, g, vnoi, gnodarr, cost, costrev, radius, vbase, heap, state, pathedge, nodearrchar, &danelims) );
	 if( danelims <= 2 * redbound )
	    da = FALSE;

      }
#endif
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &pathedge);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &costrev);
   SCIPfreeBufferArray(scip, &radius);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &nodearrchar);

   return SCIP_OKAY;
}

/** basic reduction package for the SAP (@todo under construction) */
static
SCIP_RETCODE reduceSap(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data structure */
   SCIP_Real*            fixed,              /**< pointer to store the offset value */
   int                   minelims            /**< minimal number of edges to be eliminated in order to reiterate reductions */
   )
{
   PATH*   vnoi;
   PATH*   path;
   SCIP_Real timelimit;
   GRAPH*  g = *graph;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    nodearrint;
   int*    nodearrint2;
   int     e;
   int     nnodes;
   int     sdnelims;
   int     rptnelims;
   int     degtnelims;
   int     redbound;
   char    sd = !TRUE;
   char    rpt = TRUE;

   nnodes = g->knots;
   redbound = MAX(nnodes / 1000, minelims);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &path, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodearrint2, nnodes) );

   /* @todo change .stp file format for SAP! */
   for( e = 0; e < g->edges; e++ )
      if( SCIPisEQ(scip, g->cost[e], 20000.0) )
	 g->cost[e] = FARAWAY;

   /* @todo: add additional tests */
   while( (sd || rpt) && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( sd )
      {
         SCIP_CALL( sdsp_sap_reduction(scip, g, vnoi, path, heap, state, vbase, nodearrint, nodearrint2, &sdnelims, 300) );
         if( sdnelims <= redbound )
            sd = FALSE;
	 sd = FALSE;
      }

      SCIP_CALL( degree_test_sap(scip, g, fixed, &degtnelims) );

      if( rpt )
      {
         SCIP_CALL( rptReduction(scip, g, fixed, &rptnelims) );
         if( rptnelims <= redbound )
            rpt = FALSE;
      }

      SCIP_CALL( degree_test_sap(scip, g, fixed, &degtnelims) );
   }


   SCIPfreeBufferArray(scip, &vnoi);
   SCIPfreeBufferArray(scip, &path);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &nodearrint);
   SCIPfreeBufferArray(scip, &nodearrint2);

   return SCIP_OKAY;
}


/** reduces the graph */
SCIP_RETCODE reduce(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph structure */
   SCIP_Real*            offset,             /**< pointer to store offset generated by reductions */
   int                   level,              /**< reduction level 0: none, 1: basic, 2: advanced */
   int                   minelims            /**< minimal amount of reductions to reiterate reduction methods */
   )
{
   int stp_type;

   assert((*graph)      != NULL);
   assert((*graph)->fixedges == NULL);
   assert(level  >= 0 && level <= 2);
   assert(minelims >= 0);
   assert((*graph)->layers == 1);

   *offset = 0.0;

   stp_type = (*graph)->stp_type;

   /* initialise ancestor list for each edge */
   SCIP_CALL( graph_init_history(scip, (*graph)) );

   /* initialise shortest path algorithms */
   SCIP_CALL( graph_path_init(scip, (*graph)) );

   level0(scip, (*graph));

   /* if no reduction methods available, return */
   if( (*graph)->stp_type == STP_DEG_CONS )
      return SCIP_OKAY;

   if( level == 1 )
   {
      if( stp_type == STP_PRIZE_COLLECTING || stp_type == STP_ROOTED_PRIZE_COLLECTING )
         SCIP_CALL( reducePc(scip, (graph), offset, minelims, FALSE) );
      else if( stp_type == STP_MAX_NODE_WEIGHT )
         SCIP_CALL( reduceMwcs(scip, (graph), offset, minelims, FALSE) );
      else if( stp_type == STP_HOP_CONS )
         SCIP_CALL( reduceHc(scip, (graph), offset, minelims) );
      else if( stp_type == STP_DIRECTED || stp_type == STP_NODE_WEIGHTS )
         SCIP_CALL( reduceSap(scip, (graph), offset, minelims) );
      else
         SCIP_CALL( reduceStp(scip, (graph), offset, minelims, FALSE) );
   }
   else if( level == 2 )
   {
      if( stp_type == STP_PRIZE_COLLECTING || stp_type == STP_ROOTED_PRIZE_COLLECTING )
         SCIP_CALL( reducePc(scip, (graph), offset, minelims, TRUE) );
      else if( stp_type == STP_MAX_NODE_WEIGHT )
         SCIP_CALL( reduceMwcs(scip, (graph), offset, minelims, TRUE) );
      else if( stp_type == STP_HOP_CONS )
         SCIP_CALL( reduceHc(scip, (graph), offset, minelims) );
      else if( stp_type == STP_DIRECTED || stp_type == STP_NODE_WEIGHTS )
         SCIP_CALL( reduceSap(scip, (graph), offset, minelims) );
      else
         SCIP_CALL( reduceStp(scip, (graph), offset, minelims, TRUE) );
   }
   SCIPdebugMessage("offset : %f \n", *offset);

   assert(graph_valid(*graph));

   graph_path_exit(scip, (*graph));

   return SCIP_OKAY;
}
