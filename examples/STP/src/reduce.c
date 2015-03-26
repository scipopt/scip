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
/*#define PRINT_TMP_PRESOL */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "heur_tm.h"
#include "portab.h"
#include "misc_stp.h"
#include "scip/scip.h"
#include "probdata_stp.h"

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
   SCIP* scip,
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
	       SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
	    }
	    //printf("degcontract: %d %d\n", i1, i);
            SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
	    count++;

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
		  //printf("degcontract: %d %d\n", i2, i);
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  count++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if (Is_term(g->term[i1]) && Is_term(g->term[i2]))
               {
                  if (SCIPisLT(scip, g->cost[e1], g->cost[e2]))
                  {
                     *fixed += g->cost[e1];
		     //printf("degcontract: %d %d\n", i1, i);
		     SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
                     SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
                  }
                  else
                  {
                     *fixed += g->cost[e2];
		     SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e2]) );
		     //printf("degcontract: %d %d\n", i2, i);
                     SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  }
                  count++;

                  break;
               }
               if (Is_term(g->term[i1]) && !Is_term(g->term[i2]) && LE(g->cost[e1], g->cost[e2]))
               {
                  *fixed += g->cost[e1];
		  SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
                  //printf("degcontract: %d %d\n", i1, i);
                  SCIP_CALL( graph_knot_contract(scip, g, i1, i) );
                  count++;

                  break;
               }
               if (Is_term(g->term[i2]) && !Is_term(g->term[i1]) && LE(g->cost[e2], g->cost[e1]))
               {
		  SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e2]) );
                  *fixed += g->cost[e2];
		  //printf("degcontract: %d %d\n", i2, i);
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  count++;

                  break;
               }
               done = FALSE;
            }
            while(FALSE);

            if (done
               && (((i1 < i) && (g->grad[i1] < 3))
                  || ((i2 < i) && (g->grad[i2] < 3))))
               rerun = TRUE;
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", count);
   assert(graph_valid(g));

   return count;
}

/* iterate NV and SL test while at least minelims many contractions are being performed */
static
SCIP_RETCODE nvsl_reduction(
   SCIP* scip,
   GRAPH*  g,
   PATH*   vnoi,
   SCIP_Real* fixed,
   int* heap,
   int* state,
   int* vbase,
   int* nelims,
   int minelims
   )
{
   int elims;
   int nvelims;
   int slelims;
   int degelims;
   int totalelims;
   //IDX* curr1;
   assert(g != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(vnoi != NULL);
   assert(minelims >= 0);

   *nelims = 0;
   totalelims = 0;
   /*
     curr1 = g->ancestors[1482];
     while( curr1 != NULL )
     {
     printf("1482 edge: %d\n", curr1->index);
     curr1 = curr1->parent;

     }
   */
   do
   {

      elims = 0;

      /* NV-reduction */
      SCIP_CALL( nv_reduction(scip, g, vnoi, fixed, heap, state, vbase, &nvelims) );
      elims += nvelims;

      SCIPdebugMessage("NV-reduction (in NVSL): %d \n", nvelims);
      /*
        printf("aftnv 1482 is edge: %d->%d:\n", g->tail[1482], g->head[1482] );
        curr1 = g->ancestors[1482];
        while( curr1 != NULL )
        {
        printf("aftnv 1482 edge: %d\n", curr1->index);
        curr1 = curr1->parent;
        //  if( i++ > 30 )
	// assert(0);
        }
      */
      /* SL-reduction */
      SCIP_CALL( sl_reduction(scip, g, vnoi, fixed, heap, state, vbase, &slelims) );
      elims += slelims;

      SCIPdebugMessage("SL-reduction (in NVSL): %d \n", slelims);

      /*printf("aftsl 1482 is edge: %d->%d:\n", g->tail[1482], g->head[1482] );
        curr1 = g->ancestors[1482];
        while( curr1 != NULL )
        {
        printf("aftsl 1482 edge: %d\n", curr1->index);
        curr1 = curr1->parent;
        }
      */
      /* trivial reduction */
      if( elims > 0 )
         degelims = degree_test(scip, g, fixed);
      else
         degelims = 0;

      elims += degelims;
      /*
        printf("aftdeg 1482 is edge: %d->%d:\n", g->tail[1482], g->head[1482] );
        curr1 = g->ancestors[1482];
        while( curr1 != NULL )
        {
        printf("aftdeg 1482 edge: %d\n", curr1->index);
        curr1 = curr1->parent;
        }
      */
      SCIPdebugMessage("Degree Test-reduction (in NVSL): %d \n", degelims);

      totalelims += elims;
      //printf("elims %d\n ", elims);
   }while( elims > minelims );

   *nelims = totalelims;
   return SCIP_OKAY;
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
SCIP_RETCODE bound_reduce(
   SCIP*  scip,
   GRAPH* graph,
   PATH* vnoi,
   double* cost,
   double* radius,
   double* costrev,
   int* heap,
   int* state,
   int* vbase,
   int* nelims,
   int fixed
   )
{
   SCIP_HEUR** heurs;
   SCIP_HEURDATA* tmheurdata;
   SCIP_Real  radiisum;
   SCIP_Real  obj;
   int* perm;
   int* result;
   int* starts;
   int* radedges;
   int e;
   int etemp;
   int radedge;
   int k;
   int nterms;
   int v;
   int r;
   int runs;
   int skip;
   int nnodes;
   int nedges;
   int nheurs;
   int nskips;
   int best_start = 0;
   unsigned int seed = 0;
   char* blocked;
   *nelims = 0;
   nedges = graph->edges;
   nnodes = graph->knots;

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &result, nedges) );
   /* debug, TODO delete */
   SCIP_CALL( SCIPallocBufferArray(scip, &radedges, nnodes) );
   /* debug, TODO delete */
   SCIP_CALL( SCIPallocBufferArray(scip, &blocked, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &starts, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &perm, nnodes) );
   e = 0;
   nterms = 0;
   for( k = 0; k < nnodes; k++ )
   {
      perm[k] = k;
      graph->mark[k] = (graph->grad[k] > 0);
      if( graph->mark[k] )
      {
         e++;
	 if( Is_term(graph->term[k]) )
            nterms++;
      }
   }
   runs = MIN(e, 50);
   r = 0;
   assert(graph->source[0] >= 0);
   if( graph->mark[graph->source[0]] )
      starts[r++] = graph->source[0];

   SCIPpermuteIntArray(perm, 0, nnodes, &(seed));

   /* use terminals (randomly permutated) as starting points for TM heuristic */
   for( k = 0; k < nnodes; k++ )
   {
      if( r >= runs || r >= nterms )
         break;

      if( Is_term(graph->term[perm[k]]) && graph->mark[perm[k]] && perm[k] != graph->source[0] )
      {
         starts[r++] = perm[k];
      }
   }

   /* still empty slots in start array? */
   if( r < runs )
   {
      /* fill empty slots randomly */
      for( k = 0; k < nnodes; k++ )
      {
         if( r >= runs )
            break;
         if( !Is_term(graph->term[perm[k]]) && graph->mark[perm[k]] )
         {
            //printf(" perm2[%d): %d \n", r, perm[k]);
            starts[r++] = perm[k];
         }
      }
   }

   SCIPfreeBufferArray(scip, &perm);

   for( e = 0; e < nedges; e++ )
   {
      blocked[e] = FALSE;
      blocked[e + 1] = FALSE;
      result[e] = UNKNOWN;
      result[e + 1] = UNKNOWN;
      if( graph->ieat[e] == EAT_FREE )
      {
         assert(graph->oeat[e] == EAT_FREE);
	 assert(graph->oeat[flipedge(e)] == EAT_FREE);
	 assert(graph->ieat[flipedge(e)] == EAT_FREE);
         cost[e] = FARAWAY;
	 costrev[e] = FARAWAY;
      }
      else
      {
	 assert(graph->oeat[e] != EAT_FREE);
	 assert(graph->oeat[flipedge(e)] != EAT_FREE);
	 assert(graph->ieat[flipedge(e)] != EAT_FREE);
	 cost[e] = graph->cost[e];
	 costrev[e] = graph->cost[flipedge(e)];
      }
   }

   voronoi_radius(scip, graph, vnoi, radius, cost, costrev, vbase, heap, state, radedges);

   SCIPsortRealInt(radius, radedges, nnodes);
   radiisum = 0.0;
   e = 0;

   for( e = 0; e < nterms - 2; e++ )
   {
      assert( SCIPisGT(scip, FARAWAY, radius[e]) );
      radiisum += radius[e];
      radedge = radedges[e];
      assert(radedge != UNKNOWN);
      blocked[radedge] = TRUE;
      blocked[flipedge(radedge)] = TRUE;
      k = graph->tail[radedge];
      assert(vbase[k] != vbase[graph->head[radedge]]);
      while( k != vbase[k] )
      {
         assert(vnoi[k].edge >= 0);
         blocked[vnoi[k].edge] = TRUE;
         blocked[flipedge(vnoi[k].edge)] = TRUE;
         k = graph->tail[vnoi[k].edge];
      }
      // printf("rad: %f\n", radius[e]);
   }

   /* get TM heuristic data */
   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);
   for( k = 0; k < nheurs; k++ )
      if( strcmp(SCIPheurGetName(heurs[k]), "TM") == 0 )
         break;
   assert(k < nheurs);
   tmheurdata = SCIPheurGetData(heurs[k]);

   SCIP_CALL( do_layer(scip, tmheurdata, graph, starts, &best_start, result, runs, graph->source[0], cost, costrev, 0.0) );
   obj = fixed;
   for( e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         obj += graph->cost[e];
   //printf("obj: %f\n", obj);

   nskips = 0;
   for( k = 0; k < nnodes; k++ )
   {
      e = graph->outbeg[k];
      while( e != EAT_LAST )
      {
         assert(e >= 0);
         //printf("e: %d->%d\n", graph->tail[e], graph->head[e]);
         if( SCIPisGT(scip, graph->cost[e] + vnoi[graph->head[e]].dist + vnoi[graph->tail[e]].dist + radiisum, obj) )
         {
            //printf("e: %d->%d outbeg: %d \n", graph->tail[e], graph->head[e], graph->oeat[e]);
            /* debug, TODO delete */
            skip = FALSE;
            if( blocked[e] )
            {
               skip = TRUE;
               printf("blocked \n");
            }
            if( !skip )
            {
               v = graph->head[e];
	       assert(vbase[v] >= 0);
               while( v != vbase[v] )
               {
                  if( vnoi[v].edge < 0 )
                  {
                     skip = TRUE;
                     printf("edge misses! \n");
                     break;
                  }
                  v = graph->tail[vnoi[v].edge];
               }
               if( !skip )
               {
                  v = graph->tail[e];
		  assert(vbase[v] >= 0);
                  while( v != vbase[v] )
                  {
                     if( vnoi[v].edge < 0 )
                     {
                        skip = TRUE;
                        printf("edge misses! \n");
                        break;
                     }
                     v = graph->tail[vnoi[v].edge];
                  }
               }
            }
            if( skip )
            {
               nskips++;
               e = graph->oeat[e];
            }
            else
            {
               (*nelims)++;
               SCIPindexListNodeFree(scip, &((graph->ancestors)[e]));
               SCIPindexListNodeFree(scip, &((graph->ancestors)[Edge_anti(e)]));
               etemp = graph->oeat[e];
               graph_edge_del(graph, e);
               e = etemp;
            }
         }
         else
         {
            e = graph->oeat[e];
         }
      }
   }

   printf("nelims in bound reduce: %d, nskips: %d \n", *nelims, nskips);
   SCIPfreeBufferArray(scip, &radedges);
   SCIPfreeBufferArray(scip, &result);
   SCIPfreeBufferArray(scip, &starts);
   SCIPfreeBufferArray(scip, &blocked);

   assert(graph_valid(graph));
   return SCIP_OKAY;
}
#if 0
static
int bound_test(
   SCIP*  scip,
   GRAPH* graph,
   int* elimins
   )
{
   PATH**      heurpath;
   PATH**      path;
   PATH*       pathtoterm;
   SCIP_Real*  edgecost;
   SCIP_Real*  edgecostrev;
   SCIP_Real*  distance;
   SCIP_Real*  radius;
   SCIP_Real** termdist;
   int*        result;
   int*        vregion;
   int*        radiushops;
   int*        heap;
   int*        state;
   int*        pred;
   int*        terms;
   int         source;
   int         termcount;
   int         i;
   int         e;
   int         nnodes;
   int         nedges;

   SCIP_Real   obj;

   int     closeterms[3] = {-1, -1, -1};
   double  closetermsdist[3] = {FARAWAY, FARAWAY, FARAWAY};
   int     closetermshops[3] = {-1, -1, -1};
   double  lowerbound = FARAWAY;
   int     hopsbound = 0;

   assert(scip != NULL);
   assert(graph != NULL);
   (*elimins) = 0;
   nnodes = graph->knots;
   nedges = graph->edges;

   path = malloc((size_t)nnodes * sizeof(PATH*));
   SCIP_CALL(SCIPallocBufferArray(scip, &heurpath, nnodes));

   assert(path != NULL);
   assert(heurpath != NULL);

   pathtoterm = malloc((size_t)nnodes * sizeof(PATH));

   assert(pathtoterm != NULL);

   edgecost = malloc((size_t)nedges * sizeof(double));
   edgecostrev = malloc((size_t)nedges * sizeof(double));

   assert(edgecost != NULL);
   assert(edgecostrev != NULL);

   radius = malloc((size_t)nnodes * sizeof(double));
   distance = malloc((size_t)nnodes * sizeof(double));
   termdist = malloc((size_t)nnodes * sizeof(double*));

   assert(distance != NULL);
   assert(radius != NULL);
   assert(termdist != NULL);

   vregion = malloc((size_t)nnodes * sizeof(int));
   radiushops = malloc((size_t)nnodes * sizeof(int));

   assert(vregion != NULL);
   assert(radiushops != NULL);

   heap  = malloc((size_t)nnodes * sizeof(int));
   state = malloc((size_t)nnodes * sizeof(int));

   assert(heap != NULL);
   assert(state != NULL);

   pred = malloc((size_t)nedges * sizeof(int));
   terms = malloc((size_t)graph->terms * sizeof(int));

   assert(pred != NULL);
   assert(terms != NULL);

   result = malloc((size_t)nedges * sizeof(int));

   assert(result != NULL);

   for( i = 0; i < nnodes; i++ )
   {
      path[i] = NULL;
      heurpath[i] = NULL;
      termdist[i] = NULL;
   }

   for( e = 0; e < nedges; e++ )
   {
      if( graph->cost[e] != FARAWAY )
         edgecost[e] = 1;
      else
         edgecost[e] = graph->cost[e];

      if( graph->cost[Edge_anti(e)] != FARAWAY )
         edgecostrev[e] = 1;
      else
         edgecostrev[e] = graph->cost[Edge_anti(e)];
   }

   /* Not currently in use */
#if 0
   /* storing the current type of the stp */
   temptype = graph->stp_type;
   graph->stp_type = STP_DIRECTED;

   SCIP_CALL( SCIPtmHeur(scip, graph, heurpath, edgecost, edgecostrev, result) );

   /* resetting the stp type */
   graph->stp_type = temptype;
#endif

   obj = 0.0;

   for( e = 0; e < nedges; e++ )
      if( result[e] == CONNECT )
         obj += edgecost[e];

   termcount = 0;
   for(i = 0; i < nnodes; i++)
   {
      if( Is_term(graph->term[i]) )
      {
         terms[termcount] = i;
         termcount++;
      }
      graph->mark[i] = (graph->grad[i] > 0);
      termdist[i] = malloc((size_t)nnodes * sizeof(double));
   }

   source = graph->source[0];

   voronoi_hop(graph, edgecost, distance, radius, pathtoterm, vregion, heap, state, pred, radiushops);

   SCIPsortRealInt(radius, radiushops, nnodes);

   calculate_distances(graph, path, edgecost, BSP_MODE);

   printf("Hop Limit: %d\n", graph->hoplimit);

   /* test each node i */
   for( i = 0; i < nnodes; i++ )
   {
      if( Is_term(graph->term[i]) || vregion[i] == source )
         continue;

      if( graph->grad[i] == 0 )
         continue;

      if( vregion[i] < 0 )
      {
         while( graph->inpbeg[i] != EAT_LAST )
         {
            e = graph->inpbeg[i];
            //if( Is_term(graph->term[graph->tail[e]]) )
            //printf("Found terminal: %d\n", graph->tail[e]);
            graph_edge_del(graph, e);
            SCIPindexListNodeFree(scip, &(graph->ancestors[e]));
            graph->ancestors[e] = NULL;
            (*elimins)++;
         }
         continue;
      }

      get_close_terms(path, closetermsdist, closetermshops, closeterms, vregion, terms, termcount, i);

      /* computing the lower bound for node i */
      lowerbound = compute_node_lb(radius, closetermsdist, closetermshops, closeterms, radiushops, termcount,
         graph->grad[i], graph->source[0], graph->stp_type, &hopsbound);

      //if( i % 1000 == 0 )
      printf("node: %d, lowerbound: %f, grad: %d, vregion: %d\n", i, lowerbound, graph->grad[i], vregion[i]);

      if( GT(lowerbound, graph->hoplimit) )
      {
         while( graph->inpbeg[i] != EAT_LAST )
         {
            e = graph->inpbeg[i];
            //if( Is_term(graph->term[graph->tail[e]]) )
            //printf("Found terminal: %d\n", graph->tail[e]);
            graph_edge_del(graph, e);
            SCIPindexListNodeFree(scip, &(graph->ancestors[e]));
            graph->ancestors[e] = NULL;
            (*elimins)++;
         }

         voronoi_hop(graph, edgecost, distance, radius, pathtoterm, vregion, heap, state, pred, radiushops);

         SCIPsortRealInt(radius, radiushops, nnodes);

         calculate_distances(graph, path, edgecost, BSP_MODE);

         continue;
      }

#if 0
      /* computing the lowerbound for the inclusion of a single edge */
      /* Resetting the lowerbound */
      lowerbound -= closetermsdist[1];
      if( graph->grad[i] >= 3 )
      {
         lowerbound -= closetermsdist[2];
         lowerbound += radius[termcount - 4];
      }


      for( j = graph->inpbeg[i]; j != EAT_LAST; j = graph->ieat[j] )
      {
         k = graph->tail[j];
         if( LT(edgecost[j], edgecostrev[j]) )
            tempcost = edgecost[j] + path[vregion[k]][k].dist;
         else
            tempcost = edgecostrev[j] + path[vregion[k]][k].dist;
         lowerbound += tempcost;

         //printf("Edge - lowerbound: %f, bestbound: %f\n", lowerbound, bestbound);
         if( GT(lowerbound, graph->hoplimit) )
         {
            graph_edge_del(graph, j);
            SCIPindexListNodeFree(scip, &(graph->ancestors[j]));
            graph->ancestors[j] = NULL;
            (*elimins)++;

            /* computing the voronoi regions inward to a node */
            voronoi_hop(graph, edgecost, distance, radius, pathtoterm, vregion, heap, state, pred, radiushops);

            /* sorting the radius values */
            SCIPsortRealInt(radius, radiushops, graph->knots);

            /* computing the shortest paths from each terminal to every other node */
            //calculate_distances(g, path, graph->cost, FSP_MODE);

            continue;
         }

         lowerbound -= tempcost;
      }
#endif
   }

   free(result);

   free(terms);
   free(pred);
   free(state);
   free(heap);

   free(vregion);

   for( i = 0; i < nnodes; i++ )
      free(termdist[i]);

   free(termdist);

   free(distance);
   free(radius);

   free(edgecostrev);
   free(edgecost);

   free(pathtoterm);

   for( i = 0; i < nnodes; i++ )
      SCIPfreeBufferArrayNull(scip, &heurpath[i]);

   SCIPfreeBufferArray(scip, &heurpath);

   for( i = 0; i < nnodes; i++ )
      free(path[i]);

   free(path);

   return 1;
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
#endif
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

#if 0
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

   degree_test(scip, g, &fixed);
   if( 0 )
   {
      //bound_test(scip, g);

      degree_test(scip, g, &fixed);

      tt_deletion(g);

      /*czv_reduction(g, &fixed);*/

      tt_aggregation(scip, g, &fixed);

      lle_reduction(g);

      degree_test(scip, g, &fixed);

      assert(graph_valid(g));
   }
   SCIPdebugMessage("Reduction Level 1: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}
#if 0
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

      /*  if (bd3_reduction(g))
          rerun = TRUE;
      */
      if (degree_test(scip, g, &fixed) > 0)
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

   degree_test(scip, g, &fixed);

   // sd_reduction(g);

   degree_test(scip, g, &fixed);
#if 0
   nsv_reduction(g, &fixed);

   bd3_reduction(g);
#endif
   // sd_reduction(g);

   degree_test(scip, g, &fixed);

   //sd_reduction(g);

   degree_test(scip, g, &fixed);

   SCIPdebugMessage("Reduction Level 3: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}
#endif

static
SCIP_RETCODE level4(
   SCIP* scip,
   GRAPH** graph,
   SCIP_Real* fixed,
   int minelims
   )
{
   PATH* vnoi;
   SCIP_Real timelimit;
   GRAPH* g = *graph;
   double*  sddist;
   double*  sdtrans;
   double*  sdrand;
   double* cost;
   double* random;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    knotexamined;
   int     i;
   int     nelims;
   int     nnodes;
   int     nedges;
   int     runnum;
   int     sdnelims;
   int     lenelims;
   int     bd3nelims;
   int     nsvnelims;
   int     nvslnelims;
   int     brednelims;
   int     degtnelims;
   int     reductbound;

   char    le = TRUE;
   char    sd = TRUE;
   char    bd3 = TRUE;
   char    nsv = TRUE;
   char    nvsl = TRUE;
   char    bred;
   char    rerun = TRUE;

#ifdef PRINT_TMP_PRESOL
   SCIP_PROBDATA* probdata;
   const char*   probname;
   char   presolvetempfile[SCIP_MAXSTRLEN];
   double offset;
#endif

#ifdef PRINT_TMP_PRESOL
   probdata = SCIPgetProbData(scip);
   probname = SCIPgetProbName(scip);
   (void)SCIPsnprintf(presolvetempfile, SCIP_MAXSTRLEN, "presol/%s-presolve-tmp.stp", probname);
   offset = SCIPprobdataGetOffset(scip);
#endif
   assert(scip != NULL);
   assert(g != NULL);
   assert(minelims >= 0);

   nnodes = g->knots;
   nedges = g->edges;

   if( SCIPisLE(scip, (double) g->terms / (double) g->knots, 0.03 ) )
      bred = TRUE;
   else
      bred = FALSE;

#ifdef PRINT_TMP_PRESOL
   SCIPprobdataSetGraph(probdata, g);
   SCIPprobdataSetOffset(probdata, offset + *fixed);
   /* Writing the problem to a temporary file */
   SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
   /* get timelimit parameter*/
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &knotexamined, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sddist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdtrans, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdrand, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &random, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );

   for( i = 0; i < g->knots; i++ )
      knotexamined[i] = -1;

   runnum = 0;

   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = MAX(g->knots / 500, minelims);
   //printf("BOUND: %d \n", reductbound);

   degree_test(scip, g, fixed);

   if( bred )
   {
      SCIP_CALL( bound_reduce(scip, g, vnoi, cost, sddist, random,  heap, state, vbase, &brednelims, *fixed) );
      //printf("bound reduction1: %d \n \n", brednelims);
   }

   while( rerun && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      sdnelims = 0;
      bd3nelims = 0;
      nsvnelims = 0;
      nvslnelims = 0;
      degtnelims = 0;
      lenelims = 0;

      if( nvsl )
      {
         SCIP_CALL( nvsl_reduction(scip, g, vnoi, fixed, heap, state, vbase, &nvslnelims, reductbound) );

         if( nvslnelims <= 0.3 * reductbound )
            nvsl = FALSE;

         //printf("nvsl: %d \n", nvslnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( le )
      {
         SCIP_CALL( ledge_reduction(scip, g, vnoi, heap, state, vbase, &lenelims) );

         if( lenelims <= 0.5 * reductbound )
            le = FALSE;

         //printf("lenelims: %d \n", lenelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd )
      {
	 nelims = 0;
         for( i = 0; i < 6; i++ )
         {
            SCIP_CALL( sd_reduction(scip, g, sddist, sdtrans, sdrand, cost, random, heap, state, knotexamined, &nelims, runnum) );
            runnum++;
	    sdnelims += nelims;
         }

#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + *fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
         //printf("sdnelims: %d \n", sdnelims);
         if( sdnelims <= reductbound )
            sd = FALSE;
      }

      degtnelims += degree_test(scip, g, fixed);

#ifdef PRINT_TMP_PRESOL
      SCIPprobdataSetGraph(probdata, g);
      SCIPprobdataSetOffset(probdata, offset + *fixed);
      /* Writing the problem to a temporary file */
      SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( nsv )
      {
         SCIP_CALL( nsv_reduction(scip, g, cost, fixed, &nsvnelims) );
         if( nsvnelims <= reductbound )
            nsv = FALSE;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + *fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
      }

      if( bd3 )
      {
         SCIP_CALL( bd3_reduction(scip, g, sddist, sdtrans, heap, state, &bd3nelims) );
         if( bd3nelims <= reductbound )
            bd3 = FALSE;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + *fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
      }

      degtnelims += degree_test(scip, g, fixed);

      if( (sdnelims + bd3nelims + nsvnelims + nvslnelims + degtnelims + lenelims) <= reductbound )
         rerun = FALSE;
   }

   /*
     brednelims = le_reduction(g);
     printf("le_reduction: %d \n\n", brednelims);
     brednelims =  tt_deletion(g);
     printf("t_reduction: %d \n\n", brednelims);
     brednelims = tt_aggregation(
     g,
     fixed);
     printf("ttag_reduction: %d \n\n", brednelims);
   */
   SCIPdebugMessage("Reduction Level 4: Fixed Cost = %.12e\n", *fixed);
   /*printf("Total Fixed: %f\n", fixed);*/

   /* free memory */
   SCIPfreeBufferArray(scip, &sddist);
   SCIPfreeBufferArray(scip, &sdtrans);
   SCIPfreeBufferArray(scip, &sdrand);
   SCIPfreeBufferArray(scip, &knotexamined);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &random);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);

   return SCIP_OKAY;
}
#if 1
static double levelm1(
   SCIP* scip,
   GRAPH* g
   )
{
   double fixed = 0.0;
   int degnelims;
   SCIP_CALL( degree_test_dir(scip, g, &fixed, &degnelims) );
   printf("level1 nelims: %d \n", degnelims);
   return fixed;
}
#endif
static double levelm4(
   SCIP* scip,
   GRAPH* g
   )
{
   SCIP_Real timelimit;
   double fixed   = 0.0;
   int    rerun   = TRUE;
   int    i;
   int    nelims = 0;
   int    redbound;
   double*  sddist;
   double*  sdtrans;
   double*  sdrand;
#if 1
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
   int*    knotexamined;
   int*     outterms;
   int     bd3nelims;
   int     runnum = 0;
   int     sdnelims;
   int     nvnelims;
   int     degnelims;
   char    sd = TRUE;
   char    nsv = TRUE;
   char    bd3 = TRUE;
   char    timebreak = FALSE;


#ifdef PRINT_TMP_PRESOL
   SCIP_PROBDATA* probdata;
   const char*   probname;
   char   presolvetempfile[SCIP_MAXSTRLEN];
   double offset;
#endif

   assert(scip != NULL);

#ifdef PRINT_TMP_PRESOL
   probdata = SCIPgetProbData(scip);
   probname = SCIPgetProbName(scip);
   (void)SCIPsnprintf(presolvetempfile, SCIP_MAXSTRLEN, "presol/%s-presolve-tmp.stp", probname);
   offset = SCIPprobdataGetOffset(scip);
#endif

   assert(g != NULL);
   redbound = MAX(g->knots / 500, 8);
   /* redbound = 0;
    */
   printf("redbound: %d \n", redbound );
   heap        = malloc((size_t)g->knots * sizeof(int));
   state       = malloc((size_t)g->knots * sizeof(int));
   knotexamined = malloc((size_t)g->knots * sizeof(int));
   sddist      = malloc((size_t)g->knots * sizeof(double));
   sdtrans     = malloc((size_t)g->knots * sizeof(double));
   sdrand      = malloc((size_t)g->knots * sizeof(double));
#if 1
   sd_indist   = malloc((size_t)g->knots * sizeof(double*));
   sd_intran   = malloc((size_t)g->knots * sizeof(double*));
   sd_outdist  = malloc((size_t)g->knots * sizeof(double*));
   sd_outtran  = malloc((size_t)g->knots * sizeof(double*));
#endif
   cost        = malloc((size_t)g->edges * sizeof(double));
   random        = malloc((size_t)g->edges * sizeof(double));
   outterms    = malloc((size_t)g->knots * sizeof(int));

   for( i = 0; i < g->knots; i++ )
      knotexamined[i] = -1;

#if 1
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

#ifdef PRINT_TMP_PRESOL
   SCIPprobdataSetGraph(probdata, g);
   SCIPprobdataSetOffset(probdata, offset + fixed);
   /* Writing the problem to a temporary file */
   SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
   if( g->stp_type == STP_HOP_CONS )
   {
#if 0
      do
      {
         printf("Bound test\n");
         bound_test(scip, g, &nelims);
         printf("Num elimins: %d\n", nelims);

         degree_test_dir(g, &fixed);

         nelims += nv_reduction_optimal(g, &fixed, runnum);

         nelims += sd_reduction_dir(g, sd_indist, sd_intran, sd_outdist, sd_outtran, cost, heap, state, outterms);
         printf("Num elimins: %d\n", nelims);
      } while(nelims > 0);
#endif

      rerun = FALSE;
   }
   else
   {
      SCIP_CALL( degree_test_dir(scip, g, &fixed, &degnelims) );

      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   }

#ifdef PRINT_TMP_PRESOL
   SCIPprobdataSetOffset(probdata, offset + fixed);
   SCIPprobdataSetGraph(probdata, g);
   /* Writing the problem to a temporary file */
   SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif

   while(rerun && !SCIPisStopped(scip))
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      sdnelims = 0;
      nvnelims = 0;
      degnelims = 0;
      bd3nelims = 0;
      rerun = FALSE;
      if( sd )
      {
         sd = FALSE;
         for(i = 0; i < 2; i++)
         {
#if 0
            if( g->stp_type == STP_HOP_CONS )
               nelims = sd_reduction_dir(g, sd_indist, sd_intran, sd_outdist, sd_outtran, cost, heap, state, outterms);
            else
#endif
               SCIP_CALL( sd_reduction(scip, g, sddist, sdtrans, sdrand, cost, random, heap, state, knotexamined, &nelims, runnum) );
            runnum++;
            sdnelims += nelims;
            printf("SD Reduction %d: %d\n", i, nelims);

            if( SCIPgetTotalTime(scip) > timelimit )
            {
               timebreak = TRUE;
               break;
            }


#ifdef PRINT_TMP_PRESOL
            SCIPprobdataSetOffset(probdata, offset + fixed);
            SCIPprobdataSetGraph(probdata, g);
            /* Writing the problem to a temporary file */
            SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
            if( nelims > redbound )
               sd = TRUE;
            else
               break;
         }
      }

      if( timebreak )
         break;

      SCIP_CALL( degree_test_dir(scip, g, &fixed, &i) );
      degnelims += i;

#ifdef PRINT_TMP_PRESOL
      /* Writing the problem to a temporary file */
      SCIPprobdataSetOffset(probdata, offset + fixed);
      SCIPprobdataSetGraph(probdata, g);
      SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( nsv )
      {
         nsv = FALSE;
         if( 1 || g->stp_type != STP_HOP_CONS )
         {

            for( i = 0; i < 4; i++ )
            {
               SCIP_CALL( nv_reduction_optimal(scip, g, &fixed, &nelims, runnum) );
               runnum++;
	       nvnelims += nelims;
               printf("NV Reduction %d: %d\n", i, nvnelims);

               if( SCIPgetTotalTime(scip) > timelimit )
               {
                  timebreak = TRUE;
                  break;
               }

               if( nelims > redbound )
               {
                  nsv = TRUE;
               }
               else
                  break;
            }
#ifdef PRINT_TMP_PRESOL
            /* Writing the problem to a temporary file */
            SCIPprobdataSetOffset(probdata, offset + fixed);
            SCIPprobdataSetGraph(probdata, g);
            SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
         }
      }

      if( timebreak )
         break;

      if( bd3 )
      {
         SCIP_CALL( bd3_reduction(scip, g, sddist, sdtrans, heap, state, &bd3nelims) );
         if( bd3nelims <= redbound )
            bd3 = FALSE;

         printf("BD3 Reduction: %d\n", bd3nelims);

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
      }

      SCIP_CALL( degree_test_dir(scip, g, &fixed, &i) );
      degnelims += i;

#ifdef PRINT_TMP_PRESOL
      /* Writing the problem to a temporary file */
      SCIPprobdataSetOffset(probdata, offset + fixed);
      SCIPprobdataSetGraph(probdata, g);
      SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
      if( degnelims + sdnelims + nvnelims + bd3nelims > redbound )
	 rerun = TRUE;
   }
   SCIPdebugMessage("Reduction Level 4: Fixed Cost = %.12e\n", fixed);

#if 1
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
#if 1
   free(sd_indist);
   free(sd_intran);
   free(sd_outdist);
   free(sd_outtran);
#endif
   free(knotexamined);
   free(heap);
   free(state);
   free(cost);
   free(random);
   free(outterms);

   return(fixed);
}


#if 0
static double level5(
   GRAPH* g)
{
   double fixed   = 0.0;
   int    rerun   = TRUE;
   int    i;
   int    maxruns = 20;

   assert(g != NULL);

   degree_test(scip, g, &fixed);
   tt_aggregation(g, &fixed);

   while(rerun && (--maxruns > 0))
   {
      rerun = FALSE;

      for(i = 0; i < 4; i++)
         /* if (sd_reduction(g))
            rerun = TRUE;
         */
         if (degree_test(scip, g, &fixed) > 0)
            rerun = TRUE;

      if (bd3_reduction(g))
         rerun = TRUE;
   }
   tt_aggregation(g, &fixed);
   degree_test(scip, g, &fixed);

   SCIPdebugMessage("Reduction Level 5: Fixed Cost = %.12e\n",
      fixed);

   return(fixed);
}
#endif

static
SCIP_RETCODE level5(
   SCIP* scip,
   GRAPH** graph,
   SCIP_Real* fixed
   )
{
   PATH* vnoi;
   SCIP_Real timelimit;
   GRAPH* g = *graph;
   double*  sddist;
   double*  sdtrans;
   double*  sdrand;
   double* cost;
   double* random;
   int*    heap;
   int*    state;
   int*    vbase;
   int*    knotexamined;
   int     i;
   int     nelims;
   int     nnodes;
   int     nedges;
   int     runnum;
   int     sdnelims;
   int     lenelims;
   int     bd3nelims;
   int     nsvnelims;
   int     nvslnelims;
   int     brednelims;
   int     degtnelims;
   int     reductbound;

   char    le = TRUE;
   char    sd = TRUE;
   char    bd3 = TRUE;
   char    nsv = TRUE;
   char    nvsl = TRUE;
   char    bred;
   char    rerun = TRUE;

#ifdef PRINT_TMP_PRESOL
   SCIP_PROBDATA* probdata;
   const char*   probname;
   char   presolvetempfile[SCIP_MAXSTRLEN];
   double offset;
#endif

#ifdef PRINT_TMP_PRESOL
   probdata = SCIPgetProbData(scip);
   probname = SCIPgetProbName(scip);
   (void)SCIPsnprintf(presolvetempfile, SCIP_MAXSTRLEN, "presol/%s-presolve-tmp.stp", probname);
   offset = SCIPprobdataGetOffset(scip);
#endif
   assert(scip != NULL);
   assert(g != NULL);

   nnodes = g->knots;
   nedges = g->edges;

   if( SCIPisLE(scip, (double) g->terms / (double) g->knots, 0.03 ) )
      bred = TRUE;
   else
      bred = FALSE;

#ifdef PRINT_TMP_PRESOL
   SCIPprobdataSetGraph(probdata, g);
   SCIPprobdataSetOffset(probdata, offset + *fixed);
   /* Writing the problem to a temporary file */
   SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
   /* get timelimit parameter*/
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &heap, nnodes + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &knotexamined, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sddist, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdtrans, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdrand, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &random, nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vbase, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vnoi, nnodes) );

   runnum = 0;

   /* define minimal number of edge/node eliminations for a reduction test to be continued */
   reductbound = 0;
   //printf("BOUND: %d \n", reductbound);

   degree_test(scip, g, fixed);

   if( bred )
   {
      SCIP_CALL( bound_reduce(scip, g, vnoi, cost, sddist, random,  heap, state, vbase, &brednelims, *fixed) );
      //printf("bound reduction1: %d \n \n", brednelims);
   }

   while( rerun && !SCIPisStopped(scip) )
   {
      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      sdnelims = 0;
      bd3nelims = 0;
      nsvnelims = 0;
      nvslnelims = 0;
      degtnelims = 0;
      lenelims = 0;

      if( nvsl )
      {
         SCIP_CALL( nvsl_reduction(scip, g, vnoi, fixed, heap, state, vbase, &nvslnelims, reductbound) );

         if( nvslnelims <= 0.3 * reductbound )
            nvsl = FALSE;

         //printf("nvsl: %d \n", nvslnelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( le )
      {
         SCIP_CALL( ledge_reduction(scip, g, vnoi, heap, state, vbase, &lenelims) );

         if( lenelims <= 0.5 * reductbound )
            le = FALSE;

         //printf("lenelims: %d \n", lenelims);
         if( SCIPgetTotalTime(scip) > timelimit )
            break;
      }

      if( sd )
      {
	 nelims = 0;
         for( i = 0; i < 6; i++ )
         {
            SCIP_CALL( sd_reduction(scip, g, sddist, sdtrans, sdrand, cost, random, heap, state, knotexamined, &nelims, runnum) );
            runnum++;
	    sdnelims += nelims;
         }

#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + *fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif

         if( SCIPgetTotalTime(scip) > timelimit )
            break;
         //printf("sdnelims: %d \n", sdnelims);
         if( sdnelims <= reductbound )
            sd = FALSE;
      }

      degtnelims += degree_test(scip, g, fixed);

#ifdef PRINT_TMP_PRESOL
      SCIPprobdataSetGraph(probdata, g);
      SCIPprobdataSetOffset(probdata, offset + *fixed);
      /* Writing the problem to a temporary file */
      SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif

      if( SCIPgetTotalTime(scip) > timelimit )
         break;

      if( nsv )
      {
         SCIP_CALL( nsv_reduction(scip, g, cost, fixed, &nsvnelims) );
         if( nsvnelims <= reductbound )
            nsv = FALSE;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + *fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
      }

      if( bd3 )
      {
         SCIP_CALL( bd3_reduction(scip, g, sddist, sdtrans, heap, state, &bd3nelims) );
         if( bd3nelims <= reductbound )
            bd3 = FALSE;

         if( SCIPgetTotalTime(scip) > timelimit )
            break;

#ifdef PRINT_TMP_PRESOL
         SCIPprobdataSetGraph(probdata, g);
         SCIPprobdataSetOffset(probdata, offset + *fixed);
         /* Writing the problem to a temporary file */
         SCIP_CALL( SCIPwriteOrigProblem(scip, presolvetempfile, NULL, FALSE) );
#endif
      }

      degtnelims += degree_test(scip, g, fixed);

      if( (sdnelims + bd3nelims + nsvnelims + nvslnelims + degtnelims + lenelims) <= reductbound )
         rerun = FALSE;
   }

   /*
     brednelims = le_reduction(g);
     printf("le_reduction: %d \n\n", brednelims);
     brednelims =  tt_deletion(g);
     printf("t_reduction: %d \n\n", brednelims);
     brednelims = tt_aggregation(
     g,
     fixed);
     printf("ttag_reduction: %d \n\n", brednelims);
   */
   SCIPdebugMessage("Reduction Level 4: Fixed Cost = %.12e\n", *fixed);
   /*printf("Total Fixed: %f\n", fixed);*/

   /* free memory */
   SCIPfreeBufferArray(scip, &sddist);
   SCIPfreeBufferArray(scip, &sdtrans);
   SCIPfreeBufferArray(scip, &sdrand);
   SCIPfreeBufferArray(scip, &knotexamined);
   SCIPfreeBufferArray(scip, &heap);
   SCIPfreeBufferArray(scip, &state);
   SCIPfreeBufferArray(scip, &cost);
   SCIPfreeBufferArray(scip, &random);
   SCIPfreeBufferArray(scip, &vbase);
   SCIPfreeBufferArray(scip, &vnoi);

   return SCIP_OKAY;
}

   SCIP_RETCODE reduce(
      SCIP*  scip,
      GRAPH** graph,
      SCIP_Real* offset,
      int    level,
      int    minelims
      )
   {
      assert((*graph)      != NULL);
      assert((*graph)->fixedges == NULL);
      assert(level  >= 0 || level == -4);
      assert(minelims >= 0);

      *offset = 0.0;

      assert((*graph)->layers == 1);

      graph_init_history(scip, (*graph), &((*graph)->orgtail), &((*graph)->orghead), &((*graph)->ancestors));
#if 0
      for( i = 0; i < (*graph)->edges; i++ )
      {
         printf("%d->%d ancestor: %d->%d \n", (*graph)->tail[i], (*graph)->head[i], (*graph)->tail[((*graph)->ancestors[i])->index], (*graph)->head[((*graph)->ancestors[i])->index] );
         assert(((*graph)->ancestors[i])->parent == NULL );
      }
      printf("level: %d \n", level);
#endif

      if( 0 && (*graph)->stp_type != STP_UNDIRECTED )
         return SCIP_OKAY;

      if( (*graph)->stp_type == STP_DEG_CONS || (*graph)->stp_type == STP_GRID || (*graph)->stp_type == STP_OBSTACLES_GRID )
         return SCIP_OKAY;

      /* initialise shortest path algorithms */
      graph_path_init((*graph));

      if( (*graph)->stp_type != STP_UNDIRECTED && (*graph)->stp_type != STP_GRID && (*graph)->stp_type != STP_OBSTACLES_GRID )
         level = level * (-1);

      if( level == 1 )
         *offset = level1(scip, (*graph));
#if 0
      if (level == 2)
         fixed = level1(scip, (*graph)) + level2(scip, (*graph));

      if (level == 3)
         fixed = level3((*graph));

#endif
      if( level == -1 )
         *offset = levelm1(scip, (*graph));

      if( level == 4 )
         SCIP_CALL( level4(scip, (graph), offset, minelims) );

      if( level == 5 )
         SCIP_CALL( level5(scip, (graph), offset) );

      if( level == -4 )
         *offset = levelm4(scip, (*graph));

      graph_path_exit((*graph));

      return SCIP_OKAY;
   }
