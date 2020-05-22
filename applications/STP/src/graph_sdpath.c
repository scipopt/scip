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

/**@file   graph_sdpath.c
 * @brief  Special distance path and walk algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various special distance path (and walk) algorithms for
 * Steiner tree and related problems.
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <assert.h>
#include "portab.h"
#include "graph.h"



inline static SCIP_Real sdwalk_getdistnewEdge(
   const int*            prevedges,          /**< previous edges per node */
   const int*            nprevedges,         /**< number of previous edges per node */
   const SCIP_Real*      cost,               /**< cost */
   const SCIP_Real*      dist,               /**< distance */
   int                   k,                  /**< previous node  */
   int                   e,                  /**< current outgoing edge  */
   int                   maxnprevs           /**< maximum number of previous  */
)
{
   const int nprevs = nprevedges[k];
   SCIP_Real dist_e;

   /* ancestor list not full? */
   if( nprevs != maxnprevs + 1 )
   {
      int i;
      const int e2 = e / 2;
      assert(nprevs <= maxnprevs);

      /* check whether m is contained in ancestor list */
      for( i = 0; i < nprevs; i++ )
      {
         const int prevedge = prevedges[maxnprevs * k + i];

         if( e2 == prevedge )
            break;
      }

      /* e2 in list? */
      if( i != nprevs )
      {
         assert(e2 == prevedges[maxnprevs * k + i]);
         dist_e = dist[k];
      }
      else
         dist_e = dist[k] + cost[e];
   }
   else
   {
      dist_e = dist[k] + cost[e];
   }

   return dist_e;
}


inline static SCIP_Real sdwalk_getdistnewPrize(
   const int*            prevNPterms,        /**< previous np terminals per node */
   const int*            nprevNPterms,       /**< number of previous np terminals per node */
   const int*            termmark,           /**< terminal mark */
   const STP_Bool*       visited,            /**< visited */
   const SCIP_Real*      prize,              /**< prize */
   int                   k,                  /**< current node  */
   int                   m,                  /**< next node  */
   SCIP_Real             distnew,            /**< distance of m */
   int                   maxnprevs           /**< maximum number of previous  */
)
{
   SCIP_Real distnewP = distnew;

   assert(termmark[m] == 1 || termmark[m] == 2 );

   if( termmark[m] == 2 || !visited[m] )
   {
      distnewP = MAX(0.0, distnewP - prize[m]);
   }
   else
   {
      const int nprevs = nprevNPterms[k];

      /* ancestor list not full? */
      if( nprevs != maxnprevs + 1 )
      {
         int i;
         assert(nprevs <= maxnprevs);

         /* check whether m is contained in ancestor list */
         for( i = 0; i < nprevs; i++ )
         {
            const int prevterm = prevNPterms[maxnprevs * k + i];

            if( m == prevterm )
               break;
         }

         /* m not in list? */
         if( i == nprevs )
            distnewP = MAX(0.0, distnewP - prize[m]);
      }
   }

   return distnewP;
}



inline static SCIP_Bool sdwalk_conflict(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   const int*            prevterms,          /**< previous terminals */
   const int*            nprevterms,         /**< number of previous terminals */
   const SCIP_Bool       nodevisited         /**< visited node already? */
   )
{
   const int nprevs = nprevterms[prednode];
   SCIP_Bool conflict = FALSE;

   assert(Is_term(g->term[node]));

   if( !nodevisited )
      return FALSE;

   if( nprevs > maxnprevs )
   {
      assert(nprevs == maxnprevs + 1);
      return TRUE;
   }

   for( int i = 0; i < nprevs; i++ )
   {
      const int prevterm = prevterms[maxnprevs * prednode + i];
      assert(prevterm >= 0);

      if( prevterm == node )
      {
         conflict = TRUE;
         break;
      }
   }

   return conflict;
}

inline static void sdwalk_update(
   const GRAPH*          g,                  /**< graph data structure */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms          /**< number of previous terminals */
   )
{
   const int predsize = nprevterms[prednode];
   const SCIP_Bool isterm = Is_term(g->term[node]);

   assert(predsize <= maxnprevs + 1);

   if( predsize == maxnprevs + 1 || (isterm && predsize == maxnprevs) )
   {
      nprevterms[node] = maxnprevs + 1;
   }
   else
   {
#ifndef NDEBUG
      for( int j = 0; j < predsize; j++ )
         assert(prevterms[maxnprevs * prednode + j] != node);
#endif

      for( int i = 0; i < predsize; i++ )
         prevterms[maxnprevs * node + i] = prevterms[maxnprevs * prednode + i];

      nprevterms[node] = predsize;

      if( isterm )
      {
         assert(predsize < maxnprevs);
         prevterms[maxnprevs * node + predsize] = node;
         nprevterms[node]++;
      }

      assert(nprevterms[node] <= maxnprevs);
   }
}

inline
static void sdwalk_updateCopy(
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   int*                  prev,               /**< previous data elements */
   int*                  nprev               /**< number of previous data elements */
   )
{
   const int predsize = nprev[prednode];

   assert(predsize <= maxnprevs);

   /* copy data from predecesseor */
   for( int i = 0; i < predsize; i++ )
      prev[maxnprevs * node + i] = prev[maxnprevs * prednode + i];

   nprev[node] = predsize;
}

static void sdwalk_update2(
   const int*            termmark,           /**< terminal mark */
   int                   node,               /**< the node to be updated */
   int                   prednode,           /**< the predecessor node */
   int                   edge,               /**< the edge to be updated */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Bool             clear,              /**< clear arrays */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  prevNPterms,        /**< previous non-proper terminals */
   int*                  nprevNPterms,       /**< number of previous non-proper terminals */
   int*                  prevedges,          /**< previous edges */
   int*                  nprevedges          /**< number of previous edges */
   )
{
   int predsize = nprevterms[prednode];

   /*** 1. proper terminals ***/

   /* not enough space? */
   if( predsize == maxnprevs + 1 || (termmark[node] == 2 && predsize == maxnprevs) )
   {
      nprevterms[node] = maxnprevs + 1;
   }
   else
   {
#ifndef NDEBUG
      for( int j = 0; j < predsize; j++ )
         assert(prevterms[maxnprevs * prednode + j] != node);
#endif

      sdwalk_updateCopy(node, prednode, maxnprevs, prevterms, nprevterms);

      if( termmark[node] == 2 )
      {
         assert(predsize < maxnprevs);
         prevterms[maxnprevs * node + predsize] = node;
         nprevterms[node]++;
      }

      assert(nprevterms[node] <= maxnprevs);
   }


   /*** 2. edges ***/

   if( clear )
   {
      nprevNPterms[node] = 0;
      nprevedges[node] = 0;
      return;
   }

   predsize = nprevedges[prednode];

   if( predsize >= maxnprevs )
   {
      assert(predsize == maxnprevs || predsize == maxnprevs + 1);

      nprevedges[node] = maxnprevs + 1;
      nprevNPterms[node] = maxnprevs + 1;
      return;
   }
   assert(predsize < maxnprevs);

   sdwalk_updateCopy(node, prednode, maxnprevs, prevedges, nprevedges);

   prevedges[maxnprevs * node + predsize] = edge / 2;
   nprevedges[node]++;

   assert(nprevedges[node] <= maxnprevs);


   /*** 3. non-proper terminals ***/

   predsize = nprevNPterms[prednode];

   if( predsize == maxnprevs + 1 || (termmark[node] == 1 && predsize == maxnprevs) )
   {
      nprevNPterms[node] = maxnprevs + 1;
   }
   else
   {
      sdwalk_updateCopy(node, prednode, maxnprevs, prevNPterms, nprevNPterms);

      if( termmark[node] == 1 )
      {
         assert(predsize < maxnprevs);
         prevNPterms[maxnprevs * node + predsize] = node;
         nprevNPterms[node]++;
      }

      assert(nprevNPterms[node] <= maxnprevs);
   }
}

inline static void sdwalk_reset(
   int                   nvisits,            /**< number of visited nodes */
   const int*            visitlist,          /**< stores all visited nodes */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   STP_Bool*             visited             /**< stores whether a node has been visited */
)
{
   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      assert(node >= 0);

      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
   }
}



inline static void correctXwalk(
   SCIP* scip,
   int* RESTRICT heap,
   int* RESTRICT state,
   int* RESTRICT count,    /* pointer to store the number of elements on the heap */
   SCIP_Real* RESTRICT pathdist,
   int    l,
   SCIP_Real newcost
   )
{
   int    t;
   int    c;
   int    j;

   pathdist[l] = newcost;

   if (state[l] == UNKNOWN)
   {
      heap[++(*count)] = l;
      state[l]      = (*count);
   }

   /* Heap shift up
    */
   j = state[l];
   c = j / 2;

   while( (j > 1) && pathdist[heap[c]] > pathdist[heap[j]] )
   {
      t              = heap[c];
      heap[c]        = heap[j];
      heap[j]        = t;
      state[heap[j]] = j;
      state[heap[c]] = c;
      j              = c;
      c              = j / 2;
   }
}


/*
 * Interface methods
 */




/** limited Dijkstra, stopping at terminals */
void graph_sdStar(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             with_zero_edges,    /**< telling name */
   int                   star_root,          /**< root of the start */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int*                  star_base,          /**< star base node, must be initially set to SDSTAR_BASE_UNSET */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Bool*            success             /**< will be set to TRUE iff at least one edge can be deleted */
   )
{
   int nchecks;
   int nstarhits;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const int star_degree = range_csr[star_root].end - range_csr[star_root].start;
   SCIP_Real distlimit;
   /* NOTE: with zero edges case is already covered with state[k] = UNKNOWN if k == star_base[k] */
   const SCIP_Real eps = graph_pc_isPcMw(g) ? 0.0 : SCIPepsilon(scip);

   assert(dcsr && g && dist && visitlist && nvisits && visited && dheap && success);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->mark[star_root] && star_degree >= 1);
   assert(dheap->size == 0);
   assert(edgelimit >= 1);

   *nvisits = 0;
   *success = FALSE;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(dist[k] == FARAWAY);
      assert(star_base[k] == SDSTAR_BASE_UNSET);
      assert(state[k] == UNKNOWN);
   }
#endif

   distlimit = 0.0;
   dist[star_root] = 0.0;
   state[star_root] = CONNECT;
   visitlist[(*nvisits)++] = star_root;

   for( int e = range_csr[star_root].start, end = range_csr[star_root].end; e < end; e++ )
   {
      const int m = head_csr[e];

      assert(g->mark[m]);
      assert(!visited[m]);

      visitlist[(*nvisits)++] = m;
      visited[m] = TRUE;
      dist[m] = cost_csr[e];
      star_base[m] = m;

      /*  add epsilon to make sure that m is removed from the heap last in case of equality */
      graph_heap_correct(m, cost_csr[e] + eps, dheap);

      if( cost_csr[e] > distlimit )
         distlimit = cost_csr[e];
   }


   nchecks = 0;
   nstarhits = 0;

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != star_root);
      assert(state[k] == CONNECT);
      assert(LE(dist[k], distlimit));

      if( with_zero_edges && k == star_base[k] )
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         assert(g->mark[m] && star_base[k] >= 0);

         if( state[m] != CONNECT )
         {
            const SCIP_Real distnew = dist[k] + cost_csr[e];

            if( GT(distnew, distlimit) )
               continue;

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               if( star_base[m] == m )
                  nstarhits++;

               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m);
            }
            else if( EQ(distnew, dist[m]) && star_base[m] == m )
            {
               if( with_zero_edges && star_base[k] == star_base[m] )
                  continue;

               assert(visited[m]);
               nstarhits++;

               assert(star_base[m] != star_base[k]);

               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m);
            }

            /* all star nodes hit already? */
            if( nstarhits == star_degree )
            {
               nchecks = edgelimit + 1;
               break;
            }
         }
         nchecks++;
      }
   }

  *success = (nstarhits > 0);
}


/** limited Dijkstra with node bias */
SCIP_RETCODE graph_sdStarBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   star_root,          /**< root of the start */
   int*                  star_base,          /**< star base node, must be initially set to SDSTAR_BASE_UNSET */
   DIJK*                 dijkdata,           /**< Dijkstra data */
   SCIP_Bool*            success             /**< will be set to TRUE iff at least one edge can be deleted */
   )
{
   int nchecks;
   int nstarhits;
   const int nnodes = graph_get_nNodes(g);
   SCIP_Real* restrict dist = dijkdata->node_distance;
   int* restrict visitlist = dijkdata->visitlist;
   STP_Bool* restrict visited = dijkdata->node_visited;
   int* node_preds;
   DHEAP* dheap = dijkdata->dheap;
   const SCIP_Real* const nodebias = dijkdata->node_bias;
   const int* const nodebias_source = dijkdata->node_biassource;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   const RANGE* const RESTRICT range_csr = dcsr->range;
   const int* const RESTRICT head_csr = dcsr->head;
   const SCIP_Real* const RESTRICT cost_csr = dcsr->cost;
   const int star_degree = range_csr[star_root].end - range_csr[star_root].start;
   int nvisits = 0;
   SCIP_Real distlimit;
   const int edgelimit = dijkdata->edgelimit;
   /* NOTE: with zero edges case is already covered with state[k] = UNKNOWN if k == star_base[k] */
   const SCIP_Real eps = graph_pc_isPcMw(g) ? 0.0 : 2.0 * SCIPepsilon(scip);

   assert(dcsr && dist && visitlist && visited && dheap && success);
   assert(nodebias && nodebias_source);
   assert(!g->extended);
   assert(g->mark[star_root] && star_degree >= 1);
   assert(dheap->size == 0);
   assert(edgelimit >= 1);

   nvisits = 0;
   *success = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &node_preds, nnodes) );

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      assert(dist[k] == FARAWAY);
      assert(star_base[k] == SDSTAR_BASE_UNSET);
      assert(state[k] == UNKNOWN);
      node_preds[k] = UNKNOWN;
   }
#endif

   distlimit = 0.0;
   dist[star_root] = 0.0;
   state[star_root] = CONNECT;
   visitlist[(nvisits)++] = star_root;

   for( int e = range_csr[star_root].start, end = range_csr[star_root].end; e < end; e++ )
   {
      const int m = head_csr[e];

      assert(g->mark[m]);
      assert(!visited[m]);

      visitlist[(nvisits)++] = m;
      visited[m] = TRUE;
      dist[m] = cost_csr[e];
      star_base[m] = m;
      node_preds[m] = star_root;

      /*  add epsilon to make sure that m is removed from the heap last in case of equality */
      graph_heap_correct(m, cost_csr[e] + eps, dheap);

      if( cost_csr[e] > distlimit )
         distlimit = cost_csr[e];
   }

   nchecks = 0;
   nstarhits = 0;

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labeled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;
      const int k_pred = node_preds[k];

      assert(k != star_root);
      assert(k_pred >= 0 && k_pred < nnodes);
      assert(state[k] == CONNECT);
      assert(LE(dist[k], distlimit));

      if( k == star_base[k] )
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];
         assert(g->mark[m] && star_base[k] >= 0);

         if( state[m] != CONNECT )
         {
            const int source = nodebias_source[k];
            const SCIP_Bool useBias = (source != m && source != k_pred);
            const SCIP_Real bias = (useBias)? MIN(cost_csr[e], nodebias[k]) : 0.0;
            const SCIP_Real distnew = dist[k] + cost_csr[e] - MIN(dist[k], bias);

            if( GT(distnew, distlimit) )
               continue;

            if( LT(distnew, dist[m]) )
            {
               if( !visited[m] )
               {
                  visitlist[(nvisits)++] = m;
                  visited[m] = TRUE;
               }

               if( star_base[m] == m )
                  nstarhits++;

               node_preds[m] = k;
               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m && m != star_root);
            }
            else if( EQ(distnew, dist[m]) && star_base[m] == m )
            {
               if( star_base[k] == star_base[m] )
                  continue;

               assert(visited[m]);
               nstarhits++;

               assert(star_base[m] != star_base[k]);

               node_preds[m] = k;
               dist[m] = distnew;
               star_base[m] = star_base[k];
               graph_heap_correct(m, distnew, dheap);

               assert(star_base[m] != m && m != star_root);
            }

            /* all star nodes hit already? */
            if( nstarhits == star_degree )
            {
               nchecks = edgelimit + 1;
               break;
            }
         }
         nchecks++;
      }
   }

  dijkdata->nvisits = nvisits;
  *success = (nstarhits > 0);

  SCIPfreeBufferArray(scip, &node_preds);

  return SCIP_OKAY;
}


/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalks(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise)*/
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost[e] - g->prize[m], 0.0);
            correctXwalk(scip, heap, state, &count, dist, m, newcost);
         }
         else
         {
            correctXwalk(scip, heap, state, &count, dist, m, cost[e]);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( (state[m] != CONNECT) && g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}



/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalks_csr(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise)*/
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   RANGE* const RESTRICT range_csr = dcsr->range;
   int* const RESTRICT head_csr = dcsr->head;
   SCIP_Real* const RESTRICT cost_csr = dcsr->cost;

   assert(dcsr && g && dist && visitlist && nvisits && visited && dheap);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->grad[start] != 0 && g->grad[end] != 0);
   assert(g->mark[start] && g->mark[end]);
   assert(dheap->size == 0);

   *nvisits = 0;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
      assert(state[k] == UNKNOWN);
#endif

   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   for( int e = range_csr[start].start; e < range_csr[start].end; e++ )
   {
      const int m = head_csr[e];
      assert(g->mark[m]);

      if( SCIPisLE(scip, cost_csr[e], distlimit) && m != end )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost_csr[e] - g->prize[m], 0.0);

            dist[m] = newcost;
            graph_heap_correct(m, newcost, dheap);
         }
         else
         {
            dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT && m != start )
         {
            SCIP_Real distnew = dist[k] + cost_csr[e];

            assert(g->mark[m]);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
         nchecks++;
      }
   }

   return success;
}


/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksTriangle(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const int*            stateprev,          /**< state of previous run or NULL */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            prizeoffset,        /**< array for storing prize offset or NULL */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;
   int* const state = dheap->position;
   DCSR* const dcsr = g->dcsr_storage;
   RANGE* const RESTRICT range_csr = dcsr->range;
   int* const RESTRICT head_csr = dcsr->head;
   SCIP_Real* const RESTRICT cost_csr = dcsr->cost;

   assert(dcsr && g && dist && visitlist && visited && dheap);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(g->grad[start] != 0 && g->grad[end] != 0);
   assert(g->mark[start] && g->mark[end]);
   assert(dheap->size == 0);

   *nvisits = 0;

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
      assert(state[k] == UNKNOWN);
#endif

   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   for( int e = range_csr[start].start; e < range_csr[start].end; e++ )
   {
      const int m = head_csr[e];
      assert(g->mark[m]);

      if( SCIPisLE(scip, cost_csr[e], distlimit) && m != end )
      {
         assert(!visited[m]);

         if( stateprev && stateprev[m] == CONNECT )
            continue;

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
         {
            const SCIP_Real newcost = MAX(cost_csr[e] - g->prize[m], 0.0);

            dist[m] = newcost;
            graph_heap_correct(m, newcost, dheap);

            if( prizeoffset )
            {
               if( g->prize[m] > cost_csr[e] )
               {
                  prizeoffset[m] = cost_csr[e];
                  assert(SCIPisZero(scip, newcost));
               }
               else
                  prizeoffset[m] = g->prize[m];
            }
         }
         else
         {
            dist[m] = cost_csr[e];
            graph_heap_correct(m, cost_csr[e], dheap);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }

   while( dheap->size > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = graph_heap_deleteMinReturnNode(dheap);
      const int k_start = range_csr[k].start;
      const int k_end = range_csr[k].end;

      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = k_start; e < k_end; e++ )
      {
         const int m = head_csr[e];

         if( state[m] != CONNECT )
         {
            SCIP_Real distnew;

            assert(m != start);

            if( stateprev && stateprev[m] == CONNECT )
                continue;

            distnew = dist[k] + cost_csr[e];

            assert(g->mark[m]);

            if( distnew > distlimit )
               continue;

            if( termmark[m] != 0 )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               if( prizeoffset && termmark[m] != 0 )
               {
                  const SCIP_Real distnew0 = dist[k] + cost_csr[e];

                  if( g->prize[m] > distnew0 )
                  {
                     prizeoffset[m] = distnew0;
                     assert(SCIPisZero(scip, distnew));
                  }
                  else
                     prizeoffset[m] = g->prize[m];
               }

               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               dist[m] = distnew;
               graph_heap_correct(m, distnew, dheap);
            }
         }
         nchecks++;
      }
   }

   return success;
}

/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksExt(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;
         sdwalk_update(g, m, start, maxnprevs, prevterms, nprevterms);

         if( Is_term(g->term[m]) )
         {
            const SCIP_Real newcost = MAX(cost[e] - g->prize[m], 0.0);
            correctXwalk(scip, heap, state, &count, dist, m, newcost);
         }
         else
         {
            correctXwalk(scip, heap, state, &count, dist, m, cost[e]);
         }

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   assert(nprevterms[start] == 0);

   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            assert(state[m] != CONNECT);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( Is_term(g->term[m]) )
               distnew = MAX(distnew - g->prize[m], 0.0);

            if( distnew < dist[m] )
            {
               const SCIP_Bool mvisited = visited[m];
               if( !mvisited )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               if( Is_term(g->term[m]) && sdwalk_conflict(g, m, k, maxnprevs, prevterms, nprevterms, mvisited) )
                  continue;

               sdwalk_update(g, m, k, maxnprevs, prevterms, nprevterms);
               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}



/** modified Dijkstra along walks for PcMw, returns special distance between start and end */
SCIP_Bool graph_sdWalksExt2(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   SCIP_Real             distlimit,          /**< distance limit of the search */
   int                   start,              /**< start */
   int                   end,                /**< end */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   int                   maxnprevs,          /**< maximum number of previous terminals to save */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  prevterms,          /**< previous terminals */
   int*                  nprevterms,         /**< number of previous terminals */
   int*                  prevNPterms,        /**< previous non-proper terminals */
   int*                  nprevNPterms,       /**< number of previous non-proper terminals */
   int*                  prevedges,          /**< previous edges */
   int*                  nprevedges,         /**< number of previous edges */
   int*                  heap,               /**< array representing a heap */
   int*                  state,              /**< array to indicate whether a node has been scanned */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited             /**< stores whether a node has been visited */
   )
{
   int count;
   int nchecks;
   SCIP_Bool success = FALSE;
   const int edgelimit1 = edgelimit / 2;

   assert(g != NULL);
   assert(heap != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);

   *nvisits = 0;

   if( g->grad[start] == 0 || g->grad[end] == 0 )
      return FALSE;

   assert(g->mark[start] && g->mark[end]);

   count = 0;
   nchecks = 0;
   dist[start] = 0.0;
   state[start] = CONNECT;
   visitlist[(*nvisits)++] = start;

   g->mark[start] = FALSE;
   g->mark[end] = FALSE;

   for( int e = g->outbeg[start]; e != EAT_LAST; e = g->oeat[e] )
   {
      const int m = g->head[e];

      if( g->mark[m] && SCIPisLE(scip, cost[e], distlimit) )
      {
         SCIP_Real distnew = cost[e];

         assert(!visited[m]);

         visitlist[(*nvisits)++] = m;
         visited[m] = TRUE;

         if( termmark[m] != 0 )
            distnew = MAX(distnew - g->prize[m], 0.0);

         sdwalk_update2(termmark, m, start, e, maxnprevs, SCIPisZero(scip, distnew),
               prevterms, nprevterms, prevNPterms, nprevNPterms, prevedges, nprevedges);
         correctXwalk(scip, heap, state, &count, dist, m, distnew);

         if( ++nchecks > edgelimit1 )
            break;
      }
   }
   assert(nprevterms[start] == 0);

   g->mark[end] = TRUE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(k != end && k != start);
      assert(SCIPisLE(scip, dist[k], distlimit));

      state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( g->mark[m] )
         {
            SCIP_Real distnew = sdwalk_getdistnewEdge(prevedges, nprevedges, cost, dist, k, e, maxnprevs);

            assert(state[m] != CONNECT);

            if( SCIPisGT(scip, distnew, distlimit) )
               continue;

            if( termmark[m] != 0 )
               distnew = sdwalk_getdistnewPrize(prevNPterms, nprevNPterms, termmark, visited, g->prize, k, m, distnew, maxnprevs);

            if( SCIPisLT(scip, distnew, dist[m]) )
            {
               const SCIP_Bool mvisited = visited[m];
               if( !mvisited )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( m == end )
               {
                  nchecks = edgelimit + 1;
                  success = TRUE;
                  break;
               }

               /* continue if m is proper terminals and is on the walk to k */
               if( termmark[m] == 2 && sdwalk_conflict(g, m, k, maxnprevs, prevterms, nprevterms, mvisited) )
                  continue;

               sdwalk_update2(termmark, m, k, e, maxnprevs, SCIPisZero(scip, distnew),
                     prevterms, nprevterms, prevNPterms, nprevNPterms, prevedges, nprevedges);
               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;
   return success;
}


/** modified Dijkstra along walks for PcMw */
SCIP_Bool graph_sdWalksConnected(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   const int*            termmark,           /**< terminal mark (2 for proper terminal, 1 for non-proper terminal, 0 otherwise) */
   const SCIP_Real*      cost,               /**< edge costs */
   const STP_Bool*       endpoint,           /**< stores whether search should be ended at vertex */
   int                   start,              /**< start vertex */
   int                   edgelimit,          /**< maximum number of edges to consider during execution */
   SCIP_Real*            dist,               /**< distances array, initially set to FARAWAY */
   int*                  visitlist,          /**< stores all visited nodes */
   int*                  nvisits,            /**< number of visited nodes */
   STP_Bool*             visited,            /**< stores whether a node has been visited */
   SCIP_Bool             resetarrays         /**< should arrays be reset? */
   )
{
   int* const heap = g->path_heap;
   int* const state = g->path_state;
   int count;
   int nchecks;
   const SCIP_Real prize = g->prize[start];

   assert(heap != NULL);
   assert(state != NULL);
   assert(dist != NULL);
   assert(cost != NULL);
   assert(visitlist != NULL);
   assert(nvisits != NULL);
   assert(visited != NULL);
   assert(graph_pc_isPcMw(g));
   assert(!g->extended);
   assert(Is_term(g->term[start]));
   assert(g->grad[start] > 0);
   assert(g->mark[start]);

#ifndef NDEBUG
   for( int k = 0; k < g->knots; k++ )
   {
      assert(state[k] == UNKNOWN);
      assert(visited[k] == FALSE);
      assert(dist[k] == FARAWAY);
   }
#endif

   *nvisits = 0;
   nchecks = 0;
   count = 1;
   heap[count] = start;
   state[start] = count;
   dist[start] = 0.0;
   visitlist[(*nvisits)++] = start;
   g->mark[start] = FALSE;

   while( count > 0 && nchecks <= edgelimit )
   {
      /* get nearest labelled node */
      const int k = nearestX(heap, state, &count, dist);
      assert(SCIPisLE(scip, dist[k], prize));

      if( termmark[k] == 2 )
         state[k] = CONNECT;
      else
         state[k] = UNKNOWN;

      /* correct incident nodes */
      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int m = g->head[e];

         if( (state[m] != CONNECT) && g->mark[m] )
         {
            SCIP_Real distnew = dist[k] + cost[e];

            if( SCIPisGT(scip, distnew, prize) )
               continue;

            if( termmark[m] != 0 )
               distnew -= g->prize[m];

            if( distnew < dist[m] )
            {
               if( !visited[m] )
               {
                  visitlist[(*nvisits)++] = m;
                  visited[m] = TRUE;
               }

               /* finished already? */
               if( endpoint != NULL && endpoint[m] )
               {
                  g->mark[start] = TRUE;
                  if( resetarrays )
                     sdwalk_reset(*nvisits, visitlist, dist, state, visited);

                  return TRUE;
               }

               correctXwalk(scip, heap, state, &count, dist, m, distnew);
            }
         }
         nchecks++;
      }
   }

   g->mark[start] = TRUE;

   if( resetarrays )
      sdwalk_reset(*nvisits, visitlist, dist, state, visited);

   return FALSE;
}
