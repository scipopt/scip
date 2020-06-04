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

/**@file   reduce_sd.c
 * @brief  special distance (bottleneck distance) reduction methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various special distance (aka bottleneck distance) based reduction methods
 * for Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 *
 */

//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "misc_stp.h"
#include "probdata_stp.h"
#include "scip/scip.h"
#include "portab.h"


#define STP_BD_MAXDEGREE 4
#define STP_BD_MAXDNEDGES 6
#define STP_SDWALK_MAXNPREVS 8

/** initializes data needed for SD star tests */
static
SCIP_RETCODE sdStarInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   int                   edgelimit,          /**< limit */
   DIJK**                dijkdata,           /**< data to be allocated */
   int*RESTRICT*         star_base,          /**< data to be allocated */
   SCIP_Bool*RESTRICT*   edge_deletable      /**< data to be allocated */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   int* RESTRICT star_base_d;
   SCIP_Bool* RESTRICT edge_deletable_d;

   assert(!graph_pc_isPcMw(g) || !g->extended);

   SCIP_CALL( graph_dijkLimited_init(scip, g, dijkdata) );
   (*dijkdata)->edgelimit = edgelimit;

#ifndef NDEBUG
   {
      SCIP_Real* dist = (*dijkdata)->node_distance;
      STP_Bool* visited = (*dijkdata)->node_visited;
      assert(graph_heap_isClean((*dijkdata)->dheap));

      for( int i = 0; i < nnodes; i++ )
      {
         assert(!visited[i]);
         assert(EQ(dist[i], FARAWAY));
      }
   }
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, star_base, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, edge_deletable, nedges / 2) );

   edge_deletable_d = *edge_deletable;
   for( int e = 0; e < nedges / 2; e++ )
      edge_deletable_d[e] = FALSE;

   star_base_d = *star_base;
   for( int i = 0; i < nnodes; i++ )
      star_base_d[i] = SDSTAR_BASE_UNSET;

   return SCIP_OKAY;
}


/** finalizes SD star tests */
static
void sdStarFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   DIJK**                dijkdata,           /**< data to be freed */
   int*RESTRICT*         star_base,          /**< data to be freed */
   SCIP_Bool*RESTRICT*   edge_deletable      /**< data to be freed */
)
{
#ifndef NDEBUG
   const int nedges = graph_get_nEdges(g);
   DCSR* dcsr = g->dcsr_storage;
   SCIP_Bool* RESTRICT edge_deletable_d = *edge_deletable;

   for( int e = 0; e < nedges / 2; e++ )
   {
      if( edge_deletable_d[e] )
         assert(dcsr->id2csredge[e * 2] == -1);
      else if( g->oeat[e * 2] != EAT_FREE )
         assert(dcsr->id2csredge[e * 2] != -1 || !g->mark[g->tail[e * 2]] || !g->mark[g->head[e * 2]]);
   }
#endif

   graph_edge_delBlocked(scip, g, *edge_deletable, TRUE);
   SCIPfreeBufferArray(scip, edge_deletable);
   SCIPfreeBufferArray(scip, star_base);

   graph_dijkLimited_free(scip, dijkdata);
}


/** resets data needed for SD star tests */
static inline
void sdStarReset(
   int                   nnodes,
   int                   nvisits,
   const int*            visitlist,
   int* RESTRICT         star_base,
   SCIP_Real* RESTRICT   dist,
   STP_Bool* RESTRICT    visited,
   DHEAP* RESTRICT       dheap               /**< Dijkstra heap */
)
{
   int* const state = dheap->position;

   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
      star_base[node] = SDSTAR_BASE_UNSET;
   }

   graph_heap_clean(FALSE, dheap);

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      assert(star_base[k] == SDSTAR_BASE_UNSET);
      assert(visited[k] == FALSE);
      assert(state[k] == UNKNOWN);
      assert(dist[k] == FARAWAY);
   }
#endif
}


/** checks node */
static
SCIP_RETCODE sdStarBiasedProcessNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node,               /**< node to process */
   const int*            edgestate,          /**< state array or NULL */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   GRAPH*                g,                  /**< graph data structure */
   DIJK*                 dijkdata,           /**< data */
   int* RESTRICT         star_base,          /**< data to be freed */
   SCIP_Bool *RESTRICT   edge_deletable,     /**< data to be freed */
   int*                  nelims              /**< point to store number of deleted edges */
)
{
   const int nnodes = graph_get_nNodes(g);
   SCIP_Bool runloop = TRUE;
   DCSR *dcsr;
   RANGE *RESTRICT range_csr;
   int *RESTRICT head_csr;
   int *RESTRICT edgeid_csr;
   const SCIP_Bool checkstate = (edgestate != NULL);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   edgeid_csr = dcsr->edgeid;

   assert(g->mark[node]);

   while( runloop )
   {
      SCIP_Bool success;
      const int start = range_csr[node].start;

      /* not more than one edge? */
      if( range_csr[node].end - start <= 1 )
         break;

      runloop = FALSE;

      /* do the actual star run */
      SCIP_CALL( graph_sdStarBiased(scip, g, sdprofit, node, star_base, dijkdata, &success) );

      if( success )
      {
         int enext;

         /* check all star nodes (neighbors of i) */
         for( int e = start; e < range_csr[node].end; e = enext )
         {
            const int starnode = head_csr[e];
            const int starbase = star_base[starnode];
            assert(star_base[starnode] >= 0);
            assert(SCIPisLE(scip, dijkdata->node_distance[starnode], dcsr->cost[e]));
            assert(star_base[starnode] == starnode || star_base[starnode] >= 0);

            enext = e + 1;

            if( checkstate )
            {
               const int orgedge = edgeid_csr[e];
               if( edgestate[orgedge] == EDGE_BLOCKED )
                  continue;
            }

            /* shorter path to current star node found? */
            if( starnode != starbase )
            {
               assert(star_base[starbase] != SDSTAR_BASE_UNSET);

               /* path still valid? todo really necessary? */
               if( star_base[starbase] != SDSTAR_BASE_KILLED )
               {
                  star_base[starnode] = SDSTAR_BASE_KILLED;
                  edge_deletable[edgeid_csr[e] / 2] = TRUE;
                  graph_dcsr_deleteEdgeBi(scip, dcsr, e);

                  (*nelims)++;
                  enext--;
               }
               else
               {
                  runloop = TRUE;
               }
            }
         } /* traverse star nodes */
      } /* if success */

      sdStarReset(nnodes, dijkdata->nvisits, dijkdata->visitlist, star_base, dijkdata->node_distance, dijkdata->node_visited, dijkdata->dheap);
   }

   return SCIP_OKAY;
}

/** resets data needed for SD walk tests */
static inline
void sdwalkReset(
   int                   nnodes,
   int                   nvisits,
   const int*            visitlist,
   SCIP_Real* RESTRICT   dist,
   int* RESTRICT         state,
   STP_Bool* RESTRICT    visited
)
{
   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
   }

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      assert(visited[k] == FALSE);
      assert(state[k] == UNKNOWN);
      assert(dist[k] == FARAWAY);
   }
#endif
}

static inline
void sdwalkResetExt(
   int                   nnodes,
   int                   nvisits,
   const  int*           visitlist,
   SCIP_Real* RESTRICT   dist,
   int* RESTRICT         nprevterms,
   int* RESTRICT         state,
   STP_Bool* RESTRICT    visited
)
{
   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
      nprevterms[node] = 0;
   }

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      assert(visited[k] == FALSE);
      assert(state[k] == UNKNOWN);
      assert(dist[k] == FARAWAY);
      assert(nprevterms[k] == 0);
   }
#endif
}

static inline
void sdwalkResetExt2(
   int                   nnodes,
   int                   nvisits,
   const  int*           visitlist,
   SCIP_Real*            dist,
   int*                  nprevterms,
   int*                  nprevNPterms,
   int*                  nprevedges,
   int*                  state,
   STP_Bool*             visited
)
{
   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      visited[node] = FALSE;
      dist[node] = FARAWAY;
      state[node] = UNKNOWN;
      nprevterms[node] = 0;
      nprevNPterms[node] = 0;
      nprevedges[node] = 0;
   }

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      assert(visited[k] == FALSE);
      assert(state[k] == UNKNOWN);
      assert(dist[k] == FARAWAY);
      assert(nprevterms[k] == 0);
      assert(nprevNPterms[k] == 0);
      assert(nprevedges[k] == 0);
   }
#endif
}


/* can edge be deleted in SD test in case of equality? If so, 'forbidden' array is adapted */
static
SCIP_Bool sddeltable(
   SCIP*                 scip,
   GRAPH*                g,
   PATH*                 path,
   int*                  vbase,
   int*                  forbidden,
   int                   tpos,
   int                   hpos,
   int                   tail,
   int                   head,
   int                   edge,
   int                   nnodes
   )
{
   int e;

   assert(g != NULL);
   assert(path != NULL);
   assert(scip != NULL);
   assert(vbase != NULL);
   assert(forbidden != NULL);

   assert(tpos >= 0);
   assert(hpos >= 0);
   assert(tail >= 0);
   assert(head >= 0);
   assert(edge >= 0);
   assert(nnodes >= 0);

   /* edge between non-terminals */
   if( !Is_term(g->term[tail]) && !Is_term(g->term[head]) )
      return TRUE;

   /* has edge been marked as forbidden? */
   if( forbidden[edge] )
      return FALSE;

   /* edge between a terminal and a non terminal */
   if( !Is_term(g->term[tail]) || !Is_term(g->term[head]) )
   {

      /* check whether edge is used in shortest path */

      if( !Is_term(g->term[tail]) && Is_term(g->term[head]) )
      {
         e = path[tail + tpos * nnodes].edge;

         if( g->ieat[e] == EAT_FREE )
            return TRUE;

         assert(g->head[e] == tail);

         if( g->tail[e] == head )
            return FALSE;
      }
      else if( Is_term(g->term[tail]) && !Is_term(g->term[head]) )
      {
         e = path[head + hpos * nnodes].edge;

         if( g->ieat[e] == EAT_FREE )
            return TRUE;

         assert(g->head[e] == head);

         if( g->tail[e] == tail )
            return FALSE;
      }
   }

   /* update forbidden edges todo check bottleneck distance between terminals, and don't delete if distance is higher than ecost */

   if( Is_term(g->term[head]) )
   {
      SCIP_Real ecost = g->cost[edge];
      for( e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(e >= 0);

         if( SCIPisEQ(scip, g->cost[e], ecost) )
         {
            if( !(forbidden[e]) && Is_term(g->term[g->head[e]]) )
            {
               forbidden[e] = TRUE;
               forbidden[flipedge(e)] = TRUE;
            }
         }
      }
   }

   if( Is_term(g->term[tail]) )
   {
      SCIP_Real ecost = g->cost[edge];
      for( e = g->outbeg[head]; e != EAT_LAST; e = g->oeat[e] )
      {
         assert(e >= 0);

         if( SCIPisEQ(scip, g->cost[e], ecost) )
         {
            if( !(forbidden[e]) && Is_term(g->term[g->head[e]]) )
            {
               forbidden[e] = TRUE;
               forbidden[flipedge(e)] = TRUE;
            }
         }
      }
   }

   return TRUE;
}


static
int getcloseterms(
   SCIP*                 scip,
   const PATH*           vnoi,
   SCIP_Real*            termdist,
   SCIP_Real             ecost,
   const int*            vbase,
   int*                  neighbterms,
   int                   i,
   int                   nnodes
   )
{
   int k;
   int pos = i;
   int nnterms = 0;

   assert(scip != NULL);
   assert(vnoi != NULL);
   assert(vbase != NULL);
   assert(termdist != NULL);
   assert(neighbterms != NULL);

   for( k = 0; k < 4; k++ )
   {
      if( SCIPisLT(scip, vnoi[pos].dist, ecost) )
      {
         neighbterms[nnterms] = vbase[pos];
         termdist[nnterms++] = vnoi[pos].dist;
      }
      else
      {
         break;
      }
      pos += nnodes;
   }

   return nnterms;
}

static
int getcloseterms2term(
   SCIP*                 scip,
   const GRAPH*          g,
   const GRAPH*          netgraph,
   SCIP_Real*            termdist,
   SCIP_Real             ecost,
   int*                  neighbterms,
   int*                  nodesid,
   int*                  nodesorg,
   int                   i
)
{
   int nnterms = 0;

   for( int k = 0; k < 4; k++ )
      neighbterms[k] = UNKNOWN;

   /* get three nearest terms */
   for( int ne = netgraph->outbeg[nodesid[i]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
   {
      if( SCIPisLT(scip, netgraph->cost[ne], ecost) )
      {
         const SCIP_Real necost = netgraph->cost[ne];
         int j = nodesorg[netgraph->head[ne]];

         assert(Is_term(g->term[j]));
         assert(j != i);

         if( nnterms < 4 )
            nnterms++;
         for( int k = 0; k < 4; k++ )
         {
            if( neighbterms[k] == UNKNOWN || SCIPisGT(scip, termdist[k], necost) )
            {
               for( int l = 3; l > k; l-- )
               {
                  neighbterms[l] = neighbterms[l - 1];
                  termdist[l] = termdist[l - 1];
               }
               neighbterms[k] = j;
               termdist[k] = necost;
               break;
            }
         }
      }
   }

   return nnterms;
}

static
int getlecloseterms(
   SCIP*                 scip,
   PATH*                 vnoi,
   SCIP_Real*            termdist,
   SCIP_Real             ecost,
   int*                  vbase,
   int*                  neighbterms,
   int                   i,
   int                   nnodes
   )
{
   int k;
   int pos = i;
   int nnterms = 0;

   assert(scip != NULL);
   assert(vnoi != NULL);
   assert(vbase != NULL);
   assert(termdist != NULL);
   assert(neighbterms != NULL);

   for( k = 0; k < 4; k++ )
   {
      if( SCIPisLE(scip, vnoi[pos].dist, ecost) )
      {
         neighbterms[nnterms] = vbase[pos];
         termdist[nnterms++] = vnoi[pos].dist;
      }
      else
      {
         break;
      }
      pos += nnodes;
   }

   return nnterms;
}


/** bias DCSR costs for PCMW */
static
void pcBiasCostsDCSR(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Bool             biasreversed,       /**< bias reversed */
   SCIP_Real*            costbiased,         /**< biased edge cost */
   SCIP_Real*            mincost_in          /**< vertex distances */
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   int* head_csr;
   SCIP_Real* cost_csr;
   const int nnodes = g->knots;

   assert(g && scip && g->dcsr_storage);
   assert(graph_pc_isPcMw(g) && !g->extended);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   cost_csr = dcsr->cost;

   assert(g->edges >= dcsr->nedges);

   BMScopyMemoryArray(costbiased, cost_csr, dcsr->nedges);

   /* compute minimum incident edge cost per vertex */

   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && g->mark[k] )
         mincost_in[k] = g->prize[k];
      else
         mincost_in[k] = 0.0;
   }

   for( int k = 0; k < nnodes; k++ )
   {
      if( g->mark[k] )
      {
         const int end = range_csr[k].end;

         for( int e = range_csr[k].start; e < end; e++ )
         {
            const int head = head_csr[e];

            if( Is_term(g->term[head]) )
            {
               assert(g->mark[head]);

               assert(LT(costbiased[e], FARAWAY) && costbiased[e] > 0.0);

               if( costbiased[e] < mincost_in[head] )
                  mincost_in[head] = costbiased[e];
            }
         }
      }
   }

   /* change edge costs */

   if( biasreversed )
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( g->mark[k] && Is_term(g->term[k]) )
         {
            const int end = range_csr[k].end;

            for( int e = range_csr[k].start; e < end; e++ )
            {
               assert(g->mark[head_csr[e]]);
               costbiased[e] -= mincost_in[k];

               assert(!SCIPisNegative(scip, costbiased[e]));

            }
         }
      }
   }
   else
   {
      for( int k = 0; k < nnodes; k++ )
      {
         if( g->mark[k] )
         {
            const int end = range_csr[k].end;

            for( int e = range_csr[k].start; e < end; e++ )
            {
               const int head = head_csr[e];

               if( Is_term(g->term[head]) )
               {
                  assert(g->mark[head]);
                  costbiased[e] -= mincost_in[head];

                  assert(!SCIPisNegative(scip, costbiased[e]));
               }
            }
         }
      }
   }
}


/** get special distance to a single edge */
static
SCIP_Real getSd(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   GRAPH*                netgraph,           /**< auxiliary graph structure */
   PATH*                 mst,                /**< MST structure */
   PATH*                 vnoi,               /**< path structure */
   SCIP_Real*            mstsdist,           /**< MST distance in aux-graph */
   SCIP_Real             sd_initial,         /**< initial sd or -1.0 */
   int*                  vbase,              /**< bases for nearest terminals */
   int*                  nodesid,            /**< nodes identification array */
   int                   i,                  /**< first vertex */
   int                   i2,                 /**< second vertex */
   int                   limit               /**< limit for incident edges to consider */
   )
{
   SCIP_Real sd;
   SCIP_Real max;
   SCIP_Real dist;
   int l;
   int e;
   int j;
   int k;
   int ne;
   int tj;
   int tk;
   int nj;
   int nk;
   int nnodes;
   int nnterms1;
   int nnterms2;
   SCIP_Real termdist1[4];
   SCIP_Real termdist2[4];
   int neighbterms1[4];
   int neighbterms2[4];

   assert(scip != NULL);
   assert(g != NULL);
   assert(netgraph != NULL);
   assert(mst != NULL);
   assert(mstsdist != NULL);

   nnodes = g->knots;
   l = 0;

   if( sd_initial >= 0.0 )
      sd = sd_initial;
   else
      sd = FARAWAY;

   /* compare restricted sd with edge cost (if existing) */
   for( e = g->outbeg[i]; (l++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
         if( g->cost[e] < sd )
            sd = g->cost[e];
         break;
      }
   }

   /* is i a terminal? If not, get three closest terminals of distance smaller sd */
   if( Is_term(g->term[i]) )
   {
      nnterms1 = 1;
      termdist1[0] = 0.0;
      neighbterms1[0] = i;
   }
   else
   {
      nnterms1 = getcloseterms(scip, vnoi, termdist1, sd, vbase, neighbterms1, i, nnodes);

      if( nnterms1 == 0 )
         return sd;
   }

   /* is i2 a terminal? If not, get three closest terminals of distance smaller sd */
   if( Is_term(g->term[i2]) )
   {
      nnterms2 = 1;
      termdist2[0] = 0.0;
      neighbterms2[0] = i2;
   }
   else
   {
      /* get closest terminals of distance smaller sd */
      nnterms2 = getcloseterms(scip, vnoi, termdist2, sd, vbase, neighbterms2, i2, nnodes);

      if( nnterms2 == 0 )
         return sd;
   }

   for( j = 0; j < nnterms1; j++ )
   {
      tj = neighbterms1[j];
      assert(tj >= 0);
      for( k = 0; k < nnterms2; k++ )
      {
         tk = neighbterms2[k];
         assert(Is_term(g->term[tk]));
         assert(Is_term(g->term[tj]));
         assert(tk >= 0);

         if( SCIPisGT(scip, termdist1[j], termdist2[k]) )
            max = termdist1[j];
         else
            max = termdist2[k];

         if( tj == tk )
         {
            if( SCIPisLT(scip, max, sd) )
               sd = max;
         }
         else
         {
            /* get sd between (terminals) tj and tk */
            nj = nodesid[tj];
            nk = nodesid[tk];
            assert(nj != nk);

            l = nj;
            dist = 0.0;
            mstsdist[l] = 0.0;

            while( l != 0 )
            {
               ne = mst[l].edge;

               assert(netgraph->head[ne] == l);
               l = netgraph->tail[ne];
               if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                  dist = netgraph->cost[ne];

               mstsdist[l] = dist;
               if( l == nk )
                  break;
            }

            if( l == nk )
            {
               l = 0;
            }
            else
            {
               l = nk;
               dist = 0.0;
            }
            while( l != 0 )
            {
               ne = mst[l].edge;
               l = netgraph->tail[ne];
               if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                  dist = netgraph->cost[ne];

               if( mstsdist[l] >= 0 )
               {
                  if( SCIPisGT(scip, mstsdist[l], dist) )
                     dist = mstsdist[l];
                  break;
               }
               if( l == 0 )
                  assert( nj == 0);
            }

            /* restore */
            l = nj;
            mstsdist[l] = -1.0;
            while( l != 0 )
            {
               ne = mst[l].edge;
               l = netgraph->tail[ne];
               mstsdist[l] = -1.0;
               if( l == nk )
                  break;
            }

            assert(SCIPisGT(scip, dist, 0.0));
            if( SCIPisGT(scip, dist, max) )
               max = dist;
            if( SCIPisLT(scip, max, sd) )
               sd = max;
         }
      } /* k < nnterms2 */
   } /* j < nnterms1 */
   return sd;
}


/** gets special distance to a pair of nodes */
static
SCIP_Real sdGetSd(
   const GRAPH*          g,                  /**< graph structure */
   int                   i,                  /**< first vertex */
   int                   i2,                 /**< second vertex */
   SCIP_Real             sd_upper,           /**< upper bound on special distance that is accepted (can be FARAWAY) */
   SCIP_Real             sd_sufficient,      /**< bound below which to terminate (can be 0.0) */
   SD*                   sddata              /**< SD */
   )
{
   SCIP_Real termdist1[4];
   SCIP_Real termdist2[4];
   int neighbterms1[4];
   int neighbterms2[4];
   SDGRAPH* sdgraph;
   SCIP_Real sd = sd_upper;
   int nnterms1;
   int nnterms2;
   const SCIP_Bool terminateEarly = GT(sd_sufficient, 0.0);

   assert(g && sddata);
   assert(i != i2);

   sdgraph = sddata->sdgraph;
   assert(sdgraph);

   /* get closest terminals of distance strictly smaller than 'sd' */
   reduce_tpathsGet4CloseTerms(g, sddata->terminalpaths, i, sd, neighbterms1, termdist1, &nnterms1);
   if( nnterms1 == 0 )
      return sd;

   reduce_tpathsGet4CloseTerms(g, sddata->terminalpaths, i2, sd, neighbterms2, termdist2, &nnterms2);
   if( nnterms2 == 0 )
      return sd;

   for( int j = 0; j < nnterms1; j++ )
   {
      const int tj = neighbterms1[j];
      assert(tj >= 0);

      for( int k = 0; k < nnterms2; k++ )
      {
         SCIP_Real sd_jk = MAX(termdist1[j], termdist2[k]);
         const int tk = neighbterms2[k];

         assert(Is_term(g->term[tk]));
         assert(Is_term(g->term[tj]));
         assert(tk >= 0);

         if( tj != tk)
         {
            /* get terminal-to-terminal special distance */
            const SCIP_Real sdt_jk = reduce_sdgraphGetSd(tj, tk, sdgraph);

            if( sdt_jk > sd_jk )
               sd_jk = sdt_jk;
         }

         if( sd_jk < sd )
            sd = sd_jk;

         if( terminateEarly && LT(sd, sd_sufficient) )
         {
            return sd;
         }
      } /* k < nnterms2 */
   } /* j < nnterms1 */

   return sd;
}


/** is node of degree 3 pseudo-deletable? */
static inline
SCIP_Bool isPseudoDeletableDeg3(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      sdsK3,              /**< (unordered) special distances of K3 */
   const int*            edgesK3,            /**< (unordered) edges of K3 */
   SCIP_Real             costK3,             /**< total edge cost */
   SCIP_Bool             allowEquality       /**< allow equality? */
)
{
   assert(scip && g && sdsK3 && edgesK3);
   assert(costK3 >= 0.0);

   if( allowEquality )
   {
      if( SCIPisLE(scip, sdsK3[0] + sdsK3[1], costK3) )
         return TRUE;

      if( SCIPisLE(scip, sdsK3[0] + sdsK3[2], costK3) )
         return TRUE;

      if( SCIPisLE(scip, sdsK3[1] + sdsK3[2], costK3) )
         return TRUE;
   }
   else
   {
      if( SCIPisLT(scip, sdsK3[0] + sdsK3[1], costK3) )
         return TRUE;

      if( SCIPisLT(scip, sdsK3[0] + sdsK3[2], costK3) )
         return TRUE;

      if( SCIPisLT(scip, sdsK3[1] + sdsK3[2], costK3) )
         return TRUE;
   }

   if( graph_pseudoAncestors_edgesInConflict(scip, g, edgesK3, 3) )
      return TRUE;

   return FALSE;
}


/** is given graph not part of at least one optimal solution? */
static
SCIP_Bool isPseudoDeletable(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const GRAPH*          auxg,               /**< auxiliary graph structure */
   const SCIP_Real*      ecost,              /**< adjacent edges costs */
   const int*            edges,              /**< adjacent edges of node */
   int                   degree              /**< degree of the node to be checked */
)
{
   SCIP_Bool success = TRUE;
   const SCIP_Real costsum = ecost[0] + ecost[1] + ecost[2] + ecost[3];
   SCIP_Real mstcost = 0.0;
   PATH mst[STP_BD_MAXDEGREE];

   for( int l = 0; l < STP_BD_MAXDEGREE; l++ )
      mst[l].dist = UNKNOWN;

   assert(graph_isMarked(auxg));
   assert(degree == 4);

   /* compute MST on all neighbors */
   graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);

#ifndef NDEBUG
   for( int l = 1; l < STP_BD_MAXDEGREE; l++ )
      assert(mst[l].dist != UNKNOWN);
#endif

   for( int l = 1; l < STP_BD_MAXDEGREE; l++ )
      mstcost += mst[l].dist;

   if( SCIPisGT(scip, mstcost, costsum) )
   {
      success = FALSE;

      if( graph_pseudoAncestors_edgesInConflict(scip, g, edges, degree) )
         success = TRUE;
   }

   /* check all 3-subsets of neighbors */
   for( int k = 0; k < 4 && success; k++ )
   {
      const int auxvertex1 = (k + 1) % 4;
      const int auxvertex2 = (k + 2) % 4;
      const int auxvertex3 = (k + 3) % 4;
      const int edgesK3[] = { edges[auxvertex1], edges[auxvertex2], edges[auxvertex3] };
      SCIP_Real sdK3[3];
      int sdcount = 0;

      assert(auxvertex1 != k && auxvertex2 != k && auxvertex3 != k);

      for( int e = auxg->outbeg[auxvertex1]; e != EAT_LAST; e = auxg->oeat[e] )
      {
         if( auxg->head[e] != k )
         {
            assert(sdcount < 2);
            sdK3[sdcount++] = auxg->cost[e];
         }
      }

      assert(sdcount == 2);

      for( int e = auxg->outbeg[auxvertex2]; e != EAT_LAST; e = auxg->oeat[e] )
      {
         if( auxg->head[e] == auxvertex3 )
         {
            sdK3[sdcount++] = auxg->cost[e];
            break;
         }
      }

      assert(sdcount == 3);
      success = isPseudoDeletableDeg3(scip, g, sdK3, edgesK3, costsum - ecost[k], TRUE);
   }

   return success;
}


/** initializes data */
static inline
SCIP_RETCODE sdCliqueInitData(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   int                   centernode,         /**< center node or - 1 */
   const GRAPH*          cliquegraph,        /**< clique graph */
   const int*            cliqueNodeMap,      /**< maps clique graph vertices to original ones */
   DIJK*                 dijkdata,           /**< data for repeated path computations */
   SDCLIQUE*             sdclique            /**< clique */
)
{
   SCIP_Real* sds_buffer;
   int* cliquenodes;
   const int nnodes_cliquegraph = graph_get_nNodes(cliquegraph);
   const int* const nodemark = cliquegraph->mark;
   int nnodes_clique = 0;

   /* NOTE: might be too big, but slightly better for cache reuse */
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquenodes, nnodes_cliquegraph) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sds_buffer, cliquegraph->edges / 2) );

   for( int i = 0; i < nnodes_cliquegraph; ++i )
   {
      if( nodemark[i] )
         cliquenodes[nnodes_clique++] = cliqueNodeMap[i];
   }

#ifndef NDEBUG
   for( int i = nnodes_clique; i < nnodes_cliquegraph; ++i )
      cliquenodes[i] = -1;
#endif

   sdclique->centernode = centernode;
   sdclique->cliqueToNodeMap = cliqueNodeMap;
   sdclique->ncliquenodes = nnodes_clique;
   sdclique->sds = sds_buffer;
   sdclique->cliquenodes = cliquenodes;
   sdclique->dijkdata = dijkdata;

   return SCIP_OKAY;
}


/** frees data */
static inline
void sdCliqueFreeData(
   SCIP*                 scip,               /**< SCIP */
   SDCLIQUE*             sdclique            /**< clique */
)
{
   SCIPfreeBufferArray(scip, &(sdclique->sds));
   SCIPfreeBufferArray(scip, &(sdclique->cliquenodes));
}


/** tries to improve SDs of clique-graph
 *  by using the star SD clique algorithm */
static inline
SCIP_RETCODE sdCliqueUpdateGraphWithStarWalks(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   GRAPH*                cliquegraph,        /**< clique graph */
   SDCLIQUE*             sdclique            /**< clique */
)
{
   const SCIP_Real* sds_buffer;
   const int* const nodemark = cliquegraph->mark;
   const int nnodes_cliquegraph = graph_get_nNodes(cliquegraph);
   int nsds = 0;

   SCIP_CALL( graph_sdComputeCliqueStar(scip, g, sdprofit, sdclique) );
   sds_buffer = sdclique->sds;

   for( int k1 = 0; k1 < nnodes_cliquegraph; k1++ )
   {
      if( !nodemark[k1] )
         continue;

      for( int e = cliquegraph->outbeg[k1]; e != EAT_LAST; e = cliquegraph->oeat[e] )
      {
         const int k2 = cliquegraph->head[e];

         if( !nodemark[k2] )
            continue;

         if( k2 > k1 )
         {
#ifndef NDEBUG
            const int* cliqueNodeMap = sdclique->cliqueToNodeMap;
            const int v1 = cliqueNodeMap[k1];
            const int v2 = cliqueNodeMap[k2];
            assert(0 <= v1 && v1 < g->knots);

            assert(0 <= v2 && v2 < g->knots);
            assert(v1 != v2);
#endif

            assert(GT(sds_buffer[nsds], 0.0));
            assert(LE(sds_buffer[nsds], cliquegraph->cost[e]));
            assert(EQ(cliquegraph->cost[e], cliquegraph->cost[flipedge(e)]));

            if( sds_buffer[nsds] < cliquegraph->cost[e] )
            {
               cliquegraph->cost[e] = sds_buffer[nsds];
               cliquegraph->cost[flipedge(e)] = sds_buffer[nsds];
            }

            nsds++;
            assert(nsds <= cliquegraph->edges / 2);
         }
      }
   }

   return SCIP_OKAY;
}


/** get SDs between all pairs of marked vertices of given clique graph by
 *  using terminal-to-terminal special distances */
static
void sdGetSdsCliqueTermWalks(
   const GRAPH*          g,                  /**< the graph */
   SD* RESTRICT          sddata,             /**< SD */
   GRAPH* RESTRICT       cliquegraph,        /**< clique graph to be filled */
   SDCLIQUE* RESTRICT    sdclique            /**< clique */
)
{
   const int nnodes_clique = graph_get_nNodes(cliquegraph);
   const SCIP_Real maxtreecost = reduce_sdgraphGetMaxCost(sddata->sdgraph);
   const int* const nodemark = cliquegraph->mark;
   const int* cliqueNodeMap = sdclique->cliqueToNodeMap;
   SCIP_Real* RESTRICT sds_buffer = sdclique->sds;
   int nsds = 0;

   assert(cliqueNodeMap && sddata);
   assert(2 <= nnodes_clique);

   for( int k1 = 0; k1 < nnodes_clique; k1++ )
   {
      int v1;

      if( !nodemark[k1] )
         continue;

      v1 = cliqueNodeMap[k1];
      assert(0 <= v1 && v1 < g->knots);

      for( int e = cliquegraph->outbeg[k1]; e != EAT_LAST; e = cliquegraph->oeat[e] )
      {
         const int k2 = cliquegraph->head[e];

         if( !nodemark[k2] )
            continue;

         if( k2 > k1 )
         {
            const int v2 = cliqueNodeMap[k2];
            assert(0 <= v2 && v2 < g->knots);
            assert(v1 != v2);

            sds_buffer[nsds] = sdGetSd(g, v1, v2, maxtreecost, 0.0, sddata);
            cliquegraph->cost[e] = sds_buffer[nsds];
            cliquegraph->cost[flipedge(e)] = sds_buffer[nsds];

            nsds++;
            assert(nsds <= cliquegraph->edges / 2);
         }
      }
   }

#ifndef NDEBUG
   for( int k = nsds; k < cliquegraph->edges / 2; k++ )
      sds_buffer[k] = FARAWAY;
#endif
}



/*
 * Interface methods
 */


/** initializes SD structure */
SCIP_RETCODE reduce_sdInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const :(
                                                  terrible design */
   SD**                  sd                  /**< to initialize */
)
{
   SD* s;
   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, sd) );
   s = *sd;

   s->isBiased = FALSE;
   s->sdprofit = NULL;
   s->sdneighbors = NULL;
   SCIP_CALL( reduce_tpathsInit(scip, g, &(s->terminalpaths)) );
   SCIP_CALL( reduce_sdgraphInit(scip, g, &(s->sdgraph)) );
   reduce_sdgraphInitOrderedMstCosts(s->sdgraph);

   return SCIP_OKAY;
}


/** initializes biased SD structure */
SCIP_RETCODE reduce_sdInitBiased(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const :(
                                                  terrible design */
   SD**                  sd                  /**< to initialize */
)
{
   SD* s;
   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, sd) );
   s = *sd;

   s->isBiased = TRUE;
   s->sdneighbors = NULL;
   SCIP_CALL( reduce_sdprofitInit(scip, g, &(s->sdprofit)) );
   SCIP_CALL( reduce_tpathsInit(scip, g, &(s->terminalpaths)) );
   SCIP_CALL( reduce_sdgraphInitBiased(scip, g, s->sdprofit, &(s->sdgraph)) );
   reduce_sdgraphInitOrderedMstCosts(s->sdgraph);

   return SCIP_OKAY;
}


/** initializes biased neighbor SD structure */
SCIP_RETCODE reduce_sdInitBiasedNeighbor(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const :(
                                                  terrible design */
   SD**                  sd                  /**< to initialize */
)
{
   SD* s;
   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, sd) );
   s = *sd;

#if 1
   s->isBiased = TRUE;
   SCIP_CALL( reduce_sdprofitInit(scip, g, &(s->sdprofit)) );
   SCIP_CALL( reduce_tpathsInit(scip, g, &(s->terminalpaths)) );
   SCIP_CALL( reduce_sdgraphInitBiased(scip, g, s->sdprofit, &(s->sdgraph)) );
   reduce_sdgraphInitOrderedMstCosts(s->sdgraph);
#else
   s->isBiased = FALSE;
     s->sdprofit = NULL;
     SCIP_CALL( reduce_tpathsInit(scip, g, &(s->terminalpaths)) );
     SCIP_CALL( reduce_sdgraphInit(scip, g, &(s->sdgraph)) );
     reduce_sdgraphInitOrderedMstCosts(s->sdgraph);
#endif

   reduce_sdneighborUpdate(scip, g, s);

   return SCIP_OKAY;
}


/** get SDs between all pairs of marked vertices of given clique graph */
SCIP_RETCODE reduce_sdGetSdsCliquegraph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< the graph */
   int                   centernode,         /**< center node or - 1 */
   const int*            cliqueNodeMap,      /**< maps clique graph vertices to original ones */
   DIJK*                 dijkdata,           /**< data for repeated path computations */
   SD*                   sddata,             /**< SD */
   GRAPH*                cliquegraph         /**< clique graph to be filled */
)
{
   SDCLIQUE sdclique;

   assert(scip && g && cliqueNodeMap && dijkdata && sddata && cliquegraph);
   assert(cliquegraph->edges == (cliquegraph->knots) * (cliquegraph->knots - 1));

   SCIP_CALL( sdCliqueInitData(scip, g, centernode, cliquegraph, cliqueNodeMap, dijkdata, &sdclique) );

   sdGetSdsCliqueTermWalks(g, sddata, cliquegraph, &sdclique);
   SCIP_CALL( sdCliqueUpdateGraphWithStarWalks(scip, g, sddata->sdprofit, cliquegraph, &sdclique) );

   sdCliqueFreeData(scip, &sdclique);

   return SCIP_OKAY;
}


/** frees SD structure */
void reduce_sdFree(
   SCIP*                 scip,               /**< SCIP */
   SD**                  sd                  /**< to free */
)
{
   SD* s;
   assert(scip && sd);

   s = *sd;
   assert(s);

   reduce_sdgraphFree(scip, &(s->sdgraph));
   reduce_tpathsFree(scip, &(s->terminalpaths));

   if( s->sdprofit )
      reduce_sdprofitFree(scip, &(s->sdprofit));

   if( s->sdneighbors )
      reduce_sdneighborFree(scip, &(s->sdneighbors));

   SCIPfreeMemory(scip, sd);
}


/** long edge implied special distance test */
SCIP_RETCODE reduce_sdImpLongEdge(
   SCIP*                 scip,               /**< SCIP data structure */
   const int*            edgestate,          /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   GRAPH*                g,                  /**< graph data structure */
   SDGRAPH*              sdgraph,            /**< SD graph data structure */
   int*                  nelims              /**< point to store number of deleted edges */
)
{
   const STP_Bool* const edges_isBlocked = reduce_sdgraphGetMstHalfMark(sdgraph);
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool checkstate = (edgestate != NULL);
   const SCIP_Real maxcost = reduce_sdgraphGetMaxCost(sdgraph);

   assert(scip && nelims);

   for( int k = 0; k < nnodes; k++ )
   {
      int e = g->outbeg[k];
      while( e != EAT_LAST )
      {
         assert(e >= 0);
         if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
         {
            e = g->oeat[e];
            continue;
         }

         if( SCIPisGE(scip, g->cost[e], maxcost) && !edges_isBlocked[e / 2] )
         {
            const int enext = g->oeat[e];
            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;

            e = enext;
         }
         else
         {
            e = g->oeat[e];
         }
      }
   }

   /* graph might have become disconnected */
   if( *nelims > 0 )
   {
      SCIP_CALL( reduceLevel0(scip, g) );
   }

   return SCIP_OKAY;
}


/** Special distance test */
SCIP_RETCODE reduce_sd(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   SCIP_Real*            edgepreds,          /**< array to store edge predecessors of auxiliary graph */
   SCIP_Real*            mstsdist,           /**< array to store mst distances in auxiliary graph */
   int*                  state,              /**< array to indicate whether a node has been scanned during SP calculation */
   int*                  vbase,              /**< Voronoi base to each vertex */
   int*                  nodesid,            /**< array to map nodes in auxiliary graph to original ones */
   int*                  nodesorg,           /**< array to map terminals of original graph vertices of auxiliary graph */
   int*                  forbidden,          /**< array to mark whether an edge may be eliminated */
   int*                  nelims,             /**< point to store number of deleted edges */
   SCIP_Bool             nodereplacing,      /**< should node replacement (by edges) be performed? */
   int*                  edgestate           /**< array to store status of (directed) edge (for propagation, can otherwise be set to NULL) */
   )
{
   GRAPH* netgraph;
   PATH* mst;
   SCIP_Real termdist1[4];
   SCIP_Real termdist2[4];
   SCIP_Real ecost;
   SCIP_Real dist;
   int neighbterms1[4];
   int neighbterms2[4];
   int e;
   int i;
   int j;
   int k;
   int l;
   int i2;
   int tj;
   int tk;
   int ne;
   int nj;
   int nk;
   int id1;
   int id2;
   int enext;
   int nnodes;
   int nterms;
   int nedges;
   int nnterms1;
   int nnterms2;
   int maxnedges;
   SCIP_Bool checkstate = (edgestate != NULL);

   assert(g != NULL);
   assert(scip != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(nelims != NULL);
   assert(nodesid != NULL);
   assert(mstsdist != NULL);
   assert(nodesorg != NULL);
   assert(edgepreds != NULL);
   assert(forbidden != NULL);

   nnodes = g->knots;
   nterms = g->terms;
   nedges = g->edges;
   *nelims = 0;
   maxnedges = MIN(nedges, (nterms - 1) * nterms);

   /* only one terminal left? */
   if( nterms == 1 )
      return SCIP_OKAY;

   /* compute nearest four terminals to all non-terminals */
   graph_get4nextTermsPaths(g, g->cost, g->cost, vnoi, vbase, state);

   /* construct auxiliary graph to compute paths between terminals */

   /* initialize the new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1) );

   j = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && g->grad[k] > 0 )
      {
         graph_knot_add(netgraph, -1);
         nodesid[k] = j;
         mstsdist[j] = -1.0;
         netgraph->mark[j] = TRUE;
         nodesorg[j++] = k;
      }
      else
      {
         nodesid[k] = UNKNOWN;
      }
   }

   for( k = 0; k < nedges; k++ )
   {
      forbidden[k] = FALSE;
      edgepreds[k] = -1;
   }

   assert(netgraph->knots == j);
   assert(netgraph->knots == nterms);

   for( k = 0; k < nnodes; k++ )
   {
      if( g->grad[k] == 0 )
         continue;

      i = vbase[k];
      id1 = nodesid[i];

      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         if( i != vbase[g->head[e]] )
         {
            i2 = vbase[g->head[e]];
            id2 = nodesid[i2];

            assert(id1 >= 0);
            assert(id2 >= 0);
            assert(Is_term(g->term[i]));
            assert(Is_term(g->term[i2]));

            for( ne = netgraph->outbeg[id1]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
               if( netgraph->head[ne] == id2 )
                  break;

            /* cost of the edge in the auxiliary graph */
            ecost = g->cost[e] + vnoi[g->head[e]].dist + vnoi[g->tail[e]].dist;

            /* does edge already exist? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(netgraph->tail[ne] == id1);
               assert(netgraph->head[ne] == id2);

               /* is the new edge better than the existing one? */
               if( SCIPisGT(scip, netgraph->cost[ne], ecost) )
               {
                  netgraph->cost[ne]            = ecost;
                  netgraph->cost[flipedge(ne)]  = ecost;

                  edgepreds[ne] = e;
                  edgepreds[flipedge(ne)] = flipedge(e);

                  assert(ne <= maxnedges);
               }
            }
            else
            {
               edgepreds[netgraph->edges] = e;
               edgepreds[netgraph->edges + 1] = flipedge(e);

               graph_edge_add(scip, netgraph, id1, id2, ecost, ecost);

               assert(netgraph->edges <= maxnedges);
            }
         }
      }
   }

   /* compute MST on netgraph */
   graph_knot_chg(netgraph, 0, 0);
   netgraph->source = 0;

   SCIP_CALL( graph_path_init(scip, netgraph) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, nterms) );

   graph_path_exec(scip, netgraph, MST_MODE, 0, netgraph->cost, mst);

   /* mark (original) edges of MST */
   for( k = 1; k < netgraph->knots; k++ )
   {
      assert(mst[k].edge != -1);
      assert((int) edgepreds[mst[k].edge] != -1);
      assert((int) edgepreds[flipedge(mst[k].edge)] != -1);

      e = (int) edgepreds[mst[k].edge];

      assert(vbase[g->tail[e]] == nodesorg[k] || vbase[g->head[e]] == nodesorg[k]);

      if( Is_term(g->tail[e]) && Is_term(g->head[e]) )
      {
         forbidden[e] = TRUE;
         forbidden[flipedge(e)] = TRUE;
      }
   }

   /* traverse all edges */
   for( i = 0; i < nnodes; i++ )
   {
      if( g->grad[i] <= 0 )
         continue;

      nnterms1 = 1;
      if( Is_term(g->term[i]) )
      {
         termdist1[0] = 0.0;
         neighbterms1[0] = i;
      }

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         e = enext;
         i2 = g->head[e];
         enext = g->oeat[e];

         if( i2 < i || !g->mark[i2] )
            continue;

         ecost = g->cost[e];

         /* is i a terminal? If not, get three closest terminals of distance <= ecost */
         if( !Is_term(g->term[i]) )
         {
            nnterms1 = getlecloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);

            if( nnterms1 == 0 )
               continue;
         }

         if( Is_term(g->term[i2]) )
         {
            nnterms2 = 1;
            termdist2[0] = 0.0;
            neighbterms2[0] = i2;
         }
         else
         {
            /* get closest terminals of distance <= ecost */
            nnterms2 = getlecloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);

            if( nnterms2 == 0 )
               continue;
         }

         for( j = 0; j < nnterms1; j++ )
         {
            /* has edge already been deleted? */
            if( g->oeat[e] == EAT_FREE )
               break;

            tj = neighbterms1[j];

            assert(tj >= 0);

            for( k = 0; k < nnterms2; k++ )
            {
               tk = neighbterms2[k];

               assert(tk >= 0);
               assert(Is_term(g->term[tk]));
               assert(Is_term(g->term[tj]));

               if( tj == tk )
               {
                  if( SCIPisGE(scip, termdist1[j], termdist2[k] ) )
                     dist = termdist1[j];
                  else
                     dist = termdist2[k];

                  assert(SCIPisGE(scip, ecost, dist));

                  if( SCIPisEQ(scip, dist, ecost) )
                     if( !sddeltable(scip, g, vnoi, vbase, forbidden, j, k, i, i2, e, nnodes ) )
                        continue;

                  if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
                     continue;

                  graph_edge_del(scip, g, e, TRUE);
                  (*nelims)++;
                  break;
               }
               else
               {
                  /* get sd between (terminals) tj and tk */
                  nj = nodesid[tj];
                  nk = nodesid[tk];

                  assert(nj != nk);

                  l = nj;
                  dist = 0.0;
                  mstsdist[l] = 0.0;

                  /* go down to the root */
                  while( l != 0 )
                  {
                     ne = mst[l].edge;

                     assert(netgraph->head[ne] == l);

                     l = netgraph->tail[ne];

                     if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                        dist = netgraph->cost[ne];

                     mstsdist[l] = dist;

                     if( l == nk )
                        break;
                  }

                  if( l == nk )
                  {
                     l = 0;
                  }
                  else
                  {
                     l = nk;
                     dist = 0.0;
                  }

                  while( l != 0 )
                  {
                     ne = mst[l].edge;
                     l = netgraph->tail[ne];

                     if( SCIPisGT(scip, netgraph->cost[ne], dist) )
                        dist = netgraph->cost[ne];

                     if( mstsdist[l] >= 0 )
                     {
                        if( SCIPisGT(scip, mstsdist[l], dist) )
                           dist = mstsdist[l];
                        break;
                     }
                     if( l == 0 )
                        assert(nj == 0);
                  }

                  /* restore */
                  l = nj;
                  mstsdist[l] = -1.0;

                  while( l != 0 )
                  {
                     ne = mst[l].edge;
                     l = netgraph->tail[ne];
                     mstsdist[l] = -1.0;
                     if( l == nk )
                        break;
                  }

                  assert(SCIPisGT(scip, dist, 0.0));

                  if( SCIPisLT(scip, dist, termdist1[j]) )
                     dist = termdist1[j];

                  if( SCIPisLT(scip, dist, termdist2[k]) )
                     dist = termdist2[k];

                  if( SCIPisGE(scip, ecost, dist) )
                  {
                     if( SCIPisEQ(scip, ecost, dist) )
                        if( !(sddeltable(scip, g, vnoi, vbase, forbidden, j, k, i, i2, e, nnodes)) )
                           continue;

                     assert(SCIPisGE(scip, ecost, termdist1[j]));
                     assert(SCIPisGE(scip, ecost, termdist2[k]));

                     if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
                        continue;

                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               } /* tj != tk (else) */
            } /* k < nnterms2 */
         } /* j < nnterms1 */

      } /* while( enext != EAT_LAST ) */
   }

   if( nodereplacing )
   {
      SCIP_CALL( reduce_bd34WithSd(scip, g, netgraph, mst, vnoi, mstsdist, vbase, nodesid, nelims) );
   }

   /* free memory*/
   SCIPfreeBufferArray(scip, &mst);
   graph_path_exit(scip, netgraph);
   graph_free(scip, &netgraph, TRUE);

   return SCIP_OKAY;
}


/** implied-profit special distance test */
SCIP_RETCODE reduce_sdBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   SD*                   sdistance,          /**< special distances storage */
   GRAPH*                g,                  /**< graph structure */
   int*                  nelims              /**< number of eliminations */
)
{
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Real maxmstcost = reduce_sdgraphGetMaxCost(sdistance->sdgraph);

   assert(LT(maxmstcost, FARAWAY));

   graph_mark(g);

   SCIPdebugMessage("Starting SD biased... \n");

   /* traverse all edges */
   for( int i = 0; i < nnodes; i++ )
   {
      int enext;

      if( !g->mark[i] )
         continue;

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         SCIP_Real sd;
         const int e = enext;
         const int i2 = g->head[e];
         const SCIP_Real ecost = g->cost[e];

         enext = g->oeat[e];

         if( i2 < i || !g->mark[i2] )
            continue;

         sd = sdGetSd(g, i, i2, maxmstcost, ecost, sdistance);

         // todo LT
         if( SCIPisLT(scip, sd, ecost) )
         {
#ifdef SCIP_DEBUG
            SCIPdebugMessage("SD biased deletes (sd=%f):  ", sd);
            graph_edge_printInfo(g, e);
#endif

            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;

            break;
         }
      }
   }

   return SCIP_OKAY;
}



/** SD test for PC */
SCIP_RETCODE reduce_sdPc(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   PATH*                 vnoi,               /**< Voronoi data structure */
   int*                  heap,               /**< heap array */
   int*                  state,              /**< array to store state of a node during Voronoi computation*/
   int*                  vbase,              /**< Voronoi base to each node */
   int*                  nodesid,            /**< array */
   int*                  nodesorg,           /**< array */
   int*                  nelims              /**< pointer to store number of eliminated edges */
   )
{
   GRAPH* netgraph;
   SCIP_Real termdist1[4];
   SCIP_Real termdist2[4];
   int neighbterms1[4];
   int neighbterms2[4];

   int j;
   int maxnedges;
   const int root = g->source;
   const int nnodes = g->knots;
   const int nterms = g->terms;
   const int nedges = g->edges;

   assert(g != NULL);
   assert(heap != NULL);
   assert(vnoi != NULL);
   assert(state != NULL);
   assert(vbase != NULL);
   assert(scip  != NULL);
   assert(nelims != NULL);
   assert(nodesid != NULL);
   assert(nodesorg != NULL);
   assert(!g->extended);

   *nelims = 0;

   if( nterms <= 1 )
      return SCIP_OKAY;
   else
   {
      const SCIP_Longint longedges = (SCIP_Longint) nedges;
      const SCIP_Longint longtermsq = (SCIP_Longint) (nterms - 1) * (SCIP_Longint) nterms;

      if( longedges <= longtermsq )
         maxnedges = nedges;
      else
         maxnedges = ((nterms - 1) * nterms);
   }

   /* compute nearest four terminals to each non-terminal */
   graph_get4nextTermsPaths(g, g->cost, g->cost, vnoi, vbase, state);

   /*
    * construct auxiliary graph to compute paths between terminals
    */

   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1) );

   for( int k = 0; k < 4; k++ )
   {
      termdist1[k] = FARAWAY;
      termdist2[k] = FARAWAY;
   }

   j = 0;
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         graph_knot_add(netgraph, -1);
         nodesid[k] = j;
         nodesorg[j++] = k;
      }
      else
      {
         nodesid[k] = UNKNOWN;
      }
   }

   assert(netgraph->knots == j);
   assert(netgraph->knots == nterms);

   /* insert Voronoi boundary paths as edges into netgraph */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !g->mark[k] )
         continue;

      for( int e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         int i;
         if( !g->mark[g->head[e]] )
            continue;
         i = vbase[k];
         assert(i != UNKNOWN);
         if( i != vbase[g->head[e]] )
         {
            SCIP_Real ecost;
            int ne;
            const int i2 = vbase[g->head[e]];

            assert(i2 != UNKNOWN);
            assert(Is_term(g->term[i]));
            assert(Is_term(g->term[i2]));
            assert(nodesid[i] >= 0);
            assert(nodesid[i2] >= 0);

            for( ne = netgraph->outbeg[nodesid[i]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
               if( netgraph->head[ne] == nodesid[i2] )
                  break;

            ecost = g->cost[e] + vnoi[g->head[e]].dist + vnoi[g->tail[e]].dist;

            /* edge exists? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(netgraph->head[ne] == nodesid[i2]);
               assert(netgraph->tail[ne] == nodesid[i]);

               if( SCIPisGT(scip, netgraph->cost[ne], ecost) )
               {
                  netgraph->cost[ne]            = ecost;
                  netgraph->cost[flipedge(ne)] = ecost;
                  assert(ne <= maxnedges);
               }
            }
            else
            {
               graph_edge_add(scip, netgraph, nodesid[i], nodesid[i2], ecost, ecost);
               assert(netgraph->edges <= maxnedges);
            }
         }
      }
   }

   /* compute four close terminals to each terminal */
   SCIP_CALL( graph_get4nextTTerms(scip, g, g->cost, vnoi, vbase, heap, state) );

   /* traverse all edges */
   for( int i = 0; i < nnodes; i++ )
   {
      int enext;
      if( !g->mark[i] )
         continue;

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         SCIP_Real ecost;
         int e = enext;
         int nnterms1;
         int nnterms2;
         const int i2 = g->head[e];
         enext = g->oeat[e];

         if( i2 < i || Is_term(g->term[i2]) || !g->mark[i2] )
            continue;
         ecost = g->cost[e];

         /* @todo: fix */
#if 1
         if( Is_term(g->term[i]) )
            nnterms1 = getcloseterms2term(scip, g, netgraph, termdist1, ecost, neighbterms1, nodesid, nodesorg, i);
         else
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);
#else
         if( Is_term(g->term[i]) )
            nnterms1 = 0;
         else
            nnterms1 = getcloseterms(scip, vnoi, termdist1, ecost, vbase, neighbterms1, i, nnodes);
#endif

         if( nnterms1 == 0 )
            continue;

         /* @todo: fix */
#if 1
         if( Is_term(g->term[i2]) )
            nnterms2 = getcloseterms2term(scip, g, netgraph, termdist2, ecost, neighbterms2, nodesid, nodesorg, i2);
         else
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
#else
         if( Is_term(g->term[i2]) )
            nnterms2 = 0;
         else
            nnterms2 = getcloseterms(scip, vnoi, termdist2, ecost, vbase, neighbterms2, i2, nnodes);
#endif

         if( nnterms2 == 0 )
            continue;

         /* todo: mark nearest terminals! */
         for( j = 0; j < nnterms1; j++ )
         {
            int tj;

            /* has edge already been deleted? */
            if( g->oeat[e] == EAT_FREE )
               break;

            tj = neighbterms1[j];

            assert(tj >= 0);
            assert(Is_term(g->term[tj]));

            for( int k = 0; k < nnterms2; k++ )
            {
               const int tk = neighbterms2[k];

               assert(tk >= 0);
               assert(Is_term(g->term[tk]));

               if( tj == tk )
               {
                  if( SCIPisGT(scip, ecost, termdist1[j] + termdist2[k] - g->prize[tj]) || tj == root )
                  {
                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               }
               else
               {
                  SCIP_Real necost = FARAWAY;
                  int e2;
                  int pos;

                  /* get distance between terminals */
                  for( e2 = netgraph->outbeg[nodesid[tj]]; e2 != EAT_LAST; e2 = netgraph->oeat[e2] )
                  {
                     if( netgraph->head[e2] == nodesid[tk] )
                     {
                        necost = netgraph->cost[e2];
                        break;
                     }
                  }
                  pos = tj;
                  for( int l = 0; l < 4; l++ )
                  {
                     if( vbase[pos] == UNKNOWN )
                        break;
                     if( vbase[pos] == tk && SCIPisLT(scip, vnoi[pos].dist, necost) )
                     {
                        necost = vnoi[pos].dist;
                        break;
                     }
                     pos += nnodes;
                  }

                  if( SCIPisGT(scip, ecost, necost)
                     && SCIPisGT(scip, ecost, necost + termdist1[j] - g->prize[tj])
                     && SCIPisGT(scip, ecost, necost + termdist2[k] - g->prize[tk])
                     && SCIPisGT(scip, ecost, necost + termdist1[j] + termdist2[k] - g->prize[tj] - g->prize[tk]) )
                  {
                     SCIPdebugMessage("SDPC delete: %d %d (%d)\n", g->tail[e], g->head[e], e);
                     graph_edge_del(scip, g, e, TRUE);
                     (*nelims)++;
                     break;
                  }
               }
            }
         }
      }
   }

   SCIPdebugMessage("SDPC eliminations: %d \n", *nelims);
   graph_free(scip, &netgraph, TRUE);

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/** get SD to a single edge*/
SCIP_RETCODE reduce_getSd(
   SCIP* scip,
   GRAPH* g,
   PATH*  pathtail,
   PATH*  pathhead,
   SCIP_Real*    sdist,
   SCIP_Real     distlimit,
   int*    heap,
   int*    statetail,
   int*    statehead,
   int*    memlbltail,
   int*    memlblhead,
   int     i,
   int     i2,
   int     limit,
   SCIP_Bool pc,
   SCIP_Bool mw
   )
{
   SCIP_Real sd;
   SCIP_Real dist;
   int k;
   int l;
   int e;
   int nnodes;
   int nlbltail;
   int nlblhead;
   const SCIP_Bool pcmw = (pc || mw);

   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(sdist != NULL);

   nnodes = g->knots;

   /* start from tail */
   graph_sdPaths(g, pathtail, g->cost, distlimit, heap, statetail, memlbltail, &nlbltail, i, i2, limit);

   /* test whether edge e can be eliminated */
   graph_sdPaths(g, pathhead, g->cost, distlimit, heap, statehead, memlblhead, &nlblhead, i2, i, limit);

   sd = FARAWAY;

   /* get restore state and path of tail and head */
   for( k = 0; k < nlbltail; k++ )
   {
      l = memlbltail[k];
      assert(statetail[l]     != UNKNOWN);
      if( statehead[l] != UNKNOWN )
      {
         assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
         assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));
         if( Is_term(g->term[l]) )
         {
            dist = 0.0;
            if( SCIPisLT(scip, dist, pathhead[l].dist) )
               dist = pathhead[l].dist;
            if( SCIPisLT(scip, dist, pathtail[l].dist) )
               dist = pathtail[l].dist;
            if( pcmw && SCIPisLT(scip, dist, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
               dist = pathhead[l].dist + pathtail[l].dist - g->prize[l];
            if( SCIPisGT(scip, sd, dist) )
               sd = dist;
         }
         else
         {
            if( mw && l != i && l != i2 )
               assert(SCIPisLE(scip, g->prize[l], 0.0));
            if( mw && SCIPisLT(scip, g->prize[l], 0.0) )
               dist = pathhead[l].dist + pathtail[l].dist + g->prize[l];
            else
               dist = pathhead[l].dist + pathtail[l].dist;
            if( SCIPisGT(scip, sd, dist) )
               sd = dist;
         }
      }

      statetail[l]     = UNKNOWN;
      pathtail[l].dist = FARAWAY;
      pathtail[l].edge = UNKNOWN;
   }
   /* restore state and path of tail and head */
   for( k = 0; k < nlblhead; k++ )
   {
      l = memlblhead[k];
      statehead[l]     = UNKNOWN;
      pathhead[l].dist = FARAWAY;
      pathhead[l].edge = UNKNOWN;
   }


   for( k = 0; k < nnodes; k++ )
   {
      assert(statetail[k]     == UNKNOWN);
      assert(pathtail[k].dist == FARAWAY);
      assert(pathtail[k].edge == UNKNOWN);
      assert(statehead[k]     == UNKNOWN);
      assert(pathhead[k].dist == FARAWAY);
      assert(pathhead[k].edge == UNKNOWN);
   }


   l = 0;
   /* compare restricted sd with edge cost (if existing) */
   for( e = g->outbeg[i]; (l++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
         if( mw )
            sd = 0.0;
         else if( SCIPisGT(scip, sd, g->cost[e]) )
            sd = g->cost[e];
         break;
      }
   }

   *sdist = sd;
   return SCIP_OKAY;
}


/** get SD to a single edge*/
SCIP_RETCODE reduce_getSdPcMw(
   SCIP* scip,
   const GRAPH* g,
   PATH*  pathtail,
   PATH*  pathhead,
   SCIP_Real*    sdist,
   SCIP_Real     distlimit,
   int*    heap,
   int*    statetail,
   int*    statehead,
   int*    memlbltail,
   int*    memlblhead,
   int*    pathmaxnodetail,
   int*    pathmaxnodehead,
   int     i,
   int     i2,
   int     limit
   )
{
   SCIP_Real sd;
   int nlbltail;
   int nlblhead;
   const int nnodes = g->knots;
   const SCIP_Bool mw = g->stp_type == STP_MWCSP;

   assert((g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG) || mw);
   assert(g != NULL);
   assert(scip  != NULL);
   assert(pathtail != NULL);
   assert(pathhead != NULL);
   assert(heap != NULL);
   assert(statetail != NULL);
   assert(statehead != NULL);
   assert(memlbltail != NULL);
   assert(memlblhead != NULL);
   assert(pathmaxnodetail != NULL);
   assert(pathmaxnodehead != NULL);
   assert(sdist != NULL);

   graph_path_PcMwSd(scip, g, pathtail, g->cost, distlimit, pathmaxnodetail, heap, statetail, NULL, memlbltail, &nlbltail, i, i2, limit);
   graph_path_PcMwSd(scip, g, pathhead, g->cost, distlimit, pathmaxnodehead, heap, statehead, statetail, memlblhead, &nlblhead, i2, i, limit);

   sd = FARAWAY;

   /* get restore state and path of tail and head */
   for( int k = 0; k < nlbltail; k++ )
   {
      const int l = memlbltail[k];
      assert(statetail[l] != UNKNOWN);

      if( statehead[l] != UNKNOWN )
      {
         SCIP_Real dist = FARAWAY;
         const int tailmaxterm = pathmaxnodetail[l];
         const int headmaxterm = pathmaxnodehead[l];

         assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
         assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));
         assert(tailmaxterm != i && headmaxterm != i);
         assert(tailmaxterm != i2 && headmaxterm != i2);

         /* any terminal on the path? */
         if( tailmaxterm >= 0 || headmaxterm >= 0 )
         {
            if( tailmaxterm == headmaxterm )
            {
               assert(tailmaxterm == l);
               assert(SCIPisPositive(scip, g->prize[tailmaxterm]));

               dist = miscstp_maxReal((SCIP_Real []) {
                      pathhead[headmaxterm].dist,
                      pathtail[tailmaxterm].dist,
                      pathhead[l].dist + pathtail[l].dist - g->prize[l]
                     }, 3);
               SCIPdebugMessage("sd1 %f \n", dist);
            }
            else if( tailmaxterm >= 0 && headmaxterm >= 0 )
            {
               const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;
               const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;

               assert(tailmaxterm != headmaxterm);
               assert(!SCIPisNegative(scip, distl2tailmax));
               assert(!SCIPisNegative(scip, distl2headmax));
               assert(SCIPisPositive(scip, g->prize[tailmaxterm]) && SCIPisPositive(scip, g->prize[headmaxterm]));

               dist = miscstp_maxReal((SCIP_Real []) {
                      pathhead[headmaxterm].dist,
                      pathtail[tailmaxterm].dist,
                      distl2tailmax + distl2headmax,
                      distl2tailmax + pathhead[l].dist - g->prize[headmaxterm],
                      distl2headmax + pathtail[l].dist - g->prize[tailmaxterm],
                      pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm] - g->prize[headmaxterm]
                     }, 6);
               SCIPdebugMessage("sd2 %f \n", dist);
            }
            else if( tailmaxterm >= 0 )
            {
               const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;

               assert(headmaxterm < 0);
               assert(SCIPisPositive(scip, g->prize[tailmaxterm]));

               dist = miscstp_maxReal((SCIP_Real []) {
                      pathtail[tailmaxterm].dist,
                      distl2tailmax + pathhead[l].dist,
                      pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm]
                     }, 3);
               SCIPdebugMessage("sd3 %f \n", dist);
            }
            else if( headmaxterm >= 0 )
            {
               const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;

               assert(tailmaxterm < 0);
               assert(SCIPisPositive(scip, g->prize[headmaxterm]));

               dist = miscstp_maxReal((SCIP_Real []) {
                      pathhead[headmaxterm].dist,
                      distl2headmax + pathtail[l].dist,
                      pathhead[l].dist + pathtail[l].dist - g->prize[headmaxterm]
                     }, 3);
               SCIPdebugMessage("sd4 %f \n", dist);
            }
         }
         else
         {
            dist = pathhead[l].dist + pathtail[l].dist;
         }

         if( dist < sd )
            sd = dist;

         if( Is_term(g->term[l]) )
         {
            dist = miscstp_maxReal((SCIP_Real []) {
                   pathhead[l].dist,
                   pathtail[l].dist,
                   pathhead[l].dist + pathtail[l].dist - g->prize[l]
                  }, 3);
            if( dist < sd )
               sd = dist;
         }
      }
   }

   /* restore state and path of tail and head */

   for( int k = 0; k < nlbltail; k++ )
   {
      const int l = memlbltail[k];
      statetail[l] = UNKNOWN;
      pathtail[l].dist = FARAWAY;
      pathtail[l].edge = UNKNOWN;
      pathmaxnodetail[l] = -1;
   }

   for( int k = 0; k < nlblhead; k++ )
   {
      const int l = memlblhead[k];
      statehead[l] = UNKNOWN;
      pathhead[l].dist = FARAWAY;
      pathhead[l].edge = UNKNOWN;
      pathmaxnodehead[l] = -1;
   }

   for( int k = 0; k < nnodes; k++ )
   {
      assert(statetail[k]     == UNKNOWN);
      assert(pathtail[k].dist == FARAWAY);
      assert(pathtail[k].edge == UNKNOWN);
      assert(statehead[k]     == UNKNOWN);
      assert(pathhead[k].dist == FARAWAY);
      assert(pathhead[k].edge == UNKNOWN);
      assert(pathmaxnodehead[k] == -1);
      assert(pathmaxnodetail[k] == -1);
   }

   /* compare restricted sd with edge cost (if existing) */
   for( int e = g->outbeg[i], count = 0; (count++ <= limit) && (e != EAT_LAST); e = g->oeat[e] )
   {
      if( g->head[e] == i2 )
      {
         if( mw )
            sd = 0.0;
         else if( sd > g->cost[e] )
            sd = g->cost[e];
         break;
      }
   }

   *sdist = sd;

   return SCIP_OKAY;
}



/*  longest edge reduction test from T. Polzin's "Algorithms for the Steiner problem in networks" (Lemma 20)
 *
 *  *** DEPRECATED! *** */
SCIP_RETCODE reduce_ledge(
   SCIP*   scip,
   GRAPH*  g,
   PATH*   vnoi,
   int* heap,
   int* state,
   int* vbase,
   int* nelims,
   int* edgestate
)
{
   GRAPH* netgraph;
   PATH* mst;
   SCIP_Real cost;
   SCIP_Real maxcost;
   int v1;
   int v2;
   int k;
   int e;
   int ne;
   int nedges;
   int nnodes;
   int nterms;
   int maxnedges;
   int netnnodes;
   int* nodesid;
   int* edgeorg;
   SCIP_Bool checkstate = (edgestate != NULL);
   STP_Bool* blocked;

   assert(g != NULL);
   assert(vnoi != NULL);
   assert(heap != NULL);
   assert(state != NULL);
   assert(vbase != NULL);

   *nelims = 0;
   nedges = g->edges;
   nnodes = g->knots;
   assert(graph_valid(scip, g));

   nterms = 0;
   for( k = 0; k < nnodes; k++ )
   {
      g->mark[k] = (g->grad[k] > 0);
      if( Is_term(g->term[k]) && g->mark[k] )
         nterms++;
   }

   if( nterms <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &blocked, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edgeorg, nedges / 2) );

   graph_voronoiTerms(g, g->cost, vnoi, vbase, state);

   if( nedges >= (nterms - 1) * nterms )
      maxnedges = (nterms - 1) * nterms;
   else
      maxnedges = nedges;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodesid, nnodes) );

   /* initialize the new graph */
   SCIP_CALL( graph_init(scip, &netgraph, nterms, maxnedges, 1) );

   e = 0;
   for( k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && g->grad[k] > 0 )
      {
         if( e == 0 )
            graph_knot_add(netgraph, 0);
         else
            graph_knot_add(netgraph, -1);
         nodesid[k] = e++;
      }
      else
      {
         nodesid[k] = UNKNOWN;
      }
   }

   netnnodes = netgraph->knots;
   assert(netnnodes == e);
   assert(netnnodes == nterms);

   for( e = 0; e < nedges / 2; e++ )
   {
      blocked[e] = FALSE;
      edgeorg[e] = UNKNOWN;
   }

   for( k = 0; k < nnodes; k++ )
   {
      for( e = g->outbeg[k]; e != EAT_LAST; e = g->oeat[e] )
      {
         v1 = vbase[k];
         assert(k == g->tail[e]);

         if( v1 != vbase[g->head[e]] )
         {
            v2 = vbase[g->head[e]];
            assert(Is_term(g->term[v1]));
            assert(Is_term(g->term[v2]));
            assert(nodesid[v1] >= 0);
            assert(nodesid[v2] >= 0);

            for( ne = netgraph->outbeg[nodesid[v1]]; ne != EAT_LAST; ne = netgraph->oeat[ne] )
               if( netgraph->head[ne] == nodesid[v2] )
                  break;

            cost = g->cost[e] + vnoi[g->head[e]].dist + vnoi[g->tail[e]].dist;
            /* edge exists? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(netgraph->head[ne] == nodesid[v2]);
               assert(netgraph->tail[ne] == nodesid[v1]);
               if( SCIPisGT(scip, netgraph->cost[ne], cost) )
               {
                  netgraph->cost[ne]            = cost;
                  netgraph->cost[Edge_anti(ne)] = cost;
                  edgeorg[ne / 2] = e;
                  assert(ne <= maxnedges);
               }
            }
            else
            {
               edgeorg[netgraph->edges / 2] = e;
               graph_edge_add(scip, netgraph, nodesid[v1], nodesid[v2], cost, cost);
               assert(netgraph->edges <= maxnedges);
            }
         }
      }
   }
   netgraph->source = 0;

   assert(graph_valid(scip, netgraph));

   for( k = 0; k < netnnodes; k++ )
      netgraph->mark[k] = TRUE;

   /* compute a MST on netgraph */
   SCIP_CALL( SCIPallocBufferArray(scip, &mst, netnnodes) );
   SCIP_CALL( graph_path_init(scip, netgraph) );
   graph_path_exec(scip, netgraph, MST_MODE, 0, netgraph->cost, mst);

   maxcost = -1;
   assert(mst[0].edge == -1);

   for( k = 1; k < netnnodes; k++ )
   {
      assert(netgraph->path_state[k] == CONNECT);
      e = mst[k].edge;
      assert(e >= 0);
      cost = netgraph->cost[e];
      if( SCIPisGT(scip, cost, maxcost) )
         maxcost = cost;

      ne = edgeorg[e / 2];
      blocked[ne / 2] = TRUE;
      for( v1 = g->head[ne]; v1 != vbase[v1]; v1 = g->tail[vnoi[v1].edge] )
         blocked[vnoi[v1].edge / 2] = TRUE;

      for( v1 = g->tail[ne]; v1 != vbase[v1]; v1 = g->tail[vnoi[v1].edge] )
         blocked[vnoi[v1].edge / 2] = TRUE;
      assert(e != EAT_LAST);
   }

   for( k = 0; k < nnodes; k++ )
   {
      e = g->outbeg[k];
      while( e != EAT_LAST )
      {
         assert(e >= 0);
         if( checkstate && (edgestate[e] == EDGE_BLOCKED) )
         {
            e = g->oeat[e];
            continue;
         }

         if( SCIPisGE(scip, g->cost[e], maxcost) && !blocked[e / 2] )
         {
            (*nelims)++;
            v1 = g->oeat[e];
            graph_edge_del(scip, g, e, TRUE);
            e = v1;
         }
         else
         {
            e = g->oeat[e];
         }
      }
   }

   /* graph might have become disconnected */
   if( *nelims > 0 )
   {
      SCIP_CALL( reduceLevel0(scip, g) );
   }

   /* free netgraph and  MST data structure */
   graph_path_exit(scip, netgraph);
   graph_free(scip, &netgraph, TRUE);
   SCIPfreeBufferArray(scip, &mst);
   SCIPfreeBufferArray(scip, &nodesid);
   SCIPfreeBufferArray(scip, &edgeorg);
   SCIPfreeBufferArray(scip, &blocked);

   assert(graph_valid(scip, g));
   return SCIP_OKAY;
}



/** Non-Terminal-Set test*/
SCIP_RETCODE reduce_nts(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   PATH*                 pathtail,           /**< array for internal use */
   PATH*                 pathhead,           /**< array for internal use */
   int*                  heap,               /**< array for internal use */
   int*                  statetail,          /**< array for internal use */
   int*                  statehead,          /**< array for internal use */
   int*                  memlbltail,         /**< array for internal use */
   int*                  memlblhead,         /**< array for internal use */
   int*                  nelims,             /**< point to return number of eliminations */
   int                   limit               /**< limit for edges to consider for each vertex */
   )
{
   PATH mst[STP_BD_MAXDEGREE];
   SCIP_Real ecost[STP_BD_MAXDEGREE];
   int outedge[STP_BD_MAXDEGREE];
   int adjvert[STP_BD_MAXDEGREE];
   GRAPH* auxg;
   int* pathmaxnodetail = NULL;
   int* pathmaxnodehead = NULL;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool pc = g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG;

   SCIPdebugMessage("NTS-Reduction: ");

   assert(scip != NULL);
   assert(g != NULL);
   assert(heap != NULL);
   assert(nelims != NULL);
   assert(limit > 0);

   // todo extra method, also called from bd34

   /* initialize new graph for nts tests */
   SCIP_CALL( graph_init(scip, &auxg, STP_BD_MAXDEGREE, 2 * STP_BD_MAXDNEDGES, 1) );

   for( int k = 0; k < STP_BD_MAXDEGREE; k++ )
   {
      graph_knot_add(auxg, -1);
      auxg->mark[k] = TRUE;
   }

   for( int k = 0; k < STP_BD_MAXDEGREE; k++ )
      for( int k2 = STP_BD_MAXDEGREE - 1; k2 >= k + 1; k2-- )
         graph_edge_add(scip, auxg, k, k2, 1.0, 1.0);

   /* init graph for mst computation */
   SCIP_CALL( graph_path_init(scip, auxg) );

   *nelims = 0;

   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodetail, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodehead, nnodes) );
   }
   else
   {
      for( int  i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
      if( pc )
      {
         pathmaxnodetail[i] = -1;
         pathmaxnodehead[i] = -1;
      }
   }

   for( int degree = 4; degree <= STP_BD_MAXDEGREE; degree++ )
   {
      for( int edge = 0; edge < nedges; edge += 2 )
      {
         /* edge not deleted? */
         if( g->oeat[edge] != EAT_FREE )
         {
            SCIP_Real csum;
            SCIP_Real mstcost;
            const SCIP_Real edgecost = g->cost[edge];
            int edgecount;
            const int tail = g->tail[edge];
            const int head = g->head[edge];
            SCIP_Bool success = TRUE;
            SCIP_Bool duplicate;

            assert(edge >= 0);
            assert(tail >= 0 && head >= 0);

            // todo
            if( g->grad[tail] != 3 || g->grad[head] != 3 )
               continue;

            if( Is_term(g->term[tail]) || Is_term(g->term[head]) )
               continue;

            edgecount = 0;
            for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == head )
                  continue;
               outedge[edgecount] = e;
               ecost[edgecount] = g->cost[e];
               adjvert[edgecount++] = g->head[e];
            }

            for( int e = g->outbeg[head]; e != EAT_LAST; e = g->oeat[e] )
            {
               if( g->head[e] == tail )
                  continue;
               outedge[edgecount] = e;
               ecost[edgecount] = g->cost[e];
               adjvert[edgecount++] = g->head[e];
            }

            assert(edgecount == degree);

            duplicate = FALSE;

            /* check for shared neighbor */
            for( int i = 0; i < degree - 1; i++ )
               for( int j = i + 1; j < degree; j++ )
                  if( adjvert[i] == adjvert[j] && adjvert[j] >= 0 )
                  {
                     adjvert[j] = -adjvert[j] - 1;
                     duplicate = TRUE;
                  }

            // todo reshuffle
            if( duplicate )
            {
               continue;
            }

            assert(g->mark[head] && g->mark[tail]);

            g->mark[head] = FALSE;
            g->mark[tail] = FALSE;

            csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];

            /* compute sd values and check 2-subsets of neighbors */
            for( int k = 0; k < degree - 1 && success; k++ )
            {
               for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  const int k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     SCIP_Real s1 = -1.0;

                     if( pc )
                        SCIP_CALL( reduce_getSdPcMw(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                                    pathmaxnodetail, pathmaxnodehead, adjvert[k], adjvert[k2], limit) );
                     else
                        SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                              adjvert[k], adjvert[k2], limit, pc, FALSE) );

                     assert(s1 >= 0);

                     if( g->tail[outedge[k]] != g->tail[outedge[k2]] )
                     {
                        const SCIP_Real innercost = ecost[k] + ecost[k2] + edgecost;

                        if( SCIPisGT(scip, s1, innercost) )
                        {
          //                 success = FALSE;
            //               break;
                        }
                     }

                     auxg->cost[e] = s1;
                     auxg->cost[flipedge(e)] = s1;
                  }
               }
            }

            g->mark[head] = TRUE;
            g->mark[tail] = TRUE;

            for( int l = 0; l < 4; l++ )
               mst[l].dist = UNKNOWN;

            /* compute MST on all neighbors */
            graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
            mstcost = 0.0;

            /* sum cost of (root == 0) */
            for( int l = 1; l < 4; l++ )
            {
               assert(mst[l].dist != UNKNOWN);
               mstcost += mst[l].dist;
            }

            success = (success && SCIPisGE(scip, csum, mstcost));


            if( success && 0)
            {
               /* compute MST on all 3-subsets of all neighbors */
               for( int k = 0; k < degree; k++ )
               {
                  assert(auxg->mark[k]);

                  auxg->mark[k] = FALSE;
                  mstcost = 0.0;

                  if( k != 0 )
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);
                     for( int l = 1; l < 4; l++ )
                        if( auxg->mark[l] )
                           mstcost += mst[l].dist;
                  }
                  else
                  {
                     graph_path_exec(scip, auxg, MST_MODE, 1, auxg->cost, mst);
                     for( int l = 2; l < 4; l++ )
                        mstcost += mst[l].dist;
                  }

                  auxg->mark[k] = TRUE;

                  if( SCIPisLT(scip, csum - ecost[k], mstcost) )
                  {
                     success = FALSE;
                     break;
                  }
               }
            }

            if( success )
            {
             //  graph_edge_del(scip, g, edge, TRUE);

               (*nelims)++;
            }
         } /* check single edge */
      } /* for all edges */
   } /* for all degrees */


   /* free memory */

   SCIPfreeBufferArrayNull(scip, &pathmaxnodehead);
   SCIPfreeBufferArrayNull(scip, &pathmaxnodetail);

   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/** SDC test for the SAP using a limited version of Dijkstra's algorithm from both endpoints of an arc */
SCIP_RETCODE reduce_sdspSap(
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
   SCIP_Real* costrev;
   int i;
   int k;
   int l;
   int e;
   int i2;
   int enext;
   int nnodes;
   int nedges;
   int nlbltail;
   int nlblhead;

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

   nedges = g->edges;
   nnodes = g->knots;
   *nelims = 0;

   for( i = 0; i < nnodes; i++ )
   {
      g->mark[i] = (g->grad[i] > 0);
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &costrev, nedges) );

   for( e = 0; e < nedges; e++ )
      costrev[e] = g->cost[flipedge(e)];

   /* iterate through all nodes */
   for( i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
         continue;
      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {

         enext = g->oeat[e];
         i2 = g->head[e];

         assert(g->mark[i2]);

         /* start limited dijkstra from i, marking all reached vertices */
         graph_sdPaths(g, pathtail, g->cost, g->cost[e], heap, statetail, memlbltail, &nlbltail, i, i2, limit);

         /* start limited dijkstra from i2, marking all reached vertices */
         graph_sdPaths(g, pathhead, costrev, g->cost[e], heap, statehead, memlblhead, &nlblhead, i2, i, limit);

         sdist = FARAWAY;
#if 0
         if( statetail[i2] != UNKNOWN )
         {
            sdist = pathtail[i2].dist;
            assert(SCIPisGT(scip, FARAWAY, sdist));
         }
         if( statehead[i] != UNKNOWN && SCIPisGT(scip, sdist, pathhead[i].dist) )
         {
            sdist = pathhead[i].dist;
            assert(SCIPisGT(scip, FARAWAY, sdist));
         }
#endif
         /* get restore state and path of tail and head */
         for( k = 0; k < nlbltail; k++ )
         {
            l = memlbltail[k];
            assert(g->mark[l]);
            assert(statetail[l]     != UNKNOWN);
            if( statehead[l] != UNKNOWN )
            {
               assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
               assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));

               if( SCIPisGT(scip, sdist, pathhead[l].dist + pathtail[l].dist) )
                  sdist = pathhead[l].dist + pathtail[l].dist;
            }

            statetail[l]     = UNKNOWN;
            pathtail[l].dist = FARAWAY;
            pathtail[l].edge = UNKNOWN;
         }
         /* restore state and path of tail and head */
         for( k = 0; k < nlblhead; k++ )
         {
            l = memlblhead[k];
            statehead[l]     = UNKNOWN;
            pathhead[l].dist = FARAWAY;
            pathhead[l].edge = UNKNOWN;
         }

#if 1
         for( k = 0; k < nnodes; k++ )
         {
            assert(statetail[k]     == UNKNOWN);
            assert(pathtail[k].dist == FARAWAY);
            assert(pathtail[k].edge == UNKNOWN);
            assert(statehead[k]     == UNKNOWN);
            assert(pathhead[k].dist == FARAWAY);
            assert(pathhead[k].edge == UNKNOWN);
         }
#endif
         /* can edge be deleted? */
         if( SCIPisGE(scip, g->cost[e], sdist) )
         {
            if( SCIPisGE(scip, costrev[e], FARAWAY) )
               graph_edge_del(scip, g, e, TRUE);
            else
               g->cost[e] = FARAWAY;
            (*nelims)++;
         }

         e = enext;
      }
   }

   SCIPfreeBufferArray(scip, &costrev);

   return SCIP_OKAY;
}

/** SD test for PcMw using only limited Dijkstra-like walk from both endpoints of an edge */
SCIP_RETCODE reduce_sdWalk_csr(
   SCIP*                 scip,
   int                   edgelimit,
   const int*            edgestate,
   GRAPH*                g,
   int*                  termmark,
   SCIP_Real*            dist,
   int*                  visitlist,
   STP_Bool*             visited,
   DHEAP*                dheap,
   int*                  nelims
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   int* head_csr;
   int* edgeid_csr;
   SCIP_Real* cost_csr;
   SCIP_Bool* edge_deletable;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g && scip && nelims && visited && visitlist && dheap);
   assert(!g->extended);
   assert(graph_pc_isPcMw(g));

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      dist[i] = FARAWAY;
   }

   graph_heap_clean(TRUE, dheap);
   graph_init_dcsr(scip, g);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   edgeid_csr = dcsr->edgeid;
   cost_csr = dcsr->cost;

   SCIP_CALL( SCIPallocBufferArray(scip, &edge_deletable, nedges / 2) );

   for( int e = 0; e < nedges / 2; e++ )
      edge_deletable[e] = FALSE;

   assert(dcsr && range_csr && edgeid_csr && cost_csr);

   graph_pc_termMarkProper(g, termmark);

   for( int i = 0; i < nnodes; i++ )
   {
      int enext;
      const int start = range_csr[i].start;

      /* degree <= 1? */
      if( range_csr[i].end - start <= 1 )
         continue;

      /* traverse neighbours */
      for( int e = start; e < range_csr[i].end; e = enext )
      {
         SCIP_Bool success;
         const SCIP_Real ecost = cost_csr[e];
         int nvisits;
         const int i2 = head_csr[e];

         assert(g->mark[i] && g->mark[i2]);

         enext = e + 1;

         if( checkstate  )
         {
            const int orgedge = edgeid_csr[e];
            if( edgestate[orgedge] == EDGE_BLOCKED )
               continue;
         }

         success = graph_sdWalks_csr(scip, g, termmark, ecost, i, i2, edgelimit, dist, visitlist, &nvisits, dheap, visited);
         sdwalkReset(nnodes, nvisits, visitlist, dist, dheap->position, visited);
         graph_heap_clean(FALSE, dheap);

         if( success )
         {
            edge_deletable[edgeid_csr[e] / 2] = TRUE;
            graph_dcsr_deleteEdgeBi(scip, dcsr, e);

            (*nelims)++;
            enext--;

            if( range_csr[i].end - start <= 1 )
               break;
         }
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges / 2; e++ )
   {
      if( edge_deletable[e] )
         assert(dcsr->id2csredge[e * 2] == -1);
      else if( g->oeat[e * 2] != EAT_FREE )
         assert(dcsr->id2csredge[e * 2] != -1 || !g->mark[g->tail[e * 2]] || !g->mark[g->head[e * 2]]);
   }
#endif

   graph_edge_delBlocked(scip, g, edge_deletable, TRUE);

   SCIPfreeBufferArray(scip, &edge_deletable);

   graph_free_dcsr(scip, g);

   return SCIP_OKAY;
}


/** SD test */
SCIP_RETCODE reduce_sdEdgeCliqueStar(
   SCIP*                 scip,
   int                   edgelimit,
   GRAPH*                g,
   int*                  nelims
)
{
   const int nnodes = graph_get_nNodes(g);
   SDPROFIT* sdprofit;
   SCIP_Real sds = FARAWAY;
   int cliquenodes[2];
   SDCLIQUE cliquedata = { .dijkdata = NULL, .cliquenodes = cliquenodes, .ncliquenodes = 2, .sds = &sds };

   SCIP_CALL( graph_dijkLimited_init(scip, g, &(cliquedata.dijkdata)) );
   graph_dijkLimited_clean(g, (cliquedata.dijkdata));
   cliquedata.dijkdata->edgelimit = edgelimit;
   SCIP_CALL( reduce_sdprofitInit(scip, g, &sdprofit) );

   graph_mark(g);

   SCIPdebugMessage("Starting SD biased... \n");

   /* traverse all edges */
   for( int i = 0; i < nnodes; i++ )
   {
      int enext;
      if( !g->mark[i] )
         continue;

      enext = g->outbeg[i];
      while( enext != EAT_LAST )
      {
         const int e = enext;
         const int i2 = g->head[e];
         const SCIP_Real ecost = g->cost[e];

         enext = g->oeat[e];

         if( i2 < i || !g->mark[i2] )
            continue;

         sds = ecost;
         cliquenodes[0] = i;
         cliquenodes[1] = i2;
         SCIP_CALL( graph_sdComputeCliqueStar(scip, g, sdprofit, &cliquedata) );

         graph_dijkLimited_reset(g, cliquedata.dijkdata);

         // todo LT
         if( SCIPisLT(scip, sds, ecost) )
         {
#ifdef SCIP_DEBUG
            printf("SD biased deletes (sd=%f):  ", sds);
            graph_edge_printInfo(g, e);
#endif

            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;

            break;
         }
      }
   }

   reduce_sdprofitFree(scip, &sdprofit);
   graph_dijkLimited_free(scip, &(cliquedata.dijkdata));

   return SCIP_OKAY;
}


/** SD test for PcMw using limited Dijkstra-like walk from both endpoints of an edge */
SCIP_RETCODE reduce_sdWalkTriangle(
   SCIP*                 scip,
   int                   edgelimit,
   const int*            edgestate,
   GRAPH*                g,
   int*                  termmark,
   SCIP_Real*            dist,
   int*                  visitlist,
   STP_Bool*             visited,
   DHEAP*                dheap,
   int*                  nelims
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   int* head_csr;
   int* edgeid_csr;
   int* state2;
   int* visitlist2;
   SCIP_Real* dist2;
   SCIP_Real* cost_csr;
   SCIP_Real* prizeoffset;
   SCIP_Bool* edge_deletable;
   STP_Bool* visited2;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g && scip && nelims && visited && visitlist && dheap);
   assert(!g->extended);
   assert(graph_pc_isPcMw(g));

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   graph_heap_clean(TRUE, dheap);
   graph_init_dcsr(scip, g);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   edgeid_csr = dcsr->edgeid;
   cost_csr = dcsr->cost;

   SCIP_CALL( SCIPallocBufferArray(scip, &dist2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &state2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &visited2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &visitlist2, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prizeoffset, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &edge_deletable, nedges / 2) );

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      visited2[i] = FALSE;
      dist[i] = FARAWAY;
      dist2[i] = FARAWAY;
      state2[i] = UNKNOWN;
   }

   for( int e = 0; e < nedges / 2; e++ )
      edge_deletable[e] = FALSE;

   assert(dcsr && range_csr && edgeid_csr && cost_csr);

   graph_pc_termMarkProper(g, termmark);

   /* main loop */
   for( int i = 0; i < nnodes; i++ )
   {
      int enext;
      const int start = range_csr[i].start;

      /* degree <= 1? */
      if( range_csr[i].end - start <= 1 )
         continue;

      /* traverse neighbors of i */
      for( int e = start; e < range_csr[i].end; e = enext )
      {
         SCIP_Bool success;
         const SCIP_Real ecost = cost_csr[e];
         int nvisits;
         const int i2 = head_csr[e];

         assert(g->mark[i] && g->mark[i2]);

         enext = e + 1;

         /* avoid double checking */
         if( i2 < i )
            continue;

         if( checkstate  )
         {
            const int orgedge = edgeid_csr[e];
            if( edgestate[orgedge] == EDGE_BLOCKED )
               continue;
         }

         success = graph_sdWalksTriangle(scip, g, termmark, NULL, ecost, i, i2, edgelimit, NULL, dist, visitlist, &nvisits, dheap, visited);

         /* could not reach i2? */
         if( !success )
         {
            int nvisits2;
            int* const state = dheap->position;

            dheap->position = state2;
            graph_heap_clean(FALSE, dheap);

#ifndef NDEBUG
            for( int k = 0; k < nnodes; k++ )
               prizeoffset[k] = -FARAWAY;
#endif

            /* run from i2 */
            success = graph_sdWalksTriangle(scip, g, termmark, state, ecost, i2, i, edgelimit, prizeoffset, dist2, visitlist2, &nvisits2, dheap, visited2);

            /* could not reach i1? */
            if( !success )
            {
               assert(nvisits2 > 0 && visitlist2[0] == i2);

               /* maybe we can connect two walks */
               for( int k = 1; k < nvisits2; k++ )
               {
                  const int node = visitlist2[k];
                  SCIP_Real dist_combined;

                  assert(visited2[node]);
                  assert(node != i && node != i2);

                  if( !visited[node] )
                     continue;

                  dist_combined = dist[node] + dist2[node];
                  assert(dist_combined < FARAWAY);

                  if( termmark[node] != 0 )
                  {
                     dist_combined += prizeoffset[node];
                     assert(prizeoffset[node] >= 0.0);
                  }

                  if( SCIPisLE(scip, dist_combined, ecost) )
                  {
                     success = TRUE;
                     break;
                  }
               }
            }

            dheap->position = state;
            sdwalkReset(nnodes, nvisits2, visitlist2, dist2, state2, visited2);
         }

         sdwalkReset(nnodes, nvisits, visitlist, dist, dheap->position, visited);
         graph_heap_clean(FALSE, dheap);

         if( success )
         {
            edge_deletable[edgeid_csr[e] / 2] = TRUE;
            graph_dcsr_deleteEdgeBi(scip, dcsr, e);

            (*nelims)++;
            enext--;

            if( range_csr[i].end - start <= 1 )
               break;
         }
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges / 2; e++ )
   {
      if( edge_deletable[e] )
         assert(dcsr->id2csredge[e * 2] == -1);
      else if( g->oeat[e * 2] != EAT_FREE )
         assert(dcsr->id2csredge[e * 2] != -1 || !g->mark[g->tail[e * 2]] || !g->mark[g->head[e * 2]]);
   }
#endif

   graph_edge_delBlocked(scip, g, edge_deletable, TRUE);

   SCIPfreeBufferArray(scip, &edge_deletable);
   SCIPfreeBufferArray(scip, &prizeoffset);
   SCIPfreeBufferArray(scip, &visitlist2);
   SCIPfreeBufferArray(scip, &visited2);
   SCIPfreeBufferArray(scip, &state2);
   SCIPfreeBufferArray(scip, &dist2);

   graph_free_dcsr(scip, g);

   return SCIP_OKAY;
}

/** SD star test for PcMw and SPG */
SCIP_RETCODE reduce_sdStar(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edgelimit,          /**< limit */
   const int*            edgestate,          /**< state array or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            dist,               /**< vertex distances */
   int*                  star_base,          /**< array of size nnodes */
   int*                  visitlist,          /**< array of size nnodes */
   STP_Bool*             visited,            /**< array of size nnodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   int*                  nelims              /**< point to store number of deleted edges */
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   int* head_csr;
   int* edgeid_csr;
   SCIP_Bool* edge_deletable;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g && scip && nelims && visited && visitlist && dheap && star_base);
   assert(!graph_pc_isPcMw(g) || !g->extended);

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   graph_heap_clean(TRUE, dheap);
   graph_init_dcsr(scip, g);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   edgeid_csr = dcsr->edgeid;

   SCIP_CALL( SCIPallocBufferArray(scip, &edge_deletable, nedges / 2) );

   for( int e = 0; e < nedges / 2; e++ )
      edge_deletable[e] = FALSE;

   assert(dcsr && range_csr && edgeid_csr);

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      dist[i] = FARAWAY;
      star_base[i] = SDSTAR_BASE_UNSET;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      SCIP_Bool runloop;

      if( !g->mark[i] )
      {
         assert(g->grad[i] == 2 || g->grad[i] == 0 || i == g->source);
         continue;
      }

      runloop = TRUE;

      while( runloop )
      {
         SCIP_Bool success;
         int nvisits;
         const int start = range_csr[i].start;

         if( range_csr[i].end - start <= 1 )
            break;

         runloop = FALSE;

         /* do the actual star run */
         graph_sdStar(scip, g, FALSE, i, edgelimit, star_base, dist, visitlist, &nvisits, dheap, visited, &success);

         if( success )
         {
            int enext;

            /* check all star nodes (neighbors of i) */
            for( int e = start; e < range_csr[i].end; e = enext )
            {
               const int starnode = head_csr[e];
               const int starbase = star_base[starnode];
               assert(star_base[starnode] >= 0);
               assert(SCIPisLE(scip, dist[starnode], dcsr->cost[e]));
               assert(star_base[starnode] == starnode || star_base[starnode] >= 0);

               enext = e + 1;

               if( checkstate )
               {
                  const int orgedge = edgeid_csr[e];
                  if( edgestate[orgedge] == EDGE_BLOCKED )
                     continue;
               }

               /* shorter path to current star node found? */
               if( starnode != starbase )
               {
                  assert(star_base[starbase] != SDSTAR_BASE_UNSET);

                  /* path still valid? todo really necessary? */
                  if( star_base[starbase] != SDSTAR_BASE_KILLED )
                  {
                     star_base[starnode] = SDSTAR_BASE_KILLED;
                     edge_deletable[edgeid_csr[e] / 2] = TRUE;
                     graph_dcsr_deleteEdgeBi(scip, dcsr, e);

                     (*nelims)++;
                     enext--;
                  }
                  else
                  {
                     runloop = TRUE;
                  }
               }
            } /* traverse star nodes */
         } /* if success */

         sdStarReset(nnodes, nvisits, visitlist, star_base, dist, visited, dheap);
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges / 2; e++ )
   {
      if( edge_deletable[e] )
         assert(dcsr->id2csredge[e * 2] == -1);
      else if( g->oeat[e * 2] != EAT_FREE )
         assert(dcsr->id2csredge[e * 2] != -1 || !g->mark[g->tail[e * 2]] || !g->mark[g->head[e * 2]]);
   }
#endif

   graph_edge_delBlocked(scip, g, edge_deletable, TRUE);

   SCIPfreeBufferArray(scip, &edge_deletable);

   graph_free_dcsr(scip, g);

   return SCIP_OKAY;
}


/** SD star test for PcMw and SPG */
SCIP_RETCODE reduce_sdStarBiased(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edgelimit,          /**< limit */
   const int*            edgestate,          /**< state array or NULL */
   GRAPH*                g,                  /**< graph data structure */
   int*                  nelims              /**< point to store number of deleted edges */
   )
{
   DIJK* dijkdata;
   SDPROFIT* sdprofit;

   int* RESTRICT star_base;
   SCIP_Bool* RESTRICT edge_deletable;
   const int nnodes = graph_get_nNodes(g);

   assert(scip && nelims);
   assert(edgelimit > 0);

   graph_init_dcsr(scip, g);

   SCIP_CALL( sdStarInit(scip, g, edgelimit, &dijkdata, &star_base, &edge_deletable) );
   SCIP_CALL( reduce_sdprofitInit(scip, g, &sdprofit) );

   for( int i = 0; i < nnodes; i++ )
   {
      if( !g->mark[i] )
      {
         assert(g->grad[i] == 2 || g->grad[i] == 0 || i == g->source);
         continue;
      }

      SCIP_CALL( sdStarBiasedProcessNode(scip, i, edgestate, sdprofit, g, dijkdata, star_base, edge_deletable, nelims) );
   }

   reduce_sdprofitFree(scip, &sdprofit);
   sdStarFinalize(scip, g, &dijkdata, &star_base, &edge_deletable);

   graph_free_dcsr(scip, g);

   return SCIP_OKAY;
}


/** SD star test for PcMw */
SCIP_RETCODE reduce_sdStarPc2(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edgelimit,          /**< limit */
   const int*            edgestate,          /**< state array or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            dist,               /**< vertex distances */
   int*                  star_base,          /**< array of size nnodes */
   int*                  visitlist,          /**< array of size nnodes */
   STP_Bool*             visited,            /**< array of size nnodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   int*                  nelims              /**< point to store number of deleted edges */
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   SCIP_Real* cost_dcsr_org;
   SCIP_Real* cost_dcsr_biased;
   int* head_csr;
   int* edgeid_csr;
   SCIP_Bool* edge_deletable;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g && scip && nelims && visited && visitlist && dheap && star_base);
   assert(graph_pc_isPcMw(g) && !g->extended);
   assert(graph_isMarked(g));

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   graph_heap_clean(TRUE, dheap);
   graph_init_dcsr(scip, g);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   edgeid_csr = dcsr->edgeid;
   cost_dcsr_org = dcsr->cost;

   SCIP_CALL( SCIPallocBufferArray(scip, &edge_deletable, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cost_dcsr_biased, dcsr->nedges) );

   for( int e = 0; e < nedges / 2; e++ )
      edge_deletable[e] = FALSE;

   assert(dcsr && range_csr && edgeid_csr);

   pcBiasCostsDCSR(scip, g, FALSE, cost_dcsr_biased, dist);

   dcsr->cost = cost_dcsr_biased;

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      dist[i] = FARAWAY;
      star_base[i] = SDSTAR_BASE_UNSET;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      SCIP_Bool runloop;

      if( !g->mark[i] )
      {
         assert(g->grad[i] == 2 || g->grad[i] == 0 || i == g->source);
         continue;
      }

      runloop = TRUE;

      while( runloop )
      {
         SCIP_Bool success;
         int nvisits;
         const int start = range_csr[i].start;

         if( range_csr[i].end - start <= 1 )
            break;

         runloop = FALSE;

         /* do the actual star run */
         graph_sdStar(scip, g, TRUE, i, 2 * edgelimit, star_base, dist, visitlist, &nvisits, dheap, visited, &success);

         if( success )
         {
            int enext;

            /* check all star nodes (neighbors of i) */
            for( int e = start; e < range_csr[i].end; e = enext )
            {
               const int starnode = head_csr[e];
               const int starbase = star_base[starnode];
               assert(star_base[starnode] >= 0);
               assert(SCIPisLE(scip, dist[starnode], dcsr->cost[e]));
               assert(star_base[starnode] == starnode || star_base[starnode] >= 0);

               enext = e + 1;

               if( checkstate )
               {
                  const int orgedge = edgeid_csr[e];
                  if( edgestate[orgedge] == EDGE_BLOCKED )
                     continue;
               }

               /* shorter path to current star node found? */
               if( starnode != starbase )
               {
                  assert(star_base[starbase] != SDSTAR_BASE_UNSET);

                  /* path still valid?   */
                  if( star_base[starbase] != SDSTAR_BASE_KILLED || !EQ(dist[starnode], dcsr->cost[e]) )
                  {
                     star_base[starnode] = SDSTAR_BASE_KILLED;
                     edge_deletable[edgeid_csr[e] / 2] = TRUE;

                     dcsr->cost = cost_dcsr_org;
                     dcsr->cost2 = cost_dcsr_biased;
                     graph_dcsr_deleteEdgeBi(scip, dcsr, e);
                     dcsr->cost = cost_dcsr_biased;
                     dcsr->cost2 = NULL;

                     (*nelims)++;
                     enext--;
                  }
                  else
                  {
                     runloop = TRUE;
                  }
               }
            } /* traverse star nodes */
         } /* if success */

         sdStarReset(nnodes, nvisits, visitlist, star_base, dist, visited, dheap);
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges / 2; e++ )
   {
      if( edge_deletable[e] )
         assert(dcsr->id2csredge[e * 2] == -1);
      else if( g->oeat[e * 2] != EAT_FREE )
         assert(dcsr->id2csredge[e * 2] != -1 || !g->mark[g->tail[e * 2]] || !g->mark[g->head[e * 2]]);
   }
#endif

   graph_edge_delBlocked(scip, g, edge_deletable, TRUE);

   dcsr->cost = cost_dcsr_org;

   SCIPfreeBufferArray(scip, &cost_dcsr_biased);
   SCIPfreeBufferArray(scip, &edge_deletable);

   graph_free_dcsr(scip, g);

   return SCIP_OKAY;
}



/** SD star test for PcMw
 *  NOTE: deprecated */
SCIP_RETCODE reduce_sdStarPc(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   edgelimit,          /**< limit */
   const int*            edgestate,          /**< state array or NULL */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Real*            dist,               /**< vertex distances */
   int*                  star_base,          /**< array of size nnodes */
   int*                  visitlist,          /**< array of size nnodes */
   STP_Bool*             visited,            /**< array of size nnodes */
   DHEAP*                dheap,              /**< Dijkstra heap */
   int*                  nelims              /**< point to store number of deleted edges */
   )
{
   DCSR* dcsr;
   RANGE* range_csr;
   int* head_csr;
   int* edgeid_csr;
   int* star_base_out;
   SCIP_Real* cost_csr;
   SCIP_Real* costbiased_out;
   SCIP_Real* costbiased_in;
   SCIP_Bool* edge_deletable;
   const int nnodes = g->knots;
   const int nedges = g->edges;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g && scip && nelims && visited && visitlist && dheap && star_base);
   assert(!graph_pc_isPcMw(g) || !g->extended);

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   graph_heap_clean(TRUE, dheap);
   graph_init_dcsr(scip, g);

   dcsr = g->dcsr_storage;
   range_csr = dcsr->range;
   head_csr = dcsr->head;
   edgeid_csr = dcsr->edgeid;
   cost_csr = dcsr->cost;

   assert(dcsr && range_csr && edgeid_csr && cost_csr);

   SCIP_CALL( SCIPallocBufferArray(scip, &edge_deletable, nedges / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costbiased_out, dcsr->nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &costbiased_in, dcsr->nedges) );
   SCIP_CALL( SCIPallocBufferArray(scip, &star_base_out, nnodes) );

   for( int e = 0; e < nedges / 2; e++ )
      edge_deletable[e] = FALSE;

   pcBiasCostsDCSR(scip, g, FALSE, costbiased_out, dist);
   pcBiasCostsDCSR(scip, g, TRUE, costbiased_in, dist);

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      dist[i] = FARAWAY;
      star_base[i] = SDSTAR_BASE_UNSET;
      star_base_out[i] = SDSTAR_BASE_UNSET;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      SCIP_Bool runloop;

      if( !g->mark[i] )
      {
         assert(g->grad[i] == 2 || g->grad[i] == 0 || i == g->source);
         continue;
      }

      runloop = TRUE;

      while( runloop )
      {
         SCIP_Bool success;
         int nvisits;
         const int start = range_csr[i].start;

         if( range_csr[i].end - start <= 1 )
            break;

         runloop = FALSE;

         /* do two star runs */

         dcsr->cost = costbiased_out;
         graph_sdStar(scip, g, TRUE, i, edgelimit, star_base, dist, visitlist, &nvisits, dheap, visited, &success);

         if( success )
         {
            for( int e = start; e < range_csr[i].end; e++ )
               star_base_out[head_csr[e]] = star_base[head_csr[e]];
         }

         sdStarReset(nnodes, nvisits, visitlist, star_base, dist, visited, dheap);

         if( success )
         {
            dcsr->cost = costbiased_in;
            graph_sdStar(scip, g, TRUE, i, edgelimit, star_base, dist, visitlist, &nvisits, dheap, visited, &success);
         }

         if( success )
         {
            int* const star_base_in = star_base;
            int enext;

            dcsr->cost = cost_csr;
            dcsr->cost2 = costbiased_in;
            dcsr->cost3 = costbiased_out;

            /* check all star nodes (neighbors of i) */
            for( int e = start; e < range_csr[i].end; e = enext )
            {
               const int starnode = head_csr[e];
               const int starbase_in = star_base_in[starnode];
               const int starbase_out = star_base_out[starnode];

               assert(starbase_in >= 0 && starbase_out >= 0);

               assert(SCIPisLE(scip, dist[starnode], costbiased_in[e]));

               enext = e + 1;

               if( checkstate )
               {
                  const int orgedge = edgeid_csr[e];
                  if( edgestate[orgedge] == EDGE_BLOCKED )
                     continue;
               }

               /* shorter path to current star node found? */
               if( starnode != starbase_in && starnode != starbase_out )
               {
                  assert(star_base_in[starbase_in] != SDSTAR_BASE_UNSET);
                  assert(star_base_out[starbase_out] != SDSTAR_BASE_UNSET);

                  /* path still valid? */
                  if( star_base_in[starbase_in] != SDSTAR_BASE_KILLED && star_base_out[starbase_out] != SDSTAR_BASE_KILLED )
                  {
                     star_base_in[starnode] = SDSTAR_BASE_KILLED;
                     star_base_out[starnode] = SDSTAR_BASE_KILLED;

                     edge_deletable[edgeid_csr[e] / 2] = TRUE;
                     graph_dcsr_deleteEdgeBi(scip, dcsr, e);

                     (*nelims)++;
                     enext--;
                  }
                  else
                  {
                     runloop = TRUE;
                  }
               }
            } /* traverse star nodes */
         } /* if success */

         /* todo bit of an overkill */
         for( int k = 0; k < nvisits; k++ )
            star_base_out[visitlist[k]] = SDSTAR_BASE_UNSET;

         sdStarReset(nnodes, nvisits, visitlist, star_base, dist, visited, dheap);

#ifndef NDEBUG
         for( int k = 0; k < nnodes; k++ )
            assert(star_base_out[k] == SDSTAR_BASE_UNSET);
#endif
      }
   }

#ifndef NDEBUG
   for( int e = 0; e < nedges / 2; e++ )
   {
      if( edge_deletable[e] )
         assert(dcsr->id2csredge[e * 2] == -1);
      else if( g->oeat[e * 2] != EAT_FREE )
         assert(dcsr->id2csredge[e * 2] != -1 || !g->mark[g->tail[e * 2]] || !g->mark[g->head[e * 2]]);
   }
#endif

   graph_edge_delBlocked(scip, g, edge_deletable, TRUE);

   SCIPfreeBufferArray(scip, &star_base_out);
   SCIPfreeBufferArray(scip, &costbiased_in);
   SCIPfreeBufferArray(scip, &costbiased_out);
   SCIPfreeBufferArray(scip, &edge_deletable);

   dcsr->cost = cost_csr;

   graph_free_dcsr(scip, g);

   return SCIP_OKAY;
}

/** SD test for PcMw using only limited Dijkstra-like walk from both endpoints of an edge */
SCIP_RETCODE reduce_sdWalk(
   SCIP*                 scip,
   int                   edgelimit,
   const int*            edgestate,
   GRAPH*                g,
   int*                  termmark,
   SCIP_Real*            dist,
   int*                  heap,
   int*                  state,
   int*                  visitlist,
   STP_Bool*             visited,
   int*                  nelims
   )
{
   const int nnodes = g->knots;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g != NULL);
   assert(scip != NULL);
   assert(heap != NULL);
   assert(nelims != NULL);
   assert(visited != NULL);
   assert(visitlist != NULL);
   assert(!g->extended);
   assert(graph_pc_isPcMw(g));

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      state[i] = UNKNOWN;
      dist[i] = FARAWAY;
   }

   graph_pc_termMarkProper(g, termmark);

   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      if( !g->mark[i] )
         continue;

      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         SCIP_Bool success;
         const SCIP_Real ecost = g->cost[e];
         int nvisits;
         const int i2 = g->head[e];
         const int enext = g->oeat[e];

         if( !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         if( checkstate && edgestate[e] == EDGE_BLOCKED )
         {
            e = enext;
            continue;
         }

         success = graph_sdWalks(scip, g, g->cost, termmark, ecost, i2, i, edgelimit, dist, heap, state, visitlist, &nvisits, visited);
         sdwalkReset(nnodes, nvisits, visitlist, dist, state, visited);

         if( success )
         {
            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;
         }

         e = enext;
      }
   }

   return SCIP_OKAY;
}


/** SD test for PcMw using only limited Dijkstra-like walk from both endpoints of an edge */
SCIP_RETCODE reduce_sdWalkExt(
   SCIP*                 scip,
   int                   edgelimit,
   const int*            edgestate,
   GRAPH*                g,
   SCIP_Real*            dist,
   int*                  heap,
   int*                  state,
   int*                  visitlist,
   STP_Bool*             visited,
   int*                  nelims
   )
{
   int* prevterms;
   int* nprevterms;
   const int nnodes = g->knots;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g != NULL);
   assert(scip != NULL);
   assert(heap != NULL);
   assert(nelims != NULL);
   assert(visited != NULL);
   assert(visitlist != NULL);
   assert(!g->extended);
   assert(graph_pc_isPcMw(g));

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &prevterms, nnodes * STP_SDWALK_MAXNPREVS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nprevterms, nnodes) );

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      state[i] = UNKNOWN;
      dist[i] = FARAWAY;
      nprevterms[i] = 0;
   }

   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      if( !g->mark[i] )
         continue;

      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         SCIP_Bool success;
         const SCIP_Real ecost = g->cost[e];
         int nvisits;
         const int i2 = g->head[e];
         const int enext = g->oeat[e];

         /* avoid double checking */
         if( !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         if( checkstate && edgestate[e] == EDGE_BLOCKED )
         {
            e = enext;
            continue;
         }

         success = graph_sdWalksExt(scip, g, g->cost, ecost, i2, i, edgelimit, STP_SDWALK_MAXNPREVS, dist, prevterms, nprevterms, heap, state, visitlist, &nvisits, visited);
         sdwalkResetExt(nnodes, nvisits, visitlist, dist, nprevterms, state, visited);

         if( success )
         {
            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;
         }

         e = enext;
      }
   }

   SCIPfreeBufferArray(scip, &nprevterms);
   SCIPfreeBufferArray(scip, &prevterms);

   return SCIP_OKAY;
}

/** SD test for PcMw using only limited Dijkstra-like walk from both endpoints of an edge */
SCIP_RETCODE reduce_sdWalkExt2(
   SCIP*                 scip,
   int                   edgelimit,
   const int*            edgestate,
   GRAPH*                g,
   int*                  termmark,
   SCIP_Real*            dist,
   int*                  heap,
   int*                  state,
   int*                  visitlist,
   STP_Bool*             visited,
   int*                  nelims
   )
{
   int* prevterms;
   int* nprevterms;
   int* prevNPterms;
   int* nprevNPterms;
   int* prevedges;
   int* nprevedges;
   const int nnodes = g->knots;
   const SCIP_Bool checkstate = (edgestate != NULL);

   assert(g != NULL);
   assert(scip != NULL);
   assert(heap != NULL);
   assert(nelims != NULL);
   assert(visited != NULL);
   assert(visitlist != NULL);
   assert(!g->extended);
   assert(graph_pc_isPcMw(g));

   if( edgelimit <= 0 )
      return SCIP_OKAY;

   assert(0 && "triggers bug in STP-DIMACS/PCSPG-hand/HAND_SMALL_ICERM/handsi04.stp");

   SCIP_CALL( SCIPallocBufferArray(scip, &prevterms, nnodes * STP_SDWALK_MAXNPREVS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nprevterms, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prevNPterms, nnodes * STP_SDWALK_MAXNPREVS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nprevNPterms, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prevedges, nnodes * STP_SDWALK_MAXNPREVS) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nprevedges, nnodes) );

   for( int i = 0; i < nnodes; i++ )
   {
      visited[i] = FALSE;
      state[i] = UNKNOWN;
      dist[i] = FARAWAY;
      nprevterms[i] = 0;
      nprevNPterms[i] = 0;
      nprevedges[i] = 0;
   }

   graph_pc_termMarkProper(g, termmark);

   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      if( !g->mark[i] )
         continue;

      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         SCIP_Bool success;
         const SCIP_Real ecost = g->cost[e];
         int nvisits;
         const int i2 = g->head[e];
         const int enext = g->oeat[e];

         /* avoid double checking */
         if( !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         if( checkstate && edgestate[e] == EDGE_BLOCKED )
         {
            e = enext;
            continue;
         }

         success = graph_sdWalksExt2(scip, g, g->cost, termmark, ecost, i2, i, edgelimit, STP_SDWALK_MAXNPREVS, dist, prevterms, nprevterms,
               prevNPterms, nprevNPterms, prevedges, nprevedges, heap, state, visitlist, &nvisits, visited);
         sdwalkResetExt2(nnodes, nvisits, visitlist, dist, nprevterms, nprevNPterms, nprevedges, state, visited);

         if( success )
         {
            graph_edge_del(scip, g, e, TRUE);
            (*nelims)++;
         }

         e = enext;
      }
   }

   SCIPfreeBufferArray(scip, &nprevedges);
   SCIPfreeBufferArray(scip, &prevedges);
   SCIPfreeBufferArray(scip, &nprevNPterms);
   SCIPfreeBufferArray(scip, &prevNPterms);
   SCIPfreeBufferArray(scip, &nprevterms);
   SCIPfreeBufferArray(scip, &prevterms);

   return SCIP_OKAY;
}

/** SD test using only limited dijkstra from both endpoints of an edge */
SCIP_RETCODE reduce_sdsp(
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
   int                   limit,
   int*                  edgestate
)
{
   int* pathmaxnodetail = NULL;
   int* pathmaxnodehead = NULL;
   const int nnodes = g->knots;
   const SCIP_Bool pc = (g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG);
   const SCIP_Bool checkstate = (edgestate != NULL);

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
   assert(!g->extended);

   *nelims = 0;

   if( limit <= 0 )
      return SCIP_OKAY;

   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodetail, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodehead, nnodes) );

      for( int i = 0; i < nnodes; i++ )
      {
         pathmaxnodetail[i] = -1;
         pathmaxnodehead[i] = -1;
      }
   }
   else
   {
      for( int i = 0; i < nnodes; i++ )
         g->mark[i] = (g->grad[i] > 0);
   }

   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
   }

   /* iterate through all nodes */
   for( int i = 0; i < nnodes; i++ )
   {
      int e;
      if( !g->mark[i] )
         continue;

      /* traverse neighbours */
      e = g->outbeg[i];
      while( e != EAT_LAST )
      {
         const SCIP_Real ecost = g->cost[e];
         const int i2 = g->head[e];
         const int enext = g->oeat[e];
         int nlbltail;
         int nlblhead;
         SCIP_Bool deletable;

         /* avoid double checking */
         if( i2 < i || !g->mark[i2] )
         {
            e = enext;
            continue;
         }

         /* execute limited Dijkstra from both sides */

         if( pc )
         {
            graph_path_PcMwSd(scip, g, pathtail, g->cost, ecost, pathmaxnodetail, heap, statetail, NULL, memlbltail, &nlbltail, i, i2, limit);
            graph_path_PcMwSd(scip, g, pathhead, g->cost, ecost, pathmaxnodehead, heap, statehead, statetail, memlblhead, &nlblhead, i2, i, limit);
         }
         else
         {
            graph_sdPaths(g, pathtail, g->cost, ecost, heap, statetail, memlbltail, &nlbltail, i, i2, limit);
            graph_sdPaths(g, pathhead, g->cost, ecost, heap, statehead, memlblhead, &nlblhead, i2, i, limit);
         }

         deletable = FALSE;

         /* check whether edge e can be deleted and restore data structures */
         for( int k = 0; k < nlbltail && !deletable; k++ )
         {
            const int l = memlbltail[k];

            assert(g->mark[l]);
            assert(statetail[l] != UNKNOWN);

            if( statehead[l] != UNKNOWN )
            {
               assert(SCIPisGT(scip, FARAWAY, pathtail[l].dist));
               assert(SCIPisGT(scip, FARAWAY, pathhead[l].dist));

               if( pc )
               {
                  const int tailmaxterm = pathmaxnodetail[l];
                  const int headmaxterm = pathmaxnodehead[l];

                  assert(tailmaxterm != i && headmaxterm != i);
                  assert(tailmaxterm != i2 && headmaxterm != i2);

                  /* any terminal on the path? */
                  if( tailmaxterm >= 0 || headmaxterm >= 0 )
                  {
                     if( tailmaxterm == headmaxterm )
                     {
                        assert(tailmaxterm == l);
                        assert(SCIPisPositive(scip, g->prize[tailmaxterm]));
                        assert(SCIPisGE(scip, ecost, pathhead[headmaxterm].dist) && SCIPisGE(scip, ecost, pathtail[tailmaxterm].dist));

                        if( SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("delete1Term \n");
                        }
                     }
                     else if( tailmaxterm >= 0 && headmaxterm >= 0 )
                     {
                        const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;
                        const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;

                        assert(tailmaxterm != headmaxterm);
                        assert(!SCIPisNegative(scip, distl2tailmax));
                        assert(!SCIPisNegative(scip, distl2headmax));
                        assert(SCIPisPositive(scip, g->prize[tailmaxterm]) && SCIPisPositive(scip, g->prize[headmaxterm]));
                        assert(SCIPisGE(scip, ecost, pathhead[headmaxterm].dist) && SCIPisGE(scip, ecost, pathtail[tailmaxterm].dist));

                        if( SCIPisGE(scip, ecost, distl2tailmax + distl2headmax)
                              && SCIPisGE(scip, ecost, distl2tailmax + pathhead[l].dist - g->prize[headmaxterm])
                              && SCIPisGE(scip, ecost, distl2headmax + pathtail[l].dist - g->prize[tailmaxterm])
                              && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm] - g->prize[headmaxterm]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("delete2Term \n");
                        }
                     }
                     else if( tailmaxterm >= 0 )
                     {
                        const SCIP_Real distl2tailmax = pathtail[l].dist - pathtail[tailmaxterm].dist;
                        // todo consider l == term?
                        assert(headmaxterm < 0);
                        assert(SCIPisGE(scip, ecost, pathtail[tailmaxterm].dist));
                        assert(SCIPisPositive(scip, g->prize[tailmaxterm]));

                        if( SCIPisGE(scip, ecost, distl2tailmax + pathhead[l].dist)
                              && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[tailmaxterm]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("deleteHalfTerm1 \n");
                        }
                     }
                     else if( headmaxterm >= 0 )
                     {
                        const SCIP_Real distl2headmax = pathhead[l].dist - pathhead[headmaxterm].dist;
                        // todo consider l == term?
                        assert(tailmaxterm < 0);
                        assert(SCIPisGE(scip, ecost, pathhead[headmaxterm].dist));
                        assert(SCIPisPositive(scip, g->prize[headmaxterm]));

                        if( SCIPisGE(scip, ecost, distl2headmax + pathtail[l].dist)
                              && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[headmaxterm]) )
                        {
                           deletable = TRUE;
                           SCIPdebugMessage("deleteHalfTerm2 \n");
                        }
                     }
                  }
                  else if( SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist) )
                  {
                     deletable = TRUE;
                  }

                  if( Is_term(g->term[l]) )
                  {
                     if( SCIPisGE(scip, ecost, pathhead[l].dist) && SCIPisGE(scip, ecost, pathtail[l].dist)
                           && SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                        deletable = TRUE;
                  }
               }
               else
               {
                  if( Is_term(g->term[l]) )
                  {
                     if( SCIPisGE(scip, ecost, pathhead[l].dist) && SCIPisGE(scip, ecost, pathtail[l].dist) )
                     {
                        deletable = TRUE;
#if 0
                        if( pc && SCIPisLT(scip, ecost, pathhead[l].dist + pathtail[l].dist - g->prize[l]) )
                           deletable = FALSE;
#endif
                     }
                  }
                  else
                  {
                     if( SCIPisGE(scip, ecost, pathhead[l].dist + pathtail[l].dist) )
                        deletable = TRUE;
                  }
               }
            }
         }

         /* restore data */

         for( int k = 0; k < nlbltail; k++ )
         {
            const int l = memlbltail[k];
            statetail[l] = UNKNOWN;
            pathtail[l].dist = FARAWAY;
            pathtail[l].edge = UNKNOWN;
            if( pc )
               pathmaxnodetail[l] = -1;
         }

         for( int k = 0; k < nlblhead; k++ )
         {
            const int l = memlblhead[k];
            statehead[l]     = UNKNOWN;
            pathhead[l].dist = FARAWAY;
            pathhead[l].edge = UNKNOWN;
            if( pc )
               pathmaxnodehead[l] = -1;
         }

#ifndef NDEBUG
         for( int k = 0; k < nnodes; k++ )
         {
            assert(statetail[k]     == UNKNOWN);
            assert(pathtail[k].dist == FARAWAY);
            assert(pathtail[k].edge == UNKNOWN);
            assert(statehead[k]     == UNKNOWN);
            assert(pathhead[k].dist == FARAWAY);
            assert(pathhead[k].edge == UNKNOWN);
            if( pc )
            {
               assert(pathmaxnodetail[k] == -1);
               assert(pathmaxnodehead[k] == -1);
            }
         }
#endif
         /* can edge be deleted? */
         if( deletable )
         {
            if( !checkstate || (edgestate[e] != EDGE_BLOCKED) )
            {
               graph_edge_del(scip, g, e, TRUE);
               (*nelims)++;
            }
         }

         e = enext;
      }
   }

   SCIPfreeBufferArrayNull(scip, &pathmaxnodehead);
   SCIPfreeBufferArrayNull(scip, &pathmaxnodetail);

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}


/** bd_k test for given Steiner bottleneck distances */
/* *** DEPRECATED ***/
SCIP_RETCODE reduce_bd34WithSd(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   GRAPH*                netgraph,           /**< auxiliary graph structure */
   PATH*                 netmst,             /**< MST structure */
   PATH*                 vnoi,               /**< path structure */
   SCIP_Real*            mstsdist,           /**< MST distance in aux-graph */
   int*                  vbase,              /**< bases for nearest terminals */
   int*                  nodesid,            /**< nodes identification array */
   int*                  nelims              /**< number of eliminations */
   )
{
   SCIP_Real cutoffs[STP_BD_MAXDNEDGES];
   SCIP_Real sd[STP_BD_MAXDEGREE];
   SCIP_Real ecost[STP_BD_MAXDEGREE];
   int edges[STP_BD_MAXDEGREE];
   int adjvert[STP_BD_MAXDEGREE];
   GRAPH* auxg;
   const int nnodes = g->knots;

   assert(scip && g && netgraph && netmst && vnoi);
   assert(!graph_pc_isPcMw(g));

   /* build auxiliary graph */
   SCIP_CALL( graph_buildCompleteGraph(scip, &auxg, STP_BD_MAXDEGREE) );
   assert(auxg->edges == 2 * STP_BD_MAXDNEDGES);

   SCIP_CALL( graph_path_init(scip, auxg) );

   SCIPdebugMessage("BD34-SD Reduction: ");

   for( int i = 0; i < STP_BD_MAXDEGREE; i++ )
      sd[i] = 0.0;

   graph_mark(g);

   for( int degree = 3; degree <= STP_BD_MAXDEGREE; degree ++ )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         if( Is_term(g->term[i]) || g->grad[i] != degree )
         {
            continue;
         }
         else
         {
            int k = 0;
            for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
            {
               edges[k] = e;
               ecost[k] = g->cost[e];
               adjvert[k++] = g->head[e];
            }
            assert(k == degree);
         }

         assert(g->mark[i]);

         /* vertex of degree 3? */
         if( degree == 3 )
         {
            const SCIP_Real costsum = ecost[0] + ecost[1] + ecost[2];

            sd[0] = getSd(scip, g, netgraph, netmst, vnoi, mstsdist, ecost[0] + ecost[1], vbase, nodesid, adjvert[0],
                  adjvert[1], 300);
            sd[1] = getSd(scip, g, netgraph, netmst, vnoi, mstsdist, ecost[1] + ecost[2], vbase, nodesid, adjvert[1],
                  adjvert[2], 300);
            sd[2] = getSd(scip, g, netgraph, netmst, vnoi, mstsdist, ecost[2] + ecost[0], vbase, nodesid, adjvert[2],
                  adjvert[0], 300);

            if( isPseudoDeletableDeg3(scip, g, sd, edges, costsum, TRUE) )
            {
               SCIP_Bool success;

               cutoffs[0] = sd[0];
               cutoffs[1] = sd[2];
               cutoffs[2] = sd[1];

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));

               assert(success);
               assert(g->grad[i] == 0);

               SCIPdebugMessage("BD3-R Reduction: %f %f %f csum: %f\n ", sd[0], sd[1], sd[2], costsum);
               (*nelims)++;
            }
         }
         /* vertex of degree 4? */
         else if( degree == 4 )
         {
            SCIP_Bool success = TRUE;

            for( int k = 0; k < 4; k++ )
            {
               auxg->mark[k] = TRUE;
               for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  const int k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     auxg->cost[e] = getSd(scip, g, netgraph, netmst, vnoi, mstsdist, ecost[k] + ecost[k2], vbase, nodesid,
                           adjvert[k], adjvert[k2], 200);
                     auxg->cost[flipedge(e)] = auxg->cost[e];
                  }
               }
            }

            success = isPseudoDeletable(scip, g, auxg, ecost, edges, 4);

            if( success )
            {
               int edgecount = 0;
               for( int k = 0; k < 3; k++ )
               {
                  for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
                  {
                     const int k2 = auxg->head[e];
                     if( k2 > k )
                        cutoffs[edgecount++] = auxg->cost[e];
                  }
               }

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));

               if( success )
                  (*nelims)++;
            }
         }
      }
   }

   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);

   return SCIP_OKAY;
}


/* C. W. Duin and A. Volganant
 *
 * "Reduction Tests for the Steiner Problem in Graphs"
 *
 * Networks, Volume 19 (1989), 549-567
 *
 * Bottleneck Degree 3,4 Test
 */
SCIP_RETCODE reduce_bd34(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   PATH*                 pathtail,           /**< array for internal use */
   PATH*                 pathhead,           /**< array for internal use */
   int*                  heap,               /**< array for internal use */
   int*                  statetail,          /**< array for internal use */
   int*                  statehead,          /**< array for internal use */
   int*                  memlbltail,         /**< array for internal use */
   int*                  memlblhead,         /**< array for internal use */
   int*                  nelims,             /**< point to return number of eliminations */
   int                   limit,              /**< limit for edges to consider for each vertex */
   SCIP_Real*            offset              /**< offset */
   )
{
   SCIP_Real cutoffs[STP_BD_MAXDNEDGES];
   int edges[STP_BD_MAXDEGREE];
   SCIP_Real sd[STP_BD_MAXDEGREE];
   SCIP_Real ecost[STP_BD_MAXDEGREE];
   GRAPH* auxg;
   int* pathmaxnodetail = NULL;
   int* pathmaxnodehead = NULL;
   int adjvert[STP_BD_MAXDEGREE];
   const int nnodes = g->knots;
   const SCIP_Bool pc = g->stp_type == STP_PCSPG || g->stp_type == STP_RPCSPG;

   SCIPdebugMessage("BD34-Reduction: ");

   assert(scip && g && heap && nelims);
   assert(!g->extended);

   /* build auxiliary graph */
   SCIP_CALL( graph_buildCompleteGraph(scip, &auxg, STP_BD_MAXDEGREE) );
   assert(auxg->edges == 2 * STP_BD_MAXDNEDGES);

   SCIP_CALL( graph_path_init(scip, auxg) );

   *nelims = 0;

   for( int i = 0; i < STP_BD_MAXDEGREE; i++ )
      sd[i] = 0.0;

   if( pc )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodetail, nnodes) );
      SCIP_CALL( SCIPallocBufferArray(scip, &pathmaxnodehead, nnodes) );
   }

   graph_mark(g);

   for( int i = 0; i < nnodes; i++ )
   {
      statetail[i]     = UNKNOWN;
      pathtail[i].dist = FARAWAY;
      pathtail[i].edge = UNKNOWN;
      statehead[i]     = UNKNOWN;
      pathhead[i].dist = FARAWAY;
      pathhead[i].edge = UNKNOWN;
      if( pc )
      {
         pathmaxnodetail[i] = -1;
         pathmaxnodehead[i] = -1;
      }
   }

   for( int degree = 3; degree <= STP_BD_MAXDEGREE; degree++ )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         int edgecount;

         if( pc )
         {
            int pcdegree;

            if( !g->mark[i] || graph_pc_knotIsFixedTerm(g, i) )
               continue;

            if( Is_term(g->term[i]) && !graph_pc_termIsNonLeafTerm(g, i) )
               continue;

            pcdegree = graph_pc_realDegree(g, i, FALSE);

            if( pcdegree != degree )
               continue;
         }
         else
         {
            if( Is_term(g->term[i]) || g->grad[i] != degree )
               continue;
         }

         edgecount = 0;
         for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];

            if( g->mark[head] )
            {
               edges[edgecount] = e;
               ecost[edgecount] = g->cost[e];
               adjvert[edgecount++] = g->head[e];
            }
            else
            {
               assert(pc && Is_term(g->term[i]));
            }
         }

         assert(edgecount == degree);

         /* vertex of degree 3? */
         if( degree == 3 )
         {
            const SCIP_Real iprize = (pc && Is_term(g->term[i])) ? g->prize[i] : 0.0;
            const SCIP_Real csum = ecost[0] + ecost[1] + ecost[2] - iprize;

            assert(edgecount == 3);
            assert(iprize <= ecost[0] && iprize <= ecost[1] && iprize <= ecost[2]);

            if( pc )
            {
               SCIP_CALL(
                     reduce_getSdPcMw(scip, g, pathtail, pathhead, &(sd[0]), csum, heap, statetail, statehead, memlbltail, memlblhead,
                           pathmaxnodetail, pathmaxnodehead, adjvert[0], adjvert[1], limit));
               SCIP_CALL(
                     reduce_getSdPcMw(scip, g, pathtail, pathhead, &(sd[1]), csum, heap, statetail, statehead, memlbltail, memlblhead,
                           pathmaxnodetail, pathmaxnodehead, adjvert[1], adjvert[2], limit));
               SCIP_CALL(
                     reduce_getSdPcMw(scip, g, pathtail, pathhead, &(sd[2]), csum, heap, statetail, statehead, memlbltail, memlblhead,
                           pathmaxnodetail, pathmaxnodehead, adjvert[2], adjvert[0], limit));
            }
            else
            {
               SCIP_CALL(
                     reduce_getSd(scip, g, pathtail, pathhead, &(sd[0]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[0], adjvert[1], limit, pc, FALSE));
               SCIP_CALL(
                     reduce_getSd(scip, g, pathtail, pathhead, &(sd[1]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[1], adjvert[2], limit, pc, FALSE));
               SCIP_CALL(
                     reduce_getSd(scip, g, pathtail, pathhead, &(sd[2]), csum, heap, statetail, statehead, memlbltail, memlblhead, adjvert[2], adjvert[0], limit, pc, FALSE));
            }

            if( isPseudoDeletableDeg3(scip, g, sd, edges, csum, !Is_term(g->term[i])) )
            {
               SCIP_Bool success;

               cutoffs[0] = sd[0];
               cutoffs[1] = sd[2];
               cutoffs[2] = sd[1];

               if( Is_term(g->term[i]) )
               {
                  assert(offset != NULL);
                  assert(graph_pc_isPcMw(g));

                  *offset += g->prize[i];
               }

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));
               assert(success);
               assert(g->grad[i] == 0);

               SCIPdebugMessage("BD3 Reduction: %f %f %f csum: %f\n ", sd[0], sd[1], sd[2], csum);
               (*nelims)++;
            }
         }
         /* vertex of degree 4? */
         else if( degree == 4 )
         {
            SCIP_Bool success = TRUE;
            const SCIP_Real csum = ecost[0] + ecost[1] + ecost[2] + ecost[3];

            if( Is_term(g->term[i]) )
            {
               assert(pc);
               continue;
            }

            assert(edgecount == 4);

            for( int k = 0; k < 4; k++ )
            {
               auxg->mark[k] = TRUE;
               for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
               {
                  const int k2 = auxg->head[e];
                  if( k2 > k )
                  {
                     SCIP_Real s1 = -1.0;
                     if( pc )
                        SCIP_CALL( reduce_getSdPcMw(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                                    pathmaxnodetail, pathmaxnodehead, adjvert[k], adjvert[k2], limit));
                     else
                        SCIP_CALL( reduce_getSd(scip, g, pathtail, pathhead, &(s1), csum, heap, statetail, statehead, memlbltail, memlblhead,
                              adjvert[k], adjvert[k2], limit, pc, FALSE));
                     assert(s1 >= 0);
                     auxg->cost[e] = s1;
                     auxg->cost[flipedge(e)] = s1;
                  }
               }
            }

            success = isPseudoDeletable(scip, g, auxg, ecost, edges, 4);

            if( success )
            {
               edgecount = 0;

               for( int k = 0; k < 3; k++ )
               {
                  for( int e = auxg->outbeg[k]; e != EAT_LAST; e = auxg->oeat[e] )
                  {
                     const int k2 = auxg->head[e];
                     if( k2 > k )
                        cutoffs[edgecount++] = auxg->cost[e];
                  }
               }

               SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, i, &success));

               if( success )
                  (*nelims)++;
            }
         }
      }
   }

   /* free memory */

   SCIPfreeBufferArrayNull(scip, &pathmaxnodehead);
   SCIPfreeBufferArrayNull(scip, &pathmaxnodetail);

   graph_path_exit(scip, auxg);
   graph_free(scip, &auxg, TRUE);

   SCIPdebugMessage("bd34: %d nodes deleted\n", *nelims);

   assert(graph_valid(scip, g));

   return SCIP_OKAY;
}
