/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reduce_sdutil.c
 * @brief  utility methods for Steiner tree reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem special distance reduction techniques.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG
#include "reduce.h"
#include "portab.h"

#define STP_SDN_NCLOSETERMS 4
#define STP_SDN_MAXDEGREE 4
#define STP_SDN_MAXNVISITS 4


/** see reduce.h */
struct special_distance_neighbors
{
   SCIP_Real*            closeterms_distN;   /**< paths to N closest terms */
   int*                  closetermsN;        /**< close terminals */
   SCIP_Bool*            nodes_isBlocked;    /**< node is blocked */
   /** temporary arrays start: */
   DIJK*                 dijkdata;
   int*                  nodes_nhits;
   SCIP_Real*            nodes_maxdist;
   int*                  hitlist;
   /** temporary arrays end */
   int                   nnodesHit;
   int                   N;                  /**< the number of closest terms per node */
};


/*
 * local methods
 */



/** inserts new edge */
static inline
void sdgraphInsertEdge(
   SCIP*                  scip,              /**< SCIP data structure */
   int                    term1,             /**< end node 1 */
   int                    term2,             /**< end node 2 */
   SCIP_Real              edgecost,          /**< cost */
   SDGRAPH*               sdgraph            /**< the SD graph */
)
{
   SCIP_Bool success;
   const int* nodesid = reduce_sdgraphGetOrgnodesToSdMap(sdgraph);
   const int sdnode1 = nodesid[term1];
   const int sdnode2 = nodesid[term2];

   reduce_sdgraphInsertEdge(scip, sdnode1, sdnode2, edgecost, -1, NULL, sdgraph, &success);
}


/** allocates memory */
static
SCIP_RETCODE sdprofitAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SCIP_Bool             useSecond,          /**< use second profit? */
   SDPROFIT**            sdprofit            /**< SD profit */
)
{
   const int nnodes = graph_get_nNodes(g);
   SDPROFIT* sdp;

   SCIP_CALL( SCIPallocMemory(scip, sdprofit) );
   sdp = *sdprofit;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_bias), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_biassource), nnodes) );

   if( useSecond )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_bias2), nnodes) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_biassource2), nnodes) );
   }
   else
   {
      sdp->nodes_bias2 = NULL;
      sdp->nodes_biassource2 = NULL;
   }

   return SCIP_OKAY;
}


/** gets SD node bias */
static
void sdprofitBuild(
   const GRAPH*          g,                  /**< the graph */
   SDPROFIT*             sdprofit            /**< SD profit */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int pseudoroot = (graph_pc_isUnrootedPcMw(g)) ? g->source : -1;
   const SCIP_Real* const costs = g->cost;
   SCIP_Real* RESTRICT node_bias = sdprofit->nodes_bias;
   int* RESTRICT node_biassource = sdprofit->nodes_biassource;
   SCIP_Real* RESTRICT node_bias2 = sdprofit->nodes_bias2;
   int* RESTRICT node_biassource2 = sdprofit->nodes_biassource2;
   const SCIP_Bool isPcMw = graph_pc_isPcMw(g);

   assert(graph_pc_isPcMw(g) || !g->extended);
   assert(node_bias2 && node_biassource2);

   for( int k = 0; k < nnodes; k++ )
   {
      node_bias[k] = 0.0;
      node_biassource[k] = k;
      node_bias2[k] = 0.0;
      node_biassource2[k] = k;
   }

   /* compute profits for non-terminals */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && k != pseudoroot )
      {
         int minneighbor = -1;
         SCIP_Real mincost = isPcMw ? g->prize[k] : FARAWAY;
         SCIP_Real mincost2 = mincost;

         for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         {
            const int neighbor = g->tail[e];

            if( neighbor == pseudoroot )
               continue;

            if( costs[e] < mincost )
            {
               assert(!graph_pc_isPcMw(g) || !graph_pc_knotIsDummyTerm(g, neighbor));

               minneighbor = neighbor;
               mincost2 = mincost;
               mincost = costs[e];
            }
            else if( costs[e] < mincost2 )
            {
               mincost2 = costs[e];
            }
         }

         assert(minneighbor >= 0 || isPcMw || g->terms == 1);

         if( minneighbor >= 0 )
         {
            const SCIP_Real bias = mincost2 - mincost;

            assert(GE(bias, 0.0));

            if( bias > node_bias[minneighbor] )
            {
               node_bias[minneighbor] = bias;
               node_biassource[minneighbor] = k;
            }
            else if( GT(bias, node_bias2[minneighbor]) )
            {
               node_bias2[minneighbor] = bias;
               node_biassource2[minneighbor] = k;

               assert(node_biassource[minneighbor] != node_biassource2[minneighbor]);
            }
         }
      }
   }

   /* correct profits for terminals */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
      {
         continue;
      }

      if( isPcMw )
      {
         if( g->prize[k] >= node_bias[k] )
         {
            node_bias[k] = g->prize[k];
            node_biassource[k] = k;
         }
         else if( g->prize[k] >= node_bias2[k]  )
         {
            node_bias2[k] = g->prize[k];
            node_biassource2[k] = k;
         }

         continue;
      }

      node_bias[k] = FARAWAY;
      node_biassource[k] = k;
      node_bias2[k] = FARAWAY;
      node_biassource2[k] = k;
   }

#ifndef NDEBUG
   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
      {
         assert(GE(node_bias[k], node_bias2[k]));
         assert(node_biassource[k] != node_biassource2[k] || node_biassource[k] == k);
      }
   }
#endif
}


/** gets SD node bias */
static
void sdprofitBuild1stOnly(
   const GRAPH*          g,                  /**< the graph */
   const SCIP_Real*      edgecosts,          /**< edge costs (w.r.t graph 'g') */
   SDPROFIT*             sdprofit            /**< SD profit */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int pseudoroot = (graph_pc_isUnrootedPcMw(g)) ? g->source : -1;
   SCIP_Real* RESTRICT node_bias = sdprofit->nodes_bias;
   int* RESTRICT node_biassource = sdprofit->nodes_biassource;
   const SCIP_Bool isPcMw = graph_pc_isPcMw(g);

   assert(graph_pc_isPcMw(g) || !g->extended);
   assert(edgecosts);

   for( int k = 0; k < nnodes; k++ )
   {
      node_bias[k] = 0.0;
      node_biassource[k] = k;
   }

   /* compute profits for non-terminals */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) && k != pseudoroot )
      {
         int minneighbor = -1;
         SCIP_Real mincost = isPcMw ? g->prize[k] : FARAWAY;
         SCIP_Real mincost2 = mincost;

         for( int e = g->inpbeg[k]; e != EAT_LAST; e = g->ieat[e] )
         {
            const int neighbor = g->tail[e];

            if( neighbor == pseudoroot )
               continue;

            if( edgecosts[e] < mincost )
            {
               assert(!graph_pc_isPcMw(g) || !graph_pc_knotIsDummyTerm(g, neighbor));

               minneighbor = neighbor;
               mincost2 = mincost;
               mincost = edgecosts[e];
            }
            else if( edgecosts[e] < mincost2 )
            {
               mincost2 = edgecosts[e];
            }
         }

         assert(minneighbor >= 0 || isPcMw || g->terms == 1);

         if( minneighbor >= 0 )
         {
            const SCIP_Real bias = mincost2 - mincost;

            assert(GE(bias, 0.0));

            if( GT(bias, node_bias[minneighbor]) )
            {
               node_bias[minneighbor] = bias;
               node_biassource[minneighbor] = k;
            }
         }
      }
   }

   /* correct profits for terminals */
   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
      {
         continue;
      }

      if( isPcMw )
      {
         if( g->prize[k] >= node_bias[k] )
         {
            node_bias[k] = g->prize[k];
            node_biassource[k] = k;
         }

         continue;
      }

      node_bias[k] = FARAWAY;
      node_biassource[k] = k;
   }
}


/** updates implied profits by using exact bottleneck distances */
static inline
void sdprofitUpdateNode(
   const GRAPH*          g,                  /**< graph */
   int                   node,               /**< the node to update */
   int                   sourceterm,         /**< the source terminal */
   SCIP_Real             edgecost,           /**< edge cost */
   SCIP_Real             bdist,              /**< bottleneck distance */
   SCIP_Bool             blctree_isUpdated,  /**< BLC tree fresh? */
   SDPROFIT*             sdprofit            /**< the SD profit */
)
{
   const SCIP_Real profit = bdist - edgecost;

   assert(GT(profit, 0.0));
   assert(GE(edgecost, 0.0) && LE(edgecost, FARAWAY));
   assert(GE(bdist, 0.0) && LE(bdist, FARAWAY));
   assert(Is_term(g->term[sourceterm]));

#ifndef NDEBUG
   if( blctree_isUpdated &&  sourceterm == sdprofit->nodes_biassource[node] )
   {
      assert(GE(profit, sdprofit->nodes_bias[node]));
   }
   else if( blctree_isUpdated && sourceterm == sdprofit->nodes_biassource2[node] )
   {
      assert(GE(profit, sdprofit->nodes_bias2[node]));
   }
#endif

   if( GT(profit, sdprofit->nodes_bias[node]) )
   {
      sdprofit->nodes_bias[node] = profit;
      sdprofit->nodes_biassource[node] = sourceterm;
   }
   else if( GT(profit, sdprofit->nodes_bias2[node]) )
   {
      sdprofit->nodes_bias2[node] = profit;
      sdprofit->nodes_biassource2[node] = sourceterm;
   }
   else if( sdprofit->nodes_biassource[node] == sdprofit->nodes_biassource2[node] )
   {
      sdprofit->nodes_bias2[node] = profit;
      sdprofit->nodes_biassource2[node] = sourceterm;
   }
}

#if 0
/** marks terminals that can be reached from given source node */
static inline
void sdneighborMarkCloseTerms(
   const GRAPH*          g,                  /**< graph to initialize from */
   int                   sourcenode,
   int* RESTRICT         nodes_nhits,
   SCIP_Real* RESTRICT   nodes_maxdist,
   SD*                   sddata              /**< SD */
)
{
   int nnterms1;
   SCIP_Real termdist1[4];
   int neighbterms1[4];

   graph_tpathsGet4CloseTerms(g, sddata->terminalpaths, sourcenode, FARAWAY, neighbterms1, NULL, termdist1, &nnterms1);

   /* go over all close terminals of the source node  */
   for( int k = 0; k < nnterms1; k++ )
   {
      const int neighborterm = neighbterms1[k];
    //  const SCIP_Real dist = MAX(ecost, termdist1[k]); // todo also try without!
      const SCIP_Real dist = termdist1[k];
      assert(Is_term(g->term[neighborterm]));

      if( nodes_nhits[neighborterm] == 0 )
      {
         // todo save neighborterm in list
      }

      nodes_nhits[neighborterm]++;

      if( dist > nodes_maxdist[neighborterm] )
         nodes_maxdist[neighborterm] = dist;
   }
}
#endif


/** marks nodes that can be reached from given source node */
static inline
SCIP_RETCODE sdneighborMarkCloseNodes(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph to initialize from */
   int                   sourcenode,         /**< (neighbor) node to mark from */
   SD*                   sddata              /**< SD */
)
{
   const SDPROFIT* sdprofit = sddata->sdprofit;
   SDN* const sdneighbors = sddata->sdneighbors;
   DIJK* RESTRICT dijkdata = sdneighbors->dijkdata;
   int* RESTRICT nodes_nhits = sdneighbors->nodes_nhits;
   SCIP_Real* RESTRICT nodes_maxdist = sdneighbors->nodes_maxdist;
   const int* const visitlist = dijkdata->visitlist;
   SCIP_Real* RESTRICT distance = dijkdata->node_distance;
   int* RESTRICT hitlist = sdneighbors->hitlist;
   int nvisits;

   assert(dijkdata && nodes_nhits && nodes_maxdist && hitlist);

   SCIP_CALL( graph_sdCloseNodesBiased(scip, g, sdprofit, sourcenode, dijkdata) );

   nvisits = dijkdata->nvisits;
   assert(nvisits >= 0);
   assert(sdneighbors->nnodesHit >= 0);

   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      const SCIP_Real dist = distance[node];

      if( node == sourcenode )
         continue;

      assert(GT(dist, 0.0));

      if( nodes_nhits[node] == 0 )
      {
         hitlist[sdneighbors->nnodesHit++] = node;
         assert(EQ(nodes_maxdist[node], 0.0));
      }

      nodes_nhits[node]++;

      if( dist > nodes_maxdist[node] )
         nodes_maxdist[node] = dist;
   }

   graph_dijkLimited_reset(g, dijkdata);

   return SCIP_OKAY;
}


/** helper */
static
SCIP_RETCODE sdneighborUpdateInit(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph  */
   SDN*                  sdn                 /**< SD neighbors */
)
{
   const int N = sdn->N;
   SCIP_Real* const RESTRICT closeterms_distN = sdn->closeterms_distN;
   int* const RESTRICT closetermsN = sdn->closetermsN;
   const int nnodes = graph_get_nNodes(g);
   int* RESTRICT nodes_nhits;
   SCIP_Real* RESTRICT nodes_maxdist;

   assert(N >= 1);
   assert(!sdn->nodes_nhits);
   assert(!sdn->nodes_maxdist);
   assert(!sdn->hitlist);
   assert(!sdn->dijkdata);

   SCIP_CALL( SCIPallocBufferArray(scip, &(sdn->nodes_nhits), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(sdn->nodes_maxdist), nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(sdn->hitlist), nnodes) );
   nodes_nhits = sdn->nodes_nhits;
   nodes_maxdist = sdn->nodes_maxdist;

   SCIP_CALL( graph_dijkLimited_init(scip, g, &(sdn->dijkdata)) );
   sdn->dijkdata->edgelimit = STP_SDN_MAXNVISITS;

   graph_init_dcsr(scip, g);

   for( int i = 0; i < nnodes; ++i )
   {
      nodes_nhits[i] = 0;
      nodes_maxdist[i] = 0.0;
      sdn->nodes_isBlocked[i] = FALSE;
   }

   for( int i = 0; i < N * nnodes; ++i )
   {
      closeterms_distN[i] = FARAWAY;
      closetermsN[i] = UNKNOWN;
   }

   return SCIP_OKAY;
}


/** helper */
static
void sdneighborUpdateFree(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph  */
   SDN*                  sdn                 /**< SD neighbors */
)
{
   graph_free_dcsr(scip, g);
   graph_dijkLimited_free(scip, &(sdn->dijkdata));

   SCIPfreeBufferArray(scip, &(sdn->hitlist));
   SCIPfreeBufferArray(scip, &(sdn->nodes_maxdist));
   SCIPfreeBufferArray(scip, &(sdn->nodes_nhits));
}


/** updates single node */
static inline
void sdneighborUpdateNode(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   int                   node,               /**< node to update */
   int                   term,               /**< terminal */
   SDN*                  sdn,                /**< SD neighbors */
   SDGRAPH*              sdgraph,            /**< SD graph */
   int*                  nupdates            /**< number of distance updates */
)
{
   SCIP_Real* RESTRICT nodes_maxdist = sdn->nodes_maxdist;
   const SCIP_Real sd_new = nodes_maxdist[node];

   assert(node != term);
   assert(GT(sd_new, 0.0) && LT(sd_new, FARAWAY));

   if( Is_term(g->term[node]) )
   {
      /* update the SD MST */
      const SCIP_Real sd_org = reduce_sdgraphGetSd(term, node, sdgraph);
      if( LT(sd_new, sd_org) )
      {
         (*nupdates)++;
         sdn->nodes_isBlocked[term] = TRUE;
         sdgraphInsertEdge(scip, term, node, sd_new, sdgraph);
  //graph_knot_printInfo(g, term);  graph_knot_printInfo(g, i); printf("term %d->%d: new=%f vs org=%f \n", term, i, nodes_maxdist[i], sdorg);
      }
   }
   else
   {
      SCIP_Real* RESTRICT closeterms_distN = sdn->closeterms_distN;
      int* RESTRICT closetermsN = sdn->closetermsN;
      SCIP_Real sd_max = 0.0;
      int node_max = -1;
      const int N = sdn->N;
      const int nnodes = g->knots;

      assert(N >= 1);

      /* get maximum old SD */
      for( int i = 0; i < N; i++ )
      {
         const int node_i = node + i * nnodes;
         if( GT(closeterms_distN[node_i], sd_max)  )
         {
            sd_max = closeterms_distN[node_i];
            node_max = node_i;
         }
      }

      assert(GT(sd_max, 0.0));

      if( LT(sd_new, sd_max) )
      {
         (*nupdates)++;
         closeterms_distN[node_max] = sd_new;
         closetermsN[node_max] = term;
      }
   }
}


/** updates SDs by using neighbor argument */
static
SCIP_RETCODE sdneighborUpdateExec(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph to initialize from */
   SD*                   sddata,             /**< SD */
   int*                  nupdates            /**< number of distance updates */
)
{
   SDN* const sdneighbors = sddata->sdneighbors;
   SDGRAPH* const sdgraph = sddata->sdgraph;
   int* RESTRICT nodes_nhits = sdneighbors->nodes_nhits;
   SCIP_Real* RESTRICT nodes_maxdist = sdneighbors->nodes_maxdist;
   int* RESTRICT hitlist = sdneighbors->hitlist;
   const int nnodes = graph_get_nNodes(g);
   const int maxdegree = STP_SDN_MAXDEGREE;

   *nupdates = 0;

   for( int term = 0; term < nnodes; ++term )
   {
      if( Is_term(g->term[term]) )
      {
         const int degree = g->grad[term];

         if( degree > maxdegree )
            continue;

#ifndef NDEBUG
         for( int i = 0; i < nnodes; ++i )
         {
            assert(nodes_nhits[i] == 0);
            assert(EQ(nodes_maxdist[i], 0.0));
         }
#endif
         sdneighbors->nnodesHit = 0;

         /* loop over all neighbors */
         for( int e = g->outbeg[term]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int neighbor = g->head[e];
            SCIP_CALL( sdneighborMarkCloseNodes(scip, g, neighbor, sddata) );
         }

         for( int i = 0; i < sdneighbors->nnodesHit; i++ )
         {
            const int hitnode = hitlist[i];

            assert(graph_knot_isInRange(g, hitnode));
            assert(0 < nodes_nhits[hitnode] && nodes_nhits[hitnode] <= degree);

            if( nodes_nhits[hitnode] == degree && hitnode != term )
            {
               sdneighborUpdateNode(scip, g, hitnode, term, sdneighbors, sdgraph, nupdates);
            }

            nodes_nhits[hitnode] = 0;
            nodes_maxdist[hitnode] = 0.0;
         }
      }
   }

   return SCIP_OKAY;
}


/** updates SDs by using neighbor argument */
static
SCIP_RETCODE sdneighborUpdate(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph to initialize from */
   SD*                   sddata,             /**< SD */
   int*                  nupdates            /**< number of distance updates */
)
{
   SDN* sdn = sddata->sdneighbors;
   SDGRAPH* sdgraph = sddata->sdgraph;
   assert(sdgraph);

//#define STP_SDN_PRINT

#ifdef STP_SDN_PRINT
   printf("start SDN update \n");
#endif

   sdneighborUpdateInit(scip, g, sdn);
   SCIP_CALL( sdneighborUpdateExec(scip, g, sddata, nupdates) );

#ifdef STP_SDN_PRINT
   printf("\n before \n");

   for( int i = 0; i < MIN(10, g->terms - 1); i++ )
      printf("%f ", sdgraph->mstcosts[i]);

   printf("\n after \n");
#endif

   SCIP_CALL( reduce_sdgraphMstBuild(scip, g, sdgraph) );
   reduce_sdgraphMstSortCosts(sdgraph);

#ifdef STP_SDN_PRINT
   for( int i = 0; i < MIN(10, g->terms - 1); i++ )
      printf("%f ", sdgraph->mstcosts[i]);

   printf("\n");
   printf("g->terms=%d nupdates=%d \n", g->terms, *nupdates);
   //exit(1);
#endif

   sdneighborUpdateFree(scip, g, sdn);

   return SCIP_OKAY;

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

   s->hasNeigborUpdate = FALSE;
   s->isBiased = FALSE;
   s->sdprofit = NULL;
   s->blctree = NULL;
   s->sdneighbors = NULL;
   SCIP_CALL( graph_tpathsInit(scip, g, &(s->terminalpaths)) );
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

   s->hasNeigborUpdate = FALSE;
   s->isBiased = TRUE;
   s->sdneighbors = NULL;
   s->blctree = NULL;
   SCIP_CALL( reduce_sdprofitInit(scip, g, &(s->sdprofit)) );
   SCIP_CALL( graph_tpathsInitBiased(scip, s->sdprofit, g, &(s->terminalpaths)) );
   SCIP_CALL( reduce_sdgraphInitBiased(scip, g, s->sdprofit, &(s->sdgraph)) );
   reduce_sdgraphInitOrderedMstCosts(s->sdgraph);

   return SCIP_OKAY;
}



/** initializes fully biased SD structure */
SCIP_RETCODE reduce_sdInitBiasedBottleneck(
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

   s->hasNeigborUpdate = FALSE;
   s->isBiased = TRUE;
   s->sdneighbors = NULL;
   SCIP_CALL( reduce_sdprofitInit(scip, g, &(s->sdprofit)) );
   SCIP_CALL( reduce_blctreeInit(scip, g, &(s->blctree)) );
   SCIP_CALL( reduce_sdprofitUpdateFromBLC(scip, g, s->blctree, TRUE, s->sdprofit) );

   SCIP_CALL( graph_tpathsInitBiased(scip, s->sdprofit, g, &(s->terminalpaths)) );
   SCIP_CALL( reduce_sdgraphInitBiasedFromTpaths(scip, g, s->sdprofit, s->terminalpaths, &(s->sdgraph)) );

//   SCIP_CALL( reduce_sdgraphInitBiased(scip, g, s->sdprofit, &(s->sdgraph)) );

   reduce_sdgraphInitOrderedMstCosts(s->sdgraph);

   return SCIP_OKAY;
}


/** repairs SD structure for imminent edge elimination */
SCIP_RETCODE reduce_sdRepair(
   SCIP*                 scip,               /**< SCIP */
   int                   edge,               /**< edge to be deleted soon */
   SCIP_Bool             withEdgeReplacement,/**< with edge replacement? */
   GRAPH*                g,                  /**< graph NOTE: will mark the graph, thus not const :(
                                                  terrible design */
   SD*                   sd                  /**< to be repaired */
)
{
   assert(scip && g && sd);
   assert(graph_edge_isInRange(g, edge));
   assert(sd->terminalpaths);
   assert(!sd->hasNeigborUpdate);
   assert(!sd->isBiased);
   assert(!sd->sdprofit);
   assert(!sd->sdneighbors);

   SCIP_CALL( graph_tpathsRepair(scip, edge, withEdgeReplacement, g, (sd->terminalpaths)) );

   if( reduce_sdgraphEdgeIsInMst(sd->sdgraph, edge) )
   {
      const SCIP_Real edgecost = g->cost[edge];

      assert(EQ(g->cost[edge], g->cost[flipedge(edge)]));

      g->cost[edge] = FARAWAY;
      g->cost[flipedge(edge)] = FARAWAY;

      reduce_sdgraphFree(scip, &(sd->sdgraph));

      SCIP_CALL( reduce_sdgraphInit(scip, g, &(sd->sdgraph)) );
      reduce_sdgraphInitOrderedMstCosts(sd->sdgraph);

      g->cost[edge] = edgecost;
      g->cost[flipedge(edge)] = edgecost;
   }

   return SCIP_OKAY;
}


/** sets up SD repairing mechanism */
SCIP_RETCODE reduce_sdRepairSetUp(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   SD*                   sd                  /**< to be repaired */
)
{
   assert(scip && g && sd);

   SCIP_CALL(  graph_tpathsRepairSetUp(g, sd->terminalpaths) );

   return SCIP_OKAY;
}


/** adds biased neighbor SD structure */
SCIP_RETCODE reduce_sdAddNeighborSd(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   SD*                   sd                  /**< to add to */
)
{
   assert(scip && g && sd);
   assert(!sd->sdneighbors);
   assert(!sd->hasNeigborUpdate);

   SCIP_CALL( reduce_sdneighborInit(scip, g, &(sd->sdneighbors)) );

   sd->hasNeigborUpdate = TRUE;

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
   graph_tpathsFree(scip, &(s->terminalpaths));

   if( s->sdprofit )
      reduce_sdprofitFree(scip, &(s->sdprofit));

   if( s->blctree )
      reduce_blctreeFree(scip, &(s->blctree));

   if( s->sdneighbors )
      reduce_sdneighborFree(scip, &(s->sdneighbors));

   SCIPfreeMemory(scip, sd);
}



/** get blocked nodes */
const SCIP_Bool* reduce_sdneighborGetBlocked(const SDN* sdneighbors)
{
   assert(sdneighbors);


   return sdneighbors->nodes_isBlocked;
}


/** gets (up to) four close terminals to given node i;
 *  with strict upper bound on allowed distances */
void reduce_sdneighborGetCloseTerms(
   const GRAPH*          g,                  /**< graph data structure */
   const SDN*            sdneighbor,         /**< SD neighbor data structure */
   int                   node,               /**< node */
   SCIP_Real             maxdist_strict,     /**< maximum valid distance (strict) */
   int* RESTRICT         closeterms,         /**< four close terminals */
   SCIP_Real* RESTRICT   closeterms_dist,    /**< four close terminal distance */
   int* RESTRICT         ncloseterms         /**< number of close terminals found */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int N = sdneighbor->N;
   const int* const terms = sdneighbor->closetermsN;
   const SCIP_Real* const dists = sdneighbor->closeterms_distN;
   int nterms = 0;

   assert(graph_knot_isInRange(g, node));
   assert(N >= 1);

   for( int i = 0; i < N; i++ )
   {
      const int node_i = node + i * nnodes;

      if( terms[node_i] < 0 )
         break;

      if( LT(dists[node_i], maxdist_strict)  )
      {
         *ncloseterms = 1;
         closeterms_dist[nterms] = dists[node_i];
         closeterms[nterms] = terms[node_i];

         assert(graph_knot_isInRange(g, closeterms[nterms]));
         assert(GE(closeterms_dist[nterms], 0.0));

         nterms++;
      }
   }

   *ncloseterms = nterms;
}


/** initializes SDN */
SCIP_RETCODE reduce_sdneighborInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDN**                 sdn                 /**< SD neighbors */
)
{
   SDN* sdneighbors;
   const int nnodes = graph_get_nNodes(g);
   const int N = STP_SDN_NCLOSETERMS;

   assert(N >= 1);

   SCIP_CALL( SCIPallocMemory(scip, sdn) );

   sdneighbors = *sdn;
   sdneighbors->nnodesHit = -1;
   sdneighbors->N = N;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdneighbors->nodes_isBlocked), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdneighbors->closeterms_distN), N * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdneighbors->closetermsN), N * nnodes) );
   sdneighbors->hitlist = NULL;
   sdneighbors->nodes_maxdist = NULL;
   sdneighbors->nodes_nhits = NULL;
   sdneighbors->dijkdata = NULL;

   return SCIP_OKAY;
}


/** frees SDN */
void reduce_sdneighborFree(
   SCIP*                 scip,               /**< SCIP */
   SDN**                 sdn                 /**< SD neighbors */
)
{
   SCIPfreeMemoryArray(scip, &((*sdn)->closetermsN));
   SCIPfreeMemoryArray(scip, &((*sdn)->closeterms_distN));
   SCIPfreeMemoryArray(scip, &((*sdn)->nodes_isBlocked));

   SCIPfreeMemory(scip, sdn);
}


/** updates SDs by using neighbor argument
 *  NOTE: invalidates certain SD routines! */
SCIP_RETCODE reduce_sdUpdateWithSdNeighbors(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph to initialize from */
   SD*                   sddata,             /**< SD */
   int*                  nupdates            /**< number of distance updates */
)
{
   assert(scip && g && sddata && nupdates);
   assert(sddata->sdneighbors);

   sdneighborUpdate(scip, g, sddata, nupdates);

   return SCIP_OKAY;
}


/** initializes SD profit */
SCIP_RETCODE reduce_sdprofitInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDPROFIT**            sdprofit            /**< the SD profit */
)
{
   assert(scip && g);

   SCIP_CALL( sdprofitAlloc(scip, g, TRUE, sdprofit) );
   sdprofitBuild(g, *sdprofit);

   return SCIP_OKAY;
}


/** initializes SD profit */
SCIP_RETCODE reduce_sdprofitInit1stOnly(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SCIP_Real*      edgecosts,          /**< edge costs (w.r.t graph 'g') */
   SDPROFIT**            sdprofit            /**< the SD profit */
)
{
   assert(scip && g);

   SCIP_CALL( sdprofitAlloc(scip, g, FALSE, sdprofit) );
   sdprofitBuild1stOnly(g, edgecosts, *sdprofit);

   return SCIP_OKAY;
}


/** prints SD profit statistics */
void reduce_sdprofitPrintStats(
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit            /**< the SD profit */
)
{
   const int nnodes = graph_get_nNodes(g);
   int nnodes_profit = 0;
   int nnodes_nonprofit = 0;
   SCIP_Real avg = 0.0;
   const SCIP_Real* node_bias;

   assert(sdprofit);

   node_bias = sdprofit->nodes_bias;

   for( int k = 0; k < nnodes; k++ )
   {
      if( !Is_term(g->term[k]) )
      {
         if( GT(node_bias[k], 0.0) )
         {
            avg += node_bias[k];
            nnodes_profit++;
         }
         else
            nnodes_nonprofit++;

         continue;
      }
   }

   printf("Steiner nodes:  profitable=%d non-profitable=%d .... Profit: sum=%f avg=%f\n", nnodes_profit, nnodes_nonprofit,
         avg, avg / (SCIP_Real) nnodes_profit);
}


/** frees SD profit */
void reduce_sdprofitFree(
   SCIP*                 scip,               /**< SCIP */
   SDPROFIT**            sdprofit            /**< the SD profit */
)
{
   SDPROFIT* sdp;
   assert(scip && sdprofit);

   sdp = *sdprofit;
   assert(sdp);

   SCIPfreeMemoryArrayNull(scip, &(sdp->nodes_biassource2));
   SCIPfreeMemoryArrayNull(scip, &(sdp->nodes_bias2));
   SCIPfreeMemoryArray(scip, &(sdp->nodes_biassource));
   SCIPfreeMemoryArray(scip, &(sdp->nodes_bias));

   SCIPfreeMemory(scip, sdprofit);
}


/** updates implied profits by using exact bottleneck distances */
SCIP_RETCODE reduce_sdprofitUpdateFromBLC(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   const BLCTREE*        blctree,            /**< BLC tree */
   SCIP_Bool             blctree_isUpdated,  /**< BLC tree fresh? */
   SDPROFIT*             sdprofit            /**< the SD profit */
)
{
   SCIP_Real* RESTRICT candidate_bottlenecks;
   int* RESTRICT candidate_edges;
   int ncandidates;

   assert(scip && g && blctree && sdprofit);

   ncandidates = reduce_blctreeGetMstNedges(blctree);
   assert(ncandidates > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &candidate_bottlenecks, ncandidates) );
   SCIP_CALL( SCIPallocBufferArray(scip, &candidate_edges, ncandidates) );
   reduce_blctreeGetMstEdges(g, blctree, candidate_edges);
   reduce_blctreeGetMstBottlenecks(g, blctree, candidate_bottlenecks);

   for( int i = 0; i < ncandidates; ++i )
   {
      const int edge = candidate_edges[i];
      const SCIP_Real bdist = candidate_bottlenecks[i];
      const SCIP_Real edgecost = g->cost[edge];

      assert(graph_edge_isInRange(g, edge));
      assert(LE(edgecost, bdist));

      if( LT(edgecost, bdist) )
      {
         const int tail = g->tail[edge];
         const int head = g->head[edge];

         if( Is_term(g->term[tail]) )
         {
            sdprofitUpdateNode(g, head, tail, edgecost, bdist, blctree_isUpdated, sdprofit);
         }

         if( Is_term(g->term[head]) )
         {
            sdprofitUpdateNode(g, tail, head, edgecost, bdist, blctree_isUpdated, sdprofit);
         }
      }
   }

   SCIPfreeBufferArray(scip, &candidate_edges);
   SCIPfreeBufferArray(scip, &candidate_bottlenecks);

   return SCIP_OKAY;
}


/** builds implied profits by using exact bottleneck distances */
SCIP_RETCODE reduce_sdprofitBuildFromBLC(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph */
   const BLCTREE*        blctree,            /**< BLC tree */
   SCIP_Bool             blctree_isUpdated,  /**< BLC tree fresh? */
   SDPROFIT*             sdprofit            /**< the SD profit */
)
{
   assert(scip && g && blctree && sdprofit);

   sdprofitBuild(g, sdprofit);
   SCIP_CALL( reduce_sdprofitUpdateFromBLC(scip, g, blctree, blctree_isUpdated, sdprofit) );

   return SCIP_OKAY;
}
