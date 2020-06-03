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


/** see reduce.h */
struct special_distance_graph
{
   GRAPH*                distgraph;          /**< (complete) distance graph */
   PATH*                 sdmst;              /**< MST on sdgraph */
   SCIP_Real*            mstcosts;           /**< maximum MST edge costs in descending order */
   SCIP_Real*            mstsdist;           /**< helper array; each entry needs to -1.0; of size nnodesorg */
   int*                  nodemapOrgToDist;   /**< node mapping from original graph to distance graph */
   STP_Bool*             halfedge_isInMst;   /**< signifies whether edge of original graph is part of MST
                                                  NOTE: operates on edges / 2! */
   SCIP_Real             mstmaxcost;         /**< maximum edge cost */
   int                   nnodesorg;          /**< number of nodes of original graph */
   int                   nedgesorg;          /**< number of edges of original graph */
   SCIP_Bool             mstcostsReady;      /**< are the mstcosts already available? */
};


/** allocates memory */
static
SCIP_RETCODE sdgraphAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   SCIP_Real* mstsdist;
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const int nterms = g->terms;
   SDGRAPH* g_sd;

   SCIP_CALL( SCIPallocMemory(scip, sdgraph) );
   g_sd = *sdgraph;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->sdmst), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstcosts), nterms) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->mstsdist), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->nodemapOrgToDist), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(g_sd->halfedge_isInMst), nedges / 2) );

   mstsdist = g_sd->mstsdist;
   for( int i = 0; i < nnodes; i++ )
   {
      mstsdist[i] = -1.0;
   }

   g_sd->nnodesorg = nnodes;
   g_sd->nedgesorg = nedges;
   g_sd->mstcostsReady = FALSE;

   return SCIP_OKAY;
}


/** gets number of edges */
static
int distgraphGetNedges(
   const GRAPH*          g                   /**< graph to initialize from */
   )
{
   int maxnedges;
   const int nedges = graph_get_nEdges(g);
   const int nterms = g->terms;
   const SCIP_Longint terms2 = (nterms - 1) * nterms;

   if( nedges >= terms2 )
   {
      assert(terms2 <= INT_MAX);
      maxnedges = terms2;
   }
   else
   {
      maxnedges = nedges;
   }

   return maxnedges;
}


/** adds nodes to distance graph */
static
void distgraphAddNodes(
   const GRAPH*          g,                  /**< graph to initialize from */
   int* RESTRICT         distnodes_id,       /**< IDs of nodes */
   GRAPH*                distgraph           /**< distance graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   int nnodes_distgraph = 0;

   /* add the nodes */
   for( int k = 0; k < nnodes; k++ )
   {
      if( Is_term(g->term[k]) )
      {
         graph_knot_add(distgraph, STP_TERM_NONE);
         distnodes_id[k] = nnodes_distgraph++;
      }
      else
      {
         distnodes_id[k] = UNKNOWN;
      }
   }

   assert(distgraph->knots == nnodes_distgraph);
   assert(distgraph->knots == g->terms);

   graph_knot_chg(distgraph, 0, STP_TERM);
   distgraph->source = 0;
}


/** gets SD for boundary edge */
static inline
SCIP_Real distgraphGetBoundaryEdgeDist(
   int                   tail,               /**< tail */
   int                   head,               /**< head */
   int                   vbase_tail,         /**< base */
   int                   vbase_head,         /**< base */
   SCIP_Real             edgecost,           /**< cost */
   const SCIP_Real*      nodes_vdist,        /**< distance */
   const SDPROFIT*       sdprofit            /**< profit */
   )
{
   const SCIP_Real profit_tail = reduce_sdprofitGetProfit(sdprofit, tail, vbase_head, vbase_tail);
   const SCIP_Real profit_head = reduce_sdprofitGetProfit(sdprofit, head, vbase_head, vbase_tail);
   SCIP_Real distSimple = edgecost + nodes_vdist[tail] + nodes_vdist[head];
   const SCIP_Real distAll = distSimple - profit_tail - profit_head;
   const SCIP_Real distTailHead = edgecost + nodes_vdist[tail] - profit_tail;
   const SCIP_Real distHeadTail = edgecost + nodes_vdist[head] - profit_head;

   const SCIP_Real dist = miscstp_maxReal((SCIP_Real[])
               { distAll, distTailHead, distHeadTail, edgecost,
                 nodes_vdist[tail], nodes_vdist[head] },
                 6);
   return dist;
}


/** adds nodes to distance graph */
static
void distgraphAddEdges(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph to initialize from */
   const int*            distnodes_id,       /**< IDs of nodes */
   const VNOI*           vnoi,               /**< Voronoi */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   int* RESTRICT         edgeorg,            /**< IDs of edges */
   GRAPH*                distgraph           /**< distance graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const int nedges = graph_get_nEdges(g);
   const SCIP_Real* nodes_vdist = vnoi->nodes_dist;
   const int* RESTRICT nodes_vbase = vnoi->nodes_base;
   const SCIP_Bool useProfit = (sdprofit != NULL);

   for( int e = 0; e < nedges / 2; e++ )
      edgeorg[e] = UNKNOWN;

   /* add the edges */
   for( int tail = 0; tail < nnodes; tail++ )
   {
      const int vbase_tail = nodes_vbase[tail];

      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         const int vbase_head = nodes_vbase[head];

         assert(tail == g->tail[e]);

         if( vbase_tail != vbase_head )
         {
            int ne;
            const SCIP_Real distance = useProfit ?
               distgraphGetBoundaryEdgeDist(tail, head, vbase_tail, vbase_head, g->cost[e], nodes_vdist, sdprofit)
               : (g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]);

           // if( LT(distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]) )
         //   printf("distance: %f < %f \n", distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]);

            assert(LE(distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]));
            assert(Is_term(g->term[vbase_tail]));
            assert(Is_term(g->term[vbase_head]));
            assert(distnodes_id[vbase_tail] >= 0);
            assert(distnodes_id[vbase_head] >= 0);

            /* find the corresponding edge in the distance graph */
            for( ne = distgraph->outbeg[distnodes_id[vbase_tail]]; ne != EAT_LAST; ne = distgraph->oeat[ne] )
            {
               if( distgraph->head[ne] == distnodes_id[vbase_head] )
                  break;
            }

            /* edge exists already? */
            if( ne != EAT_LAST )
            {
               assert(ne >= 0);
               assert(distgraph->head[ne] == distnodes_id[vbase_head]);
               assert(distgraph->tail[ne] == distnodes_id[vbase_tail]);

               if( distgraph->cost[ne] > distance )
               {
                  distgraph->cost[ne]            = distance;
                  distgraph->cost[Edge_anti(ne)] = distance;
                  edgeorg[ne / 2] = e;
                  assert(ne <= distgraphGetNedges(g));
               }
            }
            else
            {
               edgeorg[distgraph->edges / 2] = e;
               graph_edge_add(scip, distgraph, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, distance);
               assert(distgraph->edges <= distgraphGetNedges(g));
            }
         }
      }
   }
}


/** builds distance graph */
static
SCIP_RETCODE sdgraphBuildDistgraph(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   SDGRAPH*              g_sd,               /**< the SD graph */
   VNOI**                vnoi,               /**< Voronoi */
   int**                 distedge2org        /**< array of size nedges / 2 */
)
{
   GRAPH* distgraph;
   int* RESTRICT distnodes_id = g_sd->nodemapOrgToDist;
   const int nedges = graph_get_nEdges(g);
   const int nedges_distgraph = distgraphGetNedges(g);

   assert(g->terms >= 1);

   SCIP_CALL( SCIPallocBufferArray(scip, distedge2org, nedges / 2) );

   /* build biased Voronoi diagram */
   SCIP_CALL( graph_vnoiInit(scip, g, TRUE, vnoi) );

   if( sdprofit )
      graph_vnoiComputeImplied(scip, g, sdprofit, *vnoi);
   else
      graph_vnoiCompute(scip, g, *vnoi);

   /* build distance graph from Voronoi diagram */
   SCIP_CALL( graph_init(scip, &(g_sd->distgraph), g->terms, nedges_distgraph, 1) );
   distgraph = g_sd->distgraph;
   distgraphAddNodes(g, distnodes_id, distgraph);
   distgraphAddEdges(scip, g, distnodes_id, *vnoi, sdprofit, *distedge2org, distgraph);

   assert(graph_valid(scip, distgraph));

   return SCIP_OKAY;
}


/** builds MST costs (ordered) for distance graph */
static
void sdgraphBuildMstcosts(
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   const GRAPH* const distgraph = g_sd->distgraph;
   const PATH* const mst = g_sd->sdmst;
   SCIP_Real* RESTRICT mstcosts = g_sd->mstcosts;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);

   assert(mst[0].edge == UNKNOWN);
   assert(nnodes_distgraph >= 1);

   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      assert(e >= 0);
      assert(GE(distgraph->cost[e], 0.0));

      mstcosts[k - 1] = distgraph->cost[e];
   }

   SCIPsortDownReal(mstcosts, nnodes_distgraph - 1);

   /* debug sentinel */
   if( nnodes_distgraph > 1 )
   {
      mstcosts[nnodes_distgraph - 1] = -FARAWAY;
   }
   else
   {
      mstcosts[0] = 0.0;
   }
}


/** builds MST on distance graph */
static
SCIP_RETCODE sdgraphBuildMst(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const VNOI*           vnoi,               /**< Voronoi */
   const int*            distedge2org,       /**< array of size nedges / 2 */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   PATH* RESTRICT mst = g_sd->sdmst;
   GRAPH* distgraph = g_sd->distgraph;
   STP_Bool* RESTRICT orgedges_isInMst = g_sd->halfedge_isInMst;
   const int* const nodes_vbase = vnoi->nodes_base;
   const int* const nodes_vpred = vnoi->nodes_predEdge;
   SCIP_Real maxcost;
   const int nedges = graph_get_nEdges(g);
   const int nnodes_distgraph = graph_get_nNodes(distgraph);

   for( int k = 0; k < nnodes_distgraph; k++ )
      distgraph->mark[k] = TRUE;

   for( int e = 0; e < nedges / 2; e++ )
      orgedges_isInMst[e] = FALSE;

   SCIP_CALL( graph_path_init(scip, distgraph) );
   graph_path_exec(scip, distgraph, MST_MODE, distgraph->source, distgraph->cost, mst);
   graph_path_exit(scip, distgraph);

   assert(mst[0].edge == -1);

   maxcost = 0.0;
   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      const int ne = distedge2org[e / 2];
      const SCIP_Real cost = distgraph->cost[e];

      assert(e >= 0);

      if( cost > maxcost )
         maxcost = cost;

      orgedges_isInMst[ne / 2] = TRUE;

      for( int v = g->head[ne]; v != nodes_vbase[v]; v = g->tail[nodes_vpred[v]] )
      {
         orgedges_isInMst[nodes_vpred[v] / 2] = TRUE;
      }

      for( int v = g->tail[ne]; v != nodes_vbase[v]; v = g->tail[nodes_vpred[v]] )
      {
         orgedges_isInMst[nodes_vpred[v]/ 2] = TRUE;
      }

      assert(e != EAT_LAST);
   }

   g_sd->mstmaxcost = maxcost;

   return SCIP_OKAY;
}


/** Gets special distance (e.g. bottleneck distance) from graph.
 *  Corresponds to bottleneck length of path between term1 and term2 on distance graph */
static
SCIP_Real sdgraphGetSd(
   int                    term1,             /**< terminal 1 */
   int                    term2,             /**< terminal 2 */
   SDGRAPH*               sdgraph            /**< the SD graph */
)
{
   SCIP_Real* RESTRICT mstsdist = sdgraph->mstsdist;
   const GRAPH* distgraph = sdgraph->distgraph;
   const PATH* sdmst = sdgraph->sdmst;
   const int* nodesid = sdgraph->nodemapOrgToDist;
   const int sdnode1 = nodesid[term1];
   const int sdnode2 = nodesid[term2];
   SCIP_Real sdist = 0.0;
   int tempnode = sdnode1;

   assert(sdnode1 != sdnode2);
   assert(distgraph->source == 0);

   mstsdist[tempnode] = 0.0;

   /* not at root? */
   while( tempnode != 0 )
   {
      const int ne = sdmst[tempnode].edge;

      assert(distgraph->head[ne] == tempnode);
      tempnode = distgraph->tail[ne];

      if( distgraph->cost[ne] > sdist )
         sdist = distgraph->cost[ne];

      mstsdist[tempnode] = sdist;
      if( tempnode == sdnode2 )
         break;
   }

   /* already finished? */
   if( tempnode == sdnode2 )
   {
      tempnode = 0;
   }
   else
   {
      tempnode = sdnode2;
      sdist = 0.0;
   }

   while( tempnode != 0 )
   {
      const int ne = sdmst[tempnode].edge;
      tempnode = distgraph->tail[ne];

      if( distgraph->cost[ne] > sdist )
         sdist = distgraph->cost[ne];

      /* already visited? */
      if( GE(mstsdist[tempnode], 0.0) )
      {
         if( mstsdist[tempnode] > sdist )
            sdist = mstsdist[tempnode];
         break;
      }

#ifndef NDEBUG
      assert(EQ(mstsdist[tempnode], -1.0));

      if( tempnode == 0 )
      {
         assert(sdnode1 == 0);
      }
#endif
   }

   /* restore mstsdist */
   tempnode = sdnode1;
   mstsdist[tempnode] = -1.0;
   while( tempnode != 0 )
   {
      const int ne = sdmst[tempnode].edge;
      tempnode = distgraph->tail[ne];
      mstsdist[tempnode] = -1.0;
      if( tempnode == sdnode2 )
         break;
   }

   assert(GT(sdist, 0.0));

   return sdist;
}


/** finalizes distance graph */
static
void sdgraphFinalize(
   SCIP*                 scip,               /**< SCIP data structure */
   VNOI**                vnoi,               /**> Voronoi */
   int**                 edgeorg             /**< array of size nedges / 2 */
)
{
   graph_vnoiFree(scip, vnoi);
   SCIPfreeBufferArray(scip, edgeorg);
}


/** allocates memory */
static
SCIP_RETCODE sdprofitAlloc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDPROFIT**            sdprofit            /**< SD profit */
)
{
   const int nnodes = graph_get_nNodes(g);
   SDPROFIT* sdp;

   SCIP_CALL( SCIPallocMemory(scip, sdprofit) );
   sdp = *sdprofit;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_bias), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_biassource), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_bias2), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdp->nodes_biassource2), nnodes) );

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



/*
 * Interface methods
 */

/** initializes SD profit */
SCIP_RETCODE reduce_sdprofitInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDPROFIT**            sdprofit            /**< the SD profit */
)
{
   assert(scip && g);

   SCIP_CALL( sdprofitAlloc(scip, g, sdprofit) );
   sdprofitBuild(g, *sdprofit);

   return SCIP_OKAY;
}


/** updates SDs */
SCIP_RETCODE reduce_sdneighborUpdate(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SD*                   sddata              /**< SD */
)
{
   const int nnodes = graph_get_nNodes(g);
   int* nodes_nhits;
   SCIP_Real* nodes_maxdist;

   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_nhits, nnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodes_maxdist, nnodes) );


   SCIP_Real termdist1[4];
   int neighbterms1[4];
   SDGRAPH* sdgraph = sddata->sdgraph;
   int nnterms1;
int nupdates = 0;

   assert(sdgraph);


   for( int i = 0; i < nnodes; ++i )
   {
      nodes_nhits[i] = 0;
      nodes_maxdist[i] = 0.0;
   }

   assert(scip && g);

   for( int term = 0; term < nnodes; ++term )
   {
      if( Is_term(g->term[term]) )
      {
         const int degree = g->grad[term];

         for( int i = 0; i < nnodes; ++i )
         {
            nodes_nhits[i] = 0;
            nodes_maxdist[i] = 0.0;
         }


         for( int e = g->outbeg[term]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int head = g->head[e];

            reduce_tpathsGet4CloseTerms(g, sddata->terminalpaths, head, FARAWAY, neighbterms1, termdist1, &nnterms1);

            // go over all 4 neighbors
            for( int k = 0; k < nnterms1; k++ )
            {
               const int neighborterm = neighbterms1[k];

               assert(Is_term(g->term[neighborterm]));

               nodes_nhits[neighborterm]++;

               if( termdist1[k] > nodes_maxdist[neighborterm] )
                  nodes_maxdist[neighborterm] = termdist1[k];
            }
         }






         for( int i = 0; i < nnodes; ++i )
         {
            if( nodes_nhits[i] == degree && i != term )
            {
               SCIP_Real sdorg = reduce_sdgraphGetSd(term, i, sdgraph);

               if( LT(nodes_maxdist[i], sdorg) )
               {
                  graph_knot_printInfo(g, term);
                  graph_knot_printInfo(g, i);
                  printf("term %d->%d: new=%f vs org=%f \n", term, i, nodes_maxdist[i], sdorg);
                  nupdates++;
               }
            }
         }
      }
   }

   printf("g->terms=%d nupdates=%d \n", g->terms, nupdates);


   SCIPfreeBufferArray(scip, &nodes_maxdist);
   SCIPfreeBufferArray(scip, &nodes_nhits);

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

   SCIPfreeMemoryArray(scip, &(sdp->nodes_biassource2));
   SCIPfreeMemoryArray(scip, &(sdp->nodes_bias2));
   SCIPfreeMemoryArray(scip, &(sdp->nodes_biassource));
   SCIPfreeMemoryArray(scip, &(sdp->nodes_bias));

   SCIPfreeMemory(scip, sdprofit);
}


/** initializes SD graph */
SCIP_RETCODE reduce_sdgraphInit(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   VNOI* vnoi;
   int* edgeorg;
   assert(scip && g && sdgraph);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );
   SCIP_CALL( sdgraphBuildDistgraph(scip, g, NULL, *sdgraph, &vnoi, &edgeorg) );
   SCIP_CALL( sdgraphBuildMst(scip, g, vnoi, edgeorg, *sdgraph) );

   sdgraphFinalize(scip, &vnoi, &edgeorg);

   return SCIP_OKAY;
}


/** initializes biased SD graph */
SCIP_RETCODE reduce_sdgraphInitBiased(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   VNOI* vnoi;
   int* edgeorg;
   assert(scip && g && sdgraph);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );
   SCIP_CALL( sdgraphBuildDistgraph(scip, g, sdprofit, *sdgraph, &vnoi, &edgeorg) );
   SCIP_CALL( sdgraphBuildMst(scip, g, vnoi, edgeorg, *sdgraph) );

   sdgraphFinalize(scip, &vnoi, &edgeorg);

   return SCIP_OKAY;
}


/** return maximum MST edge cost */
SCIP_Real reduce_sdgraphGetMaxCost(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(GE(sdgraph->mstmaxcost, 0.0));

   return sdgraph->mstmaxcost;
}


/** initializes SD graph */
const STP_Bool* reduce_sdgraphGetMstHalfMark(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->halfedge_isInMst);

   return sdgraph->halfedge_isInMst;
}


/** MST costs in descending order available? */
SCIP_Bool reduce_sdgraphHasOrderedMstCosts(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->mstcosts);

   return sdgraph->mstcostsReady;
}


/** Gets special distance (e.g. bottleneck distance) from distance graph.
 *  Only works if both nodes are terminals!  */
SCIP_Real reduce_sdgraphGetSd(
   int                    term1,             /**< node 1 */
   int                    term2,             /**< node 2 */
   SDGRAPH*               sdgraph            /**< the SD graph */
)
{
#ifndef NDEBUG
   assert(sdgraph);
   assert(term1 != term2);
   assert(sdgraph->mstcosts);
   assert(sdgraph->mstsdist);
   assert(0 <= term1 && term1 < sdgraph->nnodesorg);
   assert(0 <= term2 && term2 < sdgraph->nnodesorg);
   assert(sdgraph->nodemapOrgToDist[term1] != UNKNOWN);
   assert(sdgraph->nodemapOrgToDist[term2] != UNKNOWN);

   for( int i = 0; i < sdgraph->nnodesorg; i++ )
      assert(EQ(sdgraph->mstsdist[i], -1.0));
#endif

   return sdgraphGetSd(term1, term2, sdgraph);
}


/** initializes all MST costs in descending order */
void reduce_sdgraphInitOrderedMstCosts(
   SDGRAPH*              sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(sdgraph->mstcosts);

   if( !sdgraph->mstcostsReady )
   {
      sdgraphBuildMstcosts(sdgraph);
      sdgraph->mstcostsReady = TRUE;
   }

   assert(reduce_sdgraphHasOrderedMstCosts(sdgraph));
   assert(sdgraph->mstcosts[0] == sdgraph->mstmaxcost);
}


/** returns all MST costs in descending order */
const SCIP_Real* reduce_sdgraphGetOrderedMstCosts(
   const SDGRAPH*         sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(reduce_sdgraphHasOrderedMstCosts(sdgraph));
   assert(sdgraph->mstcosts[0] == sdgraph->mstmaxcost);

   return sdgraph->mstcosts;
}


/** frees SD graph */
void reduce_sdgraphFree(
   SCIP*                 scip,               /**< SCIP */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   SDGRAPH* g_sd;
   assert(scip && sdgraph);

   g_sd = *sdgraph;
   assert(g_sd);

   SCIPfreeMemoryArray(scip, &(g_sd->halfedge_isInMst));
   SCIPfreeMemoryArray(scip, &(g_sd->nodemapOrgToDist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstsdist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstcosts));
   SCIPfreeMemoryArray(scip, &(g_sd->sdmst));

   graph_free(scip, &(g_sd->distgraph), TRUE);

   SCIPfreeMemory(scip, sdgraph);
}
