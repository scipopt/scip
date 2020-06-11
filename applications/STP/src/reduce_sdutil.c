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

#define STP_SDN_NCLOSETERMS 4


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
   SCIP_Bool             edgemarkReady;      /**< edge mark available? */
};

/** see reduce.h */
struct special_distance_neighbors
{
   SCIP_Real*            closeterms_distN;   /**< paths to N closest terms */
   int*                  closetermsN;        /**< close terminals */
   SCIP_Bool*            nodes_isBlocked;    /**< node is blocked */
   int                   N;                  /**< the number of closest terms per node */
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
   g_sd->edgemarkReady = TRUE;

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


/** gets SD for boundary edge */
static inline
SCIP_Real distgraphGetBoundaryEdgeDist2(
   int                   tail,               /**< tail */
   int                   head,               /**< head */
   int                   vbase_tail,         /**< base */
   int                   vbase_head,         /**< base */
   SCIP_Real             edgecost,           /**< cost */
   SCIP_Real             dist_tail,          /**< distance */
   SCIP_Real             dist_head,          /**< distance */
   const SDPROFIT*       sdprofit            /**< profit */
   )
{
   const SCIP_Real profit_tail = reduce_sdprofitGetProfit(sdprofit, tail, vbase_head, vbase_tail);
   const SCIP_Real profit_head = reduce_sdprofitGetProfit(sdprofit, head, vbase_head, vbase_tail);
   SCIP_Real distSimple = edgecost + dist_tail + dist_head;
   const SCIP_Real distAll = distSimple - profit_tail - profit_head;
   const SCIP_Real distTailHead = edgecost + dist_tail - profit_tail;
   const SCIP_Real distHeadTail = edgecost + dist_head - profit_head;

   const SCIP_Real dist = miscstp_maxReal((SCIP_Real[])
               { distAll, distTailHead, distHeadTail, edgecost,
                   dist_tail, dist_head},
                 6);
   return dist;
}


/** gets SD for boundary edge ... choose among nearest terminals (w.r.t. implied SD) */
static inline
SCIP_Real distgraphGetBoundaryEdgeDistBest(
   const GRAPH*          g,                  /**< graph  */
   const TPATHS*         tpaths,             /**< terminal paths */
   int                   tail,               /**< tail */
   int                   head,               /**< head */
   SCIP_Real             edgecost,           /**< cost */
   const SDPROFIT*       sdprofit,           /**< profit */
   int*                  base_tail,          /**< base of tail */
   int*                  base_head           /**< base of head */
   )
{
   SCIP_Real dists_tail[3];
   SCIP_Real dists_head[3];
   int terms_tail[3];
   int terms_head[3];
   int nterms_tail;
   int nterms_head;
   SCIP_Real dist = FARAWAY;

   assert(g && tpaths && sdprofit);

   *base_tail = -1;
   *base_head = -1;

   graph_tpathsGet3CloseTerms(g, tpaths, tail, FARAWAY, terms_tail, NULL, dists_tail, &nterms_tail);
   graph_tpathsGet3CloseTerms(g, tpaths, head, FARAWAY, terms_head, NULL, dists_head, &nterms_head);

   for( int i = 0; i < nterms_tail; ++i )
   {
      for( int j = 0; j < nterms_head; ++j )
      {
         if( terms_tail[i] != terms_head[j] )
         {
            const SCIP_Real distnew =
            distgraphGetBoundaryEdgeDist2(tail, head, terms_tail[i], terms_head[j], edgecost, dists_tail[i], dists_head[j], sdprofit);

            if( distnew < dist )
            {
               *base_tail = terms_tail[i];
               *base_head = terms_head[j];
               dist = distnew;
            }
         }
      }
   }

   return dist;
}


/** inserts new edge */
static inline
void distgraphInsertEdge(
   SCIP*                  scip,               /**< SCIP data structure */
   int                    sdnode1,           /**< end node 1 */
   int                    sdnode2,           /**< end node 2 */
   SCIP_Real              edgecost,          /**< cost */
   int                    edgeid,            /**< ID or -1 */
   int* RESTRICT          edgeorg,           /**< IDs of edges or NULL */
   GRAPH*                 distgraph          /**< the SD graph */
)
{
   int ne;

#ifndef NDEBUG
   assert(sdnode1 != sdnode2);
   assert(0 <= sdnode1 && sdnode1 < distgraph->knots);
   assert(0 <= sdnode2 && sdnode2 < distgraph->knots);
   assert(GT(edgecost, 0.0));
   assert(edgeorg != NULL || edgeid == -1);
#endif

   /* find the corresponding edge in the distance graph */
   for( ne = distgraph->outbeg[sdnode1]; ne != EAT_LAST; ne = distgraph->oeat[ne] )
   {
      if( distgraph->head[ne] == sdnode2 )
         break;
   }

   /* edge exists already? */
   if( ne != EAT_LAST )
   {
      assert(ne >= 0);
      assert(distgraph->head[ne] == sdnode2);
      assert(distgraph->tail[ne] == sdnode1);

      if( distgraph->cost[ne] > edgecost )
      {
         distgraph->cost[ne]            = edgecost;
         distgraph->cost[Edge_anti(ne)] = edgecost;

         if( edgeorg != NULL )
            edgeorg[ne / 2] = edgeid;
      }
   }
   else
   {
      if( edgeorg != NULL )
         edgeorg[distgraph->edges / 2] = edgeid;

      graph_edge_add(scip, distgraph, sdnode1, sdnode2, edgecost, edgecost);
   }
}



/** adds edges to distance graph, given a Voronoi diagram */
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
            const SCIP_Real distance = useProfit ?
               distgraphGetBoundaryEdgeDist(tail, head, vbase_tail, vbase_head, g->cost[e], nodes_vdist, sdprofit)
               : (g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]);

           // if( LT(distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]) )
         //   printf("distance: %f < %f \n", distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]);

            assert(LE(distance, g->cost[e] + nodes_vdist[tail] + nodes_vdist[head]));
            assert(Is_term(g->term[vbase_tail]));
            assert(Is_term(g->term[vbase_head]));

            distgraphInsertEdge(scip, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, e,
                  edgeorg, distgraph);

         }
      }
   }
}


/** adds edges to distance graph, given terminal paths */
static
void distgraphAddEdgesFromTpaths(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph to initialize from */
   const int*            distnodes_id,       /**< IDs of nodes */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const TPATHS*         tpaths,             /**< terminal paths */
   GRAPH*                distgraph           /**< distance graph */
)
{
   const int nnodes = graph_get_nNodes(g);
   const SCIP_Bool useProfit = (sdprofit != NULL);

   /* add the edges */
   for( int tail = 0; tail < nnodes; tail++ )
   {
      SCIP_Real dist_tail;
      int vbase_tail;
      graph_tpathsGetClosestTerm(g, tpaths, tail, &vbase_tail, NULL, &dist_tail);

      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         SCIP_Real dist_head;
         int vbase_head;

         graph_tpathsGetClosestTerm(g, tpaths, head, &vbase_head, NULL, &dist_head);

         if( vbase_tail != vbase_head )
         {
            const SCIP_Real distance = useProfit ?
               distgraphGetBoundaryEdgeDist2(tail, head, vbase_tail, vbase_head, g->cost[e], dist_tail, dist_head, sdprofit)
               : (g->cost[e] + dist_tail + dist_head);

            assert(LE(distance, g->cost[e] + dist_tail + dist_head));
            assert(Is_term(g->term[vbase_tail]) && Is_term(g->term[vbase_head]));

            distgraphInsertEdge(scip, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, -1, NULL, distgraph);
         }
      }
   }
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
   const int* nodesid = sdgraph->nodemapOrgToDist;
   const int sdnode1 = nodesid[term1];
   const int sdnode2 = nodesid[term2];
   assert(sdgraph->distgraph);

   distgraphInsertEdge(scip, sdnode1, sdnode2, edgecost, -1, NULL, sdgraph->distgraph);
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


/** builds distance graph */
static
SCIP_RETCODE sdgraphBuildDistgraphFromTpaths(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const TPATHS*         tpaths,             /**< terminal paths */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   GRAPH* distgraph;
   int* RESTRICT distnodes_id = g_sd->nodemapOrgToDist;
   const int nedges_distgraph = distgraphGetNedges(g);

   assert(g->terms >= 1);

   SCIP_CALL( graph_init(scip, &(g_sd->distgraph), g->terms, nedges_distgraph, 1) );
   distgraph = g_sd->distgraph;
   distgraphAddNodes(g, distnodes_id, distgraph);
   distgraphAddEdgesFromTpaths(scip, g, distnodes_id, sdprofit, tpaths, distgraph);

   assert(graph_valid(scip, distgraph));

   return SCIP_OKAY;
}


/** updates distance graph */
static
void sdgraphUpdateDistgraphFromTpaths(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< profit or NULL */
   const TPATHS*         tpaths,             /**< terminal paths */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   GRAPH* RESTRICT distgraph = g_sd->distgraph;
   const int* distnodes_id = g_sd->nodemapOrgToDist;
   const int nnodes = graph_get_nNodes(g);
   const int edgelimit = MIN(2 * distgraph->edges, distgraph->esize);

   assert(LE(distgraph->edges, distgraph->esize));

   for( int tail = 0; tail < nnodes && distgraph->edges < edgelimit; tail++ )
   {
      for( int e = g->outbeg[tail]; e != EAT_LAST; e = g->oeat[e] )
      {
         const int head = g->head[e];
         int vbase_tail;
         int vbase_head;
         const SCIP_Real distance =
            distgraphGetBoundaryEdgeDistBest(g, tpaths, tail, head, g->cost[e], sdprofit, &vbase_tail, &vbase_head);

         if( LT(distance, FARAWAY) && LT(distance, sdgraphGetSd(vbase_tail, vbase_head, g_sd)) )
         {
            assert(vbase_tail >= 0 && vbase_head >= 0);

            //SCIPdebugMessage("add biased MST edge %d->%d (%f<%f) \n", vbase_tail, vbase_head, distance, sdgraphGetSd(vbase_tail, vbase_head, g_sd));

            distgraphInsertEdge(scip, distnodes_id[vbase_tail], distnodes_id[vbase_head], distance, -1, NULL, distgraph);

            if( distgraph->edges >= edgelimit )
               break;
         }
      }
   }
   assert(sdprofit);
}


/** builds MST costs (ordered) for distance graph */
static
void sdgraphMstSortCosts(
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
#ifndef NDEBUG
      const int e = mst[k].edge;
      assert(e >= 0);
      assert(GE(distgraph->cost[e], 0.0));
      assert(EQ(mst[k].dist, distgraph->cost[e]));
#endif

      mstcosts[k - 1] = mst[k].dist;
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


/** marks original edges corresponding to MST */
static
void sdgraphMstMarkOrgEdges(
   const GRAPH*          g,                  /**< graph to initialize from */
   const VNOI*           vnoi,               /**< Voronoi */
   const int*            distedge2org,       /**< array of size nedges / 2 */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   GRAPH* distgraph = g_sd->distgraph;
   PATH* RESTRICT mst = g_sd->sdmst;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);
   const int* const nodes_vbase = vnoi->nodes_base;
   const int* const nodes_vpred = vnoi->nodes_predEdge;
   STP_Bool* RESTRICT orgedges_isInMst = g_sd->halfedge_isInMst;

#ifndef NDEBUG
   const int nedges = graph_get_nEdges(g);

   for( int e = 0; e < nedges / 2; e++ )
      assert(!orgedges_isInMst[e]);
#endif

   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      const int ne = distedge2org[e / 2];

      assert(e >= 0);
      assert(e < nedges);

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
}

/** builds MST on distance graph */
static
SCIP_RETCODE sdgraphMstBuild(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph to initialize from */
   SDGRAPH*              g_sd                /**< the SD graph */
)
{
   PATH* RESTRICT mst = g_sd->sdmst;
   GRAPH* distgraph = g_sd->distgraph;
   SCIP_Real maxcost;
   const int nnodes_distgraph = graph_get_nNodes(distgraph);
   STP_Bool* RESTRICT orgedges_isInMst = g_sd->halfedge_isInMst;
   const int nedges = graph_get_nEdges(g);

   if( orgedges_isInMst)
   {
      for( int e = 0; e < nedges / 2; e++ )
         orgedges_isInMst[e] = FALSE;
   }

   for( int k = 0; k < nnodes_distgraph; k++ )
      distgraph->mark[k] = TRUE;

   SCIP_CALL( graph_path_init(scip, distgraph) );
   graph_path_exec(scip, distgraph, MST_MODE, distgraph->source, distgraph->cost, mst);
   graph_path_exit(scip, distgraph);

   assert(mst[0].edge == -1);

   maxcost = 0.0;
   for( int k = 1; k < nnodes_distgraph; k++ )
   {
      const int e = mst[k].edge;
      const SCIP_Real cost = distgraph->cost[e];

      if( cost > maxcost )
         maxcost = cost;
   }

   g_sd->mstmaxcost = maxcost;

   return SCIP_OKAY;
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
   int                   sourcenode,
   int* RESTRICT         nodes_nhits,
   SCIP_Real* RESTRICT   nodes_maxdist,
   DIJK*                 dijkdata,
   SD*                   sddata              /**< SD */
)
{
   const int* const visitlist = dijkdata->visitlist;
   SCIP_Real* RESTRICT distance = dijkdata->node_distance;
   int nvisits;
int todo;

   assert(dijkdata && g);

// todo try bias!
   SCIP_CALL( graph_sdCloseNodesBiased(scip, g, NULL, sourcenode, dijkdata) );

   nvisits = dijkdata->nvisits;
   assert(nvisits >= 0);

   for( int k = 0; k < nvisits; k++ )
   {
      const int node = visitlist[k];
      const SCIP_Real dist = distance[node];

      if( node == sourcenode )
         continue;

      assert(GT(dist, 0.0));

      if( nodes_nhits[node] == 0 )
      {
         // todo save node in list
      }

      nodes_nhits[node]++;

      if( dist > nodes_maxdist[node] )
         nodes_maxdist[node] = dist;
   }

   graph_dijkLimited_reset(g, dijkdata);

   return SCIP_OKAY;
}


/** updates SDs by using neighbor argument */
static
SCIP_RETCODE sdneighborUpdate(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph to initialize from */
   SDN*                  sdn,                /**< SD neighbors */
   SD*                   sddata,             /**< SD */
   int*                  nupdates            /**< number of distance updates */
)
{
   SDGRAPH* sdgraph = sddata->sdgraph;
   DIJK* dijkdata;
   int* RESTRICT nodes_nhits;
   SCIP_Real* RESTRICT nodes_maxdist;
   SCIP_Real* RESTRICT closeterms_dist = sdn->closeterms_distN;
   int* RESTRICT closeterms = sdn->closetermsN;
   const int nnodes = graph_get_nNodes(g);
   const int maxdegree = 4;

   SCIP_CALL(SCIPallocBufferArray(scip, &nodes_nhits, nnodes));
   SCIP_CALL(SCIPallocBufferArray(scip, &nodes_maxdist, nnodes));

   SCIP_CALL( graph_dijkLimited_init(scip, g, &dijkdata) );
   dijkdata->edgelimit = 100;

   graph_init_dcsr(scip, g);

   // todo extra method
   for( int i = 0; i < nnodes; ++i )
   {
      closeterms_dist[i] = FARAWAY;
      closeterms[i] = UNKNOWN;
      nodes_nhits[i] = 0;
      nodes_maxdist[i] = 0.0;
      sdn->nodes_isBlocked[i] = FALSE;
   }

   assert(sdgraph);

   *nupdates = 0;

   // todo extra method
   for( int term = 0; term < nnodes; ++term )
   {
      if( Is_term(g->term[term]) )
      {
         const int degree = g->grad[term];

         if( degree > maxdegree )
            continue;

         // todo clean-up method
         for( int i = 0; i < nnodes; ++i )
         {
            nodes_nhits[i] = 0;
            nodes_maxdist[i] = 0.0;
         }

         /* loop over all neighbors */
         for( int e = g->outbeg[term]; e != EAT_LAST; e = g->oeat[e] )
         {
            const int neighbor = g->head[e];

            SCIP_CALL( sdneighborMarkCloseNodes(scip, g, neighbor, nodes_nhits, nodes_maxdist, dijkdata, sddata) );
         }

         // todo loop over list
         for( int i = 0; i < nnodes; ++i )
         {
            assert(nodes_nhits[i] <= degree);

            if( nodes_nhits[i] == degree && i != term )
            {
               const SCIP_Real sd_new = nodes_maxdist[i];
               assert(GT(sd_new, 0.0) && LT(sd_new, FARAWAY));

               if( Is_term(g->term[i]) )
               {
                  /* update the SD MST */
                  const SCIP_Real sd_org = reduce_sdgraphGetSd(term, i, sdgraph);
                  if( LT(sd_new, sd_org) )
                  {
                     (*nupdates)++;
                     sdn->nodes_isBlocked[term] = TRUE;
                     sdgraphInsertEdge(scip, term, i, sd_new, sdgraph);

     //     graph_knot_printInfo(g, term);  graph_knot_printInfo(g, i); printf("term %d->%d: new=%f vs org=%f \n", term, i, nodes_maxdist[i], sdorg);
                  }
               }
               else
               {
                  int todo; //  todo proper update
                  /* update the terminal paths */
                  if( LT(sd_new, closeterms_dist[i]) )
                  {
                     (*nupdates)++;
                     closeterms_dist[i] = sd_new;
                     closeterms[i] = term;
                  }
               }
            }
         }

      }
   }


   printf("\n before \n");

   for( int i = 0; i < MIN(10, g->terms - 1); i++ )
   {
      printf("%f ", sdgraph->mstcosts[i]);
   }
   printf("\n after \n");

   SCIP_CALL( sdgraphMstBuild(scip, g, sdgraph) );
   sdgraphMstSortCosts(sdgraph);


   for( int i = 0; i < MIN(10, g->terms - 1); i++ )
   {
      printf("%f ", sdgraph->mstcosts[i]);

   }
   printf("\n");

   printf("g->terms=%d nupdates=%d \n", g->terms, *nupdates);
  // exit(1);

   graph_free_dcsr(scip, g);

   graph_dijkLimited_free(scip, &dijkdata);

   SCIPfreeBufferArray(scip, &nodes_maxdist);
   SCIPfreeBufferArray(scip, &nodes_nhits);

   return SCIP_OKAY;

}


/*
 * Interface methods
 */

/** get blocked nodes */
const SCIP_Bool* reduce_sdneighborGetBlocked(const SDN* sdneighbors)
{
   assert(sdneighbors);


   return sdneighbors->nodes_isBlocked;
}


/** gets (up to) four close terminals to given node i;
 *  with strict upper bound on allowed distances */
void reduce_sdneighborGetCloseTerms(
   const GRAPH*          g,                  /**< graph */
   const SDN*            sdneighbor,         /**<  */
   int                   node,               /**< node */
   SCIP_Real             maxdist_strict,     /**< maximum valid distance (strict) */
   int* RESTRICT         closeterms,         /**< four close terminals */
   SCIP_Real* RESTRICT   closeterms_dist,    /**< four close terminal distance */
   int* RESTRICT         ncloseterms         /**< number of close terminals found */
)
{
   const int* terms = sdneighbor->closetermsN;
   const SCIP_Real* dists = sdneighbor->closeterms_distN;

   assert(graph_knot_isInRange(g, node));

   *ncloseterms = 0;

   if( terms[node] >= 0 && LT(dists[node], maxdist_strict)  )
   {
      *ncloseterms = 1;
      closeterms_dist[0] = dists[node];
      closeterms[0] = terms[node];

      assert(graph_knot_isInRange(g, closeterms[0]));
      assert(GE(closeterms_dist[0], 0.0));
   }

#if 0
   for( int i = 0; i < sdneighbor->N; i++ )
   {

   }
#endif

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
   sdneighbors->N = N;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdneighbors->nodes_isBlocked), nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdneighbors->closeterms_distN), N * nnodes) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(sdneighbors->closetermsN), N * nnodes) );

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

   sdneighborUpdate(scip, g, sddata->sdneighbors, sddata, nupdates);

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

   SCIP_CALL( sdprofitAlloc(scip, g, sdprofit) );
   sdprofitBuild(g, *sdprofit);

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
   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );
   sdgraphMstMarkOrgEdges(g, vnoi, edgeorg, *sdgraph);

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
   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );
   sdgraphMstMarkOrgEdges(g, vnoi, edgeorg, *sdgraph);

   sdgraphFinalize(scip, &vnoi, &edgeorg);

   return SCIP_OKAY;
}


/** initializes biased SD graph from given terminal paths */
SCIP_RETCODE reduce_sdgraphInitBiasedFromTpaths(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph to initialize from */
   const SDPROFIT*       sdprofit,           /**< SD profit */
   const TPATHS*         tpaths,             /**< terminal paths */
   SDGRAPH**             sdgraph             /**< the SD graph */
)
{
   assert(scip && g && sdgraph && tpaths);

   SCIP_CALL( sdgraphAlloc(scip, g, sdgraph) );
   SCIP_CALL( sdgraphBuildDistgraphFromTpaths(scip, g, sdprofit, tpaths, *sdgraph) );
   SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );

   if( sdprofit )
   {
      sdgraphUpdateDistgraphFromTpaths(scip, g, sdprofit, tpaths, *sdgraph);
      SCIP_CALL( sdgraphMstBuild(scip, g, *sdgraph) );
   }

   /* NOTE probably we never need that...for extending reductions we anyway should only take biased paths */
   (*sdgraph)->edgemarkReady = FALSE;
   SCIPfreeMemoryArray(scip, &(*sdgraph)->halfedge_isInMst);

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


/** returns edge mark */
const STP_Bool* reduce_sdgraphGetMstHalfMark(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   assert(reduce_sdgraphHasMstHalfMark(sdgraph));

   return sdgraph->halfedge_isInMst;
}


/** has edge mark? */
SCIP_Bool reduce_sdgraphHasMstHalfMark(
   const SDGRAPH*        sdgraph             /**< the SD graph */
)
{
   assert(sdgraph);
   if( sdgraph->edgemarkReady )
   {
      assert(sdgraph->halfedge_isInMst);
      return TRUE;
   }

   assert(!sdgraph->halfedge_isInMst);

   return FALSE;
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
      sdgraphMstSortCosts(sdgraph);
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

   SCIPfreeMemoryArrayNull(scip, &(g_sd->halfedge_isInMst));
   SCIPfreeMemoryArray(scip, &(g_sd->nodemapOrgToDist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstsdist));
   SCIPfreeMemoryArray(scip, &(g_sd->mstcosts));
   SCIPfreeMemoryArray(scip, &(g_sd->sdmst));

   graph_free(scip, &(g_sd->distgraph), TRUE);

   SCIPfreeMemory(scip, sdgraph);
}
