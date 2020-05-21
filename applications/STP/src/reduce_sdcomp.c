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

/**@file   reduce_sdcomp.c
 * @brief  special distance (bottleneck distance) component reduction methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various special distance (aka bottleneck distance) based component reduction methods
 * for Steiner problems.
 *
 * A list of all interface methods can be found in reduce.h.
 *
 *
 */

//#define SCIP_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "scip/scip.h"
#include "portab.h"


#define STP_BDKIMP_MAXDEGREE 5
#define STP_BDKIMP_MAXNEDGES 10


/** BD_k storage */
typedef struct bottleneck_distance_storage
{
   SD*                   sdistance;          /**< special distance storage */
   STAR*                 star;               /**< star structure for neighborhood of node */
   GRAPH*                cliquegraph;        /**< complete graph on adjacent vertices */
   PATH*                 clique_mst;          /**< MST on cliquegraph */
   int*                  node_outedges;      /**< for node: outgoing edges (size STP_BDKIMP_MAXNEDGES) */
   int*                  node_neighbors;     /**< for node: adjacent vertices (size STP_BDKIMP_MAXDEGREE) */
   const int*            star_outedges;      /**< for star: outgoing edges */
   int                   star_degree;        /**< degree of star */
} BDK;


/** initializes data for bdk test */
static
SCIP_RETCODE bdkInit(
   SCIP*                 scip,               /**< SCIP data structure */
   SD*                   sdistance,          /**< special distance storage */
   BDK**                 bdk                 /**< storage */
)
{
   BDK* bdk_d;
   GRAPH* cliquegraph;

   SCIP_CALL( SCIPallocMemory(scip, bdk) );
   bdk_d = *bdk;

   bdk_d->sdistance = sdistance;
   bdk_d->star_degree = -1;
   bdk_d->star_outedges = NULL;
   SCIP_CALL( reduce_starInit(scip, STP_BDKIMP_MAXDEGREE, &(bdk_d->star)) );

   SCIP_CALL( graph_buildCompleteGraph(scip, &cliquegraph, STP_BDKIMP_MAXDEGREE) );
   SCIP_CALL( graph_path_init(scip, cliquegraph) );
   assert(cliquegraph->edges == 2 * STP_BDKIMP_MAXNEDGES);
   bdk_d->cliquegraph = cliquegraph;
   for( int i = 0; i < STP_BDKIMP_MAXDEGREE; i++ )
      cliquegraph->mark[i] = TRUE;

   SCIP_CALL( SCIPallocMemoryArray(scip, &(bdk_d->clique_mst), STP_BDKIMP_MAXDEGREE) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(bdk_d->node_outedges), STP_BDKIMP_MAXNEDGES) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(bdk_d->node_neighbors), STP_BDKIMP_MAXDEGREE) );

   return SCIP_OKAY;
}

/** frees data for bdk test */
static
void bdkFree(
   SCIP*                 scip,               /**< SCIP data structure */
   BDK**                 bdk                 /**< storage */
)
{
   BDK* bdk_d = *bdk;
   GRAPH* cliquegraph = bdk_d->cliquegraph;

   SCIPfreeMemoryArray(scip, &(bdk_d->node_neighbors));
   SCIPfreeMemoryArray(scip, &(bdk_d->node_outedges));
   SCIPfreeMemoryArray(scip, &(bdk_d->clique_mst));

   graph_path_exit(scip, cliquegraph);
   graph_free(scip, &cliquegraph, TRUE);

   reduce_starFree(scip, &(bdk_d->star));

   SCIPfreeMemory(scip, bdk);
}


/** gets neighborhood information for bdk test */
static inline
void bdkGetNeighborhood(
   const GRAPH*          g,                 /**< graph data structure */
   int                   i,                 /**< the node */
   BDK*                  bdk                /**< storage */
)
{
   int* RESTRICT edges = bdk->node_outedges;
   int* RESTRICT adjverts = bdk->node_neighbors;
   int k = 0;

   for( int e = g->outbeg[i]; e != EAT_LAST; e = g->oeat[e] )
   {
      edges[k] = e;
      adjverts[k++] = g->head[e];
   }
}

/** gets SDs bdk test; stored in cliquegraph */
static inline
void bdkGetCliqueSds(
   const GRAPH*          g,                 /**< graph data structure */
   int                   node,              /**< the node */
   BDK*                  bdk                /**< storage */
)
{
   GRAPH* cliquegraph = bdk->cliquegraph;
   const int degree = g->grad[node];
   int* nodemark = cliquegraph->mark;

   assert(degree >= 3);
   assert(STP_BDKIMP_MAXDEGREE == cliquegraph->knots);
   assert(nodemark[0] && nodemark[1] && nodemark[2]);

   for( int k = 3; k < degree; k++ )
      nodemark[k] = TRUE;

   for( int k = degree; k < STP_BDKIMP_MAXDEGREE; k++ )
      nodemark[k] = FALSE;

   reduce_sdGetSdsCliquegraph(g, bdk->node_neighbors, bdk->sdistance, cliquegraph);
}


/** gets SDs bdk test; stored in cliquegraph */
static inline
void bdkGetCutoffs(
   const GRAPH*          g,                 /**< graph data structure */
   const BDK*            bdk,               /**< storage */
   int                   node,              /**< the node */
   SCIP_Real*            cutoffs            /**< cutoffs */
)
{
   const GRAPH* cliquegraph = bdk->cliquegraph;
   const int degree = g->grad[node];
   const int* const nodemark = cliquegraph->mark;
   int edgecount = 0;

  // printf("go degree=%d \n", degree);

   for( int k = 0; k < STP_BDKIMP_MAXDEGREE - 1; k++ )
   {
      if( !nodemark[k] )
         continue;

      for( int e = cliquegraph->outbeg[k]; e != EAT_LAST; e = cliquegraph->oeat[e] )
      {
         const int k2 = cliquegraph->head[e];

         if( !nodemark[k2] )
            continue;

         if( k2 > k )
         {
         //   printf("%d, %d \n", k, k2);

            cutoffs[edgecount++] = cliquegraph->cost[e];
         }
      }
   }
}


/** return cost of star */
static inline
SCIP_Real bdkStarGetCost(
   int                   i,                 /**< the node */
   const GRAPH*          g,                 /**< graph data structure */
   const BDK*            bdk                /**< storage */
   )
{
   const int* const star_edges = bdk->star_outedges;
   SCIP_Real costsum = 0.0;
   const int star_degree = bdk->star_degree;

   for( int j = 0; j < star_degree; j++ )
   {
      const int outedge = star_edges[j];
      assert(outedge >= 0 && g->tail[outedge] == i);
      costsum += g->cost[outedge];
   }

   return costsum;
}

#if 0
/** can star be replaced by SD MST? */
static inline
SCIP_Bool bdkStarIsSdMstReplacable(
   SCIP*                 scip,              /**< SCIP data structure */
   int                   i,                 /**< the node */
   const GRAPH*          g,                 /**< graph data structure */
   BDK*                  bdk                /**< storage */
)
{
   PATH mst[STP_BD_MAXDEGREE];

   for( int l = 0; l < STP_BD_MAXDEGREE; l++ )
      mst[l].dist = UNKNOWN;


   /* compute MST  */
   graph_path_exec(scip, auxg, MST_MODE, 0, auxg->cost, mst);

#ifndef NDEBUG
   for( int l = 1; l < STP_BD_MAXDEGREE; l++ )
      assert(mst[l].dist != UNKNOWN);
#endif

   for( int l = 1; l < STP_BD_MAXDEGREE; l++ )
      mstcost += mst[l].dist;
   // combine costs...
}
#endif


/** can star be replaced by SD MST? */
static inline
SCIP_Bool bdkStarIsSdReplacable(
   SCIP*                 scip,              /**< SCIP data structure */
   int                   i,                 /**< the node */
   const GRAPH*          g,                 /**< graph data structure */
   BDK*                  bdk                /**< storage */
)
{
   const SCIP_Real* const sdtreecosts = reduce_sdgraphGetOrderedMstCosts(bdk->sdistance->sdgraph);
   SCIP_Real treecost = 0.0;
   SCIP_Real costsum = bdkStarGetCost(i, g, bdk);
   const int star_degree = bdk->star_degree;

   for( int j = 0; j < star_degree - 1; j++ )
   {
      assert(GE(sdtreecosts[j], 0.0));
      treecost += sdtreecosts[j];
   }

   /* NOTE: sd can be always equal to costsum, because in the case it cannot contain the whole star! */
   if( SCIPisLE(scip, treecost, costsum) )
      return TRUE;


   return FALSE;
}


/** can vertex of degree 3 be replaced? */
static inline
SCIP_Bool bdkStarIsReplacableDeg3(
   SCIP*                 scip,              /**< SCIP data structure */
   int                   i,                 /**< the node */
   const GRAPH*          g,                 /**< graph data structure */
   BDK*                  bdk                /**< storage */
)
{
   const SCIP_Real* const maxcosts = reduce_sdgraphGetOrderedMstCosts(bdk->sdistance->sdgraph);
   const SCIP_Real* const gCost = g->cost;
   const int* const star_edges = bdk->star_outedges;
   const SCIP_Real costsum = gCost[star_edges[0]] + gCost[star_edges[1]] + gCost[star_edges[2]];

   assert(bdk->star_degree == 3);
   assert(GE(maxcosts[0], 0.0) && GE(maxcosts[1], 0.0));

   /* NOTE: sd can be always equal to costsum, because in the case it cannot contain the whole star! */
   if( SCIPisLE(scip, maxcosts[0] + maxcosts[1], costsum) )
      return TRUE;

   // todo here get the sds and combine
   // new method

   if( graph_pseudoAncestors_edgesInConflict(scip, g, bdk->star_outedges, 3) )
      return TRUE;

   return FALSE;
}


/** can star of degree 4 or greater be replaced? */
static inline
SCIP_Bool bdkStarIsReplacableDegGe4(
   SCIP*                 scip,              /**< SCIP data structure */
   int                   i,                 /**< the node */
   const GRAPH*          g,                 /**< graph data structure */
   BDK*                  bdk                /**< storage */
)
{
   const int star_degree = bdk->star_degree;

   assert(g->terms >= star_degree);
   assert(4 <= star_degree && star_degree <= STP_BDKIMP_MAXDEGREE);

   if( bdkStarIsSdReplacable(scip, i, g, bdk) )
      return TRUE;

   if( graph_pseudoAncestors_edgesInConflict(scip, g, bdk->star_outedges, star_degree) )
      return TRUE;

   return FALSE;
}


/** does bdk test for vertex of degree 3
 *  NOTE: one could also use DegGe4 instead, but this method is slightly more efficient */
static inline
SCIP_RETCODE bdkTryDeg3(
   SCIP*                 scip,              /**< SCIP data structure */
   int                   node,              /**< the node */
   GRAPH*                g,                 /**< graph data structure */
   BDK*                  bdk,               /**< storage */
   int*                  nelims             /**< number of eliminations */
)
{
   SCIP_Real sd[3];

   assert(g->grad[node] == 3);
   assert(g->terms >= 3);

   bdk->star_degree = 3;
   bdk->star_outedges = bdk->node_outedges;

   if( bdkStarIsReplacableDeg3(scip, node, g, bdk) )
   {
      SCIP_Real cutoffs[3];
      SCIP_Bool success;
      bdkGetCutoffs(g, bdk, node, cutoffs);

      SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, node, &success));

      assert(success);
      assert(g->grad[node] == 0);

      SCIPdebugMessage("BD3-implied reduction of node %d with SDs: %f %f %f \n ",node, sd[0], sd[1], sd[2]);
      (*nelims)++;
   }

   return SCIP_OKAY;
}


/** does bdk test for vertex of degree 4 or more */
static inline
SCIP_RETCODE bdkTryDegGe4(
   SCIP*                 scip,              /**< SCIP data structure */
   int                   node,              /**< the node */
   GRAPH*                g,                 /**< graph data structure */
   BDK*                  bdk,               /**< storage */
   int*                  nelims             /**< number of eliminations */
)
{
   STAR* const star = bdk->star;
   SCIP_Bool isPseudoDeletable = TRUE;

   assert(4 <= g->grad[node] && g->grad[node] <= STP_BDKIMP_MAXDEGREE);

   reduce_starReset(g, node, star);

   /* check all stars of degree >= 3, as long as they can be ruled-out */
   while ( isPseudoDeletable )
   {
      bdk->star_outedges = reduce_starGetNext(bdk->star, &(bdk->star_degree));

      if( bdk->star_degree == 3 )
      {
         isPseudoDeletable = bdkStarIsReplacableDeg3(scip, node, g, bdk);
      }
      else
      {
         isPseudoDeletable = bdkStarIsReplacableDegGe4(scip, node, g, bdk);
      }

      if( reduce_starAllAreChecked(star) )
         break;
   }

   if( isPseudoDeletable )
   {
      SCIP_Real cutoffs[STP_BDKIMP_MAXNEDGES];
      bdkGetCutoffs(g, bdk, node, cutoffs);

      SCIP_CALL(graph_knot_delPseudo(scip, g, g->cost, cutoffs, NULL, node, &isPseudoDeletable));

      if( isPseudoDeletable )
      {
         SCIPdebugMessage("BDk-implied reduction of node %d \n ", node);
         (*nelims)++;
      }
   }

   return SCIP_OKAY;
}



/*
 * Interface methods
 */


/** bd_k test without given Steiner bottleneck distances */
SCIP_RETCODE reduce_bdk(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph structure */
   int*                  nelims              /**< number of eliminations */
   )
{
   SD* sdistance;

   SCIP_CALL( reduce_sdInit(scip, g, &sdistance) );

   SCIP_CALL( reduce_bdkWithSd(scip, sdistance, g, nelims) );

   reduce_sdFree(scip, &sdistance);

   return SCIP_OKAY;
}


/** bd_k test for given Steiner bottleneck distances */
SCIP_RETCODE reduce_bdkWithSd(
   SCIP*                 scip,               /**< SCIP data structure */
   SD*                   sdistance,          /**< special distances storage */
   GRAPH*                g,                  /**< graph structure */
   int*                  nelims              /**< number of eliminations */
   )
{
   BDK* bdk;
   const int nnodes = graph_get_nNodes(g);
   const int maxdegree = MIN(g->terms, STP_BDKIMP_MAXDEGREE);
   const int nelims_initial = *nelims;

   assert(scip && sdistance && nelims);
   assert(nelims_initial >= 0);
   assert(!graph_pc_isPcMw(g));

   /* NOTE: in the case of g->terms < 3 the method does not work properly, and the case is easy enough to ignore it */
   if( g->terms < 3  )
      return SCIP_OKAY;

   SCIP_CALL( bdkInit(scip, sdistance, &bdk) );
   SCIPdebugMessage("starting BDK-SD Reduction: ");
   graph_mark(g);

   for( int degree = 3; degree <= maxdegree; degree ++ )
   {
      for( int i = 0; i < nnodes; i++ )
      {
         if( g->grad[i] != degree || Is_term(g->term[i]) )
            continue;

         bdkGetNeighborhood(g, i, bdk);
         bdkGetCliqueSds(g, i, bdk);

         if( degree == 3 )
         {
            SCIP_CALL( bdkTryDeg3(scip, i, g, bdk, nelims) );
         }
         else
         {
            SCIP_CALL( bdkTryDegGe4(scip, i, g, bdk, nelims) );
         }
      }
   }

   bdkFree(scip, &bdk);

   if( *nelims > nelims_initial  )
   {
      SCIP_CALL( reduceLevel0(scip, g) );
   }

   return SCIP_OKAY;
}
