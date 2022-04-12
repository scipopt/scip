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

/**@file   extreduce_extmstbiased.c
 * @brief  extended-reduction specific biased bottleneck distance MST algorithms for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements MST algorithms for extended reduction techniques for Steiner problems, using
 * the implied bottleneck Steiner distance
 * Furthermore, one can check for tree bottlenecks.
 *
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
//#define STP_DEBUG_EXT

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"





/**@name Local methods
 *
 * @{
 */



/** Can we (peripherally) rule out simple star of degree 3?
 *  Works by using biased distance */
static
SCIP_Bool mst3StarNeighborRuleOut(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real       dists_def[3],       /**< distances between adjacent nodes */
   const SCIP_Real       dists_bias[3],      /**< biased distances between adjacent nodes */
   SCIP_Real             starcost,           /**< cost of the star */
   SCIP_Bool             allowEquality       /**< allow equality? */
)
{
   if( allowEquality )
   {
      if( SCIPisLE(scip, dists_def[0] + dists_bias[1], starcost) )
         return TRUE;
      if( SCIPisLE(scip, dists_def[1] + dists_bias[0], starcost) )
         return TRUE;

      if( SCIPisLE(scip, dists_def[0] + dists_bias[2], starcost) )
         return TRUE;
      if( SCIPisLE(scip, dists_def[2] + dists_bias[0], starcost) )
         return TRUE;

      if( SCIPisLE(scip, dists_def[1] + dists_bias[2], starcost) )
         return TRUE;
      if( SCIPisLE(scip, dists_def[2] + dists_bias[1], starcost) )
         return TRUE;
   }
   else
   {
      if( SCIPisLT(scip, dists_def[0] + dists_bias[1], starcost) )
         return TRUE;
      if( SCIPisLT(scip, dists_def[1] + dists_bias[0], starcost) )
         return TRUE;

      if( SCIPisLT(scip, dists_def[0] + dists_bias[2], starcost) )
         return TRUE;
      if( SCIPisLT(scip, dists_def[2] + dists_bias[0], starcost) )
         return TRUE;

      if( SCIPisLT(scip, dists_def[1] + dists_bias[2], starcost) )
         return TRUE;
      if( SCIPisLT(scip, dists_def[2] + dists_bias[1], starcost) )
         return TRUE;
   }

   return FALSE;
}


/** gets biased SDs between leaves */
static inline
void mst3LeafTreeGetSds(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Real             sds_def[3],         /**< array to store the SDs */
   SCIP_Real             sds_bias[3]         /**< array to store the SDs */
   )
{
   const int* const leaves = extdata->tree_leaves;
   const REDDATA* const reddata = extdata->reddata;
   const MLDISTS* const sds_vertical = reddata->sds_vertical;
   const MLDISTS* const sdsbias_vertical = reddata->sdsbias_vertical;
   const int* const extstack_data = extdata->extstack_data;
   const int stackpos = extStackGetPosition(extdata);
   const int topedges_start = extStackGetTopOutEdgesStart(extdata, stackpos);
   const int topedges_end = extStackGetTopOutEdgesEnd(extdata, stackpos);
   const int ntopleaves = topedges_end - topedges_start;
   const int topedge1 = extstack_data[topedges_start];
   int topsibling1 = graph->head[topedge1];

   assert(1 <= ntopleaves && ntopleaves <= 2);

   if( ntopleaves == 2 )
   {
      const MLDISTS* const sds_horizontal = reddata->sds_horizontal;
      const MLDISTS* const sdsbias_horizontal = reddata->sdsbias_horizontal;
      const int topedge2 = extstack_data[topedges_start + 1];
      int topsibling2 = graph->head[topedge2];

      assert(graph->tail[topedge1] == graph->tail[topedge2]);
      assert(topsibling1 != topsibling2);
      assert(topsibling1 != leaves[0] && topsibling2 != leaves[0]);

      if( topsibling1 == leaves[2]  )
      {
         SWAP_INTS(topsibling1, topsibling2);
      }

      assert(topsibling1 == leaves[1] && topsibling2 == leaves[2]);

      /* SDs from root to siblings */
      sds_def[0] = extreduce_mldistsTopTargetDists(sds_vertical, topsibling1)[0];
      sds_bias[0] = extreduce_mldistsTopTargetDists(sdsbias_vertical, topsibling1)[0];
      sds_def[1] = extreduce_mldistsTopTargetDists(sds_vertical, topsibling2)[0];
      sds_bias[1] = extreduce_mldistsTopTargetDists(sdsbias_vertical, topsibling2)[0];

      /* SDs between siblings */
      sds_def[2] = extreduce_mldistsTopTargetDist(sds_horizontal, topsibling1, topsibling2);
      sds_bias[2] = extreduce_mldistsTopTargetDist(sdsbias_horizontal, topsibling1, topsibling2);
   }
   else
   {
      DISTDATA* const distdata_default = extdata->distdata;
      DISTDATA* const distdata_biased = extdata->distdata_biased;
      assert(topsibling1 == leaves[2]);

      sds_def[0] = extreduce_distDataGetSdDouble(scip, graph, leaves[0], leaves[1], distdata_default);
      sds_bias[0] = extreduce_distDataGetSdDouble(scip, graph, leaves[0], leaves[1], distdata_biased);

      sds_def[1] = extreduce_mldistsTopTargetDists(sds_vertical, topsibling1)[0];
      sds_bias[1] = extreduce_mldistsTopTargetDists(sdsbias_vertical, topsibling1)[0];

      sds_def[2] = extreduce_mldistsTopTargetDists(sds_vertical, topsibling1)[1];
      sds_bias[2] = extreduce_mldistsTopTargetDists(sdsbias_vertical, topsibling1)[1];
   }
}


/**@} */

/**@name Interface methods
 *
 * @{
 */



/** can current 3-leaf tree be ruled-out? */
SCIP_Bool extreduce_mstbiased3LeafTreeRuleOut(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             tree_cost,          /**< tree cost */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Bool ruledOut = FALSE;
   SCIP_Real dists_def[3];
   SCIP_Real dists_bias[3];

   assert(scip && graph && extdata);
   assert(GE(tree_cost, 0.0));
   assert(extdata->tree_nleaves == 3);

   mst3LeafTreeGetSds(scip, graph, extdata, dists_def, dists_bias);

#ifndef NDEBUG
   {
      DISTDATA* distdata_default = extdata->distdata;
      DISTDATA* distdata_biased = extdata->distdata_biased;
      const int* const leaves = extdata->tree_leaves;

      assert(EQ(dists_def[0], extreduce_distDataGetSdDouble(scip, graph, leaves[0], leaves[1], distdata_default)));
      assert(EQ(dists_def[1], extreduce_distDataGetSdDouble(scip, graph, leaves[0], leaves[2], distdata_default)));
      assert(EQ(dists_def[2], extreduce_distDataGetSdDouble(scip, graph, leaves[1], leaves[2], distdata_default)));

      assert(EQ(dists_bias[0], extreduce_distDataGetSdDouble(scip, graph, leaves[0], leaves[1], distdata_biased)));
      assert(EQ(dists_bias[1], extreduce_distDataGetSdDouble(scip, graph, leaves[0], leaves[2], distdata_biased)));
      assert(EQ(dists_bias[2], extreduce_distDataGetSdDouble(scip, graph, leaves[1], leaves[2], distdata_biased)));
   }
#endif

   if( mst3StarNeighborRuleOut(scip, dists_def, dists_bias, tree_cost, FALSE) )
   {
      SCIPdebugMessage("biased 3-leaf MST rule-out \n");
      ruledOut = TRUE;
   }

   return ruledOut;
}


/** checks node of degree 3 for possible pseudo-elimination by using bias bottleneck Steiner distance */
SCIP_RETCODE extreduce_mstbiasedCheck3NodeSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   node,               /**< the node */
   DISTDATA*             distdata_default,   /**< data for distance computations */
   DISTDATA*             distdata_biased,    /**< data for distance computations */
   SCIP_Bool*            isPseudoDeletable   /**< is node pseudo-deletable? */
)
{
   SCIP_Real dists_def[3];
   SCIP_Real dists_bias[3];
   int neighbors[3];
   int edgecount = 0;
   SCIP_Real starcost = 0.0;

   assert(scip && graph && distdata_default && distdata_biased && isPseudoDeletable);
   assert(distdata_biased->sdistdata && distdata_biased->sdistdata->sdprofit);
   assert(*isPseudoDeletable == FALSE);
   assert(graph_knot_isInRange(graph, node));
   assert(!Is_term(graph->term[node]));
   assert(graph->grad[node] == 3);
   assert(!graph_pc_isPcMw(graph));
   assert(EQ(reduce_sdprofitGetProfit(distdata_biased->sdistdata->sdprofit, node, -1, -1), 0.0));

   for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
   {
      assert(EQ(graph->cost[e], graph->cost[flipedge(e)]));

      neighbors[edgecount++] = graph->head[e];
      starcost += graph->cost[e];
   }

   assert(edgecount == 3);

   dists_def[0] = extreduce_distDataGetSdDouble(scip, graph, neighbors[0], neighbors[1], distdata_default);
   dists_def[1] = extreduce_distDataGetSdDouble(scip, graph, neighbors[0], neighbors[2], distdata_default);
   dists_def[2] = extreduce_distDataGetSdDouble(scip, graph, neighbors[1], neighbors[2], distdata_default);

   dists_bias[0] = extreduce_distDataGetSdDouble(scip, graph, neighbors[0], neighbors[1], distdata_biased);
   dists_bias[1] = extreduce_distDataGetSdDouble(scip, graph, neighbors[0], neighbors[2], distdata_biased);
   dists_bias[2] = extreduce_distDataGetSdDouble(scip, graph, neighbors[1], neighbors[2], distdata_biased);

   *isPseudoDeletable = mst3StarNeighborRuleOut(scip, dists_def, dists_bias, starcost, TRUE);

   return SCIP_OKAY;
}

/**@} */
