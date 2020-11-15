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
   const SCIP_Real       dists_def[3],
   const SCIP_Real       dists_bias[3],
   SCIP_Real             starcost,           /**< cost of the star */
   int                   node,               /**< the node */
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



/**@} */

/**@name Interface methods
 *
 * @{
 */



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

   *isPseudoDeletable = mst3StarNeighborRuleOut(scip, dists_def, dists_bias, starcost, node, TRUE);

   return SCIP_OKAY;
}
