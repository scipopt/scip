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

/**@file   extreduce_extspg.c
 * @brief  extended-reduction specific SPG algorithms for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements special distance Steiner tree algorithms for extended reduction techniques for Steiner problems.
 * Allows one to efficiently compute and store special distance (SD) Steiner trees between the leaves of extension tree.
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



/** ruled out possible by using interior point on path between 'pathstart' and 'pathend'? */
static inline
SCIP_Bool spg4VerticesRuleOut(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   node_pathstart,     /**< */
   int                   node_pathend,       /**< */
   int                   node_other,         /**< */
   SCIP_Real             tree_cost,          /**< tree cost */
   int*                  pathnodes,          /**< buffer */
   EXTDATA*              extdata             /**< extension data */
)
{
   DISTDATA* const distdata = extdata->distdata;
   const int* const tree_deg = extdata->tree_deg;
   int npathnodes;
   const SCIP_Real pathcost
      = extreduce_distDataGetSp(scip, graph, node_pathstart, node_pathend, pathnodes, &npathnodes, distdata);

   SCIP_Bool isRuledOut = FALSE;

   assert(tree_deg[node_pathstart] > 0);
   assert(tree_deg[node_pathend] > 0);
   assert(tree_deg[node_other] > 0);

   for( int i = 0; i < npathnodes; i++ )
   {
      SCIP_Real dist2;
      const int node_middle = pathnodes[i];

      assert(GE(pathcost, 0.0));

      if( tree_deg[node_middle] > 0 )
         continue;

      dist2 = extreduce_distDataGetSdDouble(scip, graph, node_middle, node_other, distdata);

      if( dist2 < -0.5 )
      {
         assert(EQ(dist2, -1.0));
         assert(graph_pc_isPc(graph));

         continue;
      }

      if( LT(pathcost + dist2, tree_cost) )
      {
         SCIPdebugMessage("STEINER TREE 3-leaf tree rule-out with %f < %f \n", pathcost + dist2, tree_cost);

         isRuledOut = TRUE;
         break;
      }
   }

   return isRuledOut;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */


/** can current 3-leaf tree be ruled-out? */
SCIP_Bool extreduce_spg3LeafTreeRuleOut(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             tree_cost,          /**< tree cost */
   EXTDATA*              extdata             /**< extension data */
)
{
   int* pathnodes;
   const int* const leaves = extdata->tree_leaves;

   assert(scip && graph && extdata);
   assert(extdata->tree_nleaves == 3);
   assert(extInitialCompIsEdge(extdata));
   assert(GE(tree_cost, 0.0));

   SCIPdebugMessage("try 3-leaf reduction by SD SPG \n");

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &pathnodes, graph->knots) );

   if( spg4VerticesRuleOut(scip, graph, leaves[0], leaves[1], leaves[2], tree_cost, pathnodes, extdata) )
   {
      SCIPfreeBufferArray(scip, &pathnodes);
      return TRUE;
   }

   if( spg4VerticesRuleOut(scip, graph, leaves[0], leaves[2], leaves[1], tree_cost, pathnodes, extdata) )
   {
      SCIPfreeBufferArray(scip, &pathnodes);
      return TRUE;
   }

   if( spg4VerticesRuleOut(scip, graph, leaves[1], leaves[2], leaves[0], tree_cost, pathnodes, extdata) )
   {
      SCIPfreeBufferArray(scip, &pathnodes);
      return TRUE;
   }

   SCIPfreeBufferArray(scip, &pathnodes);

   return FALSE;
}


/** check node for possible pseudo-elimination by using simple Steiner tree  */
SCIP_RETCODE extreduce_spgCheckNodeSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   node,               /**< the node */
   DISTDATA*             distdata,           /**< data for distance computations */
   SCIP_Bool*            isPseudoDeletable   /**< is node pseudo-deletable? */
)
{
   int neighbors[3];
   int edgecount = 0;
   SCIP_Real starcost = 0.0;
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   assert(scip && graph && distdata && isPseudoDeletable);
   assert(*isPseudoDeletable == FALSE);
   assert(!Is_term(graph->term[node]));

   for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
   {
      neighbors[edgecount++] = graph->head[e];
      starcost += graph->cost[e];
   }

   assert(edgecount == 3);

   for( int i = 0; i < 3 && !(*isPseudoDeletable); i++ )
   {
      const int n0 = neighbors[i];
      const int n1 = neighbors[(i + 1) % 3];
      const int n2 = neighbors[(i + 2) % 3];
      edgecount = 0;

      assert(n1 != n0 && n2 != n0);

      for( int e = graph->outbeg[n0]; e != EAT_LAST; e = graph->oeat[e] )
      {
         const int head = graph->head[e];
         SCIP_Real c0;
         SCIP_Real c1;
         SCIP_Real c2;

         if( head == node )
            continue;

         if( head == n1 || head == n2 )
            continue;

         if( isPc && Is_pseudoTerm(graph->term[head]) )
            continue;

         if( edgecount++ > 8 )
            break;

         c0 = graph->cost[e];
         c1 = extreduce_distDataGetSdDouble(scip, graph, head, n1, distdata);

         if( isPc && c1 < -0.5 )
         {
            assert(EQ(c1, -1.0));
            continue;
         }

         if( GT(c0 + c1, starcost) )
            continue;

         c2 = extreduce_distDataGetSdDouble(scip, graph, head, n2, distdata);
       //  c1 = (head == n1) ? 0.0 : extreduce_distDataGetSdDouble(scip, graph, head, n1, distdata);
       //  c2 = (head == n2) ? 0.0 : extreduce_distDataGetSdDouble(scip, graph, head, n2, distdata);

         if( isPc && c2 < -0.5 )
         {
            assert(EQ(c2, -1.0));
            continue;
         }

         assert(GE(c1, 0.0) && GE(c2, 0.0));

         if( LE(c0 + c1 + c2, starcost) )
         {
            SCIPdebugMessage("STEINER TREE pseudo-delete %d with %f %f \n", node, c0 + c1 + c2, starcost);
            *isPseudoDeletable = TRUE;
            break;
         }
      }
   }

   return SCIP_OKAY;
}
