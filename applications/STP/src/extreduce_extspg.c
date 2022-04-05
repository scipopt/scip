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


#define STP_EXTSPG_MAXNCHECKS_EQ 8
#define STP_EXTSPG_MAXNCHECKS 4


/**@name Local methods
 *
 * @{
 */



/** ruled out possible by using interior point on path between 'pathstart' and 'pathend'? */
static inline
SCIP_Bool spg4VerticesRuleOut(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   node_pathstart,     /**< start of path */
   int                   node_pathend,       /**< end of path */
   int                   node_other,         /**< other node */
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

   if( GE(pathcost, tree_cost) )
   {
      return FALSE;
   }

   SCIPdebugMessage("%d---%d and %d: number of path nodes: %d (pathcost=%f, tree_cost=%f)\n",
         node_pathstart, node_pathend, node_other, npathnodes, pathcost, tree_cost);

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



/** Can we (peripherally) rule out simple star of degree 3?
 *  Works by using simple Steiner tree  */
static
SCIP_Bool spg3StarNeighborRuleOut(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             starcost,           /**< cost of the star */
   int                   node,               /**< the node */
   SCIP_Bool             allowEquality,      /**< allow equality? */
   const int             neighbors[3],       /**< the neighbors */
   DISTDATA*             distdata            /**< data for distance computations */
)
{
   SCIP_Bool isPseudoDeletable = FALSE;
   const SCIP_Bool isPc = graph_pc_isPc(graph);
   const int maxnchecks = allowEquality ? STP_EXTSPG_MAXNCHECKS_EQ : STP_EXTSPG_MAXNCHECKS;

   for( int i = 0; i < 3 && !isPseudoDeletable; i++ )
   {
      const int n0 = neighbors[i];
      const int n1 = neighbors[(i + 1) % 3];
      const int n2 = neighbors[(i + 2) % 3];
      int edgecount = 0;

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

         if( edgecount++ > maxnchecks )
            break;

         c0 = graph->cost[e];
         c1 = extreduce_distDataGetSdDouble(scip, graph, head, n1, distdata);

         if( isPc && c1 < -0.5 )
         {
            assert(EQ(c1, -1.0));
            continue;
         }

         if( allowEquality ? GT(c0 + c1, starcost) : GE(c0 + c1, starcost) )
            continue;

         c2 = extreduce_distDataGetSdDouble(scip, graph, head, n2, distdata);

         if( isPc && c2 < -0.5 )
         {
            assert(EQ(c2, -1.0));
            continue;
         }

         assert(GE(c1, 0.0) && GE(c2, 0.0));

         if( allowEquality ? LE(c0 + c1 + c2, starcost) : LT(c0 + c1 + c2, starcost) )
         {
            SCIPdebugMessage("STEINER TREE pseudo-delete %d with %f %f \n", node, c0 + c1 + c2, starcost);
            isPseudoDeletable = TRUE;
            break;
         }
      }
   }

   return isPseudoDeletable;
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



/** checks component for possible pseudo-elimination by using simple Steiner tree */
SCIP_RETCODE extreduce_spgCheck3ComponentSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int                   node,               /**< the node */
   const EXTCOMP*        extcomp,            /**< component to be checked */
   SCIP_Bool             allowEquality,      /**< allow equality? */
   DISTDATA*             distdata,           /**< data for distance computations */
   SCIP_Bool*            isPseudoDeletable   /**< is component pseudo-deletable? */
)
{
   int neighbors[3];
   SCIP_Real starcost = 0.0;

   assert(scip && graph && distdata && isPseudoDeletable);
   assert(*isPseudoDeletable == FALSE);
   assert(graph_knot_isInRange(graph, node));
   assert(!Is_term(graph->term[node]) || graph_pc_isPc(graph));
   /* NOTE: because the component has been reverted... */
   assert(extcomp->nextleaves == 1);
   assert(extcomp->ncompedges == 3);
   assert(graph->head[extcomp->compedges[0]] == node);
   assert(graph->tail[extcomp->compedges[0]] == extcomp->comproot);

   neighbors[0] = extcomp->comproot;
   for( int i = 0; i < 2; i++ )
   {
      const int leaf = extcomp->extleaves[i];
      assert(graph_knot_isInRange(graph, leaf));
      assert(leaf != neighbors[0]);

      neighbors[i + 1] = leaf;
   }

   for( int i = 0; i < 3; i++ )
   {
      const int edge = extcomp->compedges[i];
      assert(graph_edge_isInRange(graph, edge));
      assert(EQ(graph->cost[edge], graph->cost[flipedge(edge)]));

      starcost += graph->cost[edge];
   }

   if( graph_pc_isPc(graph) )
   {
	  assert(GE(graph->prize[node], 0.0));
	  starcost -= graph->prize[node];

	  assert(GE(starcost, 0.0));
   }

   *isPseudoDeletable = spg3StarNeighborRuleOut(scip, graph, starcost, node, allowEquality, neighbors, distdata);

   return SCIP_OKAY;
}



/** checks node of degree 3 for possible pseudo-elimination by using simple Steiner tree  */
SCIP_RETCODE extreduce_spgCheck3NodeSimple(
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

   assert(scip && graph && distdata && isPseudoDeletable);
   assert(*isPseudoDeletable == FALSE);
   assert(graph_knot_isInRange(graph, node));
   assert(!Is_term(graph->term[node]));
   assert(graph->grad[node] == 3);

   for( int e = graph->outbeg[node]; e != EAT_LAST; e = graph->oeat[e] )
   {
      assert(EQ(graph->cost[e], graph->cost[flipedge(e)]));

      neighbors[edgecount++] = graph->head[e];
      starcost += graph->cost[e];
   }

   assert(edgecount == 3);

   *isPseudoDeletable = spg3StarNeighborRuleOut(scip, graph, starcost, node, TRUE, neighbors, distdata);

   return SCIP_OKAY;
}

/**@} */
