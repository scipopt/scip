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

/**@file   extreduce_bottleneck.c
 * @brief  extended-reduction specific tree bottleneck algorithms for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file implements tree bottleneck algorithms for extended reduction techniques for Steiner problems.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
// #define SCIP_DEBUG
//#define STP_DEBUG_EXT

#include <stdio.h>
#include <assert.h>
#include "graph.h"
#include "portab.h"
#include "extreduce.h"
#include "extreducedefs.h"



/**@name Local methods
 *
 * @{
 */


/** computes the tree bottleneck between vertices in the current tree,
 * for which vertex_pathmarked root path has been marked already */
static
SCIP_Real bottleneckGetDist(
   const GRAPH*          graph,              /**< graph data structure */
   const EXTDATA*        extdata,            /**< extension data */
#ifndef NDEBUG
   int                   vertex_pathmarked,  /**< vertex with marked rootpath */
#endif
   int                   vertex_unmarked     /**< second vertex */
   )
{
   const SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const SCIP_Real* const parentEdgeCost = extdata->tree_parentEdgeCost;
   const int* const parentNode = extdata->tree_parentNode;
   int* const tree_deg = extdata->tree_deg;
   SCIP_Real bottleneck;
   const int tree_root = extdata->tree_root;
   int currentNode;

   assert(bottleneckDist_node && parentEdgeCost && parentNode);
   assert(bottleneckDist_node[vertex_pathmarked] == -1.0 || vertex_pathmarked == tree_root);
   assert(bottleneckDist_node[vertex_unmarked] == -1.0 || vertex_unmarked == tree_root || tree_deg[vertex_unmarked] > 1);
   assert(bottleneckDist_node[tree_root] >= 0.0);
   assert(vertex_pathmarked != vertex_unmarked);

   /* go down from vertex_unmarked up to lowest common ancestor with vertex_pathmarked  */
   bottleneck = 0.0;

   if( vertex_unmarked == tree_root )
   {
      currentNode = vertex_unmarked;
   }
   else
   {
      SCIP_Real bottleneck_local = 0.0;
      const SCIP_Bool isPc = graph_pc_isPc(graph);

      assert(parentNode[vertex_unmarked] >= 0);

      for( currentNode = vertex_unmarked; bottleneckDist_node[currentNode] < -0.5; currentNode = parentNode[currentNode] )
      {
         assert(tree_deg[currentNode] >= 0 && parentEdgeCost[currentNode] >= 0.0);
         assert(bottleneckDist_node[currentNode] == -1.0);
         assert(currentNode != vertex_pathmarked);

         if( tree_deg[currentNode] == 2 )
         {
            bottleneck_local += parentEdgeCost[currentNode];
            if( isPc && Is_term(graph->term[currentNode]) )
            {
               assert(graph_pc_termIsNonLeafTerm(graph, currentNode) && graph->prize[currentNode] > 0.0);
               bottleneck_local -= graph->prize[currentNode];
            }
         }
         else
            bottleneck_local = parentEdgeCost[currentNode];

         if( bottleneck < bottleneck_local )
            bottleneck = bottleneck_local;

         assert(parentNode[currentNode] >= 0 && parentNode[currentNode] != vertex_unmarked);
      }
   }

   bottleneck = MAX(bottleneck, bottleneckDist_node[currentNode]);

   return bottleneck;
}


/** helper */
static inline
void bottleneckMarkEqualityPath(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   path_start,         /**< vertex to start from */
   int                   path_end,           /**< vertex to end at */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Bool* const edges_isEqForbidden = extdata->sdeq_edgesIsForbidden;
   const int* const parentNode = extdata->tree_parentNode;

   assert(edges_isEqForbidden);
   assert(path_start != path_end);
   assert(graph_knot_isInRange(graph, path_start));
   assert(graph_knot_isInRange(graph, path_end));

   for( int currentNode = path_start; currentNode != path_end; currentNode = parentNode[currentNode] )
   {
      int e;
      const int parent = parentNode[currentNode];

      assert(graph_knot_isInRange(graph, parent));

      for( e = graph->outbeg[parent]; e != EAT_LAST; e = graph->oeat[e] )
      {
         if( graph->head[e] == currentNode )
         {
            assert(EQ(graph->cost[e], extdata->tree_parentEdgeCost[currentNode]));
            if( !edges_isEqForbidden[e / 2] )
            {
#ifdef SCIP_DEBUG
               SCIPdebugMessage("forbid equality edge: ");
               graph_edge_printInfo(graph, e);
#endif
               edges_isEqForbidden[e / 2] = TRUE;
               extdata->sdeq_hasForbiddenEdges = TRUE;
               StpVecPushBack(scip, extdata->sdeq_resetStack, e / 2);
            }
            break;
         }
      }

      assert(e != EAT_LAST);
   }
}


/** markes bottleneck edges used for equality rule-out */
static inline
void bottleneckMarkEqualityEdges(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   SCIP_Real             dist_eq,            /**< distance that was used for equality rule-out */
   int                   vertex_pathmarked,  /**< vertex with marked rootpath */
   int                   vertex_unmarked,    /**< second vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const SCIP_Real* const parentEdgeCost = extdata->tree_parentEdgeCost;
   const int* const parentNode = extdata->tree_parentNode;
   int* const tree_deg = extdata->tree_deg;
   SCIP_Real bottleneck_local;
   const int tree_root = extdata->tree_root;
   int ancestor = UNKNOWN;
   int bottleneck_start;
   const SCIP_Bool isPc = graph_pc_isPc(graph);

   assert(bottleneckDist_node && parentEdgeCost && parentNode);
   assert(bottleneckDist_node[vertex_pathmarked] == -1.0 || vertex_pathmarked == tree_root);
   assert(bottleneckDist_node[vertex_unmarked] == -1.0 || vertex_unmarked == tree_root || tree_deg[vertex_unmarked] > 1);
   assert(bottleneckDist_node[tree_root] >= 0.0);
   assert(vertex_pathmarked != vertex_unmarked);

   /* 1. go down from vertex_unmarked to lowest common ancestor with vertex_pathmarked */

   if( vertex_unmarked == tree_root )
   {
      ancestor = vertex_unmarked;
   }
   else
   {
      int currentNode;
      assert(parentNode[vertex_unmarked] >= 0);
      bottleneck_start = UNKNOWN;
      bottleneck_local = 0.0;

      for( currentNode = vertex_unmarked; bottleneckDist_node[currentNode] < -0.5; currentNode = parentNode[currentNode] )
      {
         assert(tree_deg[currentNode] >= 0 && parentEdgeCost[currentNode] >= 0.0);
         assert(EQ(bottleneckDist_node[currentNode], -1.0));
         assert(currentNode != vertex_pathmarked);

         if( tree_deg[currentNode] == 2 )
         {
            bottleneck_local += parentEdgeCost[currentNode];
            if( isPc && Is_term(graph->term[currentNode]) )
            {
               assert(graph_pc_termIsNonLeafTerm(graph, currentNode) && graph->prize[currentNode] > 0.0);
               bottleneck_local -= graph->prize[currentNode];
            }
         }
         else
         {
            bottleneck_start = currentNode;
            bottleneck_local = parentEdgeCost[currentNode];
         }

         if( EQ(bottleneck_local, dist_eq) )
         {
            assert(parentNode[currentNode] >= 0);
            bottleneckMarkEqualityPath(scip, graph, bottleneck_start, parentNode[currentNode], extdata);

            return;
         }

         assert(parentNode[currentNode] >= 0 && parentNode[currentNode] != vertex_unmarked);
      }

      ancestor = currentNode;
      assert(GE(bottleneckDist_node[ancestor], 0.0));
   }


   /* 2. go down from vertex_marked to ancestor */

   assert(parentNode[vertex_pathmarked] >= 0);
   assert(ancestor != UNKNOWN);
   bottleneck_start = UNKNOWN;
   bottleneck_local = 0.0;

   for( int currentNode = vertex_pathmarked; currentNode != ancestor; currentNode = parentNode[currentNode] )
   {
      assert(tree_deg[currentNode] >= 0 && parentEdgeCost[currentNode] >= 0.0);
      assert(currentNode != vertex_unmarked);

      if( tree_deg[currentNode] == 2 )
      {
         bottleneck_local += parentEdgeCost[currentNode];
         if( isPc && Is_term(graph->term[currentNode]) )
         {
            assert(graph_pc_termIsNonLeafTerm(graph, currentNode) && graph->prize[currentNode] > 0.0);
            bottleneck_local -= graph->prize[currentNode];
         }
      }
      else
      {
         bottleneck_start = currentNode;
         bottleneck_local = parentEdgeCost[currentNode];
      }

      if( EQ(bottleneck_local, dist_eq) )
      {
         // todo not quite sure whether this is correct
         if( bottleneck_start == UNKNOWN  )
         {
            assert(extIsAtInitialGenStar(extdata));
            bottleneck_start = vertex_pathmarked;
         }

         assert(parentNode[currentNode] >= 0);
         bottleneckMarkEqualityPath(scip, graph, bottleneck_start, parentNode[currentNode], extdata);

         return;
      }

      assert(parentNode[currentNode] >= 0 && parentNode[currentNode] != vertex_unmarked);
   }

   assert(0 && "should never arrive here!");
}


/** markes single bottleneck edge used for equality rule-out */
static inline
void bottleneckMarkEqualityEdge(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   int                   edge,               /**< the edge to mark */
   EXTDATA*              extdata             /**< extension data */
)
{
   SCIP_Bool* const edges_isEqForbidden = extdata->sdeq_edgesIsForbidden;

   assert(edges_isEqForbidden);
   assert(graph_edge_isInRange(g, edge));

   if( !edges_isEqForbidden[edge / 2] )
   {
      edges_isEqForbidden[edge / 2] = TRUE;
      extdata->sdeq_hasForbiddenEdges = TRUE;
      StpVecPushBack(scip, extdata->sdeq_resetStack, edge / 2);
   }
}


/** helper to check the case of quality */
static inline
SCIP_Bool bottleneckIsEqualityDominated(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          g,                  /**< graph data structure */
   SCIP_Real             dist_eq,            /**< critical distance */
   int                   edge_forbidden,     /**< forbidden edge */
   int                   vertex1,            /**< first vertex */
   int                   vertex2,            /**< second vertex */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Real sd_eq = extreduce_distDataGetSdDoubleForbiddenEq(scip, g, dist_eq,
         edge_forbidden, vertex1, vertex2, extdata);

   if( sd_eq < -0.5 )
      return FALSE;

   assert(GE(sd_eq, dist_eq));

   if( LE(sd_eq, dist_eq) )
   {
      assert(EQ(sd_eq, dist_eq));
      return TRUE;
   }

   return FALSE;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */



/** marks bottleneck array on path to tree root */
void extreduce_bottleneckMarkRootPath(
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex,             /**< vertex to start from */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const SCIP_Real* const parentEdgeCost = extdata->tree_parentEdgeCost;
   const int* const parentNode = extdata->tree_parentNode;
   const int* const tree_deg = extdata->tree_deg;
   const int tree_root = extdata->tree_root;

   assert(bottleneckDist_node && parentEdgeCost && parentNode && tree_deg);
   assert(vertex >= 0 && vertex < graph->knots);
   assert(bottleneckDist_node[vertex] == -1.0);
   assert(bottleneckDist_node[tree_root] == -1.0);

   if( vertex == tree_root )
   {
      bottleneckDist_node[vertex] = 0.0;
   }
   else
   {
      /* go down from vertex */

      SCIP_Real bottleneck = 0.0;
      SCIP_Real bottleneck_local = 0.0;
      int childNode = vertex;
      int currentNode = parentNode[vertex];
      const SCIP_Bool isPc = graph_pc_isPc(graph);

      assert(currentNode != -1);
      assert(!extInitialCompIsEdge(extdata) || tree_deg[childNode] == 1);
      assert(childNode == extdata->tree_starcenter || extInitialCompIsGenStar(extdata) || tree_deg[childNode] == 1);

      while( currentNode != -1 )
      {
         assert(currentNode >= 0 && tree_deg[currentNode] >= 0);
         assert(parentEdgeCost[childNode] >= 0.0 && bottleneckDist_node[currentNode] == -1.0);
         assert(currentNode != vertex);
         assert(!isPc || !graph_pc_knotIsDummyTerm(graph, currentNode));

         if( tree_deg[childNode] == 2 )
         {
            bottleneck_local += parentEdgeCost[childNode];
            if( isPc && Is_term(graph->term[childNode]) )
            {
               assert(graph_pc_termIsNonLeafTerm(graph, childNode) && graph->prize[childNode] > 0.0);
               bottleneck_local -= graph->prize[childNode];
            }
         }
         else
            bottleneck_local = parentEdgeCost[childNode];

         if( bottleneck < bottleneck_local )
            bottleneck = bottleneck_local;

         bottleneckDist_node[currentNode] = bottleneck;
         childNode = currentNode;
         currentNode = parentNode[currentNode];
      }

      assert(childNode == tree_root);
   }
}


/** unmarks bottleneck array on path to tree root */
void extreduce_bottleneckUnmarkRootPath(
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex,             /**< vertex to start from */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real* const bottleneckDist_node = extdata->tree_bottleneckDistNode;
   const int* const parentNode = extdata->tree_parentNode;
   const int tree_root = extdata->tree_root;

   assert(extdata && bottleneckDist_node && parentNode);
   assert(bottleneckDist_node[vertex] == -1.0 || vertex == tree_root);
   assert(bottleneckDist_node[tree_root] >= 0.0);

   if( vertex == tree_root )
   {
      bottleneckDist_node[vertex] = -1.0;
      assert(parentNode[vertex] == -1);
   }
   else
   {
      assert(parentNode[vertex] >= 0);
   }

   /* go down from vertex and reset bottleneckDist_node */
   for( int currentNode = parentNode[vertex]; currentNode != -1; currentNode = parentNode[currentNode]  )
   {
      assert(currentNode >= 0);
      assert(extdata->tree_deg[currentNode] >= 0);
      assert(bottleneckDist_node[currentNode] >= 0.0);

      bottleneckDist_node[currentNode] = -1.0;
   }

   assert(bottleneckDist_node[tree_root] == -1.0);
}



/** Does a special distance approximation dominate the tree bottleneck distance between
 *  vertex_pathmarked and vertex_unmarked in the current tree.
 *  NOTE: makes additional checks in case of equality */
SCIP_Bool extreduce_bottleneckIsDominated(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex_pathmarked,  /**< vertex for which bottleneck path to root has been marked */
   int                   vertex_unmarked,    /**< second vertex */
   SCIP_Real             specialDist,        /**< best computed special distance approximation (-1.0 if unknown) */
   int                   edge_forbidden,     /**< forbidden edge */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real bottleneckDist;
   const SCIP_Bool hasSpecialDist = extSdIsNonTrivial(specialDist);

   assert(graph_knot_isInRange(graph, vertex_pathmarked));
   assert(graph_knot_isInRange(graph, vertex_unmarked));

   if( !hasSpecialDist || vertex_pathmarked == vertex_unmarked )
   {
      return FALSE;
   }

#ifndef NDEBUG
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_pathmarked, vertex_unmarked);
#else
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_unmarked);
#endif

   SCIPdebugMessage("domination test %d->%d: sd=%f bottleneck=%f \n", vertex_pathmarked, vertex_unmarked, specialDist, bottleneckDist);

   if( LT(specialDist, bottleneckDist) )
   {
      return TRUE;
   }
   else if( LE(specialDist, bottleneckDist) )
   {
#ifdef EXT_DOUBLESD_ALWAYS
      assert(EQ(extreduce_extGetSdDouble(scip, graph, vertex_pathmarked, vertex_unmarked, extdata), specialDist));
#else
      assert(LE(extreduce_extGetSdDouble(scip, graph, vertex_pathmarked, vertex_unmarked, extdata), specialDist));
#endif

      if( bottleneckIsEqualityDominated(scip, graph, specialDist, edge_forbidden,
         vertex_pathmarked, vertex_unmarked, extdata) )
      {
         SCIPdebugMessage("...ruled out with equality! \n");
         bottleneckMarkEqualityEdges(scip, graph, specialDist, vertex_pathmarked, vertex_unmarked, extdata);

         return TRUE;
      }
   }

   return FALSE;
}


/** As above, but for biased special distance */
SCIP_Bool extreduce_bottleneckIsDominatedBiased(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   vertex_pathmarked,  /**< vertex for which bottleneck path to root has been marked */
   int                   vertex_unmarked,    /**< second vertex */
   SCIP_Real             specialDistBiased,  /**< best computed special distance approximation (-1.0 if unknown) */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real bottleneckDist;
   const SCIP_Bool hasSpecialDist = extSdIsNonTrivial(specialDistBiased);

   assert(graph_knot_isInRange(graph, vertex_pathmarked));
   assert(graph_knot_isInRange(graph, vertex_unmarked));

   if( !hasSpecialDist || vertex_pathmarked == vertex_unmarked )
   {
      return FALSE;
   }

#ifndef NDEBUG
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_pathmarked, vertex_unmarked);
#else
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_unmarked);
#endif

   SCIPdebugMessage("biased domination test %d->%d: sd=%f bottleneck=%f \n", vertex_pathmarked, vertex_unmarked, specialDistBiased, bottleneckDist);

   if( LT(specialDistBiased, bottleneckDist) )
   {
      return TRUE;
   }

   return FALSE;
}


/** Does a special distance approximation dominate the tree bottleneck distance of
 *  extension edge (i.e. its edge cost) or bottleneck distance between vertex_pathmarked
 *  and vertex_unmarked in the current tree.
 *  NOTE: makes additional checks in case of equality */
SCIP_Bool extreduce_bottleneckWithExtedgeIsDominated(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   extedge,            /**< edge along which we want to extend the tree */
   int                   vertex_pathmarked,  /**< vertex for which bottleneck path to root has been marked */
   int                   vertex_unmarked,    /**< second vertex */
   SCIP_Real             specialDist,        /**< best computed special distance approximation (-1.0 if unknown) */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real bottleneckDist;
   const SCIP_Bool hasSpecialDist = extSdIsNonTrivial(specialDist);

   assert(graph_edge_isInRange(graph, extedge));
   assert(vertex_pathmarked == graph->tail[extedge]);

   if( !hasSpecialDist )
      return FALSE;

   if( LT(specialDist, graph->cost[extedge]) )
   {
      return TRUE;
   }
   else if( LE(specialDist, graph->cost[extedge]) )
   {
      const int vertex1 = graph->head[extedge];

      if( bottleneckIsEqualityDominated(scip, graph, specialDist, extedge,
         vertex1, vertex_unmarked, extdata) )
      {
         bottleneckMarkEqualityEdge(scip, graph, extedge, extdata);
         SCIPdebugMessage("...ruled out with equality by single edge ! \n");

         return TRUE;
      }
   }

   if( vertex_pathmarked == vertex_unmarked )
      return FALSE;

#ifndef NDEBUG
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_pathmarked, vertex_unmarked);
#else
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_unmarked);
#endif

   SCIPdebugMessage("extedge domination test %d->%d: sd=%f bottleneck=%f \n", vertex_pathmarked, vertex_unmarked, specialDist, bottleneckDist);

   if( LT(specialDist, bottleneckDist) )
   {
      return TRUE;
   }
   else if( LE(specialDist, bottleneckDist) )
   {
      const int vertex1 = graph->head[extedge];

      assert(vertex1 != vertex_unmarked);
      assert(vertex1 != vertex_pathmarked);
#ifdef EXT_DOUBLESD_ALWAYS
      assert(EQ(extreduce_extGetSdDouble(scip, graph, vertex1, vertex_unmarked, extdata), specialDist));
#else
      assert(LE(extreduce_extGetSdDouble(scip, graph, vertex1, vertex_unmarked, extdata), specialDist));
#endif

      if( bottleneckIsEqualityDominated(scip, graph, specialDist, extedge,
         vertex1, vertex_unmarked, extdata) )
      {
         bottleneckMarkEqualityEdges(scip, graph, specialDist, vertex_pathmarked, vertex_unmarked, extdata);
         SCIPdebugMessage("...ruled out with equality! \n");

         return TRUE;
      }
   }

   return FALSE;
}


/** as above, but for biased special distance */
SCIP_Bool extreduce_bottleneckWithExtedgeIsDominatedBiased(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   extedge,            /**< edge along which we want to extend the tree */
   int                   vertex_pathmarked,  /**< vertex for which bottleneck path to root has been marked */
   int                   vertex_unmarked,    /**< second vertex */
   SCIP_Real             specialDistBiased,  /**< best computed biased special distance approximation (-1.0 if unknown) */
   EXTDATA*              extdata             /**< extension data */
   )
{
   SCIP_Real bottleneckDist;
   const SCIP_Bool hasSpecialDist = extSdIsNonTrivial(specialDistBiased);

   assert(graph_edge_isInRange(graph, extedge));
   assert(vertex_pathmarked == graph->tail[extedge]);

   if( !hasSpecialDist )
      return FALSE;

   if( LT(specialDistBiased, graph->cost[extedge]) )
      return TRUE;

   if( vertex_pathmarked == vertex_unmarked )
      return FALSE;

#ifndef NDEBUG
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_pathmarked, vertex_unmarked);
#else
   bottleneckDist = bottleneckGetDist(graph, extdata, vertex_unmarked);
#endif

   SCIPdebugMessage("extedge biased domination test %d->%d: sd=%f bottleneck=%f \n", vertex_pathmarked, vertex_unmarked, specialDistBiased, bottleneckDist);

   if( LT(specialDistBiased, bottleneckDist) )
      return TRUE;

   return FALSE;
}


/** Does a special distance approximation dominate the tree bottleneck distance between
 *  vertex_pathmarked and vertex_unmarked in the current tree?
 *  NOTE: makes additional checks in case of equality */
SCIP_Bool extreduce_bottleneckToSiblingIsDominated(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   extedge,            /**< edge for extension */
   int                   edge2sibling,       /**< edge to sibling of extedge head */
   SCIP_Real             specialDist,        /**< best computed special distance approximation (FARAWAY if unknown) */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Bool hasSpecialDist = LT(specialDist, FARAWAY);

   assert(specialDist >= 0.0);
   assert(extedge >= 0 && edge2sibling >= 0);
   assert(extedge != edge2sibling);
   assert(graph->tail[extedge] == graph->tail[edge2sibling]);

   if( !hasSpecialDist )
   {
      return FALSE;
   }
   else
   {
      const SCIP_Real* const edgecost = graph->cost;

      assert(GE(specialDist, 0.0));

      if( LT(specialDist, edgecost[edge2sibling]) )
         return TRUE;

      if( LT(specialDist, edgecost[extedge]) )
         return TRUE;

      if( LE(specialDist, edgecost[edge2sibling]) )
      {
         const int vertex1 = graph->head[edge2sibling];
         const int vertex2 = graph->head[extedge];

         if( bottleneckIsEqualityDominated(scip, graph, specialDist, edge2sibling,
            vertex1, vertex2, extdata) )
         {
            bottleneckMarkEqualityEdge(scip, graph, edge2sibling, extdata);
            SCIPdebugMessage("...ruled out edge1 with equality! \n");

            return TRUE;
         }
      }

      if( LE(specialDist, edgecost[extedge]) )
      {
         const int vertex1 = graph->head[edge2sibling];
         const int vertex2 = graph->head[extedge];

         if( bottleneckIsEqualityDominated(scip, graph, specialDist, extedge,
            vertex1, vertex2, extdata) )
         {
            bottleneckMarkEqualityEdge(scip, graph, extedge, extdata);
            SCIPdebugMessage("...ruled out edge2 with equality! \n");

            return TRUE;
         }
      }
   }

   return FALSE;
}


/** as above, but for biased special distance */
SCIP_Bool extreduce_bottleneckToSiblingIsDominatedBiased(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   extedge,            /**< edge for extension */
   int                   edge2sibling,       /**< edge to sibling of extedge head */
   SCIP_Real             specialDistBiased,  /**< best computed special distance approximation (FARAWAY if unknown) */
   EXTDATA*              extdata             /**< extension data */
)
{
   const SCIP_Bool hasSpecialDist = LT(specialDistBiased, FARAWAY);

   assert(GE(specialDistBiased, 0.0));
   assert(extedge >= 0 && edge2sibling >= 0);
   assert(extedge != edge2sibling);
   assert(graph->tail[extedge] == graph->tail[edge2sibling]);

   if( !hasSpecialDist )
   {
      return FALSE;
   }
   else
   {
      const SCIP_Real* const edgecost = graph->cost;

      if( LT(specialDistBiased, edgecost[edge2sibling]) )
         return TRUE;

      if( LT(specialDistBiased, edgecost[extedge]) )
         return TRUE;
   }

   return FALSE;
}


/** checks tree bottleneck distances to non-leaves of the tree that were marked before */
void extreduce_bottleneckCheckNonLeaves_pc(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2neighbor,      /**< the edge from the tree to the neighbor */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            ruledOut            /**< could the extension be ruled out */
)
{
   const PCDATA* const pcdata = extdata->pcdata;
   const int* const pcSdCands = pcdata->pcSdCands;
   const int* const tree_deg = extdata->tree_deg;
   const int nPcSdCands = pcdata->nPcSdCands;
   const int neighbor = graph->head[edge2neighbor];
   const int neighbor_base = graph->tail[edge2neighbor];

   assert(pcSdCands);
   assert(ruledOut);
   assert(!(*ruledOut));
   assert(nPcSdCands >= 0);

   /* also check non-leaves */
   for( int c = 0; c < nPcSdCands; c++ )
   {
      SCIP_Real specialDist;
      const int cand = pcSdCands[c];

      assert(cand >= 0 && cand < graph->knots);

      /* leaf, or not contained? */
      if( tree_deg[cand] <= 1 )
         continue;

      specialDist = extreduce_extGetSd(scip, graph, neighbor, cand, extdata);

      if( extreduce_bottleneckWithExtedgeIsDominated(scip, graph, edge2neighbor, neighbor_base, cand, specialDist, extdata) )
      {
         SCIPdebugMessage("---non-leaf bottleneck rule-out---\n");
         *ruledOut = TRUE;

         return;
      }
   }
}


/** checks tree bottleneck distances to non-leaves of the tree */
void extreduce_bottleneckCheckNonLeaves(
   SCIP*                 scip,               /**< SCIP */
   const GRAPH*          graph,              /**< graph data structure */
   int                   edge2neighbor,      /**< the edge from the tree to the neighbor */
   EXTDATA*              extdata,            /**< extension data */
   SCIP_Bool*            ruledOut            /**< could the extension be ruled out */
)
{
   const int* const innerNodes = extdata->tree_innerNodes;
   const int nInnerNodes = extdata->tree_ninnerNodes;
   const int neighbor = graph->head[edge2neighbor];
   const int neighbor_base = graph->tail[edge2neighbor];

   assert(ruledOut);
   assert(!(*ruledOut));

   /* also check non-leaves */
   for( int i = 0; i < nInnerNodes; i++ )
   {
      SCIP_Real specialDist;
      const int node = innerNodes[i];
      assert(graph_knot_isInRange(graph, node));
      assert(extdata->tree_deg[node] > 1);
      assert(node != neighbor_base);

      specialDist = extreduce_extGetSd(scip, graph, neighbor, node, extdata);

      if( extreduce_bottleneckWithExtedgeIsDominated(scip, graph, edge2neighbor, neighbor_base, node, specialDist, extdata) )
      {
         SCIPdebugMessage("---non-leaf bottleneck rule-out---\n");
         *ruledOut = TRUE;
         return;
      }
   }
}

/**@} */
