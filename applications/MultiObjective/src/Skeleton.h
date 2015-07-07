/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Skeleton.h
 * @brief  Weight space polyhedron
 * @author Timo Strunk
 *
 * This class represents the lifted weight space polyhedron.  It supplies weights for the solver to test.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef CLASS_WEIGHTGRAPH
#define CLASS_WEIGHTGRAPH

#include <map>
#include <set>
#include <vector>
#include <queue>

#undef GCC_VERSION
#include "lemon/list_graph.h"

#include "scip/scip.h"
#include "scip/def.h"

class WeightSpaceVertex;

/** the 1-skeleton of the weight space polyhedron */
class Skeleton
{
 public:
   /** constructor */
   Skeleton(
      SCIP*              scip                /**< SCIP solver */
      );

   /** destructor */
   ~Skeleton();

   /** initialize the polyhedron with the first solution */
   void init(
      const std::vector<SCIP_Real>*     cost_vector,                       /**< cost vector of first solution */
      std::vector< const std::vector<SCIP_Real>* >*   cost_rays = NULL     /**< list of known unbounded cost rays */
      );

   /** whether there is an untested weight left */
   bool hasNextWeight();

   /** get the next untested weight */
   const std::vector<SCIP_Real> * nextWeight();

   /** returns true and updates the polyhedron if cost vector is a new nondominated point */
   bool isExtremal(
      const std::vector<SCIP_Real>*     cost_vector         /**< potential new nondominated point */
      );

   /** like isExtremal but check all vertices for obsolecity (not just the last returned one)*/
   bool isExtremalThorough(
      const std::vector<SCIP_Real>*     cost_vector         /**< potential new nondominated point */
      );

   /** adds a weight space constraint after finding a primal ray with unbounded weighted objective */
   void addPrimalRay(
      const std::vector<SCIP_Real>*     cost_ray            /**< cost vector of the unbounded primal ray */
      );

   /** like addPrimalRay but check all vertices for obsolecity (not just the last returned one)*/
   void addPrimalRayThorough(
      const std::vector<SCIP_Real>*     cost_ray            /**< cost vector of the unbounded primal ray */
      );

   /** add multiple rays with unbounded weighted objective */
   void addPrimalRays(
      const std::vector< const std::vector<SCIP_Real>* >*   cost_rays /**< rays to add */
      );

   /** get number of vertices added in last isExtremal call*/
   int getNNewVertices() const;

   /** get number of vertices processed in last isExtremal call*/
   int getNProcessedVertices() const;

 private:
   SCIP*                                          scip_;                   /**< SCIP solver */
   std::set<lemon::ListGraph::Node>               untested_nodes_;         /**< nodes not passed to solver yet */
   lemon::ListGraph                               graph_;                  /**< the graph structure */
   lemon::ListGraph::NodeMap<WeightSpaceVertex*>  vertex_map_;             /**< map from nodes to vertices */
   lemon::ListGraph::Node                         last_returned_node_;     /**< last tested node */
   std::vector<WeightSpaceVertex*>                vertices_;               /**< list of all generated vertices */
   int                                            n_new_nodes_;            /**< number of vertices added in last call*/
   int                                            n_proc_nodes_;           /**< number of vertices processed in
									    * last call to isExtremal() */
   std::vector< const std::vector<SCIP_Real>* >   facets_;                 /**< all facets of the polyhedron */

   /* data structures for temporary use in update step*/
   const std::vector<SCIP_Real>*                  new_facet_;         /**< new facet of weight space polyhedron */
   std::vector<WeightSpaceVertex*>*               new_vertices_;      /**< new generated vertices */
   std::queue<lemon::ListGraph::Node>*            unscanned_nodes_;   /**< nodes left to scan for obsolecity */
   std::set<lemon::ListGraph::Node>*              obsolete_nodes_;    /**< nodes identified as obsolete */
   std::vector<lemon::ListGraph::Edge>*           cut_edges_;         /**< edges from obsolete to nonobsolete nodes */

   /** create all facets defining the inital weight space polyhedron */
   void createInitialFacets(
      const std::vector<SCIP_Real>*     cost_vector         /**< cost vector of first solution */
      );

   /** create corner vertex of initial weight space polyhedron */
   WeightSpaceVertex* createCorner(
      int                               index               /**< index where weight is 1 */
      );

   /** wether the new solution makes a given weight space vertex obsolete */
   bool isMakingObsolete(
      const std::vector<SCIP_Real>*     cost_vector,        /**< cost vector of a solution */
      const WeightSpaceVertex*          vertex,             /**< vertex that might be obsolete */
      bool                              strict=false        /**< no tolerance for slight obsolecity */
      );

   /** returns node made obsolete by facet or INVALID */
   lemon::ListGraph::Node findObsoleteNode(
      const std::vector<SCIP_Real>* facet
      );

   /** updates the polyhedron with the new facet */
   void addFacet(
      const std::vector<SCIP_Real>*     facet               /**< new facet */
      );

   /** tests all the neighbours of an obsolete node for obsolecity and then removes the node */
   void scanNode(
      lemon::ListGraph::Node            obs_node            /**< a skeleton node marked as obsolete */
      );

   /** apply changes calculated by add facet */
   void updateGraph();

   /** calculate new vertices from obsolete vertices and add them to the graph */
   void createNewVertices();

   /** calculate and add edges between all pairs of combinatorially adjacent new vertices */
   void createNewEdges();

   /** creates a new node between an obsolete and a non obsolete node */
   void makeIntermediateVertex(
      lemon::ListGraph::Edge            cut_edge            /**< an edge between an obsolete and a non-obsolete node */
      );

   /** special method dealing with obsolete nodes that are also corners of the weight space */
   void updateCorner(
      WeightSpaceVertex*    obsolete_vertex     /**< an obsolete vertex that is a corner */
      );

   /** adds a graph node corresponding to new_vertex */
   lemon::ListGraph::Node addNode(
      WeightSpaceVertex*    new_vertex,         /**< new vertex not yet represented in the graph */
      bool                  mark_untested=true  /**< set true if weight should be tested at some point */
      );

   /** returns true if data is combinatorically consistent */
   bool graphIsValid() const;

   /** returns facet vector corresponding to point */
   const std::vector<SCIP_Real>* createFacetFromCost(
       const std::vector<SCIP_Real>*    cost_vector    /**< cost vector of a solution */
       ) const;

   /** returns facet vector corresponding to an unbounded cost ray */
   const std::vector<SCIP_Real>* createFacetFromRay(
      const std::vector<SCIP_Real>*        cost_ray      /**< cost vector of a primal ray */
      ) const;


};

#endif
