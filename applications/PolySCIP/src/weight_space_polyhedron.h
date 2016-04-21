/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** @file   weight_space_polyhedron.h
 *  @brief  The (partial) weight space polyhedron
 *  @author Sebastian Schenker
 *  @author Timo Strunk
 *
 * This class represents the (partial) weight space polyhedron P =
 * {(w,a) \in \Lambda \times R : w \cdot y >= a \forall y \in Y'}
 * where Y' is the set of non-dominated points computed so far and
 * \Lambda is the set of normalized weights
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED 
#define POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED 

#include <memory> // std::shared_ptr
#include <queue>
#include <tuple>
#include <unordered_map>
#include <utility> // std::pair
#include <vector>

#include "lemon/list_graph.h"

#include "polyscip.h"
#include "weight_space_facet.h"
#include "weight_space_vertex.h"

namespace polyscip {

  /** @brief 1-skeleton of the (partial) weight space polyhedron. */
  class WeightSpacePolyhedron {
  public:
    /** Container used to store the weight space vertices; needs to support: push_back */
    using VertexContainer = std::vector< std::shared_ptr<WeightSpaceVertex> >;

    /** Container used to store the weight space facets; needs to support; push_back */
    using FacetContainer = std::vector< std::shared_ptr<WeightSpaceFacet> >;

    /** Container used to store the untested weight space vertices;
	Needs to support: empty(),
    */
    using UntestedVertexContainer = std::queue< std::shared_ptr<WeightSpaceVertex> >; 

    using RayInfoType = std::tuple< std::shared_ptr<RayType>, std::shared_ptr<WeightType> >;
    
    /** @brief Creates the skeleton of the initial (partial) weight
	space polyhedron P = {(w,a) \in \Lambda \times R : w \cdot y^1
	>= a}
	@param num_objs number of objectives of given problem
	@param point first computed (weakly non-dominated) point
	@param point_weighted_obj_val weighted objective value of first point
	@param unit_weight_used indicates whether first point was computed by using 
	a unit weight; if unit_weight_used.first is true, then unit_weight_used.second containts 
	index of 1 in used unit weight
    */
    WeightSpacePolyhedron(unsigned num_objs,
			  const Polyscip::PointType& point, 
			  PolySCIP::ValueType point_weighted_obj_val,
			  std::pair<bool,unsigned> unit_weight_used);

    /** @brief Destructor */
    ~WeightSpacePolyhedron();
  
    /** @brief Checks whether there is an untested weight space vertex
     *  @return true if there is an untested weight space vertex; false otherwise
     */
    bool hasUntestedVertex() const;

    /** @brief Returns an untested weight space vertex
     *  @return an untested weight space vertex
     */
    std::shared_ptr<WeightSpaceVertex> getUntestedVertex();

    /** @brief Incorporates a newly found non-dominated point into the
     * (partial) weight space polyhedron 
     *  @param old_vertex the vertex (yielding weight and weight
     * objective value) that was considered in last computation 
     *  @param new_point the newly found non-dominated point that was 
     * computed by considering the weight and weighted objective 
     * value given by old_vertex
     */
    void addNondomPoint(std::shared_ptr<WeightSpaceVertex> old_vertex, 
			std::shared_ptr<PointType> new_point);

    /** @brief Incorporates an newly found unbounded non-dominated ray
     * into the (partial) weight space polyhedron
     *  @param old_vertex the vertex (yielding weight and weight
     * objective value) that was considered in last computation 
     *  @param new_ray the newly found non-dominated ray that was 
     * computed by considering the weight and weighted objective 
     * value given by old_vertex
     */
    void addNondomRay(std::shared_ptr<WeightSpaceVertex old_vertex,
		      std::shared_ptr<RayType> new_ray);

  private:
    using Graph = lemon::ListGraph;
    using Node = Graph::Node;
    using NodeMap = Graph::NodeMap< std::shared_ptr<WeightSpaceVertex> >;
    using VertexMap = std::unordered_map<std::shared_ptr<WeightSpaceVertex>, Node>;

    /** Creates initial weight space boundary facets w_i >= 0 for i \in [num_objs]
     *  @param num_objs number of objectives of given problem
     */
    void createInitialBoundaryFacets(unsigned num_objs);

    /** Creates initial weight space vertices
     *  @param num_objs number of objectives of given problem
     *  @param point first computed (weakly non-dominated) point
     *	@param weighted_obj_val weighted objective value of given point
    */
    void createInitialWeightSpaceVertices(unsigned num_objs,
					  const Polyscip::PointType& point,
					  Polyscip::ValueType weighted_obj_val);

    /** Creates initial 1-skeleton of complete graph with number of
	objectives many vertices
     */
    void createInitialSkeleton() {};

    /** Establishes the vertices that need to be tested
     *  @param unit_weight_used true if first (weakly non-dominated) point was computed by 
     *  using a unit weight, otherwise false
     *  @param unit_weight_index index of 1 in unit weight if unit weight was used, otherwise 
     *  not used
     */
    void establishUntestedVertices(bool unit_weight_used,
				   unsigned unit_weight_index);


    VertexContainer vertices_;                   /**< all computed weight space vertices */
    FacetContainer facets_;                      /**< all computed weight space facets */
    UntestedVertexContainer untested_vertices_;  /**< nodes not passed to solver yet */
    Graph skeleton_;                             /**< 1-skeleton of the weight space polyhedron */
    NodeMap nodes_to_vertices_;                  /**< maps nodes to vertices */
    VertexMap vertices_to_nodes_;                /**< maps vertices to nodes */
    

  
    void updateWeightSpace(const std::vector<const std::vector<SCIP_Real>* >& cost_rays);

    /* data structures for temporary use in update step*/
    //  const std::vector<SCIP_Real>* new_facet_;      /**< new facet of weight space polyhedron */
  
    //  std::vector<WeightSpaceVertex*>* new_vertices_; /**< new generated vertices */

    //  std::queue<lemon::ListGraph::Node>* unscanned_nodes_; /**< nodes left to scan for obsolecity */
  
    //  std::set<lemon::ListGraph::Node>* obsolete_nodes_;    /**< nodes identified as obsolete */
  
    //  std::vector<lemon::ListGraph::Edge>* cut_edges_;    /**< edges from obsolete to nonobsolete nodes */
    /* end of data structures for temporary use in update step*/

    /** wether the new solution makes a given weight space vertex obsolete */
    bool isMakingObsolete(
			  const std::vector<SCIP_Real>* cost_vector,     /**< cost vector of a solution */
			  const WeightSpaceVertex* vertex,               /**< vertex that might be obsolete */
			  bool strict=false                              /**< no tolerance for slight obsolecity */
			  );

    /** returns node made obsolete by facet or INVALID */
    lemon::ListGraph::Node findObsoleteNode(const std::vector<SCIP_Real>* facet);

    /** updates the polyhedron with the new facet */
    void addFacet(
		  const std::vector<SCIP_Real>* facet            /**< new facet */
		  );

    /** tests all the neighbours of an obsolete node for obsolecity and then removes the node */
    void scanNode(
		  lemon::ListGraph::Node obs_node            /**< a skeleton node marked as obsolete */
		  );

    /** apply changes calculated by add facet */
    void updateGraph();

    /** calculate new vertices from obsolete vertices and add them to the graph */
    void createNewVertices();

    /** calculate and add edges between all pairs of combinatorially adjacent new vertices */
    void createNewEdges();

    /** creates a new node between an obsolete and a non obsolete node */
    void makeIntermediateVertex(
				lemon::ListGraph::Edge cut_edge  /**< an edge between an obsolete and a non-obsolete node */
				);

    /** special method dealing with obsolete nodes that are also corners of the weight space */
    void updateCorner(
		      WeightSpaceVertex* obsolete_vertex     /**< an obsolete vertex that is a corner */
		      );

    /** adds a graph node corresponding to new_vertex */
    lemon::ListGraph::Node addNode(
				   WeightSpaceVertex* new_vertex,  /**< new vertex not yet represented in the graph */
				   bool mark_untested=true         /**< set true if weight should be tested at some point */
				   );

    /** returns true if data is combinatorically consistent */
    bool graphIsValid() const;

    /** returns facet vector corresponding to point */
    const std::vector<SCIP_Real>* createFacetFromPoint(
						       const std::vector<SCIP_Real>* point  /**< cost vector of a solution */
						       ) const;

    /** returns facet vector corresponding to point */
    const std::vector<SCIP_Real>* createFacetFromCost(
						      const std::vector<SCIP_Real>* cost_vector      /**< cost vector of a solution */
						      ) const;

    /** returns facet vector corresponding to an unbounded cost ray */
    const std::vector<SCIP_Real>* createFacetFromRay(
						     const std::vector<SCIP_Real>* cost_ray        /**< cost vector of a primal ray */
						     ) const;

  };

}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_POLYHEDRON_H_INCLUDED 
