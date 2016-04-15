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

/**@file   Skeleton.h
 * @brief  Weight space polyhedron
 * @author Sebastian Schenker, Timo Strunk
 *
 * This class represents the lifted weight space polyhedron.  
 * It supplies weights for the solver to test.
 */

#ifndef CLASS_WEIGHTGRAPH
#define CLASS_WEIGHTGRAPH

#undef GCC_VERSION /* lemon/core.h redefines GCC_VERSION additionally to scip/def.h */
#include "lemon/list_graph.h"
#include "scip/scip.h"
#include "scip/def.h"

#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

class WeightSpaceVertex;

using Vec = std::vector<SCIP_Real>;

/** the 1-skeleton of the weight space polyhedron */
class Skeleton
{
 public:
  
  /** constructor */
  Skeleton( unsigned nObjs_);

  /** destructor */
  ~Skeleton();
  
  /** whether there is an untested weight left */
  bool hasNextWeight();

  /** get the next untested weight */
  const Vec * nextWeight();

  /** updates the polyhedron with respect to given point */
  void updatesWeightSpacePolyhedron(
				   const Vec* point  /**< potential new nondominated point */
			 );

   /* /\** like isExtremal but check all vertices for obsolecity (not just the last returned one)*\/ */
   /* bool isExtremalThorough( */
   /*    const std::vector<SCIP_Real>* cost_vector     /\**< potential new nondominated point *\/ */
   /*    ); */

  /** adds a weight space constraint after finding a primal ray with unbounded weighted objective */
  void addPrimalRay(
    const Vec* cost_ray       /**< cost vector of the unbounded primal ray */
    );

  /* /\** like addPrimalRay but check all vertices for obsolecity (not just the last returned one)*\/ */
  /* void addPrimalRayThorough( */
  /*   const std::vector<SCIP_Real>* cost_ray      /\**< cost vector of the unbounded primal ray *\/ */
  /*   ); */

  /** add multiple rays with unbounded weighted objective */
  void addPrimalRays(
    const std::vector<const std::vector<SCIP_Real>* >& cost_rays    /**< rays to add */
    );

 private:
  
  using Graph = lemon::ListGraph;
  using Node = Graph::Node;
  using vert_node_map = std::unordered_map<WeightSpaceVertex*, Node>;

  /** initialize the polyhedron with the first solution */
  void init(
	    const std::vector<SCIP_Real>* cost_vector, /**< cost vector of first solution */
	    const std::vector< const std::vector<SCIP_Real>* >& cost_rays,  /**< list of known unbounded cost rays */
	    bool nondom_point_from_unit_weight,
	    unsigned unit_weight_index
     );

  /** create all facets defining the inital weight space polyhedron */
  WeightSpaceVertex const * const  createInitialFacetsAndWeightSpaceVerts(
     const std::vector<SCIP_Real>* nondom_point,      /**< first nondom_point that was found */
     bool nondom_point_from_unit_weight,
     unsigned unit_weight_index
									  );

  void createInitialWeightSpacePolyhedron(const std::vector<SCIP_Real>* nondom_point,
					  const std::vector<const std::vector<SCIP_Real>* >& cost_rays,
					  WeightSpaceVertex const * const unit_weight_index);

  struct Facet {  //form of facet: w^T w_coeffs - a_coeff*a >= 0
    Facet(std::vector<SCIP_Real>* w_coeffs, SCIP_Real a_coeff) {
      w_coeffs_ = w_coeffs;
      a_coeff_ = a_coeff;
    };
    ~Facet() {
      delete w_coeffs_;
    };
    
    Vec* w_coeffs_;
    SCIP_Real a_coeff_; 
  };

  std::set<Node> untested_nodes_;     /**< nodes not passed to solver yet */

  Graph graph_;                             /**< the graph structure */

  unsigned nObjs_;    /**< number of objectives of underlying problem */
  
  Graph::NodeMap<WeightSpaceVertex*> nodes_to_verts_;  /**< map from nodes to vertices */

  Node last_considered_node_;     /**< last considered node */

  vert_node_map verts_to_nodes_;                  /**< map from vertices to graph nodes */

  std::vector<const WeightSpaceVertex*> vertices_;      /**< list of all generated vertices */

  std::vector< const Skeleton::Facet* > facets_; /**< all facets of the polyhedron */

  


  void updateWeightSpace(const std::vector<const std::vector<SCIP_Real>* >& cost_rays);

  /* data structures for temporary use in update step*/
  //  const std::vector<SCIP_Real>* new_facet_;      /**< new facet of weight space polyhedron */
  
  //  std::vector<WeightSpaceVertex*>* new_vertices_; /**< new generated vertices */

  //  std::queue<lemon::ListGraph::Node>* unscanned_nodes_; /**< nodes left to scan for obsolecity */
  
  //  std::set<lemon::ListGraph::Node>* obsolete_nodes_;    /**< nodes identified as obsolete */
  
  //  std::vector<lemon::ListGraph::Edge>* cut_edges_;    /**< edges from obsolete to nonobsolete nodes */
  /* end of data structures for temporary use in update step*/

  /** returns true if point is a new non-dominated point; otherwise false */
  bool isNewNondomPoint(
			const std::vector<SCIP_Real>* point, /**< potential new nondominated point */
			double comp_val = 0.0              /**< value used for checking inequality */ 
			) const;


  
  /* /\** create corner vertex of initial weight space polyhedron *\/ */
  /* WeightSpaceVertex* createCorner( */
  /*   int index                       /\**< index where weight is 1 *\/ */
  /*   ); */

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

#endif
