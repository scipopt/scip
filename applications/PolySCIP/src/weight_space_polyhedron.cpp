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

#include "weight_space_polyhedron.h"

#include <iostream> 
#include <memory> // std::shared_ptr
#include <numeric> // std::inner_product
#include <utility> // std::move, std::pair
#include <vector>

#include "polyscip.h"
#include "weight_space_facet.h"
#include "weight_space_vertex.h"

using std::inner_product;
using std::make_shared;
using std::shared_ptr;
using std::vector;

namespace polyscip {

  using PointType = Polyscip::PointType;
  using ValueType = Polyscip::ValueType;
  using WeightType = Polyscip::WeightType;
  using RayContainer = Polyscip::RayContainer;
  using FacetContainer = WeightSpaceVertex::FacetContainer;

  WeightSpacePolyhedron::WeightSpacePolyhedron(unsigned num_objs,
					       PointType point,
					       ValueType weighted_obj_val,
					       std::pair<bool,WeightType::size_type> unit_weight_info) 
    : skeleton_(),
      nodes_to_vertices_(skeleton_)
  {
    /* create initial boundary facets w_i > 0 for i \in [num_objs] */
    auto facets = FacetContainer{};
    for (auto i=0; i<num_objs; ++i)
      facets.emplace_back(make_shared<const WeightSpaceFacet>(num_objs, i));
    createInitialVertices(num_objs, std::move(point), weighted_obj_val, std::move(facets));
    createInitialSkeleton();
    if (unit_weight_info.first)
      setMarkedVertex(unit_weight_info.second);
    //    if (!computed_rays.empty()) 
    //updateVertexWeightsWithRayInfo(ray_info);
  }

  WeightSpacePolyhedron::~WeightSpacePolyhedron() {}

  bool WeightSpacePolyhedron::hasUnmarkedVertex() const {
    return !unmarked_vertices_.empty();
  }
  
  void WeightSpacePolyhedron::createInitialVertices(unsigned num_objs,
						    PointType point,
						    ValueType weighted_obj_val,
						    FacetContainer boundary_facets) {
    /* make facet: point \cdot w >= weighted_obj_val to facets_ */
    auto point_facet = make_shared<const WeightSpaceFacet>(point, weighted_obj_val);
    /* create initial weight space vertices */
    for (auto i=0; i<num_objs; ++i) {
      auto facets = FacetContainer(begin(boundary_facets), end(boundary_facets));
      facets.at(i) = point_facet; // replace i-th boundary facet 
      auto weight = make_shared<WeightType>(num_objs, 0.);
      (*weight)[i] = 1.; // set unit weight
      unmarked_vertices_.emplace_back(make_shared<WeightSpaceVertex>(std::move(facets), 
								     std::move(weight), 
								     point.at(i)));
    }
  }

  void WeightSpacePolyhedron::createInitialSkeleton() {
    for (const auto vertex : unmarked_vertices_) {
      Node new_node = skeleton_.addNode(); 
      for (const auto& pair : vertices_to_nodes_) // add edge between all previous nodes and new_node
	skeleton_.addEdge(pair.second, new_node); 
      nodes_to_vertices_[new_node] = vertex;
      vertices_to_nodes_.insert({vertex, new_node});
    }
  }

  void WeightSpacePolyhedron::setMarkedVertex(WeightType::size_type unit_weight_index) {
    for (auto it=begin(unmarked_vertices_); it!=end(unmarked_vertices_); ++it) {
      if ((*it)->hasUnitWeight(unit_weight_index)) {
	marked_vertices_.emplace_back(std::move(*it));
	unmarked_vertices_.erase(it);
	return;
      }
    }

  }
  
  // void WeightSpacePolyhedron::updateVertexWeightsWithRayInfo(const RayContainer& ray_info, 
  // 							     double zero) {
  //   for (const auto& rayPair : ray_info) {
  //     for (shared_ptr<const WeightSpaceVertex>& v : untested_vertices_) {
  // 	/* check whether vertex weight \cdot ray < zero, i.e., whether the ray 
  // 	   yields an unbounded objective value w.r.t the vertex weight */
  // 	assert (v->getWeight->size() == rayPair.first.size());
  // 	ValueType res = inner_product(begin(*(v->getWeight())),
  // 				      end(*(v->getWeight())),
  // 				      begin(rayPair.first),
  // 				      0.);
  // 	if (res < zero) { // ray yields unbounded value w.r.t vertex weight
  // 	  ;
  // 	}

  //     }
  //   }

  // }

  void WeightSpacePolyhedron::printUnmarkedVertices(bool detailed) const {
    std::cout << "UNMARKED VERTICES: \n";
    for (const auto& vertex : unmarked_vertices_)
      vertex->print(detailed);
  }

  void WeightSpacePolyhedron::printMarkedVertices(bool detailed) const {
    std::cout << "MARKED VERTICES: \n";
    for (const auto& vertex : marked_vertices_)
      vertex->print(detailed);
  }

}

