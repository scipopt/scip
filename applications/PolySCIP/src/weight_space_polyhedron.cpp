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

#include <memory>
#include <vector>

#include "polyscip.h"

using std::make_shared;
using std::vector;

namespace polyscip {

  using Polyscip::PointType;
  using Polyscip::ValueType;
  using Polyscip::WeightType;

  WeightSpacePolyhedron::WeightSpacePolyhedron(unsigned num_objs,
					       const PointType& point,
					       ValueType weighted_obj_val,
					       pair<bool,unsigned> unit_weight_used) {
    createInitialBoundaryFacets(num_objs); 
    createInitialWeightSpaceVertices(num_objs, point, weighted_obj_val);
    createInitialSkeleton();
    addVerticesToUntestedVertices
  };

  WeightSpacePolyhedron::~WeightSpacePolyhedron() {};

  bool WeightSpacePolyhedron::hasUntestedVertex() const {
    return !untested_vertices_.empty();
  };
  
  void WeightSpacePolyhedron::createInitialBoundaryFacets(unsigned num_objs) {
    for (unsigned i=0; i<num_objs; ++i) 
      facets_.push_back(make_shared<WeightSpaceFacet>(num_objs, i));
  };

  void WeightSpacePolyhedron::createInitialWeightSpaceVertices(unsigned num_objs,
							       const PointType& point,
							       ValueType weighted_obj_val) {
    /* make facet: point \cdot w >= weighted_obj_val to facets_ */
    auto new_facet = make_shared<WeightSpaceFacet>(point, weighted_obj_val);
    /* create initial weight space vertices */
    for (unsigned i=0; i<num_objs; ++i) {
      vector< shared_ptr<WeightSpaceFacet> > defining_facets(facets_.begin(), facets_.end());
      defining_facets[i] = new_facet;
      auto weight = make_shared<WeightType>(num_objs, 0.);
      (*weight)[i] = 1.; // set unit weight
      vertices_.push_back(make_shared<WeightSpaceVertex>(defining_facets, std::move(weight), 
							 point[i]));
    }
  };

  void WeightSpacePolyhedron::createInitialSkeleton() {
    for (shared_ptr<WeightSpaceVertex> v : vertices_) {
      Node new_node = skeleton_.addNode(); 
      for (const auto& pair : verts_to_nodes_)    // add edge between all previous nodes and new_node
	skeleton_.addEdge(pair.second, new_node); 
      nodes_to_vertices_[new_node] = v;
      verts_to_nodes.insert({v, new_node});
    }
  };

  void WeightSpacePolyhedron::establishUntestedVertices(bool unit_weight_used, 
							unsigned unit_weight_index) {
    for ();
  };

  };
    
    

}

