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

/**@file   weight_space_vertex.cpp
 * @brief  Weight space vertex class definitions
 * @author Sebastian Schenker
 * @author Timo Strunk
 */

#include "weight_space_vertex.h"
#include "scip/def.h"

#include <algorithm> // std::copy, std::set_intersection
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator, std::inserter
#include <memory>
#include <numeric>   // std::inner_product;

using std::cout;
using std::shared_ptr;
using std::vector;

namespace polyscip {

  WeightSpaceVertex::WeightSpaceVertex(const vector< shared_ptr<const WeightSpaceFacet> >& incident_facets,
				       shared_ptr<const WeightType> weight,
				       ValueType weighted_obj_val)
    : incident_facets_(incident_facets.begin(),incident_facets_.end()),
      weight_(weight),
      weighted_obj_val_(weighted_obj_val)
  {};

  WeightSpaceVertex::~WeightSpaceVertex() {};

  // bool WeightSpaceVertex::isAdjacent(shared_ptr<const WeightSpaceVertex> other_vertex) const {
    
  //   set<unsigned> common_facets;
  //   set<unsigned>* otherFacets = point.getFacets();

  //   std::set_intersection(facet_indices_.begin(),
  // 			  facet_indices_.end(),
  // 			  otherFacets->begin(),
  // 			  otherFacets->end(),
  // 			  std::inserter( common_facets, common_facets.begin() ) );

  //   return common_facets.size() == nobjs_-1;
  // }

  ValueType WeightSpaceVertex::getWeightedObjVal() const {
    return weighted_obj_val_;
  }

  shared_ptr<const WeightType> WeightSpaceVertex::getWeight() const {
    return weight_;
  }

  shared_ptr<const FacetContainer> WeightSpaceVertex::getFacets() const {
    return facets_;
  }

  void WeightSpaceVertex::print(bool printFacets) const {
    cout << "Weight space vertex: weighted objective value = ";
    cout << weighted_obj_val_ << " weight = [ ";
    std::ostream_iterator<ValueType> out_it(cout, " ");
    std::copy(weight_->begin(), weight_->end(), out_it);
    cout << "]";
    if (printFacets) {
      cout << "defining facets: \n";
      for (const auto& facet : facets_) 
	facet->print();
    }
  }

}
