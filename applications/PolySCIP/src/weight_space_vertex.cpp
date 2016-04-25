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

#include <algorithm> // std::copy, std::set_intersection
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator, std::inserter
#include <memory>    // std::shared_ptr;
#include <numeric>   // std::inner_product;

#include "polyscip.h"

using std::cout;
using std::shared_ptr;
using std::vector;

namespace polyscip {

  using WeightType = Polyscip::WeightType;
  using ValueType = Polyscip::ValueType;

  WeightSpaceVertex::WeightSpaceVertex(FacetContainer incident_facets,
				       shared_ptr<const WeightType> weight,
				       ValueType weighted_obj_val)
    : incident_facets_(std::move(incident_facets)),
      weight_(std::move(weight)),
      weighted_obj_val_(weighted_obj_val)
  {}

  WeightSpaceVertex::~WeightSpaceVertex() {}

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

  bool WeightSpaceVertex::hasUnitWeight(Polyscip::WeightType::size_type index) const {
    assert (index < weight_->size());
    for (auto i=0; i!=weight_->size(); ++i) {
      if (i == index) { // check whether i-th element is 1.
	if ((*weight_)[i] != 1.)
	  return false;
      }
      else { // check whether i-th element is 0.
	if ((*weight_)[i] != 0.)
	  return false;
      }
    }
    return true;
  }

  void WeightSpaceVertex::print(bool printFacets) const {
    cout << "WeightSpaceVertex:\n weight = [ ";
    std::ostream_iterator<ValueType> out_it(cout, " ");
    std::copy(weight_->cbegin(), weight_->cend(), out_it);
    cout << "]\n weighted objective value = " << weighted_obj_val_ << "\n";
    if (printFacets) {
      cout << " defining facets: \n";
      for (const auto& facet : incident_facets_)
	facet->print();
    }
  }

}
