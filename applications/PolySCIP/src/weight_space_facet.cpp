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

#include "weight_space_facet.h"

#include <iostream>

#include "polyscip.h"

using std::cout;

namespace polyscip {

  using PointType = Polyscip::PointType;
  using ValueType = Polyscip::ValueType;

  WeightSpaceFacet::WeightSpaceFacet(unsigned num_objs,
				     unsigned index) 
    : lhs_(num_objs,0.),
      rhs_(0.) 
  {
    lhs_.at(index) = 1.; 
  }

  WeightSpaceFacet::WeightSpaceFacet(PointType point,
				     ValueType weighted_obj_val) 
    : lhs_(std::move(point)),
      rhs_(weighted_obj_val)
  {}
  
  WeightSpaceFacet::~WeightSpaceFacet() {}

  void WeightSpaceFacet::print() const {
    cout << " Facet: lhs_coeffs = [ ";
    for (const auto& coeff : lhs_)
      cout << coeff << " ";
    cout << "], rhs = " << rhs_ << "\n";
  }
    
}
