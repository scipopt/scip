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
#include <vector>

#include "polyscip.h"

using std::cout;
using std::vector;


namespace polyscip {
  
  using Polyscip::PointType;
  using Polyscip::ValueType;

  // WeightSpaceFacet::WeightSpaceFacet(const CoeffContainer& coeffs, ValueType rhs) 
  //   : coeffs_(coeffs),
  //     rhs_(rhs)
  // {};

  WeightSpaceFacet::WeightSpaceFacet(unsigned num_objs,
				     unsigned index) 
    : lhs_(num_objs,0.),
      rhs_(0) {
    lhs_.at(index) = 1.; 
  }

  WeightSpaceFacet::WeightSpaceFacet(unsigned num_objs,
				     const PointType& point,
				     ValueType weighted_obj_val) 
    : lhs_(point.begin(),point.end()),
      rhs_(weighted_obj_val)
  {};
  
  WeightSpaceFacet::~WeightSpaceFacet() {};

  void WeightSpaceFacet::print() const {
    cout << "Facet: lhs_coeffs = [ ";
    for (const auto& val : coeffs_)
      cout << val << " ";
    cout << "], rhs = " << rhs_ << "\n";
  };      
    
}
