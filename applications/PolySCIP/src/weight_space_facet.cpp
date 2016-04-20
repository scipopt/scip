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

/**@file   weight_space_facet.cpp
 * @brief  Weight space facet class definitions
 * @author Sebastian Schenker
 */

#include "weight_space_facet.h"

#include <iostream>

using std::cout;

namespace polyscip {

  WeightSpaceFacet::WeightSpaceFacet(const CoeffContainer& coeffs, ValueType rhs) 
    : coeffs_(coeffs),
      rhs_(rhs)
  {};

  WeightSpaceFacet::~WeightSpaceFacet() {};

  void WeightSpaceFacet::print() const {
    cout << "Facet: lhs_coeffs = [ ";
    for (const auto& val : coeffs_)
      cout << val << " ";
    cout << "], rhs = " << rhs_ << "\n";
  };      
    
}
