/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*        This file is part of the program PolySCIP                          */
/*                                                                           */
/*    Copyright (C) 2012-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  PolySCIP is distributed under the terms of the ZIB Academic License.     */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with PolySCIP; see the file LICENCE.                               */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * @file weight_space_facet.cpp
 * @brief Implements class representing a facet of the weight space polyhedron
 * @author Sebastian Schenker
 *
 */

#include "weight_space_facet.h"

#include <algorithm>
#include <cstddef> // std::size_t
#include <iterator> // std::prev
#include <numeric>
#include <ostream>

using std::ostream;
using std::size_t;

namespace polyscip {

    /**
    *  Default constructor
    *  @param outcome Values for w_coeffs_
    *  @param wov_coeff Values for wov_coeff
    */
    WeightSpaceFacet::WeightSpaceFacet(OutcomeType outcome,
                                       ValueType wov_coeff)
            : w_coeffs_(std::move(outcome)),
              wov_coeff_(std::move(wov_coeff)) { }


    /**
    * Print function
    * @param os Output stream to print to
    */
    void WeightSpaceFacet::print(ostream &os) const {
        os << " Facet: ";
        size_t counter{0};
        std::for_each(w_coeffs_.cbegin(), std::prev(w_coeffs_.cend()),
                      [&](ValueType val){os << val << " w_" << counter++ << " + ";});
        os << w_coeffs_.back() << " w_" << counter << " >= ";
        if (wov_coeff_ == 0.)
            os << "0\n";
        else
            os << wov_coeff_ << " wov\n";
    }

}
