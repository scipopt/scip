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

/**
 * @brief Class representing a facet of the weight space polyhedron
 * @author Sebastian Schenker
 *
 * Data structure representing a facet of the (partial) weight space
 * polyhedron P={(w,a) : w \cdot y >= a \forall y \in Y_N} where Y_N is
 * the (current) set of non-dominated points. A facet (w_coeffs_, wov_coeff_) is
 * represented by coefficients 'w_coeffs_' and a right hand side 'wov_coeff_'
 * yielding an inequality of the form w_coeffs_ \cdot w >= wov_coeff_ * wov
 * where wov stands for weighted objective value
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


    WeightSpaceFacet::WeightSpaceFacet(OutcomeType outcome,
                                       ValueType wov_coeff)
            : w_coeffs_(std::move(outcome)),
              wov_coeff_(std::move(wov_coeff)) { }


    ValueType WeightSpaceFacet::getWeightedWeight(const WeightType& weight) const {
        assert (weight.size() == w_coeffs_.size());
        return std::inner_product(begin(w_coeffs_), end(w_coeffs_),
                                  begin(weight), 0.);
    }



    void WeightSpaceFacet::print(ostream &os) const {
        os << " Facet: ";
        size_t counter{0};
        std::for_each(w_coeffs_.cbegin(), std::prev(w_coeffs_.cend()),
                      [&](ValueType val){os << val << " w_" << counter++ << " + ";});
        os << w_coeffs_.back() << " w_" << counter << " >= ";
        if (wov_coeff_ == 0.0)
            os << "0\n";
        else
            os << wov_coeff_ << " wov\n";
    }

}
