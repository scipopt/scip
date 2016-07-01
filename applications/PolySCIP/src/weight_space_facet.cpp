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

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <tuple> // std::tie

using std::ostream;

namespace polyscip {


    WeightSpaceFacet::WeightSpaceFacet(OutcomeType outcome,
                                       ValueType wov_coeff)
            : w_coeffs_(outcome),
              wov_coeff_{wov_coeff} { }

    ValueType WeightSpaceFacet::getWeightedWeight(const WeightType& weight) const {
        assert (weight.size() == w_coeffs_.size());
        return std::inner_product(begin(w_coeffs_), end(w_coeffs_),
                                  begin(weight), 0.);
    }



    void WeightSpaceFacet::print(ostream &os) const {
        os << " Facet: ";
        std::size_t counter{0};
        std::for_each(w_coeffs_.cbegin(), --w_coeffs_.cend(),
                      [&](ValueType val){os << val << " w_" << counter++ << " + ";});
        os << w_coeffs_.back() << " w_" << counter << " >= ";
        if (wov_coeff_ == 0.0)
            os << "0\n";
        else
            os << wov_coeff_ << " wov\n";
    }

}
