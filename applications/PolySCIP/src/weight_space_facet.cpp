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
#include <ostream>
#include <tuple> // std::tie

using std::ostream;

namespace polyscip {

    bool operator<(const WeightSpaceFacet &facet1, const WeightSpaceFacet &facet2) {
        return std::tie(facet1.wov_coeff_, facet1.w_coeffs_) <
               std::tie(facet2.wov_coeff_, facet2.w_coeffs_);
    }

    WeightSpaceFacet::WeightSpaceFacet(unsigned num_objs,
                                       unsigned index)
            : w_coeffs_(num_objs, 0.),
              wov_coeff_(0.) {
        w_coeffs_.at(index) = 1.;
    }

    WeightSpaceFacet::WeightSpaceFacet(const OutcomeType& outcome,
                                       ValueType wov_coeff)
            : w_coeffs_(begin(outcome), end(outcome)),
              wov_coeff_{wov_coeff} { }

    void WeightSpaceFacet::print(ostream &os) const {
        os << " Facet: ";
        unsigned counter{0};
        std::for_each(w_coeffs_.cbegin(), --w_coeffs_.cend(),
                      [&](ValueType val){os << val << " w_" << counter++ << " + ";});
        os << w_coeffs_.back() << " w_" << counter << " >= ";
        if (wov_coeff_ != 0.0)
            os << "wov\n";
        else
            os << "0\n";
    }

}
