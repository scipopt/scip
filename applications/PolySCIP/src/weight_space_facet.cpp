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

#include <ostream>
#include <tuple> // std::tie
#include <vector>

#include "polyscip.h"

using std::ostream;
using std::vector;

namespace polyscip {

    using OutcomeType = Polyscip::OutcomeType;
    using ValueType = Polyscip::ValueType;

    bool operator<(const WeightSpaceFacet& facet1, const WeightSpaceFacet& facet2) {
        return std::tie(facet1.rhs_,facet1.lhs_) < std::tie(facet2.rhs_,facet2.lhs_);
    }

    WeightSpaceFacet::WeightSpaceFacet(unsigned num_objs,
                                       unsigned index)
            : lhs_(num_objs, 0.),
              rhs_{0.} {
      lhs_.at(index) = 1.;
    }

    WeightSpaceFacet::WeightSpaceFacet(const OutcomeType& point,
                                       ValueType weighted_obj_val)
            : lhs_(point.begin(), point.end()),
              rhs_{weighted_obj_val} {}

    void WeightSpaceFacet::print(ostream& os) const {
      os << " Facet: lhs_coeffs = [ ";
      for (const auto &coeff : lhs_)
        os << coeff << " ";
      os << "], rhs = " << rhs_ << "\n";
    }

}
