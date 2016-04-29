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

#include "weight_space_vertex.h"

#include <algorithm> // std::copy
#include <iterator>  // std::ostream_iterator
#include <numeric>   // std::inner_product;
#include <ostream>

#include "polyscip.h"

using std::ostream;
using std::vector;

namespace polyscip {

    using OutcomeType = Polyscip::OutcomeType;
    using ValueType = Polyscip::ValueType;
    using WeightType = Polyscip::WeightType;

    WeightSpaceVertex::WeightSpaceVertex(FacetContainer incident_facets,
                                         WeightType weight,
                                         ValueType weighted_obj_val)
            : incident_facets_(std::move(incident_facets)),
              weight_(std::move(weight)),
              weighted_obj_val_{weighted_obj_val} {}


    WeightType WeightSpaceVertex::getWeight() const {
        return weight_;
    }

    bool WeightSpaceVertex::isMadeObsolete(const OutcomeType& outcome, ValueType rhs) const {
        assert (outcome.size() == weight_.size());
        ValueType res = std::inner_product(begin(outcome),
                                           end(outcome),
                                           begin(weight_),
                                           -rhs);
        return res < 0.;
    }


    bool WeightSpaceVertex::isMadeObsolete(const OutcomeType& outcome) const {
        return isMadeObsolete(outcome, weighted_obj_val_);
    }

    bool WeightSpaceVertex::hasSameWeight(const Polyscip::WeightType& weight) {
        return weight_ == weight;
    }

    bool WeightSpaceVertex::hasUnitWeight(Polyscip::WeightType::size_type index) const {
        assert(index < weight_.size());
        for (auto i = 0; i != weight_.size(); ++i) {
            if (i == index) { // check if value one in both weights are at the same position
                if (weight_[i] != 1.)
                    return false;
            }
            else { // check if other elements of vertex weight are zero
                if (weight_[i] != 0.)
                    return false;
            }
        }
        return true;
    }


    void WeightSpaceVertex::print(ostream& os, bool printFacets) const {
        os << "WeightSpaceVertex:\n weight = [ ";
        std::ostream_iterator <ValueType> out_it(os, " ");
        std::copy(weight_.cbegin(), weight_.cend(), out_it);
        os << "]\n weighted objective value = " << weighted_obj_val_ << "\n";
        if (printFacets) {
            os << " defining facets: \n";
            for (const auto &facet : incident_facets_)
                facet->print(os);
        }
    }

}
