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

#include <algorithm> // std::copy, std::set_intersection, std::sort
#include <iterator>  // std::ostream_iterator, std::back_inserter
#include <memory> //std::shared
#include <numeric>   // std::inner_product;
#include <ostream>

#include "polyscip.h"
#include "weight_space_facet.h"
#include "weight_space_polyhedron.h"

using std::ostream;
using std::shared_ptr;
using std::sort;
using std::vector;

namespace polyscip {

    bool compare_facet_ptr(const std::shared_ptr<const WeightSpaceFacet>& f1,
                           const std::shared_ptr<const WeightSpaceFacet>& f2) {
        return *f1 < *f2;
    }

    using OutcomeType = Polyscip::OutcomeType;
    using ValueType = Polyscip::ValueType;
    using WeightType = Polyscip::WeightType;

    WeightSpaceVertex::WeightSpaceVertex(WeightSpacePolyhedron::FacetContainer incident_facets,
                                         WeightType weight,
                                         ValueType weighted_obj_val,
                                         bool sort_facets)
            : incident_facets_(std::move(incident_facets)),
              weight_(std::move(weight)),
              weighted_obj_val_{weighted_obj_val} {
        if (sort_facets) // sort facets in order to be able to use std::set_intersection in other constructor
            sort(begin(incident_facets_), end(incident_facets_), compare_facet_ptr);
    }

    ValueType WeightSpaceVertex::getWOV() const {
        return weighted_obj_val_;
    }

    WeightSpaceVertex::WeightSpaceVertex(const WeightSpaceVertex *obs,
                                         const WeightSpaceVertex *non_obs,
                                         const Polyscip::OutcomeType &outcome,
                                         bool outcome_is_ray) {
        // get intersection of facets of obs and non_obs
        std::set_intersection(obs->incident_facets_.cbegin(),
                              obs->incident_facets_.cend(),
                              non_obs->incident_facets_.cbegin(),
                              non_obs->incident_facets_.cend(),
                              std::back_inserter(incident_facets_),
                              compare_facet_ptr);
        assert(incident_facets_.size() + 1 == obs->incident_facets_.size());
        // add additional facet with respect to outcome
        auto wov_coeff = outcome_is_ray ? 0.0 : 1.0;
        auto new_facet = std::make_shared<const WeightSpaceFacet>(outcome, wov_coeff);
        auto upper_it = std::upper_bound(begin(incident_facets_),
                                         end(incident_facets_),
                                         new_facet, compare_facet_ptr);
        incident_facets_.insert(upper_it, std::move(new_facet));
        if (outcome_is_ray) {
            assert (obs->isMadeObsolete(outcome, 0.));
            assert (!non_obs->isMadeObsolete(outcome, 0.));
        }
        else {
            assert (obs->isMadeObsolete(outcome));
            assert (!non_obs->isMadeObsolete(outcome));
        }
        // compute convex combination of weights of obsolete vertex and non-obsolete vertex
        auto weight_non_obs = non_obs->getWeight();
        auto weight_obs = obs->getWeight();
        auto h = calculateCombinationValue(weight_non_obs, weight_obs, outcome);
        assert (0. < h && h < 1.0);
        if (outcome_is_ray) // shift combination towards non-obsolete vertex
            h += 1e-7;
        weight_ = calculateWeightCombination(std::move(weight_non_obs),
                                             std::move(weight_obs),
                                             h);
        // computed weighted objective value
        weighted_obj_val_ = h*non_obs->getWOV() + (1.0-h)*obs->getWOV();
    }

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

    ValueType  WeightSpaceVertex::calculateCombinationValue(const WeightType& weight1,
                                                            const WeightType& weight2,
                                                            const OutcomeType& outcome) {
        assert(weight1.size() == weight2.size());
        // h = \frac{-weight2 \cdot outcome}{weight1 \cdot outcome - weight2 \cdot outcome}
        ValueType numerator = -1.0 * inner_product(begin(weight2),
                                                   end(weight2),
                                                   begin(outcome),
                                                   0.);
        ValueType denominator = inner_product(begin(weight1),
                                              end(weight1),
                                              begin(outcome),
                                              0.)
                                + numerator;
        assert(denominator != 0.);
        return numerator / denominator;
    }


    WeightType WeightSpaceVertex::calculateWeightCombination(WeightType weight1,
                                                             WeightType weight2,
                                                             ValueType h) {
        assert (weight1.size() == weight2.size());
        // set weight1 = h*weight1
        transform(begin(weight1), end(weight1), begin(weight1),
                  [h](ValueType w){return h*w;});
        // set weight2 = (1-h)*weight2
        transform(begin(weight2), end(weight2), begin(weight2),
                  [h](ValueType w){return (1.0-h)*w;});
        // set weight1 = weight1 + weight2
        transform(begin(weight1), end(weight1), begin(weight2),
                  begin(weight1), std::plus<ValueType>());
        return weight1;
    }

    bool WeightSpaceVertex::isMadeObsolete(const OutcomeType& outcome) const {
        return isMadeObsolete(outcome, weighted_obj_val_);
    }

    bool WeightSpaceVertex::hasSameWeight(const WeightType& weight) {
        return weight_ == weight;
    }

    bool WeightSpaceVertex::hasUnitWeight(unsigned index) {
        assert(index < weight_.size());
        auto weight = WeightType(weight_.size(),0.);
        weight[index] = 1.;
        return weight_ == weight;
    }

    void WeightSpaceVertex::print(ostream& os, bool printFacets) const {
        os << "WeightSpaceVertex:\n weight = [ ";
        std::ostream_iterator <ValueType> out_it(os, " ");
        std::copy(weight_.cbegin(), weight_.cend(), out_it);
        os << "]\n wov = " << weighted_obj_val_ << "\n";
        if (printFacets) {
            os << " defining facets: \n";
            for (const auto &facet : incident_facets_)
                facet->print(os);
        }
    }

}
