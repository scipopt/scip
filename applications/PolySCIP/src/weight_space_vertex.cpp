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
 * @brief  Weight space vertex
 * @author Sebastian Schenker, Timo Strunk
 *
 * Data structure storing combinatorial and geometric information
 * about a vertex of the weight space polyhedron. A weight space
 * vertex is represented by a weight 'w' and an weighted objective
 * value 'a' and is a vertex of the (partial) weight space polyhedron
 * P = {(w,a) \in \Lambda \times R : w \cdot y >= a \forall y \in Y'}
 * where Y' is the (current) set of non-dominated points and \Lambda
 * is the set of normalized weights
 */

#include "weight_space_vertex.h"

#include <algorithm> // std::copy, std::set_intersection, std::sort
#include <cstddef>
#include <iterator>  // std::back_inserter
#include <limits>
#include <memory> //std::shared
#include <numeric>   // std::inner_product;
#include <ostream>

#include "global_functions.h" //global::print
#include "polyscip_types.h"
#include "weight_space_facet.h"
#include "weight_space_polyhedron.h"

using std::ostream;
using std::shared_ptr;
using std::sort;
using std::vector;

namespace polyscip {


    WeightSpaceVertex::WeightSpaceVertex(WeightSpacePolyhedron::FacetContainer incident_facets,
                                         WeightType weight,
                                         ValueType weighted_obj_val,
                                         bool sort_facets)
            : vertex_status_(VertexStatus::unmarked),
              incident_facets_(std::move(incident_facets)),
              weight_(std::move(weight)),
              weighted_obj_val_{weighted_obj_val} {
        if (sort_facets) // sort facets in order to be able to use std::set_intersection in other constructor
            sort(begin(incident_facets_), end(incident_facets_), WeightSpaceFacet::Compare());
    }

    ValueType WeightSpaceVertex::getCurrentWOV() const {
        return weighted_obj_val_;
    }

    WeightSpaceVertex::WeightSpaceVertex(double obs_coeff,
                                         double non_obs_coeff,
                                         const WeightSpaceVertex* obs,
                                         const WeightSpaceVertex* non_obs,
                                         const shared_ptr<const WeightSpaceFacet>& incident_facet,
                                         std::size_t wsp_dimension)
            : vertex_status_(VertexStatus::unmarked)
    {
        assert (obs != non_obs);
        // get intersection of facets of obs and non_obs
        std::set_intersection(obs->incident_facets_.cbegin(),
                              obs->incident_facets_.cend(),
                              non_obs->incident_facets_.cbegin(),
                              non_obs->incident_facets_.cend(),
                              std::back_inserter(incident_facets_),
                              WeightSpaceFacet::Compare());
        auto upper_it = std::upper_bound(begin(incident_facets_),
                                         end(incident_facets_),
                                         incident_facet, WeightSpaceFacet::Compare());
        incident_facets_.insert(upper_it, incident_facet);
        assert(incident_facets_.size() >= wsp_dimension);

        std::transform(begin(obs->weight_), end(obs->weight_),
                       begin(non_obs->weight_), std::back_inserter(weight_),
                       [obs_coeff, non_obs_coeff](ValueType obs_val, ValueType non_obs_val) {
                           return obs_coeff * obs_val - non_obs_coeff * non_obs_val;
                       });
        auto normalize_val = *(std::max_element(begin(weight_), end(weight_)));
        assert (normalize_val > 0);
        std::transform(begin(weight_), end(weight_), begin(weight_),
                       [normalize_val](const ValueType &val) { return val / normalize_val; });
        weighted_obj_val_ = obs_coeff*obs->weighted_obj_val_ - non_obs_coeff*non_obs->weighted_obj_val_; // return m_coeff * ray_minus - p_coeff * ray_plus
        weighted_obj_val_ /= normalize_val;
    }

    WeightType WeightSpaceVertex::getWeight() const {
        return weight_;
    }

    OutcomeType WeightSpaceVertex::getIncFacetsBounds(std::function<ValueType()> limit,
                                                      std::function<ValueType(const ValueType&, const ValueType&)> cmp) const {
        auto bounds = OutcomeType(weight_.size(), limit());
        for (const auto& f : incident_facets_) {
            if (f->hasNonZeroWOVCoeff()) {
                std::transform(f->coeffsBegin(),
                               f->coeffsEnd(),
                               bounds.begin(),
                               bounds.begin(),
                               cmp);
            }
        }
        return bounds;
    }

    OutcomeType WeightSpaceVertex::getIncFacetsLowerBounds() const {
        return getIncFacetsBounds([](){return std::numeric_limits<ValueType>::max();},
                                  [](const ValueType& val1, const ValueType& val2){return std::min(val1, val2);});
    }

    OutcomeType WeightSpaceVertex::getIncFacetsUpperBounds() const {
        return getIncFacetsBounds([](){return std::numeric_limits<ValueType>::lowest();},
                                  [](const ValueType& val1, const ValueType& val2){return std::max(val1, val2);});
    }

    double WeightSpaceVertex::getWeightedOutcome(const OutcomeType& outcome) const {
        assert (outcome.size() == weight_.size());
        return std::inner_product(begin(outcome),
                                  end(outcome),
                                  begin(weight_),
                                  0.);
    }

    double WeightSpaceVertex::computeSlack(const OutcomeType& outcome, bool outcome_is_ray) const {
        double slack = getWeightedOutcome(outcome);
        if (!outcome_is_ray)
            slack -= getCurrentWOV();
        return slack;
    }

    const WeightType WeightSpaceVertex::calculateWeightCombination(double h,
                                                                   const WeightType& weight1,
                                                                   const WeightType& weight2) {
        assert (weight1.size() == weight2.size());
        auto new_weight = WeightType(weight1.size(),0.);
        // set new_weight = h*weight1 + (1-h)*weight2
        transform(begin(weight1), end(weight1), begin(weight2),
                  begin(new_weight), [h](ValueType w1, ValueType w2){return h*w1 + (1.-h)*w2;});
        return new_weight;
    }

    bool WeightSpaceVertex::hasSameWeight(const WeightType& weight) const {
        return weight_ == weight;
    }

    bool WeightSpaceVertex::hasUnitWeight() const {
        for (std::size_t i=0; i<weight_.size(); ++i) {
            auto weight = WeightType(weight_.size(), 0.);
            weight[i] = 1.;
            if (hasSameWeight(weight))
                return true;
        }
        return false;
    }

    bool WeightSpaceVertex::hasZeroWeight() const {
        for (std::size_t i=0; i<weight_.size(); ++i) {
            if (weight_[i] != 0)
                return false;
        }
        return true;
    }

    void WeightSpaceVertex::print(ostream& os, bool printFacets) const {
        global::print(weight_, "WeightSpaceVertex: weight = [", "]", os);
        os << "\n wov = " << weighted_obj_val_ << "\n";
        if (printFacets) {
            os << " defining facets: \n";
            for (const auto &facet : incident_facets_)
                facet->print(os);
        }
        os << "Vertex-Status: ";
        switch (vertex_status_) {
            case VertexStatus::marked:
                os << "marked\n";
                break;
            case VertexStatus::obsolete:
                os << "obsolete\n";
                break;
            case VertexStatus::unmarked:
                os << "unmarked\n";
                break;
            case VertexStatus::special:
                os << "special\n";
                break;
        }
    }

}
