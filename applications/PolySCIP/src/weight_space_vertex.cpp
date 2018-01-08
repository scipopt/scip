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
 * @file weight_space_vertex.cpp
 * @brief Implements class representing vertex of (partial) weight space polyhedron
 * @author Sebastian Schenker
 * @author Timo Strunk
 *
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

    /**
     * Default constructor
     * @param incident_facets Incident facets of the constructed vertex
     * @param weight Lhs coefficients of the constructed vertex
     * @param weighted_obj_val Rhs coefficient of the constructed vertex
     * @param sort_facets if true, incident facets are sorted
     */
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

    /**
     * Returns weighted objective value of vertex
     * @return weighted objective value
     */
    ValueType WeightSpaceVertex::getCurrentWOV() const {
        return weighted_obj_val_;
    }


    /**
     * Constructor creating a new vertex from an obsolete and non-obsolete vertex via
     * the equality new_vertex = obs_coeff * obs + non_obs_coeff * non_obs
     * @param obs_coeff Coefficient of obsolete vertex
     * @param non_obs_coeff Coefficient of non-obsolete vertex
     * @param obs Obsolete vertex
     * @param non_obs Non-obsolete vertex
     * @param incident_facet Incident facet of new vertex
     * @param wsp_dimension Dimension of the corresponding weight space polyhedron
     */
    WeightSpaceVertex::WeightSpaceVertex(double obs_coeff,
                                         double non_obs_coeff,
                                         const WeightSpaceVertex* obs,
                                         const WeightSpaceVertex* non_obs,
                                         const std::shared_ptr<const WeightSpaceFacet>& incident_facet,
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

    /**
     * Returns associated weight vector of weight space vertex
     * @return Weight vector of vertex
     */
    WeightType WeightSpaceVertex::getWeight() const {
        return weight_;
    }

    /**
     * Outcome vector \\cdot weight vector
     * @param outcome Outcome vector
     * @return Scalarproduct of outcome and weight vector of vertex
     */
    double WeightSpaceVertex::getWeightedOutcome(const OutcomeType& outcome) const {
        assert (outcome.size() == weight_.size());
        return std::inner_product(begin(outcome),
                                  end(outcome),
                                  begin(weight_),
                                  0.);
    }

    /**
     * Computes slack
     * @param outcome Outcome vector
     * @param outcome_is_ray Indicates whether given outcome corresponds to ray
     * @return outcome \\cdot weight - weighted_obj_val if outcome corresponds to point;
     * else outcome \\cdot weight
     */
    double WeightSpaceVertex::computeSlack(const OutcomeType& outcome, bool outcome_is_ray) const {
        double slack = getWeightedOutcome(outcome);
        if (!outcome_is_ray)
            slack -= getCurrentWOV();
        return slack;
    }

    /**
     * Compute convex combination of weights
     * @param weight1 First weight vector
     * @param weight2 Second weight vector
     * @param h Coefficient for convex combination
     * @return h * weight1 + (1-h) * weight2
     */
    const WeightType WeightSpaceVertex::calculateWeightCombination(double h,
                                                                   const WeightType& weight1,
                                                                   const WeightType& weight2) {
        assert (weight1.size() == weight2.size());
        auto new_weight = WeightType(weight1.size(),0.);
        // new_weight = h*weight1 + (1-h)*weight2
        transform(begin(weight1), end(weight1), begin(weight2),
                  begin(new_weight), [h](ValueType w1, ValueType w2){return h*w1 + (1.-h)*w2;});
        return new_weight;
    }

    /**
     * Compare weight vectors
     * @param weight Weight to check against
     * @return true if given weight vector coincides with weight vector of vertex; otherwise false
     */
    bool WeightSpaceVertex::hasSameWeight(const WeightType& weight) const {
        return weight_ == weight;
    }

    /**
     * Indicates whether weight vector of vertex corresponds to some unit vector
     * @return true if weight vector of vertex corresponds to some unit vector; otherwise false
     */
    bool WeightSpaceVertex::hasUnitWeight() const {
        for (std::size_t i=0; i<weight_.size(); ++i) {
            auto weight = WeightType(weight_.size(), 0.);
            weight[i] = 1.;
            if (hasSameWeight(weight))
                return true;
        }
        return false;
    }

    /**
     * Indicates whether weight vector of vertex corresponds to zero vector
     * @return true if weight vector of vertex is zero vector; otherwise false
     */
    bool WeightSpaceVertex::hasZeroWeight() const {
        for (std::size_t i=0; i<weight_.size(); ++i) {
            if (weight_[i] != 0)
                return false;
        }
        return true;
    }

    /**
     * Print function
     * @param os Output stream to write to
     * @param printFacets Indicate whehter incident facets should be printed
     */
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
