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
 * @file weight_space_vertex.h
 * @brief Class representing a vertex of the (partial) weight space polyhedron
 * @author Sebastian Schenker
 * @author Timo Strunk
 *
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED
#define POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED

#include <cstddef>
#include <functional>
#include <iosfwd>
#include <memory> // std::shared_ptr
#include <vector>

#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_polyhedron.h"

namespace polyscip {

    /**
     * @class WeightSpaceVertex
     * @brief Class representing a vertex of the (partial) weight space polyhedron
     * @details A vertex of the (partial) weight space polyhedron is given
     * by a weight vector 'weight_' and an weighted objective value 'weighted_obj_val' and
     * is active for a certain number of inequalities of the form
     * 'weight_ \\cdot y >= weighted_obj_val' where y are non-dominated extreme points
     * of the considered problem
     */
    class WeightSpaceVertex {

    public:
        /**
         * Different statuses a weight space vertex can take
         */
        enum class VertexStatus {marked, obsolete, unmarked, special};

        /**
         * Default constructor
         * @param incident_facets Incident facets of the constructed vertex
         * @param weight Lhs coefficients of the constructed vertex
         * @param weighted_obj_val Rhs coefficient of the constructed vertex
         * @param sort_facets if true, incident facets are sorted
         */
        explicit  WeightSpaceVertex(WeightSpacePolyhedron::FacetContainer incident_facets,
                                    WeightType weight,
                                    ValueType weighted_obj_val,
                                    bool sort_facets = true);

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
        explicit WeightSpaceVertex(double obs_coeff,
                                   double non_obs_coeff,
                                   const WeightSpaceVertex* obs,
                                   const WeightSpaceVertex* non_obs,
                                   const std::shared_ptr<const WeightSpaceFacet>& incident_facet,
                                   std::size_t wsp_dimension);

        /**
         * Returns associated weight vector of weight space vertex
         * @return Weight vector of vertex
         */
        WeightType getWeight() const;

        /**
         * Returns weighted objective value of vertex
         * @return weighted objective value
         */
        ValueType getCurrentWOV() const;

        /**
         * Outcome vector \\cdot weight vector
         * @param outcome Outcome vector
         * @return Scalarproduct of outcome and weight vector of vertex
         */
        double getWeightedOutcome(const OutcomeType& outcome) const;

        /**
         * Computes slack
         * @param outcome Outcome vector
         * @param outcome_is_ray Indicates whether given outcome corresponds to ray
         * @return outcome \\cdot weight - weighted_obj_val if outcome corresponds to point;
         * else outcome \\cdot weight
         */
        double computeSlack(const OutcomeType& outcome, bool outcome_is_ray) const;


        /**
         * Indicates whether weight vector of vertex corresponds to some unit vector
         * @return true if weight vector of vertex corresponds to some unit vector; otherwise false
         */
        bool hasUnitWeight() const;

        /**
         * Indicates whether weight vector of vertex corresponds to zero vector
         * @return true if weight vector of vertex is zero vector; otherwise false
         */
        bool hasZeroWeight() const;

        /**
         * Get vertex status
         * @return Status of vertex
         */
        VertexStatus getStatus() const {return vertex_status_;};

        /**
         * Set vertex status
         * @param status New status of vertex
         */
        void setStatus(VertexStatus status) {vertex_status_ = status;};

        /**
         * Compare weight vectors
         * @param weight Weight to check against
         * @return true if given weight vector coincides with weight vector of vertex; otherwise false
         */
        bool hasSameWeight(const WeightType& weight) const;

        /**
         * Print function
         * @param os Output stream to write to
         * @param printFacets Indicate whehter incident facets should be printed
         */
        void print(std::ostream& os, bool printFacets = false) const;

    private:

        /**
         * @relates WeightSpacePolyhedron
         * @param v
         * @param w
         * @return
         */
        friend bool WeightSpacePolyhedron::areAdjacent(const WeightSpaceVertex* v,
                                                       const WeightSpaceVertex* w);

        /**
         * Compute convex combination of weights
         * @param weight1 First weight vector
         * @param weight2 Second weight vector
         * @param h Coefficient for convex combination
         * @return h * weight1 + (1-h) * weight2
         */
        static const WeightType calculateWeightCombination(double h,
                                                           const WeightType& weight1,
                                                           const WeightType& weight2);


        VertexStatus vertex_status_; ///< Corresponding status of vertex
        WeightSpacePolyhedron::FacetContainer incident_facets_; ///< Incident facets of vertex
        WeightType weight_; ///< Corresponding weight vector of vertex
        ValueType weighted_obj_val_; ///< Corresponding weighted objective value of vertex
    };

}

#endif //POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED
