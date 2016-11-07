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

    /** Vertex of the weight space polyhedron. */
    class WeightSpaceVertex {

    public:
        enum class VertexStatus {marked, obsolete, unmarked, special};

        /** Creates a vertex of the (partial) weight space polyhedron.
         * @param incident_facets Facets defining of the weight space polyhedron defining the vertex
         * @param weight Corresponding weight
         * @param weighted_obj_val Corresponding maximal weight objective val in weight space polyhedron
         * @param sort_facets if true, incident facets are sorted
         */
        explicit  WeightSpaceVertex(WeightSpacePolyhedron::FacetContainer incident_facets,
                                    WeightType weight,
                                    ValueType weighted_obj_val,
                                    bool sort_facets = true);

        explicit WeightSpaceVertex(double convCombValue,
                                   const WeightSpaceVertex* obs,
                                   const WeightSpaceVertex* non_obs,
                                   const OutcomeType& outcome,
                                   bool outcome_is_ray,
                                   std::size_t wsp_dimension) = delete;

        explicit WeightSpaceVertex(double obs_coeff,
                                   double non_obs_coeff,
                                   const WeightSpaceVertex* obs,
                                   const WeightSpaceVertex* non_obs,
                                   const std::shared_ptr<const WeightSpaceFacet>& incident_facet,
                                   std::size_t wsp_dimension);

        /** Return associated weight of weight space vertex
         * @return weight of vertex
         */
        WeightType getWeight() const;

        /** Return weighted objective value of associated vertex
         * @return weighted objective value
         */
        ValueType getCurrentWOV() const;

        double getWeightedOutcome(const OutcomeType& outcome) const;

        double computeSlack(const OutcomeType& outcome, bool outcome_is_ray) const;

        OutcomeType getIncFacetsBounds(std::function<ValueType()> limit,
                                       std::function<ValueType(const ValueType&, const ValueType&)> cmp) const;
        OutcomeType getIncFacetsLowerBounds() const;
        OutcomeType getIncFacetsUpperBounds() const;

        /** Checks whether weight of vertex corresponds to unit weight
         * @param index index of 1 in unit weight
         * @return true if weight of vertex is unit weight with 1 at index; false otherwise
         */
        bool hasUnitWeight() const;
        bool hasZeroWeight() const;

        VertexStatus getStatus() const {return vertex_status_;};
        void setStatus(VertexStatus status) {vertex_status_ = status;};

        /** Checks whether weight of vertex corresponds with given weight
         * @param weight weight to check against
         */
        bool hasSameWeight(const WeightType& weight) const;

        /** Prints weight space vertex information to output stream.
         * @param printFacets if true, then defining facets are printed
         */
        void print(std::ostream& os, bool printFacets = false) const;

    private:

        friend bool WeightSpacePolyhedron::areAdjacent(const WeightSpaceVertex* v, const WeightSpaceVertex* w);
        friend double WeightSpacePolyhedron::calculateConvexCombValue(const WeightSpaceVertex* obs,
                                                                      const WeightSpaceVertex* non_obs,
                                                                      const OutcomeType& outcome,
                                                                      bool outcome_is_ray);


        /** Returns the weight h * weight1 + (1-h) * weight2
         * @param weight1 weight of vertex
         * @param weight2 weight of another vertex
         * @param h combination coefficient
         * @return convex combination h * weight1 + (1-h) * weight2
         */
        static const WeightType calculateWeightCombination(double h,
                                                           const WeightType& weight1,
                                                           const WeightType& weight2);


        VertexStatus vertex_status_;
        /**< incident facets */
        WeightSpacePolyhedron::FacetContainer incident_facets_;
        /**< used weight */
        WeightType weight_;
        /**< corresponding weighted objective value */
        ValueType weighted_obj_val_;
    };

}

#endif //POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED
