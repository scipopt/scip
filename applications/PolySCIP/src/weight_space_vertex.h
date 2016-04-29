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

/** @brief  Class representing a weight space vertex
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

#include <iosfwd>
#include <memory> // std::shared_ptr
#include <vector>

#include "polyscip.h"
#include "weight_space_facet.h"

namespace polyscip {

    /** Vertex of the weight space polyhedron. */
    class WeightSpaceVertex {

    public:
        using FacetContainer = std::vector<std::shared_ptr<const WeightSpaceFacet>>;

        /** Creates a vertex of the (partial) weight space polyhedron.
         * @param incident_facets Facets defining of the weight space polyhedron defining the vertex
         * @param weight Corresponding weight
         * @param weighted_obj_val Corresponding maximal weight objective val in weight space polyhedron
         */
        explicit WeightSpaceVertex(FacetContainer incident_facets,
                                   Polyscip::WeightType weight,
                                   Polyscip::ValueType weighted_obj_);

        /** Checks whether given outcome makes vertex obsolete, i.e. whether
         * vertex weight \cdot outcome >= rhs
         * @param outcome point or ray in objective space
         * @param rhs value to compare
         * @return true if vertex is made obsolete, false otherwise
         */
        bool isMadeObsolete(const Polyscip::OutcomeType& outcome, Polyscip::ValueType rhs) const;

        /** Checks whether given outcome makes vertex obsolete, i.e. whether
         * vertex weight \cdot outcome >= weighted objective value of vertex
         * @param outcome point or ray in objective space
         * @return true if vertex is made obsolete, false otherwise
         */
        bool isMadeObsolete(const Polyscip::OutcomeType& outcome) const;

        /** Return associated weight of weight space vertex
         * @return weight of vertex
         */
        Polyscip::WeightType getWeight() const;

        /** Checks whether weight of vertex corresponds to unit weight
         * @param index index of 1 in unit weight
         * @return true if weight of vertex is unit weight with 1 at index; false otherwise
         */
        bool hasUnitWeight(Polyscip::WeightType::size_type index) const;

        /** Checks whether weight of vertex corresponds with given weight
         * @param weight weight to check against
         */
        bool hasSameWeight(const Polyscip::WeightType& weight);

        /** Prints weight space vertex information to output stream.
         * @param printFacets if true, then defining facets are printed
         */
        void print(std::ostream& os, bool printFacets = false) const;

    private:
        /**< incident facets */
        FacetContainer incident_facets_;
        /**< used weight */
        Polyscip::WeightType weight_;
        /**< corresponding weighted objective value */
        Polyscip::ValueType weighted_obj_val_;

    };

}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_VERTEX_H_INCLUDED
