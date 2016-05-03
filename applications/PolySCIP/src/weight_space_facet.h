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

/** @brief Class representing a facet of the weight space polyhedron
 *
 * Data structure representing a facet of the (partial) weight space
 * polyhedron P={(w,a) : w \cdot y >= a \forall y \in Y} where Y is
 * the (current) set of non-dominated points. A facet (lhs,rhs) is
 * represented by coefficients 'lhs_' and a right hand side 'rhs_'
 * yielding an inequality of the form lhs_ \cdot w >= rhs_ 
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED
#define POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED

#include <iosfwd>
#include <vector>

#include "polyscip.h"

namespace polyscip {

    /** Facet of the (partial) weight space polyhedron. */
    class WeightSpaceFacet {
    public:
        bool friend operator<(const WeightSpaceFacet& facet1,
                              const WeightSpaceFacet& facet2);

        /** Creates the facet point \cdot w >= weighted_obj_val
         *  @param point computed (weakly non-dominated) point in objective space
         *  @param weighted_obj_val weighted objective value of point
         */
        explicit WeightSpaceFacet(const Polyscip::OutcomeType& point,
                                  Polyscip::ValueType weighted_obj_val);

        /** Creates the weight space facet w_i >= 0
         *  @param num_objs number of objectives of given problem
         *  @param index index i of w_i >= 0
         */
        explicit WeightSpaceFacet(unsigned num_objs, unsigned index);

        /** Prints facet information to output stream.
         */
        void print(std::ostream& os) const;

    private:
        /**< left hand side coefficients of the facet inequality */
        std::vector<Polyscip::ValueType> lhs_;
        /**< right hand side value of the facet inequality */
        Polyscip::ValueType rhs_;
    };

}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED
