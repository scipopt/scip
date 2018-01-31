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
 * @file weight_space_facet.h
 * @brief Class representing a facet of the weight space polyhedron
 * @author Sebastian Schenker
 *
 */

#ifndef POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED
#define POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED

#include <iosfwd>
#include <memory>
#include <vector>
#include <tuple>

#include "polyscip_types.h"

namespace polyscip {

    /**
     * @class WeightSpaceFacet
     * @brief Class representing a facet of the (partial) weight space polyhedron
     * @details A facet is described by the inequality w_coeffs \\cdot w >= wov_coeff * a
     */
    class WeightSpaceFacet {
    public:

        /**
         * Simple class for less than comparison
         */
        struct Compare {
            /**
             * Callable less than comparator
             * @param f1 lhs facet
             * @param f2 rhs facet
             * @return true if lhs facet is lexicographically smaller than rhs facet
             */
            bool operator()(const std::shared_ptr<const WeightSpaceFacet> &f1,
                            const std::shared_ptr<const WeightSpaceFacet> &f2) {
                return std::tie(f1->wov_coeff_, f1->w_coeffs_) <
                       std::tie(f2->wov_coeff_, f2->w_coeffs_);
            }
        };


        /**
         *  Default constructor
         *  @param outcome Values for w_coeffs_
         *  @param wov_coeff Values for wov_coeff
         */
        explicit WeightSpaceFacet(OutcomeType outcome,
                                  ValueType wov_coeff);


        /**
         * Print function
         * @param os Output stream to print to
         */
        void print(std::ostream& os) const;

    private:
        OutcomeType w_coeffs_; ///< Corresponding lhs coefficients for facet inequality
        ValueType wov_coeff_; ///< Corresponding rhs coefficients for facet inequality
    };
}

#endif // POLYSCIP_SRC_WEIGHT_SPACE_FACET_H_INCLUDED
