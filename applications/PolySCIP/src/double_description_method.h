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

/** @brief  Double description method
 *
 * Yields v-representation from h-representation.
 */

#ifndef POLYSCIP_SRC_DOUBLE_DESCRIPTION_METHOD_H_INCLUDED
#define POLYSCIP_SRC_DOUBLE_DESCRIPTION_METHOD_H_INCLUDED

#include <iostream>
#include <map>
#include <ostream>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"

namespace polyscip {

    class VRepresentation {
    public:
        using H_RepT = std::pair<OutcomeType, ValueType>;
        using V_RepT = std::pair<WeightType, ValueType>;

        VRepresentation(SCIP* scip, const ResultContainer& bounded, const ResultContainer& unbounded);

        void computeVRep();

        std::vector<V_RepT> getVRep() const {return v_rep_;};

    private:
        using SlackMap = std::map<std::size_t, std::vector<std::size_t>>;

        /** Computes initial v-representation for the following h-representation:
         * 1) bounded \cdot (w_1,...,w_k) - a >= 0
         * 2) w_1 >= 0
         * ...
         * k+1) w_k >= 0
         */
        void computeInitialRep(const OutcomeType& bounded);

        std::vector<std::size_t> computeZeroSlackSet(const V_RepT& ray) const;

        bool rayPairIsAdjacent(std::size_t plus_index,
                               std::size_t minus_index,
                               const SlackMap& zero_slacks,
                               const std::vector<V_RepT>& current_v_rep) const;

        /** Check whether first parameter is multiple of second or third parameter
         */
        bool isMultiple(const V_RepT& ray, const V_RepT& ray2) const;

        std::vector<std::pair<std::size_t, std::size_t>> computeAdjacentPairs(const std::vector<std::size_t>& plus_indices,
                                                                              const std::vector<std::size_t>& minus_indices,
                                                                              const SlackMap& zero_slacks,
                                                                              const std::vector<V_RepT>& current_v_rep) const;

        std::vector<V_RepT> extendVRep(std::vector<V_RepT> current_v_rep, const H_RepT& new_constraint);

        V_RepT computeNewRay(const V_RepT& plus_ray, const V_RepT& minus_ray, const H_RepT& new_constraint) const;

        void normalizeVRep();
        void deleteZeroWeightRay();

        template<typename ContainerOfPairs>
        void printVRep(const ContainerOfPairs& container, std::ostream &os = std::cout) const;

        SCIP* scip_;
        std::vector<OutcomeType> bounded_;
        std::vector<OutcomeType> unbounded_;
        std::vector<H_RepT> current_h_rep_;
        std::vector<V_RepT> initial_v_rep_;
        std::vector<V_RepT> v_rep_;
    };

    template<typename ContainerOfPairs>
    void VRepresentation::printVRep(const ContainerOfPairs& container, std::ostream& os) const {
        for (const auto& elem : container) {
            os << "A_coeff = " << elem.second;
            global::print(elem.first, "Weight coeffs = ", os);
            os << "\n";
        }
    }

}
#endif //POLYSCIP_SRC_DOUBLE_DESCRIPTION_METHOD_H_INCLUDED