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
#include "weight_space_facet.h"

namespace polyscip {

    class DoubleDescription {
    public:
        using H_RepT = std::pair<OutcomeType, ValueType>;
        using V_RepT = std::pair<WeightType, ValueType>;
        using V_RepContainer = std::vector<V_RepT>;
        using H_RepContainer = std::vector<H_RepT>;

        const double kLimitForNormalization = 1e+6;

        //todo computedAdjacency!!!

        DoubleDescription(SCIP* scip, const ResultContainer& bounded, const ResultContainer& unbounded);

        void computeVRep();

        void printVRep(std::ostream &os = std::cout) const {printRep(v_rep_, os);};
        std::size_t size() const {return v_rep_.size();};

        V_RepContainer getVRep() const {return v_rep_;};
        H_RepContainer getHRep() const {return h_rep_;};

        V_RepContainer&& moveVRep() {return std::move(v_rep_);};
        H_RepContainer&& moveHRep() {return std::move(h_rep_);};

    private:
        using SlackContainer = std::vector<std::vector<std::size_t>>;

        bool shouldNormalize(const H_RepT& constraint, const V_RepContainer& current_v_rep) const;

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
                               const SlackContainer& zero_slacks,
                               const std::vector<V_RepT>& current_v_rep) const;

        /** Check whether first parameter is multiple of second or third parameter
         */
        bool isMultiple(const V_RepT& ray, const V_RepT& ray2) const;

        std::vector<std::pair<std::size_t, std::size_t>> computeAdjacentPairs(const std::vector<std::size_t>& plus_indices,
                                                                              const std::vector<std::size_t>& minus_indices,
                                                                              const SlackContainer& zero_slacks,
                                                                              const std::vector<V_RepT>& current_v_rep) const;

        std::vector<V_RepT> extendVRep(std::vector<V_RepT> current_v_rep, const H_RepT& new_constraint);

        V_RepT computeNewRay(const V_RepT& plus_ray, const V_RepT& minus_ray, const H_RepT& new_constraint) const;

        void normalizeVRep(V_RepContainer& v_rep);
        void deleteZeroWeightRay();

        template<typename ContainerOfPairs>
        void printRep(const ContainerOfPairs& container, std::ostream &os = std::cout) const;

        SCIP* scip_;
        std::vector<OutcomeType> bounded_;
        std::vector<OutcomeType> unbounded_;
        H_RepContainer h_rep_;
        V_RepContainer initial_v_rep_;
        V_RepContainer v_rep_;

    };

    template<typename ContainerOfPairs>
    void DoubleDescription::printRep(const ContainerOfPairs& container, std::ostream& os) const {
        for (const auto& elem : container) {
            os << "A_coeff = " << elem.second;
            global::print(elem.first, " Weight coeffs = ", os);
            os << "\n";
        }
    }

}
#endif //POLYSCIP_SRC_DOUBLE_DESCRIPTION_METHOD_H_INCLUDED