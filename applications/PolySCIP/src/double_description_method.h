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

#include <utility>
#include <vector>

#include "objscip/objscip.h"
#include "polyscip_types.h"

namespace polyscip {

    class V_Representation {
    public:
        using H_RepT = std::pair<OutcomeType, ValueType>;
        using V_RepT = std::pair<WeightType, ValueType>;

        V_Representation(SCIP* scip, const ResultContainer& bounded, const ResultContainer& unbounded);

        void computeV_Rep();

        std::vector<V_RepT> getV_Rep() const {return v_rep_;};

    private:
        /** Computes initial v-representation for the following h-representation:
         * 1) bounded \cdot (w_1,...,w_k) - a >= 0
         * 2) w_1 >= 0
         * ...
         * k+1) w_k >= 0
         */
        void computeInitialRep(const OutcomeType& bounded);

        std::vector<std::size_t> computeZeroSet(const std::pair<WeightType, ValueType>& ray) = delete;

        std::vector<V_RepT> extendV_Rep(std::vector<V_RepT> current_v_rep, const H_RepT& ineq);

        V_RepT getNewVertex(const V_RepT& ray_plus, const V_RepT& ray_minus, const H_RepT& ineq) const;

        SCIP* scip_;
        std::vector<OutcomeType> bounded_;
        std::vector<OutcomeType> unbounded_;
        std::vector<H_RepT> current_h_rep_;
        std::vector<V_RepT> initial_v_rep_;
        std::vector<V_RepT> v_rep_;
    };

}
#endif //POLYSCIP_SRC_DOUBLE_DESCRIPTION_METHOD_H_INCLUDED