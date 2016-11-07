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
 *  @author Sebastian Schenker
 *
 * Defines classes for double description method.
 */

#ifndef POLYSCIP_SRC_DOUBLE_DESCRIPTION_H_INCLUDED
#define POLYSCIP_SRC_DOUBLE_DESCRIPTION_H_INCLUDED

#include <algorithm>
#include <bitset>
#include <cstddef>
#include <functional> // std::reference_wrapper
#include <iostream>
#include <memory>
#include <unordered_map>
#include <ostream>
#include <tuple> // std::tie
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "PolySCIPConfig.h" // defines MAX_NO_OBJS
#include "polyscip_types.h"
#include "weight_space_facet.h"

namespace polyscip {

    namespace doubledescription {

        using H_RepT = std::pair<OutcomeType, ValueType>;
        using H_RepC = std::vector<H_RepT>;

        class V_RepT {
        public:
            using SlackMap = std::unordered_map<std::size_t, ValueType>;
            constexpr static std::size_t kMaxInitialHrepSize = 2*POLYSCIP_MAX_NO_OBJS;

            friend class DoubleDescriptionMethod;

            V_RepT(WeightType weight, ValueType wov) = delete;

            explicit V_RepT(SCIP* scip, WeightType&& weight, ValueType&& wov, const H_RepC& current_h_rep);

            explicit V_RepT(SCIP* scip,
                            const V_RepT& plus,
                            const V_RepT& minus,
                            std::size_t index_of_constraint_in_hrep,
                            const H_RepC& current_h_rep);

            inline bool operator==(const V_RepT& rhs) const;

            inline bool operator!=(const V_RepT& rhs) const;

            void print(std::ostream& os, bool withIncidentFacets, const H_RepC& h_rep) const;

            bool hasNonZeroWeight() const;

            //todo incorporate way for facets when WSP is not full-dimensional

            /** Checks whether given parameter is a subset of member zero_slack_indices_
             * @Return true if member zero_slack_indices_ is superset of given parameter; false otherwise
             */
            bool hasZeroSlackSuperSet(const std::vector<std::size_t>& indices) const = delete;
            bool hasZeroIndsSuperSet(const std::bitset<kMaxInitialHrepSize>& common_zero_inds) const;

            std::size_t getMinInfeasIndex() const;
            bool isZeroSlackIndex(std::size_t index) const;

            WeightType&& moveWeight() {return std::move(weight_);};
            //WeightType getWeight() const {return weight_;};
            ValueType getWov() const {return wov_;};

        private:
            //todo weight space phase ändert sich abhängig von Größe!! Testen warum!
            constexpr static double kNormalizingThreshold = 1e+5;
            bool shouldNormalize(double threshold) const;
            void normalize(double normalizing_val);
            ValueType getSlack(std::size_t index) const;

            void setSlacksAndMinInfeasInd(SCIP* scip, const H_RepC& h_rep);

            WeightType weight_;
            ValueType wov_;
            std::pair<bool, std::size_t> min_infeas_ind_;
            SlackMap inds_to_slacks_;
            std::bitset<kMaxInitialHrepSize> zero_slacks_;
        };

        using V_RepC = std::vector<std::shared_ptr<V_RepT>>;


        class DoubleDescriptionMethod {
        public:
            using AdjPair = std::pair<std::reference_wrapper<const V_RepT>, std::reference_wrapper<const V_RepT>>;
            using AdjPairContainer = std::vector<AdjPair>;

            explicit DoubleDescriptionMethod(SCIP *scip, std::size_t no_all_obj, const ResultContainer &bounded, const ResultContainer &unbounded);

            void computeVRep();

            void computeVRep_Var1();

            void printVRep(std::ostream &os = std::cout, bool withIncidentFacets = false) const;

            std::size_t size() const { return v_rep_.size(); };

            H_RepC getHRep() {return h_rep_;};
            V_RepC getVRep() {return v_rep_;};

            H_RepC&& moveHRep() {return std::move(h_rep_);};
            V_RepC&& moveVRep() {return std::move(v_rep_);};

        private:
            using SlackContainer = std::vector<std::vector<std::size_t>>;

            enum class VarOrder {keep_var_order, change_var_order};

            /** Computes initial v-representation for the following h-representation:
             * 1) bounded \cdot (w_1,...,w_k) - a >= 0
             * 2) w_1 >= 0
             * ...
             * k+1) w_k >= 0
             */
            //void computeInitialRep(const OutcomeType &bounded) = delete;
            V_RepC computeInitialVRep() const;

            //std::vector<std::size_t> getCommonZeroSlackIndices(const V_RepT& v, const V_RepT& w) const;

            std::bitset<V_RepT::kMaxInitialHrepSize> getCommonZeroSlackIndices(const V_RepT &v, const V_RepT &w) const;

            std::tuple<bool, VarOrder, std::size_t> minInfeasCondition(const V_RepT& r1, const V_RepT& r2) const;

            void applyInfeasCondition(const V_RepT& r1,
                                      const V_RepT& r2,
                                      const V_RepC& current_v_rep,
                                      std::size_t index,
                                      bool with_adjacency_test);

            /* see function description in DoubleDescriptionRevisited
             */
            void conditionalStoreEdge(const V_RepT& r1,
                                      const V_RepT& r2,
                                      std::size_t k,
                                      std::size_t i,
                                      const V_RepC& v_rep,
                                      bool with_adjacency_test);

            /*bool rayPairIsAdjacent(std::size_t plus_index,
                                       std::size_t minus_index,
                                       const V_RepC& current_v_rep) const;*/

            bool rayPairIsAdjacent(const V_RepT& r1,
                                   const V_RepT& r2,
                                   const V_RepC& v_rep) const;

            /** Check whether first parameter is multiple of second parameter
             */
            bool isMultiple(const V_RepT& v, const V_RepT& w) const;
            bool weightIsMultiple(SCIP* scip, double v_multiple, const V_RepT& v, const V_RepT& w) const;

            /*std::vector<std::pair<std::size_t, std::size_t>> computeAdjacentPairs(const std::vector<std::size_t>& plus_inds,
                                                                                      const std::vector<std::size_t>& minus_inds,
                                                                                      const V_RepC& current_rep) const;*/


            AdjPairContainer computeAdjacentPairs(const V_RepC& plus,
                                                  const V_RepC& minus,
                                                  const V_RepC& current_v_rep) const;

            //V_RepC extendVRep(V_RepC&& current_v_rep);
            V_RepC extendVRep(V_RepC&& current_v_rep);

            V_RepC extendVRep_Var1(V_RepC&& current_v_rep);

            /*V_RepC extendVRep(std::vector<V_RepT> current_rep,
                                      const H_RepT &constraint,
                                      std::size_t index_of_constraint_in_hrep) = delete;*/

            SCIP *scip_;
            std::size_t outcome_dimension_;
            std::size_t current_hrep_index_;
            H_RepC h_rep_;
            V_RepC v_rep_;
            std::vector< AdjPairContainer > adj_pairs_;
        };
    }
}
#endif //POLYSCIP_SRC_DOUBLE_DESCRIPTION_H_INCLUDED