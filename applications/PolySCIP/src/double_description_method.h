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
 * @file double_description_method.h
 * @brief  Double description method for transforming a polyhedron given
 * via its h-representation into its v-representation
 * @author Sebastian Schenker
 *
 * The underlying algorithm is described in the paper: "Double description method revisited" by
 * Komei Fukuda, Alain Prodon
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

        using H_RepT = std::pair<OutcomeType, ValueType>; ///< Type for element of h-representation
        using H_RepC = std::vector<H_RepT>; ///< Container for h-representations

        /**
         * @class V_RepT
         * @brief Class for element of v-representation
         */
        class V_RepT {
        public:

            using SlackMap = std::unordered_map<std::size_t, ValueType>; ///< Maps indices to slacks
            constexpr static std::size_t kMaxInitialHrepSize = 2*POLYSCIP_MAX_NO_OBJS; ///< Maximal initial size of h-representation
            friend class DoubleDescriptionMethod;

            /**
             * Constructor
             * @param scip SCIP pointer
             * @param weight Corresponding weight for the constructed vertex
             * @param wov Corresponding weighted objective value for the constructed vertex
             * @param h_rep Current h-representation
             */
            explicit V_RepT(SCIP* scip,
                            WeightType&& weight,
                            ValueType&& wov,
                            const H_RepC& h_rep);

            /**
             * Constructor
             * @param scip SCIP pointer
             * @param plus
             * @param minus
             * @param index_of_ineq
             * @param h_rep
             */
            explicit V_RepT(SCIP* scip,
                            const V_RepT& plus,
                            const V_RepT& minus,
                            std::size_t index_of_ineq,
                            const H_RepC& h_rep);

            /**
             * Equal operator
             * @param rhs v-representation object to compare with
             * @return true if objects are considered equal; otherwise false
             */
            inline bool operator==(const V_RepT& rhs) const;

            /**
             * Not equal operator
             * @param rhs v-representation object to compare with
             * @return true if objects are considered not equal; otherwise false
             */
            inline bool operator!=(const V_RepT& rhs) const;

            /**
             * Print function
             * @param os Output stream to print to
             * @param withIncidentFacets Boolean indicating whether incident facets should also be printed
             * @param h_rep
             */
            void print(std::ostream& os,
                       bool withIncidentFacets,
                       const H_RepC& h_rep) const;

            /**
             * Indicates whether corresponding weight is neither unit vector nor zero vector
             * @return true if weight is neither unit vector nor zero vector; false otherwise
             */
            bool hasNonUnitAndNonZeroWeight() const;

            /**
             * Indicates whether zero indices are subset of given zero indices
             * @param common_zero_inds Zero indices
             * @return true if zero indices are subset of given common_zero_inds
             */
            bool hasZeroIndsSuperSet(const std::bitset<kMaxInitialHrepSize>& common_zero_inds) const;

            /**
             * Get minimal infeasibility index
             * @return minimal infeasibility index
             */
            std::size_t getMinInfeasIndex() const;

            /**
             * Indicates whether given index is zero index
             * @param index Index to check
             * @return true if index is zero index; false otherwise
             */
            bool isZeroSlackIndex(std::size_t index) const;

            /**
             * Move weight
             * @return Rvalue reference to corresponding weight
             */
            WeightType&& moveWeight() {return std::move(weight_);};

            /**
             * Get weighted objective value
             * @return weighted objective value
             */
            ValueType getWov() const {return wov_;};

        private:
            constexpr static double kNormalizingThreshold = 1e+5; ///< Normalizing threshold

            /**
             * Indicates whether corresponding weight should be normalized
             * @param threshold Normalizing threshold
             * @return true if weight should be normalized; false otherwise
             */
            bool shouldNormalize(double threshold) const;

            /**
             * Normalizing weight function
             * @param normalizing_val Value by which weights are normalized
             */
            void normalize(double normalizing_val);

            /**
             * Get slack value
             * @param index Corresponding index to get slack for
             * @return Slack value
             */
            ValueType getSlack(std::size_t index) const;

            /**
             * Set slack and minimal infeasibility index
             * @param scip SCIP pointer
             * @param h_rep h-representation
             */
            void setSlacksAndMinInfeasInd(SCIP* scip,
                                          const H_RepC& h_rep);

            WeightType weight_; ///< Corresponding weight of v-representation
            ValueType wov_; ///< Corresponding weighted objective value of v-representation
            std::pair<bool, std::size_t> min_infeas_ind_; ///< Minimal infeasibility index
            SlackMap inds_to_slacks_; ///< Maps index to slack
            std::bitset<kMaxInitialHrepSize> zero_slacks_; ///< Indices which have zero slack
        };

        using V_RepC = std::vector<std::shared_ptr<V_RepT>>; ///< Container for v-representations

        /**
         * @class DoubleDescriptionMethod
         * @brief Algorithm for transforming h-representation to v-representation
         * @details Details of the underlying algorithm can be found in "Double Description Method Revisited" by K. Fukuda and A. Prodon
         */
        class DoubleDescriptionMethod {
        public:
            using AdjPair = std::pair<std::reference_wrapper<const V_RepT>, std::reference_wrapper<const V_RepT>>; ///< Pair of adjacent v-representation elements
            using AdjPairContainer = std::vector<AdjPair>; ///< Container for adjacent pairs

            /**
             * Constructor
             * @param scip SCIP pointer
             * @param no_all_obj Number of objectives
             * @param bounded Bounded results
             * @param unbounded Unbounded results
             */
            explicit DoubleDescriptionMethod(SCIP *scip,
                                             std::size_t no_all_obj,
                                             const ResultContainer &bounded,
                                             const ResultContainer &unbounded);

            /**
             * Standard double description algorithm
             */
            void computeVRep();

            /**
             * Variation 1 of double description algorithm
             */
            void computeVRep_Var1();

            /**
             * Print function
             * @param os Output stream to write to
             * @param withIncidentFacets Indicates whether corresponding facets should be printed
             */
            void printVRep(std::ostream &os = std::cout,
                           bool withIncidentFacets = false) const;

            /**
             * Get size of v-representation
             * @return Number of elements in v-representation
             */
            std::size_t size() const { return v_rep_.size(); };

            /**
             * Copy entire h-representation
             * @return Copy of container storing entire h-representation
             */
            H_RepC getHRep() {return h_rep_;};

            /**
             * Copy entire v-representation
             * @return Copy of container storing entire v-representation
             */
            V_RepC getVRep() {return v_rep_;};

            /**
             * Make entire h-representation movable
             * @return Rvalue reference to container storing entire h-representation
             */
            H_RepC&& moveHRep() {return std::move(h_rep_);};

            /**
             * Make entire v-representation movable
             * @return Rvalue reference to container storing entire v-representation
             */
            V_RepC&& moveVRep() {return std::move(v_rep_);};

        private:
            using SlackContainer = std::vector<std::vector<std::size_t>>; ///< Container for slacks

            /**
             * Variable orders
             */
            enum class VarOrder {
                keep_var_order, ///< Keep variable order
                change_var_order ///< Exchange variable order
            };

            /**
             * Computes initial v-representation of following h-representation:
             * 1) bounded \\cdot (w_1,...,w_k) -a >= 0
             * 2) w_1 >= 0
             * ...
             * k+1) w_k >= 0
             * where (w_1,...,w_k) is a weight vector
             * @return Initial v-representation
             */
            V_RepC computeInitialVRep() const;

            /**
             * Computes the intersection of zero indices for two given v-representation elements
             * @param v First element
             * @param w Second element
             * @return bitset with ones indicating common zero slacks
             */
            std::bitset<V_RepT::kMaxInitialHrepSize> getCommonZeroSlackIndices(const V_RepT &v,
                                                                               const V_RepT &w) const;

            /**
             * Check for minimal infeasibility condition
             * @param r1 First element
             * @param r2 Second element
             * @return Tuple
             */
            std::tuple<bool, VarOrder, std::size_t> minInfeasCondition(const V_RepT& r1,
                                                                       const V_RepT& r2) const;

            /**
             * Applies infeasibility condition
             * @param r1
             * @param r2
             * @param current_v_rep
             * @param index
             * @param with_adjacency_test
             */
            void applyInfeasCondition(const V_RepT& r1,
                                      const V_RepT& r2,
                                      const V_RepC& current_v_rep,
                                      std::size_t index,
                                      bool with_adjacency_test);


            /**
             * Detailed description of algorithm can be found on page 13 in "Double Description Method Revisited"
             * @param r1
             * @param r2
             * @param k
             * @param i
             * @param v_rep
             * @param with_adjacency_test
             */
            void conditionalStoreEdge(const V_RepT& r1,
                                      const V_RepT& r2,
                                      std::size_t k,
                                      std::size_t i,
                                      const V_RepC& v_rep,
                                      bool with_adjacency_test);


            /**
             * Indicates whether given elements are adjacent
             * @param ray1
             * @param ray2
             * @param v_rep
             * @return True if ray1 and ray2 are adjancent; false otherwise
             */
            bool rayPairIsAdjacent(const V_RepT& ray1,
                                   const V_RepT& ray2,
                                   const V_RepC& v_rep) const;

            /**
             * Indicates whether first element is multiple of second element
             * @param v First v-representation element
             * @param w Second v-representation element
             * @return true if v is multiple of w; false otherwise
             */
            bool isMultiple(const V_RepT& v, const V_RepT& w) const;

            /**
             * Indicates whether weight of first element is a multiple (with value v_multiple) of weight of second element
             * @param scip SCIP pointer
             * @param v_multiple Value
             * @param v First element
             * @param w Second element
             * @return true if weight of first element is a v_multiple of weight of second element; false otherwise
             */
            bool weightIsMultiple(SCIP* scip,
                                  double v_multiple,
                                  const V_RepT& v,
                                  const V_RepT& w) const;

            /**
             * Compute adjacent pairs of elements in v-representation
             * @param plus Elements in plus partition
             * @param minus Elements in minus partition
             * @param current_v_rep Current v-representation
             * @return Adjacent elements
             */
            AdjPairContainer computeAdjacentPairs(const V_RepC& plus,
                                                  const V_RepC& minus,
                                                  const V_RepC& current_v_rep) const;

            /**
             * Incorporate next hyperplane into v-representation; standard algorithm
             * @param cur_v_rep Current v-representation
             * @return Extended v-representation
             */
            V_RepC extendVRep(V_RepC&& cur_v_rep);

            /**
             * Incorporate next hyperplane into v-representation; Variation1 algorithm
             * @param current_v_rep Current v-representation
             * @return Extended v-representation
             */
            V_RepC extendVRep_Var1(V_RepC&& current_v_rep);


            SCIP *scip_; ///< SCIP pointer
            std::size_t outcome_dimension_; ///< Dimension of outcome in objective space
            std::size_t current_hrep_index_; ///< hyperplane currently considered
            H_RepC h_rep_; ///< h-representation
            V_RepC v_rep_; ///< v-representation
            std::vector< AdjPairContainer > adj_pairs_; ///< Adjacent pairs
        };
    }
}
#endif //POLYSCIP_SRC_DOUBLE_DESCRIPTION_H_INCLUDED