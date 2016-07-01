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

#ifndef POLYSCIP_SRC_POLYTOPE_REPRESENTATION_H_INCLUDED
#define POLYSCIP_SRC_POLYTOPE_REPRESENTATION_H_INCLUDED

#include <algorithm>
#include <cstddef>
#include <functional> // std::function
#include <iostream>
#include <unordered_map>
#include <ostream>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "weight_space_facet.h"

namespace polyscip {

    namespace polytoperepresentation {

        using H_RepT = std::pair<OutcomeType, ValueType>;
        using H_RepContainer = std::vector<H_RepT>;

        class V_RepT {
        public:
            friend class DoubleDescriptionMethod;
            using SlackMap = std::unordered_map<std::size_t, ValueType>;

            explicit V_RepT(WeightType weight, ValueType wov);

            explicit V_RepT(SCIP* scip, WeightType weight, ValueType wov, const H_RepContainer& current_h_rep);

            explicit V_RepT(SCIP* scip,
                            const V_RepT& v,
                            const V_RepT& w,
                            std::size_t index_of_constraint_in_hrep,
                            const H_RepContainer& current_h_rep);


            ValueType getSlack(std::size_t index) const {return inds_to_slacks_.at(index);}
            void addZeroSlackIndex(std::size_t index);
            void addSlack(std::size_t index, double val);
            void print(std::ostream& os, bool withIncidentFacets, const H_RepContainer& h_rep) const;

            //todo incorporate way for facets when WSP is not full-dimensional

            /** Checks whether given parameter is a subset of member zero_slack_indices_
             * @Return true if member zero_slack_indices_ is superset of given parameter; false otherwise
             */
            bool hasZeroSlackSuperSet(const std::vector<std::size_t>& indices) const;

            WeightType&& moveWeight() {return std::move(weight_);};
            WeightType getWeight() const {return weight_;};
            ValueType getWov() const {return wov_;};

        private:
            static constexpr double kNormalizingThreshold = 1e+4;
            void normalize();

            WeightType weight_;
            ValueType wov_;
            SlackMap inds_to_slacks_;
            std::vector<std::size_t> zero_slacks_; // indices in h_rep_ which this object fulfills with equality
        };

        using V_RepContainer = std::vector<V_RepT>;


        class DoubleDescriptionMethod {
        public:

            DoubleDescriptionMethod(SCIP *scip, const ResultContainer &bounded, const ResultContainer &unbounded);

            void computeVRep();

            void printVRep(std::ostream &os = std::cout, bool withIncidentFacets = false) const;

            std::size_t size() const { return v_rep_.size(); };

            V_RepContainer getVRep() const { return v_rep_; };

            H_RepContainer getHRep() const { return h_rep_; };

            //V_RepContainer &&moveVRep() { return std::move(v_rep_); };

            //H_RepContainer &&moveHRep() { return std::move(h_rep_); };

        private:
            using SlackContainer = std::vector<std::vector<std::size_t>>;

            /** Computes initial v-representation for the following h-representation:
             * 1) bounded \cdot (w_1,...,w_k) - a >= 0
             * 2) w_1 >= 0
             * ...
             * k+1) w_k >= 0
             */

            void computeInitialRep(const OutcomeType &bounded);

            std::vector<std::size_t> getCommonZeroSlacks(const V_RepT& v, const V_RepT& w) const;

            bool rayPairIsAdjacent(std::size_t plus_index,
                                       std::size_t minus_index,
                                       const std::vector<V_RepT> &current_v_rep) const;

            /** Check whether first parameter is multiple of second or third parameter
             */
            bool isMultiple(const V_RepT& v, const V_RepT& w) const;
            bool weightIsMultiple(SCIP* scip, double v_multiple, const V_RepT& v, const V_RepT& w) const;

            std::vector<std::pair<std::size_t, std::size_t>> computeAdjacentPairs(const std::vector<std::size_t>& plus_inds,
                                                                                      const std::vector<std::size_t>& minus_inds,
                                                                                      const std::vector<V_RepT>& current_rep) const;

            std::vector<V_RepT> extendVRep(std::vector<V_RepT> current_rep,
                                               const H_RepT &constraint,
                                               std::size_t index_of_constraint_in_hrep);

            SCIP *scip_;
            std::vector<OutcomeType> bounded_;
            std::vector<OutcomeType> unbounded_;
            H_RepContainer h_rep_;
            V_RepContainer initial_v_rep_;
            V_RepContainer v_rep_;
        };
    }
}
#endif //POLYSCIP_SRC_POLYTOPE_REPRESENTATION_H_INCLUDED