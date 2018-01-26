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
 * @file double_description_method.cpp
 * @brief Implements the double description method for transforming a polyhedron given
 * via its h-representation into its v-representation
 * @author Sebastian Schenker
 *
 */

#include "double_description_method.h"

#include <algorithm>
#include <numeric>
#include <bitset>
#include <cmath> //std::fabs
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "PolySCIPConfig.h"

using std::make_shared;
using std::size_t;
using std::vector;

namespace polyscip {

    namespace doubledescription {

        /**
         * Constructor
         * @param scip SCIP pointer
         * @param weight Corresponding weight for the constructed vertex
         * @param wov Corresponding weighted objective value for the constructed vertex
         * @param h_rep Current h-representation
         */
        V_RepT::V_RepT(SCIP* scip, WeightType&& weight, ValueType&& wov, const H_RepC& h_rep)
                : weight_(weight),
                  wov_(wov),
                  min_infeas_ind_(false, 0)
        {
            if (shouldNormalize(kNormalizingThreshold))
                normalize(kNormalizingThreshold);

            setSlacksAndMinInfeasInd(scip, h_rep);
        }

        /**
         * Constructor
         * @param scip SCIP pointer
         * @param plus
         * @param minus
         * @param index_of_ineq
         * @param h_rep
         */
        V_RepT::V_RepT(SCIP* scip,
                       const V_RepT& plus,
                       const V_RepT& minus,
                       std::size_t index_of_ineq,
                       const H_RepC& h_rep)
                : min_infeas_ind_(false, 0)
        {
            auto m_coeff = plus.getSlack(index_of_ineq);
            auto p_coeff = minus.getSlack(index_of_ineq);
            assert (SCIPisPositive(scip, m_coeff));
            assert (SCIPisNegative(scip, p_coeff));

            std::transform(minus.weight_.cbegin(), minus.weight_.cend(),
                           plus.weight_.cbegin(), std::back_inserter(weight_),
                           [m_coeff, p_coeff](ValueType m_val, ValueType p_val) {
                               return m_coeff * m_val - p_coeff * p_val;
                           });
            wov_ = m_coeff*minus.wov_ - p_coeff*plus.wov_; // return m_coeff * ray_minus - p_coeff * ray_plus

            if (shouldNormalize(kNormalizingThreshold))
                normalize(kNormalizingThreshold);

            setSlacksAndMinInfeasInd(scip, h_rep);
        }

        /**
         * Equal operator
         * @param rhs v-representation object to compare with
         * @return true if objects are considered equal; otherwise false
         */
        bool V_RepT::operator==(const V_RepT& rhs) const {
            return std::tie(weight_, wov_) == std::tie(rhs.weight_, rhs.wov_);
        }

        /**
         * Not equal operator
         * @param rhs v-representation object to compare with
         * @return true if objects are considered not equal; otherwise false
         */
        bool V_RepT::operator!=(const V_RepT& rhs) const {
            return !(operator==(rhs));
        }

        /**
         * Get slack value
         * @param index Corresponding index to get slack for
         * @return Slack value
         */
        ValueType V_RepT::getSlack(std::size_t index) const {
            assert (inds_to_slacks_.count(index) != 0);
            return inds_to_slacks_.at(index);
        }

        /**
         * Get minimal infeasibility index
         * @return minimal infeasibility index
         */
        size_t V_RepT::getMinInfeasIndex() const {
            assert (min_infeas_ind_.first);
            return min_infeas_ind_.second;
        }

        /**
         * Indicates whether given index is zero index
         * @param index Index to check
         * @return true if index is zero index; false otherwise
         */
        bool V_RepT::isZeroSlackIndex(size_t index) const {
            return zero_slacks_.test(index);
        }

        /**
         * Indicates whether corresponding weight is neither unit vector nor zero vector
         * @return true if weight is neither unit vector nor zero vector; false otherwise
         */
        bool V_RepT::hasNonUnitAndNonZeroWeight() const {
            return (std::count_if(weight_.cbegin(), weight_.cend(), [](const ValueType& val){return val != 0.;}) > 1);
        }


        /**
         * Set slack and minimal infeasibility index
         * @param scip SCIP pointer
         * @param h_rep h-representation
         */
        void V_RepT::setSlacksAndMinInfeasInd(SCIP* scip, const H_RepC& h_rep) {
            for (size_t i=0; i<h_rep.size(); ++i) {
                auto result = std::inner_product(weight_.cbegin(), weight_.cend(),
                                                 h_rep[i].first.cbegin(), -(h_rep[i].second * wov_));

                inds_to_slacks_.emplace(i, result);
                if (SCIPisZero(scip, result)) {
                    inds_to_slacks_[i] = 0.;
                    zero_slacks_.set(i, 1);
                }
                else if (SCIPisNegative(scip, result) && !min_infeas_ind_.first) {
                    min_infeas_ind_.first = true;
                    min_infeas_ind_.second = i;
                }
            }
            if (!min_infeas_ind_.first) {
                min_infeas_ind_.first = true;
                min_infeas_ind_.second = h_rep.size();
            }
        }

        /**
         * Indicates whether corresponding weight should be normalized
         * @param threshold Normalizing threshold
         * @return true if weight should be normalized; false otherwise
         */
        bool V_RepT::shouldNormalize(double threshold) const {
            if (std::fabs(wov_) > threshold)
                return true;
            else if (std::any_of(weight_.cbegin(), weight_.cend(), [threshold](ValueType w){return w > threshold;}))
                return true;
            return false;
        }


        /**
         * Normalizing weight function
         * @param normalizing_val Value by which weights are normalized
         */
        void V_RepT::normalize(double normalizing_val) {
            std::transform(begin(weight_), end(weight_), begin(weight_),
                           [normalizing_val](const ValueType &val) { return val / normalizing_val; });
            wov_ /= normalizing_val;
        }

        /**
         * Print function
         * @param os Output stream to print to
         * @param withIncidentFacets Boolean indicating whether incident facets should also be printed
         * @param h_rep
         */
        void V_RepT::print(std::ostream& os, bool withIncidentFacets, const H_RepC& h_rep) const {
            global::print(weight_, "Weight = [", "]", os);
            os << " Coeff = " << wov_ << "\n";
            os << "MinInfeasIndex = " << min_infeas_ind_.second << "\n";
            if (withIncidentFacets) {
                os << "Facets: \n";
                for (size_t i=0; i<kMaxInitialHrepSize; ++i) {
                    if (zero_slacks_[i]) {
                        global::print(h_rep[i].first, "", "", os);
                        os << " " << h_rep[i].second << "\n";
                    }
                }
                os << "\n";
            }
        }


        /**
         * Indicates whether zero indices are subset of given zero indices
         * @param common_zero_inds Zero indices
         * @return true if zero indices are subset of given common_zero_inds
         */
        bool V_RepT::hasZeroIndsSuperSet(const std::bitset<kMaxInitialHrepSize>& common_zero_inds) const {
            return (common_zero_inds & zero_slacks_).count() >= common_zero_inds.count();
        }


        /**
         * Constructor
         * @param scip SCIP pointer
         * @param no_all_obj Number of objectives
         * @param bounded Bounded results
         * @param unbounded Unbounded results
         */
        DoubleDescriptionMethod::DoubleDescriptionMethod(SCIP *scip,
                                                         size_t no_all_objs,
                                                         const ResultContainer& bounded_results,
                                                         const ResultContainer& unbounded_results)
                : scip_(scip),
                  outcome_dimension_(no_all_objs),
                  current_hrep_index_(no_all_objs),
                  adj_pairs_(no_all_objs+bounded_results.size()+unbounded_results.size(), AdjPairContainer{})
        {
            /* build up h-representation for which v-representation is to be found */
            assert (!bounded_results.empty());

            for (size_t i=0; i<outcome_dimension_; ++i) {
                auto unit_vec = OutcomeType(outcome_dimension_, 0.);
                unit_vec[i] = 1.;
                h_rep_.emplace_back(std::move(unit_vec), 0.);  // add constraint: e_i 0 >= 0 with e_i being i-th unit vector
            }

            for (const auto &bd : bounded_results) {
                h_rep_.emplace_back(bd.second, 1.); // add constraint; bd_outcome -1 >= 0
            }
            for (const auto &unbd : unbounded_results) {
                h_rep_.emplace_back(unbd.second, 0.);    // add constraint; unbd_outcome 0 >= 0
            }
        }

        /**
         * Check for minimal infeasibility condition
         * @param r1 First element
         * @param r2 Second element
         * @return Tuple
         */
        std::tuple<bool, DoubleDescriptionMethod::VarOrder, size_t> DoubleDescriptionMethod::minInfeasCondition(const V_RepT& r1, const V_RepT& r2) const {
            auto minInfeasInd1 = r1.getMinInfeasIndex();
            auto minInfeasInd2 = r2.getMinInfeasIndex();
            if (minInfeasInd1 < minInfeasInd2 && !r2.isZeroSlackIndex(minInfeasInd1)) {
                assert (current_hrep_index_ < minInfeasInd1 && minInfeasInd1 < h_rep_.size());
                return std::make_tuple(true, VarOrder::change_var_order, minInfeasInd1);
            }
            else if (minInfeasInd2 < minInfeasInd1 && !r1.isZeroSlackIndex(minInfeasInd2)) {
                assert (current_hrep_index_ < minInfeasInd2 && minInfeasInd2 < h_rep_.size());
                return std::make_tuple(true, VarOrder::keep_var_order, minInfeasInd2);
            }
            else {
                return std::make_tuple(false, VarOrder::keep_var_order, 0);
            }
        }

        /**
         * Detailed description of algorithm can be found on page 13 in "Double Description Method Revisited"
         * @param r1
         * @param r2
         * @param k
         * @param i
         * @param v_rep
         * @param with_adjacency_test
         */
        void DoubleDescriptionMethod::conditionalStoreEdge(const V_RepT& plus,
                                                           const V_RepT& minus,
                                                           size_t k,
                                                           size_t i,
                                                           const V_RepC& v_rep,
                                                           bool with_adjacency_test) {
            assert (i <= k);
            for (auto index=i+1; index<k; ++index) {
                if (plus.isZeroSlackIndex(index) && minus.isZeroSlackIndex(index))
                    return;
            }
            if (!with_adjacency_test || rayPairIsAdjacent(plus, minus, v_rep)) {
                adj_pairs_.at(k).push_back({std::cref(plus), std::cref(minus)});
            }
        }

        /**
         * Print function
         * @param os Output stream to write to
         * @param withIncidentFacets Indicates whether corresponding facets should be printed
         */
        void DoubleDescriptionMethod::printVRep(std::ostream &os, bool withIncidentFacets) const {
            for (const auto& v : v_rep_)
                v->print(os, withIncidentFacets, h_rep_);
        }

        /**
         * Standard double description algorithm
         */
        void DoubleDescriptionMethod::computeVRep() {
            auto current_v_rep = computeInitialVRep();
            ++current_hrep_index_;
            while (current_hrep_index_ < h_rep_.size()) {
                current_v_rep = extendVRep(std::move(current_v_rep));
                ++current_hrep_index_;
            }
            v_rep_ = current_v_rep;
        }

        /**
         * Applies infeasibility condition
         * @param r1
         * @param r2
         * @param current_v_rep
         * @param index
         * @param with_adjacency_test
         */
        void DoubleDescriptionMethod::applyInfeasCondition(const V_RepT& r1,
                                                           const V_RepT& r2,
                                                           const V_RepC& v_rep,
                                                           size_t index,
                                                           bool with_adjacency_test
        ) {
            auto infeasTuple = minInfeasCondition(r1, r2);
            if (std::get<0>(infeasTuple)) {
                if (std::get<1>(infeasTuple) == VarOrder::keep_var_order){
                    conditionalStoreEdge(r1, r2, std::get<2>(infeasTuple), index, v_rep, with_adjacency_test);
                }
                else {
                    conditionalStoreEdge(r2, r1, std::get<2>(infeasTuple), index, v_rep, with_adjacency_test);
                }
            }
        }

        /**
         * Variation 1 of double description algorithm
         */
        void DoubleDescriptionMethod::computeVRep_Var1() {
            auto current_v_rep = computeInitialVRep();

            for (auto r1_it=begin(current_v_rep); r1_it!=std::prev(end(current_v_rep)); ++r1_it) {
                for (auto r2_it=std::next(r1_it); r2_it!=end(current_v_rep); ++r2_it) {
                    const V_RepT& r1 = r1_it->operator*();
                    const V_RepT& r2 = r2_it->operator*();
                    applyInfeasCondition(r1, r2, current_v_rep, current_hrep_index_, true);
                }
            }

            ++current_hrep_index_;
            while (current_hrep_index_ < h_rep_.size()) {
                current_v_rep = extendVRep_Var1(std::move(current_v_rep));
                ++current_hrep_index_;
            }
            v_rep_ = current_v_rep;
        }


        /**
         * Incorporate next hyperplane into v-representation; Variation1 algorithm
         * @param current_v_rep Current v-representation
         * @return Extended v-representation
         */
        V_RepC DoubleDescriptionMethod::extendVRep_Var1(V_RepC&& current_v_rep) {

            auto extended_v_rep = V_RepC{};
            const H_RepT& h = h_rep_[current_hrep_index_];

            for (const auto& v : current_v_rep) { // partition current v-representation
                auto result = std::inner_product(v->weight_.cbegin(), v->weight_.cend(),
                                                 h.first.cbegin(), -(v->wov_ * h.second));
                if (SCIPisZero(scip_, result)) {
                    extended_v_rep.push_back(v); // element will also be in extended v-representation
                }
                else if (SCIPisPositive(scip_, result)) {
                    extended_v_rep.push_back(v);
                }
            }

            for (auto p : adj_pairs_[current_hrep_index_]) {
                auto new_ray = make_shared<V_RepT>(scip_, p.first.get(), p.second.get(), current_hrep_index_, h_rep_);
                extended_v_rep.push_back(new_ray);
                applyInfeasCondition(*new_ray, p.first.get(), extended_v_rep, current_hrep_index_, false);
            }

            for (auto r1_it=extended_v_rep.cbegin(); r1_it!=std::prev(extended_v_rep.cend()); ++r1_it) {
                for (auto r2_it=std::next(r1_it); r2_it!=extended_v_rep.cend(); ++r2_it) {
                    const V_RepT& r1 = r1_it->operator*();
                    const V_RepT& r2 = r2_it->operator*();
                    if (r1.isZeroSlackIndex(current_hrep_index_) && r2.isZeroSlackIndex(current_hrep_index_)) {
                        applyInfeasCondition(r1, r2, extended_v_rep, current_hrep_index_, true);
                    }
                }
            }

            return extended_v_rep;
        }

        /**
         * Incorporate next hyperplane into v-representation; standard algorithm
         * @param cur_v_rep Current v-representation
         * @return Extended v-representation
         */
        V_RepC DoubleDescriptionMethod::extendVRep(V_RepC&& cur_v_rep) {
            auto extended_v_rep = V_RepC{};
            auto plus = V_RepC{};
            auto minus = V_RepC{};
            const H_RepT& hrep = h_rep_[current_hrep_index_];

            for (const auto& v : cur_v_rep) {
                auto result = std::inner_product(v->weight_.cbegin(), v->weight_.cend(),
                                                 hrep.first.cbegin(), -(v->wov_ * hrep.second));
                if (SCIPisNegative(scip_, result)) {
                    minus.push_back(v);
                }
                else if (SCIPisZero(scip_, result)) {
                    extended_v_rep.push_back(std::move(v)); // element will also be in extended v-representation
                }
                else {
                    assert(SCIPisPositive(scip_, result));
                    plus.push_back(v);
                }
            }

            auto adj_pairs = computeAdjacentPairs(plus, minus, cur_v_rep);
            for (const auto& p : adj_pairs)
                extended_v_rep.push_back(make_shared<V_RepT>(scip_, p.first.get(), p.second.get(), current_hrep_index_, h_rep_));

            for (const auto& p : plus)
                extended_v_rep.push_back(std::move(p));

            return extended_v_rep;
        }


        /**
         * Compute adjacent pairs of elements in v-representation
         * @param plus Elements in plus partition
         * @param minus Elements in minus partition
         * @param current_v_rep Current v-representation
         * @return Adjacent elements
         */
        DoubleDescriptionMethod::AdjPairContainer DoubleDescriptionMethod::computeAdjacentPairs(const V_RepC& plus,
                                                                                                const V_RepC& minus,
                                                                                                const V_RepC& current_v_rep) const {
            auto adj_pairs = AdjPairContainer{};
            for (const auto& p_ptr : plus) {
                for (const auto& m_ptr : minus) {
                    assert (*p_ptr != *m_ptr);
                    if (rayPairIsAdjacent(*p_ptr, *m_ptr, current_v_rep))
                        adj_pairs.push_back({std::cref(*p_ptr), std::cref(*m_ptr)});
                }
            }
            return adj_pairs;
        }

        /**
         * Computes the intersection of zero indices for two given v-representation elements
         * @param v First element
         * @param w Second element
         * @return bitset with ones indicating common zero slacks
         */
        std::bitset<V_RepT::kMaxInitialHrepSize> DoubleDescriptionMethod::getCommonZeroSlackIndices(const V_RepT &v,
                                                                                                    const V_RepT &w) const {
            return v.zero_slacks_ & w.zero_slacks_;
        }


        /**
         * Indicates whether given elements are adjacent
         * @param ray1
         * @param ray2
         * @param v_rep
         * @return True if r1 and r2 are adjancent; false otherwise
         */
        bool DoubleDescriptionMethod::rayPairIsAdjacent(const V_RepT& ray1,
                                                        const V_RepT& ray2,
                                                        const V_RepC& v_rep) const {
            auto common_zero_inds = getCommonZeroSlackIndices(ray1, ray2);

            for (const auto& ray_sptr : v_rep) {
                if (ray1 == *ray_sptr || ray2 == *ray_sptr) {
                    continue;
                }
                else if (ray_sptr->hasZeroIndsSuperSet(common_zero_inds)) {
                    /* check whether current_rep[i] is multiple of current_rep[index1] or current_rep[index2] */
                    if (!isMultiple(*ray_sptr, ray1) && !isMultiple(*ray_sptr, ray2))
                        return false;
                }
            }
            return true;
        }


        /**
         * Indicates whether first element is multiple of second element
         * @param v First v-representation element
         * @param w Second v-representation element
         * @return true if v is multiple of w; false otherwise
         * @todo Rewrite in less complicated fashion
         */
        bool DoubleDescriptionMethod::isMultiple(const V_RepT& v, const V_RepT& w) const {
            assert (v.weight_.size() == w.weight_.size());
            auto scip_ptr = scip_; // needed for lambda functions
            if (SCIPisEQ(scip_, v.wov_, w.wov_)) { //v.wov = w.wov
                auto mismatch_pair = std::mismatch(begin(v.weight_),
                                                   end(v.weight_),
                                                   begin(w.weight_),
                                                   [scip_ptr](ValueType v_val, ValueType w_val)
                                                   {return SCIPisEQ(scip_ptr, v_val, w_val);});
                if (mismatch_pair.first == end(v.weight_)) {// v.weight = w.weight
                    return true; // multiple is 1
                }
                else { // v.weight != w.weight
                    if (SCIPisZero(scip_, v.wov_)) { // v.wov=0 && w.wov=0
                        auto v_val = *mismatch_pair.first;
                        auto w_val = *mismatch_pair.second;
                        if (SCIPisZero(scip_, v_val) || SCIPisZero(scip_, w_val)) { // implies v_val!=0 || w_val!=0
                            return false; // v cannot be multiple of w
                        }
                        else { // v.wov=w.wov=0 && v_val!=w_val && v_val!=0 && w_val!=0
                            auto multiple = w_val / v_val; // multiple * v_val = w_val
                            return weightIsMultiple(scip_ptr, multiple, v, w);
                        }
                    }
                    else { // v.wov = w.wov && v.wov!=0 && v.weight!=w.weight
                        return false;
                    }
                }
            }
            else { // v.wov != w.wov_
                if (SCIPisZero(scip_, v.wov_)) { // -> w.wov!=0 implying v == k*w only if k=0 && v=0 implying something is wrong since v.weight should != 0
                    assert (!std::all_of(begin(v.weight_), end(v.weight_), [scip_ptr](ValueType val){return SCIPisZero(scip_ptr,val);})); // assert (v.weight!=0)
                    return false;
                }
                else if (SCIPisZero(scip_, w.wov_)) { // -> v.wov!=0 implying w == k*v only if k=0 && w=0 implying something is wrong since w.weight should != 0
                    assert (!std::all_of(begin(w.weight_), end(w.weight_), [scip_ptr](ValueType val){return SCIPisZero(scip_ptr,val);})); // assert (w.weight!=0)
                    return false;
                }
                else { // v.wov!=w.wov && v.wov!=0 && w.wov!=0
                    auto multiple = w.wov_ / v.wov_; // multiple * v.wov = w.wov
                    return weightIsMultiple(scip_ptr, multiple, v, w);
                }
            }
        }

        /**
         * Indicates whether weight of first element is a multiple (with value v_multiple) of weight of second element
         * @param scip SCIP pointer
         * @param v_multiple Value
         * @param v First element
         * @param w Second element
         * @return true if weight of first element is a v_multiple of weight of second element; false otherwise
         */
        bool DoubleDescriptionMethod::weightIsMultiple(SCIP* scip, double v_multiple, const V_RepT& v, const V_RepT& w) const {
            auto mismatch_weight = std::mismatch(begin(v.weight_),
                                                 end(v.weight_),
                                                 begin(w.weight_),
                                                 [scip, v_multiple](ValueType v_val, ValueType w_val)
                                                 {return SCIPisEQ(scip, v_multiple * v_val, w_val);});
            return mismatch_weight.first == end(v.weight_);
        }


        /**
         * Computes initial v-representation of following h-representation:
         * 1) bounded \\cdot (w_1,...,w_k) -a >= 0
         * 2) w_1 >= 0
         * ...
         * k+1) w_k >= 0
         * where (w_1,...,w_k) is a weight vector
         * @return Initial v-representation
         */
        V_RepC DoubleDescriptionMethod::computeInitialVRep() const {
            // create initial v_rep
            auto init_v_rep = V_RepC{};
            init_v_rep.push_back(make_shared<V_RepT>(scip_, WeightType(outcome_dimension_, 0.), -1., h_rep_)); // add v_rep: 0 0 ... 0 -1
            for (size_t i=0; i<current_hrep_index_; ++i) {
                auto unit_vec = WeightType(outcome_dimension_, 0.);
                unit_vec[i] = 1.;
                auto wov = h_rep_[current_hrep_index_].first.at(i);
                init_v_rep.push_back(make_shared<V_RepT>(scip_, std::move(unit_vec), std::move(wov), h_rep_)); // add v_rep: e_i, vrep_wov with e_i being i-th unit vector
            }
            return init_v_rep;
        }

    }

}
