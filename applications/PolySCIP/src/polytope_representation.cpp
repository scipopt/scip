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

#include "polytope_representation.h"

#include <algorithm>
#include <bitset>
#include <cmath> //std::fabs
#include <functional>
#include <iterator>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "PolySCIPConfig.h"

using std::make_shared;
using std::pair;
using std::size_t;
using std::vector;

namespace polyscip {

    namespace polytoperepresentation {

        V_RepT::V_RepT(SCIP* scip, WeightType&& weight, ValueType&& wov, const H_RepContainer& h_rep)
                : weight_(std::move(weight)),
                  wov_(std::move(wov)),
                  min_infeas_ind_(false, 0)
        {
            if (shouldNormalize(kNormalizingThreshold))
                normalize(kNormalizingThreshold);

            setSlacksAndMinInfeasInd(scip, h_rep);
        }

        V_RepT::V_RepT(SCIP* scip,
                       const V_RepT& plus,
                       const V_RepT& minus,
                       size_t index_of_ineq,
                       const H_RepContainer& h_rep)
                : min_infeas_ind_(false, 0)
        {
            auto m_coeff = plus.getSlack(index_of_ineq);
            assert (SCIPisPositive(scip, m_coeff));
            auto p_coeff = minus.getSlack(index_of_ineq);
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

        bool V_RepT::operator==(const V_RepT& rhs) const {
            return std::tie(weight_, wov_) == std::tie(rhs.weight_, rhs.wov_);
        }

        bool V_RepT::operator!=(const V_RepT& rhs) const {
            return !(operator==(rhs));
        }

        ValueType V_RepT::getSlack(std::size_t index) const {
            assert (inds_to_slacks_.count(index) != 0);
            return inds_to_slacks_.at(index);
        }

        size_t V_RepT::getMinInfeasIndex() const {
            assert (min_infeas_ind_.first);
            return min_infeas_ind_.second;
        }

        bool V_RepT::isZeroSlackIndex(size_t index) const {
            return zero_slacks_.test(index);
        }

        void V_RepT::setSlacksAndMinInfeasInd(SCIP* scip, const H_RepContainer& h_rep) {
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

        //todo test ratio between wov_ und weight
        bool V_RepT::shouldNormalize(double threshold) const {
            if (std::fabs(wov_) > threshold)
                return true;
            else if (std::any_of(weight_.cbegin(), weight_.cend(), [threshold](ValueType w){return w > threshold;}))
                return true;
            return false;
        }

        void V_RepT::normalize(double normalizing_val) {
            std::transform(begin(weight_), end(weight_), begin(weight_),
                           [normalizing_val](const ValueType &val) { return val / normalizing_val; });
            wov_ /= normalizing_val;
        }

        void V_RepT::print(std::ostream& os, bool withIncidentFacets, const H_RepContainer& h_rep) const {
            global::print(weight_, "Weight = [", "]", os);
            os << " Coeff = " << wov_ << "\n";
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


        bool V_RepT::hasZeroIndsSuperSet(const std::bitset<kMaxInitialHrepSize>& common_zero_inds) const {
            return (common_zero_inds & zero_slacks_).count() >= common_zero_inds.count();
        }

        DoubleDescriptionMethod::DoubleDescriptionMethod(SCIP *scip, const ResultContainer& bounded_results,
                                                         const ResultContainer& unbounded_results)
                : scip_(scip),
                  adj_pairs_(bounded_results.size()+unbounded_results.size()-1, AdjPairContainer{})
        {
            /* build up h-representation for which v-representation is to be found */
            assert (!bounded_results.empty());
            outcome_dimension_ = bounded_results.front().second.size();
            adj_pairs_ind_offset_ = outcome_dimension_+1;
            for (size_t i=0; i<outcome_dimension_; ++i) {
                auto unit_vec = OutcomeType(outcome_dimension_, 0.);
                unit_vec[i] = 1.;
                h_rep_.push_back({std::move(unit_vec), 0.});  // add constraint: e_i 0 >= 0 with e_i being i-th unit vector
            }
            for (const auto &bd : bounded_results) {
                h_rep_.emplace_back(bd.second, 1.); // add constraint; bd_outcome -1 >= 0
            }
            current_hrep_index_ = outcome_dimension_;
            for (const auto &unbd : unbounded_results) {
                h_rep_.emplace_back(unbd.second, 0.);    // add constraint; unbd_outcome 0 >= 0
            }
        }

        std::pair<bool, size_t> DoubleDescriptionMethod::minInfeasCondition(const V_RepT& r1, const V_RepT& r2) const {
            auto minInfeasInd1 = r1.getMinInfeasIndex();
            auto minInfeasInd2 = r2.getMinInfeasIndex();
            if (minInfeasInd1 < minInfeasInd2 && !r2.isZeroSlackIndex(minInfeasInd1)) {
                assert (current_hrep_index_ < minInfeasInd1 && minInfeasInd1 < h_rep_.size());
                return {true, minInfeasInd1};
            }
            else {
                return {false, 0};
            }
        }

        /*void DoubleDescriptionMethod::conditionalStoreEdge(const V_RepT& r1,
                                                           const V_RepT& r2,
                                                           size_t k,
                                                           size_t i,
                                                           const V_RepContainer& v_rep) {
            assert (i <= k);
            for (auto index=i+1; index<k; ++index) {
                if (r1.isZeroSlackIndex(index) && r2.isZeroSlackIndex(index))
                    return;
            }
            if (rayPairIsAdjacent(r1, r2, v_rep)) {
                adj_pairs_.at(k-adj_pairs_ind_offset_).push_back({&r2, &r1});
            }
        }*/

        void DoubleDescriptionMethod::printVRep(std::ostream &os, bool withIncidentFacets) const {
            for (const auto& v : v_rep_)
                v->print(os, withIncidentFacets, h_rep_);
        }

        void DoubleDescriptionMethod::computeVRep() {
            auto current_v_rep = computeInitialVRep();

            ++current_hrep_index_;
            while (current_hrep_index_ < h_rep_.size()) {
                current_v_rep = extendVRep(std::move(current_v_rep));
                ++current_hrep_index_;
            }

            assert (v_rep_.empty());
            auto scip_ptr = scip_;
            std::remove_copy_if(
                    begin(current_v_rep),   // copy all elements from current_v_rep with non-zero weights to v_rep_
                    end(current_v_rep),
                    std::back_inserter(v_rep_),
                    [scip_ptr](const std::shared_ptr<V_RepT>& v) {
                        return SCIPisZero(scip_ptr, std::accumulate(begin(v->weight_),
                                                                    end(v->weight_),
                                                                    0.0));
                    });
            assert (v_rep_.size() + 1 == current_v_rep.size()); // assert that only (0,0,...,0,-1) was removed
        }

        /*void DoubleDescriptionMethod::computeVRep_Variation1() {
            auto current_v_rep = computeInitialVRep();

            for (auto r1_it=begin(current_v_rep); r1_it!=std::prev(end(current_v_rep)); ++r1_it) {
                for (auto r2_it=std::next(r1_it); r2_it!=end(current_v_rep); ++r2_it) {
                    auto minInfeasPair = minInfeasCondition(*r1_it, *r2_it);
                    if (minInfeasPair.first) {
                        conditionalStoreEdge(*r1_it, *r2_it, minInfeasPair.second, current_hrep_index_, current_v_rep);
                    }
                    else {
                        minInfeasPair = minInfeasCondition(*r2_it, *r1_it);
                        if (minInfeasPair.first)
                            conditionalStoreEdge(*r2_it, *r1_it, minInfeasPair.second, current_hrep_index_, current_v_rep);
                    }
                }
            }

            ++current_hrep_index_;
            while (current_hrep_index_ < h_rep_.size()) {
                std::cout << "current h_rep index = " << current_hrep_index_ << "\n";
                current_v_rep = extendVRep_Variation1(std::move(current_v_rep));
                ++current_hrep_index_;
            }

            assert (v_rep_.empty());
            auto scip_ptr = scip_;
            std::remove_copy_if(
                    begin(current_v_rep),   // copy all elements from current_v_rep with non-zero weights to v_rep_
                    end(current_v_rep),
                    std::back_inserter(v_rep_),
                    [scip_ptr](const V_RepT &v) {
                        return SCIPisZero(scip_ptr, std::accumulate(begin(v.weight_),
                                                                    end(v.weight_),
                                                                    0.0));
                    });
            assert (v_rep_.size() + 1 == current_v_rep.size()); // assert that only (0,0,...,0,-1) was removed
        }


        V_RepContainer DoubleDescriptionMethod::extendVRep_Variation1(V_RepContainer&& cur_v_rep) {
            auto extended_v_rep = V_RepContainer {};

            const H_RepT& constraint = h_rep_[current_hrep_index_];

            for (const auto& v_rep : cur_v_rep) { // partition current v-representation
                auto result = std::inner_product(v_rep.weight_.cbegin(), v_rep.weight_.cend(),
                                                 constraint.first.cbegin(), -(v_rep.wov_ * constraint.second));

                if (SCIPisZero(scip_, result)) {
                    extended_v_rep.push_back(v_rep); // element will also be in extended v-representation
                }
                else if (SCIPisPositive(scip_, result)) {
                    extended_v_rep.push_back(v_rep);
                }
            }

            for (auto p : adj_pairs_[current_hrep_index_-adj_pairs_ind_offset_]) {
                extended_v_rep.emplace_back(scip_, *(p.first), *(p.second), current_hrep_index_, h_rep_);
            }
            for (auto r1=extended_v_rep.cbegin(); r1!=std::prev(extended_v_rep.cend()); ++r1) {
                for (auto r2=std::next(r1); r2!=extended_v_rep.cend(); ++r2) {
                    if (r1->isZeroSlackIndex(current_hrep_index_) && r2->isZeroSlackIndex(current_hrep_index_)) {
                        auto minInfeasPair = minInfeasCondition(*r1, *r2);
                        if (minInfeasPair.first) {
                            conditionalStoreEdge(*r1, *r2, minInfeasPair.second, current_hrep_index_, extended_v_rep);
                        }
                        else {
                            auto minInfeasPair = minInfeasCondition(*r2, *r1);
                            if (minInfeasPair.first)
                                conditionalStoreEdge(*r2, *r1, minInfeasPair.second, current_hrep_index_, extended_v_rep);
                        }
                    }
                }
            }

            return extended_v_rep;
        }*/

        V_RepContainer DoubleDescriptionMethod::extendVRep(V_RepContainer&& cur_v_rep) {
            auto extended_v_rep = V_RepContainer {};
            auto plus_inds = vector<std::size_t> {};
            auto minus_inds = vector<std::size_t> {};
            const H_RepT& constraint = h_rep_[current_hrep_index_];

            for (size_t i = 0; i < cur_v_rep.size(); ++i) { // partition current v-representation
                auto result = std::inner_product(cur_v_rep[i]->weight_.cbegin(), cur_v_rep[i]->weight_.cend(),
                                                 constraint.first.cbegin(), -(cur_v_rep[i]->wov_ * constraint.second));
                if (SCIPisNegative(scip_, result)) {
                    minus_inds.push_back(i);
                }
                else if (SCIPisZero(scip_, result)) {
                    extended_v_rep.push_back(cur_v_rep[i]); // element will also be in extended v-representation
                }
                else {
                    assert(SCIPisPositive(scip_, result));
                    plus_inds.push_back(i);
                }
            }
            auto adj_pairs = computeAdjacentPairs(plus_inds, minus_inds, cur_v_rep);
            for (const auto &elem : adj_pairs)
                extended_v_rep.push_back(make_shared<V_RepT>(scip_, *cur_v_rep[elem.first], *cur_v_rep[elem.second], current_hrep_index_, h_rep_));

            for (const auto& i : plus_inds)
                extended_v_rep.push_back(cur_v_rep[i]);

            return extended_v_rep;
        }

        vector<pair<size_t, size_t>> DoubleDescriptionMethod::computeAdjacentPairs(const vector<size_t> &plus_inds,
                                                                                   const vector<size_t> &minus_inds,
                                                                                   const V_RepContainer& cur_v_rep) const {
            auto adj_pairs = vector<pair<size_t, size_t>> {};
            for (const auto& plus : plus_inds) {
                for (const auto& minus : minus_inds) {
                    assert (plus != minus);
                    if (rayPairIsAdjacent(plus, minus, cur_v_rep))
                        adj_pairs.push_back({plus, minus});
                }
            }
            return adj_pairs;
        }

        std::bitset<V_RepT::kMaxInitialHrepSize> DoubleDescriptionMethod::getCommonZeroSlackIndices(const V_RepT &v,
                                                                                                    const V_RepT &w) const {
            return v.zero_slacks_ & w.zero_slacks_;
        }

        bool DoubleDescriptionMethod::rayPairIsAdjacent(size_t index1,
                                                        size_t index2,
                                                        const V_RepContainer& cur_v_rep) const {

            auto common_zero_inds = getCommonZeroSlackIndices(*cur_v_rep[index1], *cur_v_rep[index2]);

            for (size_t i = 0; i < cur_v_rep.size(); ++i) {
                if (i == index1 || i == index2)
                    continue;
                else if (cur_v_rep[i]->hasZeroIndsSuperSet(common_zero_inds)) {
                    /* check whether current_rep[i] is multiple of current_rep[index1] or current_rep[index2] */
                    if (!isMultiple(*cur_v_rep[i], *cur_v_rep[index1]) &&
                        !isMultiple(*cur_v_rep[i], *cur_v_rep[index2]))
                        return false;
                }
            }
            return true;
        }


        /*bool DoubleDescriptionMethod::rayPairIsAdjacent(const V_RepT& ray1,
                                                        const V_RepT& ray2,
                                                        const V_RepContainer& v_rep) const {
            auto common_zero_inds = getCommonZeroSlackIndices(ray1, ray2);

            for (const auto& ray : v_rep) {
                if (ray == ray1 || ray == ray2) {
                    continue;
                }
                else if (ray.hasZeroIndsSuperSet(common_zero_inds)) {
                    *//* check whether current_rep[i] is multiple of current_rep[index1] or current_rep[index2] *//*
                    if (!isMultiple(ray, ray1) && !isMultiple(ray, ray2))
                        return false;
                }
            }
            return true;
        }*/


        //todo Check function thoroughly
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

        bool DoubleDescriptionMethod::weightIsMultiple(SCIP* scip, double v_multiple, const V_RepT& v, const V_RepT& w) const {
            auto mismatch_weight = std::mismatch(begin(v.weight_),
                                                 end(v.weight_),
                                                 begin(w.weight_),
                                                 [scip, v_multiple](ValueType v_val, ValueType w_val)
                                                 {return SCIPisEQ(scip, v_multiple * v_val, w_val);});
            return mismatch_weight.first == end(v.weight_);
        }


        V_RepContainer DoubleDescriptionMethod::computeInitialVRep() const {
            // create initial v_rep
            auto init_v_rep = V_RepContainer {};
            init_v_rep.push_back(make_shared<V_RepT>(scip_, WeightType(outcome_dimension_, 0.), -1., h_rep_)); // add v_rep: 0 0 ... 0 -1
            for (size_t i=0; i<outcome_dimension_; ++i) {
                auto unit_vec = WeightType(outcome_dimension_, 0.);
                unit_vec[i] = 1.;
                auto wov = h_rep_[current_hrep_index_].first.at(i);
                init_v_rep.push_back(make_shared<V_RepT>(scip_, std::move(unit_vec), std::move(wov), h_rep_)); // add v_rep: e_i, vrep_wov with e_i being i-th unit vector
            }
            return init_v_rep;
        }

    }

}
