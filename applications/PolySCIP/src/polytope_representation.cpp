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
#include <cmath> //std::fabs
#include <iterator>
#include <set>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"

using std::pair;
using std::size_t;
using std::vector;

namespace polyscip {

    namespace polytoperepresentation {

        V_RepT::V_RepT(WeightType weight, ValueType wov)
                : weight_(weight),
                  wov_(wov)
        {}

        V_RepT::V_RepT(SCIP* scip, WeightType weight, ValueType wov, const H_RepContainer& current_h_rep)
                : weight_(weight),
                  wov_(wov)
        {
            for (size_t i=0; i<current_h_rep.size(); ++i) {
                auto result = std::inner_product(begin(weight_), end(weight_),
                                                 begin(current_h_rep[i].first), -(current_h_rep[i].second * wov_));
                inds_to_slacks_.emplace(i, result);
                if (SCIPisZero(scip, result)) {
                    inds_to_slacks_[i] = 0.;
                    addZeroSlackIndex(i);
                }
            }
        }

        V_RepT::V_RepT(SCIP* scip,
                       const V_RepT& plus,
                       const V_RepT& minus,
                       size_t index_of_ineq,
                       const H_RepContainer& current_h_rep) {
            auto m_coeff = plus.getSlack(index_of_ineq);
            assert (SCIPisPositive(scip, m_coeff));
            auto p_coeff = minus.getSlack(index_of_ineq);
            assert (SCIPisNegative(scip, p_coeff));

            std::transform(begin(minus.weight_), end(minus.weight_),
                           begin(plus.weight_), std::back_inserter(weight_),
                           [m_coeff, p_coeff](ValueType m_val, ValueType p_val) {
                               return m_coeff * m_val - p_coeff * p_val;
                           });
            wov_ = m_coeff*minus.wov_ - p_coeff*plus.wov_; // return m_coeff * ray_minus - p_coeff * ray_plus

            auto normalizing_threshold = kNormalizingThreshold;
            if (std::any_of(weight_.cbegin(), weight_.cend(), [normalizing_threshold](ValueType w){return w > normalizing_threshold;}))
                normalize();

            for (size_t i=0; i<current_h_rep.size(); ++i) {
                auto res = m_coeff * minus.getSlack(i) - p_coeff * plus.getSlack(i);
                inds_to_slacks_.emplace(i, res);
                if (SCIPisZero(scip, res)) {
                    inds_to_slacks_[i] = 0.;
                    addZeroSlackIndex(i);
                }
            }
        }

        void V_RepT::normalize() {
            auto normalizing_val = kNormalizingThreshold;
            std::transform(begin(weight_), end(weight_), begin(weight_),
                           [normalizing_val](const ValueType &val) { return val / normalizing_val; });
            wov_ /= normalizing_val;
        }

        void V_RepT::addSlack(std::size_t index, double val) {
            auto res = inds_to_slacks_.emplace(index, val);
            assert (res.second);
        }

        void V_RepT::addZeroSlackIndex(size_t index) {
            zero_slacks_.push_back(index);
            assert (std::is_sorted(begin(zero_slacks_), end(zero_slacks_)));
        }

        void V_RepT::print(std::ostream& os, bool withIncidentFacets, const H_RepContainer& h_rep) const {
            global::print(weight_, "Weight = [", "]", os);
            os << " Coeff = " << wov_ << "\n";
            if (withIncidentFacets) {
                os << "Facets: \n";
                for (const auto& ind : zero_slacks_) {
                    global::print(h_rep[ind].first, "", "", os);
                    os << " " << h_rep[ind].second << "\n";
                }
                os << "\n";
            }
        }

        bool V_RepT::hasZeroSlackSuperSet(const std::vector<size_t>& indices) const {
            return std::includes(begin(zero_slacks_),
                                 end(zero_slacks_),
                                 begin(indices),
                                 end(indices));
        }

        DoubleDescriptionMethod::DoubleDescriptionMethod(SCIP *scip, const ResultContainer& bounded_results,
                                                         const ResultContainer& unbounded_results)
                : scip_(scip) {
            for (const auto &bd : bounded_results)
                bounded_.push_back(bd.second);
            for (const auto &unbd : unbounded_results)
                unbounded_.push_back(unbd.second);
        }

        void DoubleDescriptionMethod::printVRep(std::ostream &os, bool withIncidentFacets) const {
            for (const auto& v : v_rep_)
                v.print(os, withIncidentFacets, h_rep_);
        }

        void DoubleDescriptionMethod::computeVRep() {
            assert (!bounded_.empty());
            computeInitialRep(bounded_.front());
            auto current_v_rep = initial_v_rep_;
            for (auto bd = std::next(begin(bounded_)); bd != end(bounded_); ++bd) {
                h_rep_.emplace_back(*bd, 1.0); // add new inequality
                current_v_rep = extendVRep(std::move(current_v_rep),
                                               h_rep_.back(),
                                               h_rep_.size()-1); // extend v-representation w.r.t. new inequality
            }

            for (auto unbd = begin(unbounded_); unbd != end(unbounded_); ++unbd) {
                h_rep_.emplace_back(*unbd, 0.0); // add new inequality
                current_v_rep = extendVRep(std::move(current_v_rep),
                                               h_rep_.back(),
                                               h_rep_.size()-1); // extend v-representation w.r.t. new inequality
            }
            //normalizeVRep(current_v_rep);
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

        vector<V_RepT> DoubleDescriptionMethod::extendVRep(vector<V_RepT> current_rep,
                                                               const H_RepT& constraint,
                                                               size_t index_of_constraint_in_hrep) {
            auto extended_v_rep = vector<V_RepT> {};
            auto plus_inds = vector<std::size_t> {};
            auto minus_inds = vector<std::size_t> {};

            for (size_t i = 0; i < current_rep.size(); ++i) { // partition current v-representation
                auto result = std::inner_product(begin(current_rep[i].weight_), end(current_rep[i].weight_),
                                                 begin(constraint.first), -(current_rep[i].wov_ * constraint.second));
                current_rep[i].addSlack(index_of_constraint_in_hrep, result);
                if (SCIPisNegative(scip_, result)) {
                    minus_inds.push_back(i);
                }
                else if (SCIPisZero(scip_, result)) {
                    current_rep[i].inds_to_slacks_[index_of_constraint_in_hrep] = 0.; // set to zero instead of some small epsilon
                    current_rep[i].addZeroSlackIndex(index_of_constraint_in_hrep);
                    extended_v_rep.push_back(current_rep[i]); // element will also be in extended v-representation
                }
                else {
                    assert(SCIPisPositive(scip_, result));
                    plus_inds.push_back(i);
                    extended_v_rep.push_back(current_rep[i]); // element will also be in extended v-representation
                }
            }
            auto adj_pairs = computeAdjacentPairs(plus_inds, minus_inds, current_rep);
            for (const auto &p : adj_pairs)
                extended_v_rep.emplace_back(scip_, current_rep[p.first], current_rep[p.second], index_of_constraint_in_hrep, h_rep_);

            return extended_v_rep;
        }

        vector<pair<size_t, size_t>> DoubleDescriptionMethod::computeAdjacentPairs(const vector<size_t> &plus_inds,
                                                                                       const vector<size_t> &minus_inds,
                                                                                       const vector<V_RepT> &current_rep) const {
            auto adj_pairs = vector<pair<size_t, size_t>> {};
            for (auto plus : plus_inds) {
                for (auto minus : minus_inds) {
                    assert (plus != minus);
                    if (rayPairIsAdjacent(plus, minus, current_rep))
                        adj_pairs.push_back({plus, minus});
                }
            }
            return adj_pairs;
        }

        vector<size_t> DoubleDescriptionMethod::getCommonZeroSlacks(const V_RepT& v, const V_RepT& w) const {

            auto common_intersec = vector<size_t> {};
            std::set_intersection(begin(v.zero_slacks_),
                                  end(v.zero_slacks_),
                                  begin(w.zero_slacks_),
                                  end(w.zero_slacks_),
                                  std::back_inserter(common_intersec));
            return common_intersec;
        }


        bool DoubleDescriptionMethod::rayPairIsAdjacent(size_t index1,
                                                            size_t index2,
                                                            const vector<V_RepT>& current_rep) const {

            auto intersec = getCommonZeroSlacks(current_rep[index1], current_rep[index2]);

            for (size_t i = 0; i < current_rep.size(); ++i) {
                if (i == index1 || i == index2)
                    continue;
                else if (current_rep[i].hasZeroSlackSuperSet(intersec)) {
                    /* check whether current_rep[i] is multiple of current_rep[index1] or current_rep[index2] */
                    if (!isMultiple(current_rep[i], current_rep[index1]) &&
                        !isMultiple(current_rep[i], current_rep[index2]))
                        return false;
                }
            }
            return true;
        }

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


        void DoubleDescriptionMethod::computeInitialRep(const OutcomeType &bd_outcome) {
            auto size = bd_outcome.size();
            // create initial h_rep
            for (size_t i=0; i<size; ++i) {
                auto constraint = WeightType(size, 0.);
                constraint[i] = 1.;
                h_rep_.push_back({constraint, 0.}); // add constraint: e_i 0 >= 0 with e_i being i-th unit vector
            }
            h_rep_.emplace_back(bd_outcome, 1.); // add constraint; bd_outcome -1 >= 0
            // create initial v_rep
            for (size_t i=0; i<size; ++i) {
                initial_v_rep_.emplace_back(scip_, h_rep_[i].first, bd_outcome[i], h_rep_); // add v_rep: e_i, bd_outcome[i] with e_i being i-th unit vector
            }
            initial_v_rep_.emplace_back(scip_, WeightType(size, 0.), -1., h_rep_); // add v_rep: 0 0 ... 0 -1
        }
    }

}
