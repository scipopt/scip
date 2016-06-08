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

#include "double_description_method.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

#include "global_functions.h"
#include "objscip/objscip.h"
#include "polyscip_types.h"

using std::pair;
using std::size_t;
using std::vector;

namespace polyscip {

    DoubleDescription::DoubleDescription(SCIP *scip, const ResultContainer &bounded_results,
                                         const ResultContainer &unbounded_results)
            : scip_(scip) {
        for (const auto &bd : bounded_results)
            bounded_.push_back(bd.second);
        for (const auto &unbd : unbounded_results)
            unbounded_.push_back(unbd.second);
    }

    void DoubleDescription::computeVRep() {
        assert (!bounded_.empty());
        computeInitialRep(bounded_.front());
        auto current_v_rep = initial_v_rep_;
        for (auto bd = std::next(begin(bounded_)); bd != end(bounded_); ++bd) {
            auto new_constraint = H_RepT(*bd,1);
            if (shouldNormalize(new_constraint, current_v_rep))
                normalizeVRep(current_v_rep);
            current_v_rep = extendVRep(std::move(current_v_rep), new_constraint);
            h_rep_.push_back(new_constraint);
        }
        for (auto unbd = begin(unbounded_); unbd != end(unbounded_); ++unbd) {
            auto new_constraint = H_RepT(*unbd, 0.);
            if (shouldNormalize(new_constraint, current_v_rep))
                normalizeVRep(current_v_rep);
            current_v_rep = extendVRep(std::move(current_v_rep), new_constraint);
            h_rep_.push_back(new_constraint);
        }
        normalizeVRep(current_v_rep);
        auto scip_ptr = scip_;
        std::remove_copy_if(begin(current_v_rep),   // copy all elements from current_v_rep with non-zero weights to v_rep_
                            end(current_v_rep),
                            std::back_inserter(v_rep_),
                            [scip_ptr](const V_RepT& v){return SCIPisZero(scip_ptr, std::accumulate(begin(v.first),
                                                                                                    end(v.first),
                                                                                                    0.0));});
        assert (v_rep_.size() + 1 == current_v_rep.size());
    }

    //todo make simpler
    bool DoubleDescription::shouldNormalize(const H_RepT& constraint, const V_RepContainer& current_v_rep) const {
        if (constraint.second > kLimitForNormalization) {
            return true;
        }
        else if (std::upper_bound(begin(constraint.first), // check whether value in constraint.first > kLimitForNormalization
                                  end(constraint.first),
                                  kLimitForNormalization) != end(constraint.first)) {
            return true;
        }
        else {
            for (const auto& vrep : current_v_rep) {
                if (vrep.second > kLimitForNormalization)
                    return true;
                else if (std::upper_bound(begin(vrep.first), // check whether value in v > kLimitForNormalization
                                          end(vrep.first),
                                          kLimitForNormalization) != end(vrep.first))
                    return true;
            }
            return false;
        }
    }

    vector<DoubleDescription::V_RepT> DoubleDescription::extendVRep(vector<V_RepT> current_rep,
                                                                    const H_RepT &constraint) {

        auto extended_v_rep = vector<V_RepT> {};
        auto plus_inds = vector<std::size_t> {};
        auto minus_inds = vector<std::size_t> {};
        auto zero_slacks = SlackContainer {};

        // partition current v-representation
        for (size_t i = 0; i < current_rep.size(); ++i) {
            zero_slacks.emplace_back(computeZeroSlackSet(current_rep[i]));
            auto result = std::inner_product(begin(current_rep[i].first), end(current_rep[i].first),
                                             begin(constraint.first), -(current_rep[i].second * constraint.second));
            if (SCIPisNegative(scip_, result)) {
                minus_inds.push_back(i);
            }
            else if (SCIPisZero(scip_, result)) {
                extended_v_rep.push_back(current_rep[i]);
            }
            else {
                assert(SCIPisPositive(scip_, result));
                plus_inds.push_back(i);
                extended_v_rep.push_back(current_rep[i]);
            }
        }

        auto adj_pairs = computeAdjacentPairs(plus_inds, minus_inds, zero_slacks, current_rep);
        for (const auto& p : adj_pairs) {
            extended_v_rep.push_back(computeNewRay(current_rep[p.first],
                                                   current_rep[p.second],
                                                   constraint));
        }
        return extended_v_rep;
    }

    vector<pair<size_t, size_t>> DoubleDescription::computeAdjacentPairs(const vector<size_t>& plus_inds,
                                                                       const vector<size_t>& minus_inds,
                                                                       const SlackContainer& zero_slacks,
                                                                       const vector<V_RepT>& current_rep) const {
        auto adj_pairs = vector<pair<size_t, size_t>> {};
        for (const auto& plus_index : plus_inds) {
            for (const auto& minus_index : minus_inds) {
                if (rayPairIsAdjacent(plus_index, minus_index, zero_slacks, current_rep))
                    adj_pairs.push_back({plus_index, minus_index});
            }
        }
        return adj_pairs;
    };

    bool DoubleDescription::rayPairIsAdjacent(size_t index1,
                                            size_t index2,
                                            const SlackContainer& zero_slacks,
                                            const vector<V_RepT>& current_rep) const {
        assert (std::max(index1, index2) < zero_slacks.size());
        auto intersec = vector<size_t> {};
        std::set_intersection(begin(zero_slacks[index1]),
                              end(zero_slacks[index1]),
                              begin(zero_slacks[index2]),
                              end(zero_slacks[index2]),
                              std::back_inserter(intersec));

        for (size_t i=0; i<zero_slacks.size(); ++i) {
            if (i == index1 || i == index2 || zero_slacks[i].size() <= intersec.size())
                continue;
            auto includes = std::includes(begin(zero_slacks[i]),
                                          end(zero_slacks[i]),
                                          begin(intersec),
                                          end(intersec));

            if (includes && !isMultiple(current_rep[i], current_rep[index1]) &&
                !isMultiple(current_rep[i], current_rep[index2]))
                return false;
        }
        return true;
    }

    bool DoubleDescription::isMultiple(const V_RepT& ray, const V_RepT& ray2) const {
        assert (ray.first.size() == ray2.first.size());
        if (SCIPisZero(scip_, ray.second) && !SCIPisZero(scip_, ray2.second)) {
            return false;
        }
        else if (!SCIPisZero(scip_, ray.second) && SCIPisZero(scip_,ray2.second)) {
            return false;
        }
        else if (SCIPisEQ(scip_, ray.second, ray2.second)) {
            for (size_t i=0; i<ray.first.size(); ++i) {
                if (!SCIPisEQ(scip_, ray.first[i], ray2.first[i]))
                    return false;
            }
            return true;
        }
        else {
            // at this point: ray.second != ray2.second, ray.second != 0, ray2.second != 0
            auto multiple = ray.second / ray2.second;
            for (size_t i=0; i<ray.first.size(); ++i) {
                if (!SCIPisEQ(scip_, ray.first[i], ray2.first[i]*multiple))
                    return false;
            }
            return true;
        }
    }

    vector<size_t> DoubleDescription::computeZeroSlackSet(const V_RepT& ray) const {
        auto zeroSet = vector<size_t> {};
        for (size_t i = 0; i < h_rep_.size(); ++i) {
            auto weight_coeff = h_rep_[i].first;
            auto a_coeff = h_rep_[i].second;
            auto result = std::inner_product(begin(weight_coeff), end(weight_coeff),
                                             begin(ray.first), -a_coeff*ray.second);
            if (SCIPisZero(scip_, result))
                zeroSet.push_back(i);
        }
        return zeroSet;
    }

    void DoubleDescription::normalizeVRep(V_RepContainer& v_rep) {
        for (auto &v : v_rep) {
            auto weight_length = std::accumulate(begin(v.first), end(v.first), 0.);
            if (SCIPisPositive(scip_, weight_length)) {
                std::transform(begin(v.first), end(v.first), begin(v.first),
                               [weight_length](const ValueType &val) { return val / weight_length; });
                v.second /= weight_length;
            }
        }
    }

    DoubleDescription::V_RepT DoubleDescription::computeNewRay(const V_RepT &plus_ray, const V_RepT &minus_ray,
                                                       const H_RepT &ineq) const {
        auto size = plus_ray.first.size();
        assert (size == minus_ray.first.size());
        assert (size == ineq.first.size());
        auto m_coeff = std::inner_product(begin(ineq.first),  // m_coeff = ineq \cdot ray_plus
                                          end(ineq.first),
                                          begin(plus_ray.first),
                                          -(ineq.second * plus_ray.second));
        auto p_coeff = std::inner_product(begin(ineq.first),  // p_coeff = ineq \cdot ray_minus
                                          end(ineq.first),
                                          begin(minus_ray.first),
                                          -(ineq.second * minus_ray.second));
        auto new_weight = WeightType(size, 0.);
        std::transform(begin(minus_ray.first), end(minus_ray.first),
                       begin(plus_ray.first), begin(new_weight),
                       [m_coeff, p_coeff](ValueType m_val, ValueType p_val) {
                           return m_coeff * m_val - p_coeff * p_val;
                       });
        // return m_coeff * ray_minus - p_coeff * ray_plus
        return {new_weight, m_coeff * minus_ray.second - p_coeff * plus_ray.second};
    }


    void DoubleDescription::computeInitialRep(const OutcomeType &bd_outcome) {
        auto size = bd_outcome.size();
        initial_v_rep_.emplace_back(WeightType(size, 0), -1.);
        h_rep_.emplace_back(bd_outcome, 1.);
        for (size_t i = 0; i < bd_outcome.size(); ++i) {
            auto ray = WeightType(size, 0.);
            ray[i] = 1.;
            initial_v_rep_.push_back({ray, bd_outcome[i]});
            h_rep_.push_back({ray, 0.});
        }
    }

    void DoubleDescription::deleteZeroWeightRay() {
        std::size_t count {0};
        for (auto it = begin(v_rep_); it != end(v_rep_); ++it) {
            auto weight_length = std::accumulate(begin(it->first), end(it->first), 0.);
            if (SCIPisZero(scip_, weight_length)) {
                v_rep_.erase(it);
                ++count;
            }
        }
        assert (count == 1);
    }
}