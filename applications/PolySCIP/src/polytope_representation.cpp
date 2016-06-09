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

    namespace polytoperepresentation {

        V_RepT::V_RepT(WeightType weight, ValueType wov)
        : weight_(weight),
          wov_(wov) {
        }

        V_RepT::V_RepT(WeightType weight, ValueType wov, size_t index)
                : V_RepT(weight, wov) {
            addIncFacetInd(index);
        }

        void V_RepT::addIncFacetInd(size_t index) {
            incident_facet_indices_.push_back(index);
        }

        void V_RepT::print(std::ostream& os, bool withIncidentFacets) const {
            global::print(weight_, "Weight = ", os);
            os << "Coeff = " << wov_ << "\n";
            if (withIncidentFacets)
                global::print(incident_facet_indices_, "Incident facet indices: ", os);
        }

        DoubleDescriptionMethod::DoubleDescriptionMethod(SCIP *scip, const ResultContainer& bounded_results,
                                                         const ResultContainer& unbounded_results)
                : scip_(scip) {
            for (const auto &bd : bounded_results)
                bounded_.push_back(bd.second);
            for (const auto &unbd : unbounded_results)
                unbounded_.push_back(unbd.second);
        }

        void DoubleDescriptionMethod::computeVRep() {
            assert (!bounded_.empty());
            computeInitialRep(bounded_.front());
            auto current_v_rep = initial_v_rep_;
            for (auto bd = std::next(begin(bounded_)); bd != end(bounded_); ++bd) {
                auto new_constraint = H_RepT(*bd, 1);
                /*if (shouldNormalize(new_constraint, current_v_rep))
                    normalizeVRep(current_v_rep);*/
                current_v_rep = extendVRep(std::move(current_v_rep), new_constraint);
                h_rep_.push_back(new_constraint);
            }
            for (auto unbd = begin(unbounded_); unbd != end(unbounded_); ++unbd) {
                auto new_constraint = H_RepT(*unbd, 0.);
                /*if (shouldNormalize(new_constraint, current_v_rep))
                    normalizeVRep(current_v_rep);*/
                current_v_rep = extendVRep(std::move(current_v_rep), new_constraint);
                h_rep_.push_back(new_constraint);
            }
            normalizeVRep(current_v_rep);
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
            assert (v_rep_.size() + 1 == current_v_rep.size());
        }

        //todo make simpler
        bool DoubleDescriptionMethod::shouldNormalize(const H_RepT &constraint, const V_RepContainer &current_v_rep) const {
            if (constraint.second > kLimitForNormalization) {
                return true;
            }
            else if (std::upper_bound(
                    begin(constraint.first), // check whether value in constraint.first > kLimitForNormalization
                    end(constraint.first),
                    kLimitForNormalization) != end(constraint.first)) {
                return true;
            }
            else {
                for (const auto &v : current_v_rep) {
                    if (v.wov_ > kLimitForNormalization)
                        return true;
                    else if (std::upper_bound(begin(v.weight_), // check whether value in v > kLimitForNormalization
                                              end(v.weight_),
                                              kLimitForNormalization) != end(v.weight_))
                        return true;
                }
                return false;
            }
        }


        void DoubleDescriptionMethod::printVRep(std::ostream &os, bool withIncidentFacets) const {
            for (const auto& v : v_rep_)
                v.print(os, withIncidentFacets);
        }

        vector<V_RepT> DoubleDescriptionMethod::extendVRep(vector<V_RepT> current_rep,
                                                           const H_RepT &constraint) {

            auto extended_v_rep = vector<V_RepT> {};
            auto plus_inds = vector<std::size_t> {};
            auto minus_inds = vector<std::size_t> {};
            auto zero_slacks = SlackContainer {};

            // partition current v-representation
            for (size_t i = 0; i < current_rep.size(); ++i) {
                zero_slacks.emplace_back(computeZeroSlackSet(current_rep[i]));
                auto result = std::inner_product(begin(current_rep[i].weight_), end(current_rep[i].weight_),
                                                 begin(constraint.first), -(current_rep[i].wov_ * constraint.second));
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
            for (const auto &p : adj_pairs) {
                extended_v_rep.push_back(computeNewRay(current_rep[p.first],
                                                       current_rep[p.second],
                                                       constraint));
            }
            return extended_v_rep;
        }

        vector<pair<size_t, size_t>> DoubleDescriptionMethod::computeAdjacentPairs(const vector<size_t> &plus_inds,
                                                                             const vector<size_t> &minus_inds,
                                                                             const SlackContainer &zero_slacks,
                                                                             const vector<V_RepT> &current_rep) const {
            auto adj_pairs = vector<pair<size_t, size_t>> {};
            for (const auto &plus_index : plus_inds) {
                for (const auto &minus_index : minus_inds) {
                    if (rayPairIsAdjacent(plus_index, minus_index, zero_slacks, current_rep))
                        adj_pairs.push_back({plus_index, minus_index});
                }
            }
            return adj_pairs;
        };

        bool DoubleDescriptionMethod::rayPairIsAdjacent(size_t index1,
                                                  size_t index2,
                                                  const SlackContainer &zero_slacks,
                                                  const vector<V_RepT> &current_rep) const {
            assert (std::max(index1, index2) < zero_slacks.size());
            auto intersec = vector<size_t> {};
            std::set_intersection(begin(zero_slacks[index1]),
                                  end(zero_slacks[index1]),
                                  begin(zero_slacks[index2]),
                                  end(zero_slacks[index2]),
                                  std::back_inserter(intersec));

            for (size_t i = 0; i < zero_slacks.size(); ++i) {
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

        bool DoubleDescriptionMethod::isMultiple(const V_RepT& v, const V_RepT& w) const {
            assert (v.weight_.size() == w.weight_.size());
            if (SCIPisZero(scip_, v.wov_) && !SCIPisZero(scip_, w.wov_)) {
                return false;
            }
            else if (!SCIPisZero(scip_, v.wov_) && SCIPisZero(scip_, w.wov_)) {
                return false;
            }
            else if (SCIPisEQ(scip_, v.wov_, w.wov_)) {
                for (size_t i = 0; i < v.weight_.size(); ++i) {
                    if (!SCIPisEQ(scip_, v.weight_[i], w.weight_[i]))
                        return false;
                }
                return true;
            }
            else {
                // at this point: ray.second != ray2.second, ray.second != 0, ray2.second != 0
                auto multiple = v.wov_ / w.wov_;
                for (size_t i = 0; i < v.weight_.size(); ++i) {
                    if (!SCIPisEQ(scip_, v.weight_[i], w.weight_[i] * multiple))
                        return false;
                }
                return true;
            }
        }

        vector<size_t> DoubleDescriptionMethod::computeZeroSlackSet(const V_RepT& v) const {
            auto zeroSet = vector<size_t> {};
            for (size_t i = 0; i < h_rep_.size(); ++i) {
                auto weight_coeff = h_rep_[i].first;
                auto a_coeff = h_rep_[i].second;
                auto result = std::inner_product(begin(weight_coeff), end(weight_coeff),
                                                 begin(v.weight_), -a_coeff * v.wov_);
                if (SCIPisZero(scip_, result))
                    zeroSet.push_back(i);
            }
            return zeroSet;
        }

        void DoubleDescriptionMethod::normalizeVRep(V_RepContainer &v_rep) {
            for (auto &v : v_rep) {
                auto weight_length = std::accumulate(begin(v.weight_), end(v.weight_), 0.);
                if (SCIPisPositive(scip_, weight_length)) {
                    std::transform(begin(v.weight_), end(v.weight_), begin(v.weight_),
                                   [weight_length](const ValueType &val) { return val / weight_length; });
                    v.wov_ /= weight_length;
                }
            }
        }

        V_RepT DoubleDescriptionMethod::computeNewRay(const V_RepT &plus_ray, const V_RepT &minus_ray,
                                                                   const H_RepT &ineq) const {
            auto size = plus_ray.weight_.size();
            assert (size == minus_ray.weight_.size());
            assert (size == ineq.first.size());
            auto m_coeff = std::inner_product(begin(ineq.first),  // m_coeff = ineq \cdot ray_plus
                                              end(ineq.first),
                                              begin(plus_ray.weight_),
                                              -(ineq.second * plus_ray.wov_));
            auto p_coeff = std::inner_product(begin(ineq.first),  // p_coeff = ineq \cdot ray_minus
                                              end(ineq.first),
                                              begin(minus_ray.weight_),
                                              -(ineq.second * minus_ray.wov_));
            auto new_weight = WeightType(size, 0.);
            std::transform(begin(minus_ray.weight_), end(minus_ray.weight_),
                           begin(plus_ray.weight_), begin(new_weight),
                           [m_coeff, p_coeff](ValueType m_val, ValueType p_val) {
                               return m_coeff * m_val - p_coeff * p_val;
                           });
            // return m_coeff * ray_minus - p_coeff * ray_plus
            return {new_weight, m_coeff * minus_ray.wov_ - p_coeff * plus_ray.wov_};
        }


        void DoubleDescriptionMethod::computeInitialRep(const OutcomeType &bd_outcome) {
            auto size = bd_outcome.size();
            initial_v_rep_.emplace_back(WeightType(size, 0), -1.);
            h_rep_.emplace_back(bd_outcome, 1.);
            for (size_t i = 0; i < bd_outcome.size(); ++i) {
                auto ray = WeightType(size, 0.);
                ray[i] = 1.;
                h_rep_.push_back({ray, 0.});
                initial_v_rep_.push_back({ray, bd_outcome[i], h_rep_.size()-1});

            }
        }
    }

}
