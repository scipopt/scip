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
#include <utility>
#include <vector>

#include "objscip/objscip.h"
#include "polyscip_types.h"

using std::pair;
using std::vector;

namespace polyscip {

    V_Representation::V_Representation(SCIP* scip, const ResultContainer& bounded_results, const ResultContainer& unbounded_results)
            : scip_(scip)
    {
        for (const auto& bd : bounded_results)
            bounded_.push_back(bd.second);
        for (const auto& unbd : unbounded_results)
            unbounded_.push_back(unbd.second);
    }

    void V_Representation::computeV_Rep() {
        assert (!bounded_.empty());
        computeInitialRep(bounded_.front());
        auto current_v_rep = initial_v_rep_;
        for (auto bd_it = std::next(begin(bounded_)); bd_it!=end(bounded_); ++bd_it) {
            current_h_rep_.emplace_back({(begin(bd_it), end(bd_it)), -1.});
            current_v_rep = extendV_Rep(std::move(current_v_rep), current_h_rep_.back());
        }
        for (auto unbd_it=begin(unbounded_); unbd_it!=end(unbounded_); ++unbd_it) {
            current_h_rep_.emplace_back({(begin(unbd_it), end(unbd_it)), 0.});
            current_v_rep = extendV_Rep(std::move(current_v_rep), current_h_rep_.back());
        }
        v_rep_ = current_v_rep;
    }

    vector<V_RepT> V_Representation::extendV_Rep(vector<V_RepT> current_v_rep, const H_RepT& ineq) {
        auto extended_v_rep = vector<V_RepT> {};
        auto plus_inds = vector<std::size_t> {};
        auto minus_inds = vector<std::size_t> {};
        auto zero_inds = vector<std::size_t> {};
        // partition current v-representation
        for (std::size_t i = 0; i < current_v_rep.size(); ++i) {
            assert (current_v_rep[i].first.size() == ineq.first.size());
            auto result = std::inner_product(begin(current_v_rep[i].first), end(current_v_rep[i].first),
                                             begin(ineq.first), -current_v_rep[i].second*ineq.second);
            if (SCIPisNegative(scip_, result)) {
                minus_inds.push_back(i);
            }
            else if (SCIPisZero(scip_, result)) {
                extended_v_rep.push_back(current_v_rep[i]);
            }
            else {
                assert(SCIPisPositive(scip_, result));
                plus_inds.push_back(i);
                extended_v_rep.push_back(current_v_rep[i]);
            }
        }
        // compute zero sets for rays of plus and minus partition
        /*auto plus_inds_zero_set = vector<vector<std::size_t>> {};
        for (const auto& index : plus_inds)
            plus_inds_zero_set.push_back(computeZeroSet(current_v_rep[index]));
        auto minus_inds_zero_set = vector<vector<std::size_t>> {};
        for (const auto& index : minus_inds)
            minus_inds_zero_set.push_back(computeZeroSet(current_v_rep[index]));*/
        // compute new rays of adjacent rays of plus and minus partition
        for (const auto& p_ind : plus_inds) {
            for (const auto& m_ind : minus_inds) {
                extended_v_rep.push_back(getNewVertex(current_v_rep[p_ind],
                                                      current_v_rep[m_ind],
                                                      ineq));
            }
        }
        return extended_v_rep;
    }

    V_RepT V_Representation::getNewVertex(const V_RepT& ray_plus, const V_RepT& ray_minus, const H_RepT& ineq) const {
        auto size = ray_plus.first.size();
        assert (size == ray_minus.first.size());
        assert (size == ineq.first.size());
        // m_coeff = ineq \cdot ray_plus
        auto m_coeff = std::inner_product(begin(ineq.first), end(ineq.first),
                                          begin(ray_plus.first), -ineq.second*ray_plus.second);
        // p_coeff = ineq \cdot ray_minus
        auto p_coeff = std::inner_product(begin(ineq.first), end(ineq.first),
                                          begin(ray_minus.first), -ineq.second*ray_minus.second);
        // return m_coeff * ray_minus - p_coeff * ray_plus
        auto new_weight = WeightType(size,0.);
        std::transform(begin(ray_minus.first), end(ray_minus.first),
                       begin(ray_plus.first), begin(new_weight),
                       [m_coeff, p_coeff](ValueType m_val, ValueType p_val){return m_coeff*m_val - p_coeff*p_val;});
        return {new_weight, m_coeff*ray_minus.second - p_coeff*ray_plus.second};
    }

    vector<std::size_t> V_Representation::computeZeroSet(const pair <WeightType, ValueType> &ray) {
        auto zeroSet = vector < std::size_t > {};
        for (std::size_t i = 0; i < current_h_rep_.size(); ++i) {
            auto lhs = current_h_rep_[i].first;
            auto rhs = current_h_rep_[i].second;
            auto result = std::inner_product(begin(lhs), end(lhs),
                                             begin(lhs), -ray.second * rhs);
            if (SCIPisZero(scip_, result))
                zeroSet.push_back(i);
        }
        assert(zeroSet.size() == current_h_rep_.size() - 1);
        std::sort(begin(zeroSet), end(zeroSet));
        return zeroSet;
    }


    void V_Representation::computeInitialRep(const OutcomeType& bd_outcome) {
        auto size = bd_outcome.size();
        initial_v_rep_.emplace_back({(size, 0), -1.});
        current_h_rep_.emplace_back({(begin(bd_outcome),end(bd_outcome)), -1.});
        for (std::size_t i = 0; i < bd_outcome.size(); ++i) {
            auto ray = WeightType(size, 0.);
            ray[i] = 1.;
            initial_v_rep_.push_back({ray, bd_outcome[i]});
            current_h_rep_.push_back({ray, 0.});
        }
    }

}