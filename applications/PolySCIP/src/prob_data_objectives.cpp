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

/**
 * @brief  Problem data for objectives
 * @author Sebastian Schenker
 *
 * Implements problem data for objectives .
 */

#include "prob_data_objectives.h"

#include <algorithm> // std::transform, std::find
#include <cstddef>
#include <functional> // std::negate, std::plus
#include <numeric> // std::inner_product
#include <stack>
#include <stdexcept>
#include <string>

#include "objscip/objscip.h"
#include "polyscip_types.h"

using std::size_t;
using std::string;
using std::vector;
using polyscip::OutcomeType;
using polyscip::ValueType;
using polyscip::WeightType;

void ProbDataObjectives::addObjName(const char* name)  {
    string obj_name(name);
    if (name_to_no_.count(obj_name) != 0)
        throw std::runtime_error("Name of objective already encountered: " + obj_name);
    name_to_no_.emplace(obj_name, getNoObjs());
    name_to_nonzero_coeffs_.emplace(obj_name, vector<SCIP_VAR*>{});
    no_to_name_.push_back(obj_name);
    non_ignored_objs_.push_back(getNoObjs());
}

void ProbDataObjectives::addObjCoeff(SCIP_VAR *var, const char* obj_name, ValueType val) {
    if (val != 0) {
        if (var_to_coeffs_.count(var) == 0)
            var_to_coeffs_.emplace(var, OutcomeType(getNoObjs(), 0));
        auto obj_no = name_to_no_.at(obj_name);
        var_to_coeffs_[var].at(obj_no) = val;
        name_to_nonzero_coeffs_.at(obj_name).push_back(var);
    }
    else
        std::cout << "addObjCoeff with zero value\n";
}

ValueType ProbDataObjectives::getObjCoeff(SCIP_VAR* var, size_t obj) {
    if (var_to_coeffs_.count(var))
        return (var_to_coeffs_[var]).at(non_ignored_objs_.at(obj));
    else
        return 0.;
}

size_t ProbDataObjectives::getNumberNonzeroCoeffs(size_t obj) const {
    auto obj_name = no_to_name_.at(non_ignored_objs_.at(obj));
    return name_to_nonzero_coeffs_.at(obj_name).size();
}

ValueType ProbDataObjectives::getWeightedObjVal(SCIP_VAR* var, const WeightType& weight) {
    assert (weight.size() == non_ignored_objs_.size());
    if (var_to_coeffs_.count(var)) {
        auto& coeffs = var_to_coeffs_[var];
        return std::inner_product(begin(weight),
                                  end(weight),
                                  begin(non_ignored_objs_),
                                  0.,
                                  std::plus<ValueType>(),
                                  [&coeffs](ValueType w, size_t obj){return w*coeffs[obj];});
        /*return std::inner_product(begin(weight),
                                  end(weight),
                                  begin(var_to_coeffs_[var]),
                                  0.);*/
    }
    else {
        return 0.;
    }
}

/*ValueType ProbDataObjectives::getWeightedObjVal(SCIP_VAR* var, const WeightType& weight,
                                                const vector<size_t>& non_redundant_objs) {
    assert (non_redundant_objs.size() == weight.size());
    if (var_to_coeffs_.count(var)) {
        auto coeffs = var_to_coeffs_[var];
        return std::inner_product(begin(weight),
                                  end(weight),
                                  begin(non_redundant_objs),
                                  0.0,
                                  std::plus<ValueType>(),
                                  [&coeffs](ValueType w, size_t index) { return w * coeffs.at(index); });
    }
    else {
        return 0;
    }
}*/

vector<SCIP_VAR*> ProbDataObjectives::getNonZeroCoeffVars(size_t obj) const {
    auto obj_name = no_to_name_.at(non_ignored_objs_.at(obj));
    return name_to_nonzero_coeffs_.at(obj_name);
}

ValueType ProbDataObjectives::getObjVal(SCIP_VAR* var, size_t obj, ValueType sol_val) {
    if (var_to_coeffs_.count(var))
        return var_to_coeffs_[var].at(non_ignored_objs_.at(obj)) * sol_val;
    else
        return 0.;
}

void ProbDataObjectives::negateAllCoeffs() {
    for (auto it=begin(var_to_coeffs_); it!=end(var_to_coeffs_); ++it) {
        std::transform(begin(it->second), end(it->second),
                       begin(it->second), std::negate<ValueType>());
    }
}

void ProbDataObjectives::ignoreObjectives(size_t obj_1, size_t obj_2) {
    auto obj_ind_1 = non_ignored_objs_.at(obj_1);
    auto obj_ind_2 = non_ignored_objs_.at(obj_2);
    auto it = std::find(begin(non_ignored_objs_), end(non_ignored_objs_), obj_ind_1);
    assert (it != end(non_ignored_objs_));
    non_ignored_objs_.erase(it);
    it = std::find(begin(non_ignored_objs_), end(non_ignored_objs_), obj_ind_2);
    assert (it != end(non_ignored_objs_));
    non_ignored_objs_.erase(it);
    ignored_obj_.push(obj_ind_1);
    ignored_obj_.push(obj_ind_2);
}

void ProbDataObjectives::unignoreObjectives() {
    auto obj = ignored_obj_.top();
    ignored_obj_.pop();
    non_ignored_objs_.push_back(obj);
    obj = ignored_obj_.top();
    ignored_obj_.pop();
    non_ignored_objs_.push_back(obj);
    std::sort(begin(non_ignored_objs_), end(non_ignored_objs_));
}
