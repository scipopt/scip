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
 * @file prob_data_objectives.cpp
 * @brief  Implements class storing multiple objectives of given problem instance
 * @author Sebastian Schenker
 *
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


/**
 * Add identifier of objective
 * @param name Objective identifier
 */
void ProbDataObjectives::addObjName(const char* name)  {
    string obj_name(name);
    if (name_to_no_.count(obj_name) != 0)
        throw std::runtime_error("Name of objective already encountered: " + obj_name);
    name_to_no_.emplace(obj_name, getNoObjs());
    name_to_nonzero_coeffs_.emplace(obj_name, vector<SCIP_VAR*>{});
    no_to_name_.push_back(obj_name);
    non_ignored_objs_.push_back(getNoObjs());
}


/**
 * Add objective coefficient
 * @param var Corresponding variable
 * @param obj_name Corresponding objective identifier
 * @param val Corresponding coefficient
 */
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


/**
 * Get objective coefficient
 * @param var Corresponding variable
 * @param obj_no Corresponding number of objective
 * @todo Const qualification
 */
ValueType ProbDataObjectives::getObjCoeff(SCIP_VAR* var, size_t obj) {
    if (var_to_coeffs_.count(var))
        return (var_to_coeffs_[var]).at(non_ignored_objs_.at(obj));
    else
        return 0.;
}


/**
 * Number of non-zero coefficients of objective
 * @param obj_no Corresponding objective number
 */
size_t ProbDataObjectives::getNumberNonzeroCoeffs(size_t obj) const {
    auto obj_name = no_to_name_.at(non_ignored_objs_.at(obj));
    return name_to_nonzero_coeffs_.at(obj_name).size();
}


/**
 * Scalar product of given weight and objectives w.r.t. given variable;
 * if given variable is unknown, return 0.0 (since var can only have zero objective
 * coefficients in given problem)
 * @param var Corresponding variable
 * @param weight Corresponding weight vector
 * @todo Const qualification
 */
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
    }
    else {
        return 0.;
    }
}


/**
 * Variables corresponding to non-zero objective coefficients
 * @param obj_no Corresponding objective number
 */
vector<SCIP_VAR*> ProbDataObjectives::getNonZeroCoeffVars(size_t obj) const {
    auto obj_name = no_to_name_.at(non_ignored_objs_.at(obj));
    return name_to_nonzero_coeffs_.at(obj_name);
}


/**
 * Product of given solution value and objective coefficient w.r.t. given
 * variable and objective number
 * @param var Corresponding variable
 * @param obj_no Corresponding objective number
 * @param sol_val Corresponding solution value
 * @todo Const qualification
 */
ValueType ProbDataObjectives::getObjVal(SCIP_VAR* var, size_t obj, ValueType sol_val) {
    if (var_to_coeffs_.count(var))
        return var_to_coeffs_[var].at(non_ignored_objs_.at(obj)) * sol_val;
    else
        return 0.;
}


/**
 * Negate all objective coefficients of all variables in all objectives
 */
void ProbDataObjectives::negateAllCoeffs() {
    for (auto it=begin(var_to_coeffs_); it!=end(var_to_coeffs_); ++it) {
        std::transform(begin(it->second), end(it->second),
                       begin(it->second), std::negate<ValueType>());
    }
}

/**
 * Ignore two objectives
 * @param obj_1 First objective to ignore
 * @param obj_2 Second objective to ignore
 */
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

/**
 * Unignore latest two objectives that were ignored previously
 */
void ProbDataObjectives::unignoreObjectives() {
    auto obj = ignored_obj_.top();
    ignored_obj_.pop();
    non_ignored_objs_.push_back(obj);
    obj = ignored_obj_.top();
    ignored_obj_.pop();
    non_ignored_objs_.push_back(obj);
    std::sort(begin(non_ignored_objs_), end(non_ignored_objs_));
}
