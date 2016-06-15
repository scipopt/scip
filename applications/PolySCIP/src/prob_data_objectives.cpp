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

#include "prob_data_objectives.h"

#include <algorithm> // std::transform
#include <cstddef>
#include <functional> // std::negate, std::plus
#include <numeric> // std::inner_product
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
    name_to_no_.emplace(obj_name, getNoAllObjs());
    name_to_nonzero_coeffs_.emplace(obj_name, vector<SCIP_VAR*>{});
    no_to_name_.push_back(obj_name);
}

void ProbDataObjectives::addObjCoeff(SCIP_VAR *var, const char* obj_name, ValueType val) {
    if (val != 0) {
        if (var_to_coeffs_.count(var) == 0)
            var_to_coeffs_.emplace(var, OutcomeType(getNoAllObjs(), 0));
        auto obj_no = name_to_no_.at(obj_name);
        var_to_coeffs_[var].at(obj_no) = val;
        name_to_nonzero_coeffs_.at(obj_name).push_back(var);
    }
    else
        std::cout << "addObjCoeff with zero value\n";
}

ValueType ProbDataObjectives::getObjCoeff(SCIP_VAR* var, size_t obj_no) {
    if (var_to_coeffs_.count(var))
        return (var_to_coeffs_[var]).at(obj_no);
    else
        return 0;
}

std::size_t ProbDataObjectives::getNumberNonzeroCoeffs(std::size_t obj_index) const {
    assert (obj_index < getNoAllObjs());
    auto obj_name = no_to_name_.at(obj_index);
    return name_to_nonzero_coeffs_.at(obj_name).size();
}

/*ValueType ProbDataObjectives::getWeightedObjVal(SCIP_VAR* var, const WeightType& weight) {
    if (var_to_coeffs_.count(var)) {
        return std::inner_product(begin(weight),
                                  end(weight),
                                  begin(var_to_coeffs_[var]),
                                  0.);
    }
    else {
        return 0.0;
    }
}*/

ValueType ProbDataObjectives::getWeightedObjVal(SCIP_VAR* var, const WeightType& weight,
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
}

vector<SCIP_VAR*> ProbDataObjectives::getNonZeroCoeffVars(std::size_t obj_index) const {
    assert (obj_index < getNoAllObjs());
    auto obj_name = no_to_name_.at(obj_index);
    return name_to_nonzero_coeffs_.at(obj_name);
}

ValueType ProbDataObjectives::getObjVal(SCIP_VAR* var, size_t obj_no, ValueType sol_val) {
    if (var_to_coeffs_.count(var))
        return var_to_coeffs_[var].at(obj_no) * sol_val;
    else
        return 0;
}

void ProbDataObjectives::negateAllCoeffs() {
    for (auto it=begin(var_to_coeffs_); it!=end(var_to_coeffs_); ++it) {
        std::transform(begin(it->second), end(it->second),
                       begin(it->second), std::negate<ValueType>());
    }
}

