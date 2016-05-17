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
#include <functional> // std::negate
#include <numeric> // std::inner_product
#include <stdexcept>
#include <string>

#include "objscip/objscip.h"
#include "polyscip_types.h"

using std::string;
using polyscip::OutcomeType;
using polyscip::ValueType;
using polyscip::WeightType;

void ProbDataObjectives::addObjName(const char* name)  {
    string obj_name(name);
    if (name_to_no_.count(obj_name) != 0)
        throw std::runtime_error("Name of objective already encountered: " + obj_name);
    name_to_no_.insert({obj_name, getNObjs()});
    no_to_name_.push_back(obj_name);
}

void ProbDataObjectives::addObjCoeff(SCIP_VAR *var, const char* obj_name, ValueType val) {
    if (var_to_coeffs_.count(var) == 0)
        var_to_coeffs_.emplace(var,OutcomeType(getNObjs(),0.));
    auto obj_no = name_to_no_.at(obj_name);
    var_to_coeffs_[var].at(obj_no) = val;
}

ValueType ProbDataObjectives::getWeightedObjVal(SCIP_VAR* var, const WeightType& weight) {
    if (var_to_coeffs_.count(var)) {
        return std::inner_product(begin(weight),
                                  end(weight),
                                  begin(var_to_coeffs_[var]),
                                  0.);
    }
    else {
        return 0.0;
    }
}

ValueType ProbDataObjectives::getObjVal(SCIP_VAR* var, std::size_t obj_no, ValueType sol_val) {
    if (var_to_coeffs_.count(var))
        return var_to_coeffs_[var].at(obj_no) * sol_val;
    else
        return 0.0;
}

void ProbDataObjectives::negateAllCoeffs() {
    for (auto it=begin(var_to_coeffs_); it!=end(var_to_coeffs_); ++it) {
        std::transform(begin(it->second), end(it->second),
                       begin(it->second), std::negate<ValueType>());
    }
}

