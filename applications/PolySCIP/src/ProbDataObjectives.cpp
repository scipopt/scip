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

/**@file   ProbDataObjectives.cpp
 * @brief  Problem data for objectives
 * @author Sebastian Schenker
 *
 * Problem data storing objectives.
 */

#include "ProbDataObjectives.h"
#include "objscip/objscip.h"

#include <algorithm> // std::transform
#include <functional> // std::negate
#include <numeric> // std::inner_product
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using std::cerr;
using std::inner_product;
using std::out_of_range;
using std::string;
using std::vector;

using uint = unsigned;
using objMap = std::unordered_map<std::string, uint>;
using varMap = std::unordered_map<SCIP_VAR*, std::vector<SCIP_Real> >;

/** constructor */
ProbDataObjectives::ProbDataObjectives() {
  objNameToObjNo_ = new objMap();
  varToObjVals_ = new varMap();
}

/** destructor */
ProbDataObjectives::~ProbDataObjectives() {
  delete objNameToObjNo_;
  delete varToObjVals_;
}

/** return number of objectives given in problem */
uint ProbDataObjectives::getNObjs() {
  return objNameToObjNo_->size();
}

/** adds name to data structure holding objective names; if name is already 
    available, name is not added again and false returned; otherwise name is 
    added and true returned. */
bool ProbDataObjectives::addObjName(const string& name) {
  if (objNameToObjNo_->count(name) > 0) {
    cerr << "ERROR in addObjName: " << name << "exists already.\n";
    return false;
  }
  else {
    objNameToObjNo_->emplace(name,getNObjs());
    objNoToObjName_.push_back(name);
  }
  return true;
}

/** adds objective coefficient w.r.t. given variable and objective; 
    return true in case of success */
bool ProbDataObjectives::addObjValue(SCIP_VAR* var, const string& objName, SCIP_Real val) {
  if (varToObjVals_->count(var) == 0) 
    varToObjVals_->emplace(var, vector<SCIP_Real>(getNObjs(),0.0) );
  try {
    uint objNo = objNameToObjNo_->at(objName);
    varToObjVals_->at(var).at(objNo) = val;
  }
  catch (const out_of_range& e) {
    cerr << "ERROR in addObjValue: " << e.what() << "\n";
    return false;
  }
  return true;
}

/** returns identifier of given objective number */
std::string ProbDataObjectives::getObjectiveName(uint index) {
  return objNoToObjName_.at(index);
}

/** return the scalar product of given weight vector and objective function 
    values with respect to given variable var 
    If var is not in varToObjVals return 0 because var has only zero objective function coefficients */
SCIP_Real ProbDataObjectives::getWeightedObjVal(SCIP_VAR* var, const vector<SCIP_Real>* weight) {
  if (varToObjVals_->count(var))
    return inner_product(weight->begin(), weight->end(), (*varToObjVals_)[var].begin(), 0.0);
  else
    return 0.0;
}

/** return product of given solution value and objective coefficient w.r.t. given 
    objective number and variable */
SCIP_Real ProbDataObjectives::getObjVal(SCIP_VAR* var, unsigned objNo, SCIP_Real varSolVal) {
  if (varToObjVals_->count(var))
    return (*varToObjVals_)[var].at(objNo)*varSolVal;
  else
    return 0.0;
}

/* negate all objective coefficients of all variables */
void ProbDataObjectives::negateAllObjCoeffs() {
  for (auto it=varToObjVals_->begin(); it!=varToObjVals_->end(); ++it)
    std::transform(it->second.begin(), it->second.end(), 
                   it->second.begin(), std::negate<SCIP_Real>());
}

/** implementation of SCIP function for releasing memory */
SCIP_RETCODE ProbDataObjectives::scip_delorig(SCIP* scip) {
  return SCIP_OKAY;
}
