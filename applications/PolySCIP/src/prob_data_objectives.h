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
 * @file prob_data_objectives.h
 * @brief Class storing multiple objectives of given problem instance
 * @author Sebastian Schenker
 *
 */

#ifndef POLYSCIP_SRC_PROB_DATA_OBJECTIVES_H_INCLUDED
#define POLYSCIP_SRC_PROB_DATA_OBJECTIVES_H_INCLUDED

#include <cstdlib>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "scip/def.h"

/**
 * @class ProbDataObjectives
 * @brief Stores coefficients and basic methods for objectives of given multi-objective problem
 * @details In order to store the coefficients of our multi-objective problem we use specialisation of SCIP's ObjProbData
 */
class ProbDataObjectives : public scip::ObjProbData {
public:

    /**
     * (Virtual) Destructor
     */
    virtual ~ProbDataObjectives() {};

    /**
     * Number of objectives
     * @return Number of objectives
     */
    std::size_t getNoObjs() const {return non_ignored_objs_.size();};

    /**
     * Add identifier of objective
     * @param name Objective identifier
     */
    void addObjName(const char* name);

    /**
     * Add objective coefficient
     * @param var Corresponding variable
     * @param obj_name Corresponding objective identifier
     * @param val Corresponding coefficient
     */
    void addObjCoeff(SCIP_VAR* var,
                     const char* obj_name,
                     polyscip::ValueType val);

    /**
     * Get objective coefficient
     * @param var Corresponding variable
     * @param obj_no Corresponding number of objective
     * @todo Const qualification
     */
    polyscip::ValueType getObjCoeff(SCIP_VAR* var,
                                    std::size_t obj_no);

    /**
     * Scalar product of given weight and objectives w.r.t. given variable;
     * if given variable is unknown, return 0.0 (since var can only have zero objective
     * coefficients in given problem)
     * @param var Corresponding variable
     * @param weight Corresponding weight vector
     * @todo Const qualification
     */
    polyscip::ValueType getWeightedObjVal(SCIP_VAR* var,
                                          const polyscip::WeightType& weight);

    /**
     * Product of given solution value and objective coefficient w.r.t. given
     * variable and objective number
     * @param var Corresponding variable
     * @param obj_no Corresponding objective number
     * @param sol_val Corresponding solution value
     * @todo Const qualification
     */
    polyscip::ValueType getObjVal(SCIP_VAR* var,
                                  std::size_t obj_no,
                                  polyscip::ValueType sol_val);

    /**
     * Variables corresponding to non-zero objective coefficients
     * @param obj_no Corresponding objective number
     */
    std::vector<SCIP_VAR*> getNonZeroCoeffVars(std::size_t obj_no) const;

    /**
     * Negate all objective coefficients of all variables in all objectives
     */
    void negateAllCoeffs();

    /**
     * Number of non-zero coefficients of objective
     * @param obj_no Corresponding objective number
     */
    std::size_t getNumberNonzeroCoeffs(std::size_t obj_no) const;

    /**
     * Ignore two objectives
     * @param obj_1 First objective to ignore
     * @param obj_2 Second objective to ignore
     */
    void ignoreObjectives(std::size_t obj_1,
                          std::size_t obj_2);

    /**
     * Unignore latest two objectives that were ignored previously
     */
    void unignoreObjectives();

    /**
     * SCIP function for releasing memory
     * @param scip Corresponding SCIP pointer
     */
    virtual SCIP_RETCODE scip_delorig(SCIP* scip) {return SCIP_OKAY;};

private:

    std::unordered_map<std::string, std::size_t> name_to_no_; ///< maps objective identifier to objective number
    std::unordered_map<SCIP_VAR*, polyscip::OutcomeType> var_to_coeffs_; ///< maps SCIP variable to objective coefficients
    std::unordered_map<std::string, std::vector<SCIP_VAR*>> name_to_nonzero_coeffs_; ///< maps objective identifer to non-zero variables
    std::vector<std::string> no_to_name_; ///< maps objective number to objective identifier
    std::vector<std::size_t> non_ignored_objs_; ///< objective indices to be considered
    std::stack<std::size_t> ignored_obj_;     ///< objective indices to be ignored

};

#endif //POLYSCIP_SRC_PROB_DATA_OBJECTIVES_H_INCLUDED
