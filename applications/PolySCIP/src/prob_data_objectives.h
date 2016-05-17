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

/* @brief  Problem data for objectives
 *
 * Problem data for objectives .
 */

#ifndef POLYSCIP_SRC_PROB_DATA_OBJECTIVES_H_INCLUDED
#define POLYSCIP_SRC_PROB_DATA_OBJECTIVES_H_INCLUDED

#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

#include "objscip/objscip.h"
#include "polyscip_types.h"
#include "scip/def.h"

/** Data structure for storing the objective coefficients */
class ProbDataObjectives : public scip::ObjProbData {
public:

    virtual ~ProbDataObjectives() {};

    /** Get number of objectives of given problem
     * @return number of objectives
     */
    std::size_t getNObjs() {return name_to_no_.size();};

    /** Adds objective identifier */
    void addObjName(const char* name);

    /** Adds objective coefficient w.r.t. given variable and objective */
    void addObjCoeff(SCIP_VAR *var, const char* obj_name, polyscip::ValueType val);

    /** Returns scalar product of given weight and objectives w.r.t. given variable;
        if given variable is unknown, return 0.0 (since var can only have zero objective
        coefficients in given problem) */
    polyscip::ValueType getWeightedObjVal(SCIP_VAR* var, const polyscip::WeightType& weight);

    /** Returns product of given solution value and objective coefficient w.r.t. given
        variable and objective number */
    polyscip::ValueType getObjVal(SCIP_VAR* var, std::size_t obj_no, polyscip::ValueType sol_val);

    /** Returns identifier of given objective number */
    std::string getName(std::size_t i) {return no_to_name_.at(i);};

    /** Negates all objective coefficients of all variables */
    void negateAllCoeffs();

    /** SCIP function for releasing memory */
    virtual SCIP_RETCODE scip_delorig(SCIP* scip) {return SCIP_OKAY;};

private:
    /**< maps objective identifier to objective number */
    std::unordered_map<std::string, std::size_t> name_to_no_;

    /**< maps SCIP variable to objective coefficients */
    std::unordered_map<SCIP_VAR*, polyscip::OutcomeType> var_to_coeffs_;

    /**< maps objective number to objective identifier */
    std::vector<std::string> no_to_name_;
};

#endif //POLYSCIP_SRC_PROB_DATA_OBJECTIVES_H_INCLUDED
				      
