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

/**@file   ProbDataObjectives.h
 * @brief  Problem data for objectives 
 * @author Sebastian Schenker
 *
 * Problem data storing objectives.
 */

#ifndef GUARD_PROBDATA_OBJECTIVES
#define GUARD_PROBDATA_OBJECTIVES

#include "objscip/objscip.h"
#include "scip/def.h"

#include <unordered_map>
#include <string>
#include <vector>

/** data structure for storing the objective coefficients */
class ProbDataObjectives : public scip::ObjProbData {
 public:
  ProbDataObjectives(); /**constructor */
  virtual ~ProbDataObjectives(); /** destructor */

  /** return number of objectives given in problem */
  unsigned getNObjs(); 

  /** add identifier of objective; return true in case of success */
  bool addObjName(const std::string& name); 
  
  /** adds objective coefficient w.r.t. given variable and objective; 
      return true in case of success */
  bool addObjValue(SCIP_VAR* var, const std::string& objName, SCIP_Real val); 
  
  /** return the scalar product of given weight and objectives w.r.t. given variable; 
      if var was not encountered before, return 0.0 because var can only have zero objective 
      function coefficients in given problem file */
  SCIP_Real getWeightedObjVal(SCIP_VAR* var, const std::vector<SCIP_Real>* weight);

  /** return product of given solution value and objective coefficient w.r.t. given 
      objective number and variable */
  SCIP_Real getObjVal(SCIP_VAR* var, unsigned objNo, SCIP_Real varSolVal); 
  
  /** returns identifier of given objective number */
  std::string getObjectiveName(unsigned i);

  /** negates all objective coefficients of all variables */
  void negateAllObjCoeffs();
  
  virtual SCIP_RETCODE scip_delorig(SCIP* scip);
  
 private:
  std::unordered_map<std::string,unsigned>* objNameToObjNo_; /**< maps objective identifier 
								 to objective number */
  
  std::unordered_map<SCIP_VAR*,std::vector<SCIP_Real> >* varToObjVals_; /**< maps problem variable to
									   objective coefficients */

  std::vector<std::string> objNoToObjName_; /**< maps objective number to objective identifier */
};


#endif					  
				      
