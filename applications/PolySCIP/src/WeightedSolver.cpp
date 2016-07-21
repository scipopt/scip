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

/**@file   WeightedSolver.cpp
 * @brief  Class providing all methods for using the algorithm
 * @author Sebastian Schenker, Timo Strunk
 *
 * Abstract superclass for a multi objective solver using weighted objective functions.
 */

#include "WeightedSolver.h"

#include "lpi/lpi.h"
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"
#include "ProbDataObjectives.h"
#include "ReaderMOP.h"

#include <algorithm> //std::equal
#include <functional> //std::less_equal
#include <iostream>
#include <string>
#include <sstream>
#include <unordered_set>

using std::string;
using std::vector;


/** SCIP style constructor */
WeightedSolver::WeightedSolver()
  : scip_(NULL),
    objSense_(SCIP_OBJSENSE_MINIMIZE), // default sense is minimization
    nObjs_(0),
    nRuns_(0),
    foundNewOptimum_(false),
    foundNewUnbounded_(false),
    solution_(NULL),
    weight_(NULL),
    cost_vector_(NULL)
{
  SCIPcreate(&scip_);
  assert(scip_ != NULL);
  SCIPincludeDefaultPlugins(scip_);
  SCIPincludeObjReader(scip_, new ReaderMOP(scip_), TRUE);
}

/** destructor */
WeightedSolver::~WeightedSolver() {
  for (auto cost : cost_rays_)
    delete cost;
  for (auto ray : primal_rays_)
    delete ray;
  for (auto weight : primal_ray_weights_)
    delete weight;
  for (auto p : points_)
    delete p;
  for (auto weight : sol_weights_)
    delete weight;
  for (auto sol: solutions_)
    SCIPfreeSol(scip_, &sol);

  SCIPfree(&scip_);
}

void WeightedSolver::printVec(std::ostream& os, const string& str, 
			      const vector<SCIP_Real>* vals, bool negate, unsigned prec) {
  std::stringstream ss;
  ss.precision(prec);
  os << str << "[";
  for (auto it=vals->begin(); it!=vals->end(); ++it) {
    if (negate)
      ss << -(*it);
    else 
      ss << *it;
    if (it != vals->end()-1)
      ss << ", ";
  }
  ss << "]";
  os << ss.str();
};


void WeightedSolver::printRayInfo(std::ostream& os, bool onlyLastAdded, bool withWeights,
				  const string& prelim, unor_Set weaklyEff) {
  unsigned start = 0, 
    size = primal_rays_.size();
  assert (size == cost_rays_.size());
  assert (size == primal_ray_weights_.size());
  if (onlyLastAdded)
    start = size-1;
  os << "\n" + prelim;
  
  bool negate = (objSense_ == SCIP_OBJSENSE_MAXIMIZE) ? true : false;
  
  for (unsigned i=start; i<size; ++i) {
    if (weaklyEff.count(i)>0)
      continue;
    if (withWeights)
      printVec(os, "for weight = ", primal_ray_weights_[i]);

    printVec(os, " ObjRay = ", cost_rays_[i], negate);
    
    os << " corresp. non-zero vars in FeasSpace: ";
    for (auto it=primal_rays_[i]->begin(); it!=primal_rays_[i]->end(); ++it)
      os << it->first << "=" << it->second << "; ";
    os << "\n";
  }
}

void WeightedSolver::printSolInfo(std::ostream& os, bool onlyLastAdded, bool withWeights,
				  const string& prelim, unor_Set weaklyEff) {
  unsigned start = 0,
    size = solutions_.size();
  assert (size == points_.size());
  assert (size == sol_weights_.size());
  if (onlyLastAdded)
    start = size-1;
  os << "\n" + prelim;
  unsigned int nVars = SCIPgetNOrigVars(scip_);
  SCIP_VAR** vars = SCIPgetOrigVars(scip_);
  SCIP_Real val;
  bool negate = (objSense_ == SCIP_OBJSENSE_MAXIMIZE) ? true : false;
  for (unsigned i=start; i<size; ++i) {
    if (weaklyEff.count(i)>0)
      continue;
    if (withWeights)
      printVec(os, "for weight = ", sol_weights_[i]);

    printVec(os, " Point = ", points_[i], negate);

    os << " corresp. non-zero vars in FeasSpace: ";
    for (unsigned j=0; j<nVars; ++j) {
      val = SCIPgetSolVal(scip_, solutions_[i], vars[j]);
      if ( !SCIPisZero(scip_,val))
    	os << SCIPvarGetName(vars[j]) << "=" << val << "; ";
    }
    os << "\n";
  }
}


void WeightedSolver::printObjName(std::ostream& os, unsigned objNo) {
  os << dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_))->getObjectiveName(objNo);
}


/** read SCIP parameter settings file */
SCIP_RETCODE WeightedSolver::readParamSettings(const string& file) {
  SCIP_CALL( SCIPreadParams(scip_, file.data()) );
  return SCIP_OKAY;
}
      
/** reads problem data from file */
SCIP_RETCODE WeightedSolver::readProblem(const string& file, bool beVerbose) {
  SCIP_CALL( SCIPreadProb(scip_, file.data(), "mop") );
  nObjs_ = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_))->getNObjs();
  SCIP_Objsense sense = SCIPgetObjsense(scip_);
  if (sense == SCIP_OBJSENSE_MAXIMIZE) { // objective sense of read problem is maximization
    objSense_ = SCIP_OBJSENSE_MAXIMIZE;
    //internally we treat problem as min problem and negate objective values
    SCIPsetObjsense(scip_, SCIP_OBJSENSE_MINIMIZE); 
    dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_))->negateAllObjCoeffs();
  }
  if (beVerbose) {
    std::cout << "No of objectives: " << nObjs_ << std::endl;
    std::cout << "Objective sense: ";
    if (objSense_ == SCIP_OBJSENSE_MAXIMIZE)
      std::cout << "MAXIMIZE\n";
    else if (objSense_ == SCIP_OBJSENSE_MINIMIZE)
      std::cout << "MINIMIZE\n";
    else
      std::cout << "UNKNOWN.\n";
  }
  return SCIP_OKAY;
}

/** returns true if the last run lead to a new efficient solution */
bool WeightedSolver::foundNewOptimum() const {
  return foundNewOptimum_;
}

bool WeightedSolver::foundNewUnbounded() const {
  return foundNewUnbounded_;
}

/** returns the number of all found pareto optima */
unsigned long WeightedSolver::getNSolutions() const {
   return points_.size();
}

/** returns the number of objective functions */
unsigned int WeightedSolver::getNObjs() const {
  return nObjs_;
}

/** returns the number of weighted runs so far */
unsigned long WeightedSolver::getNRuns() const {
   return nRuns_;
}

void WeightedSolver::deleteWeaklyNondomPointsAndCorrespSols() {
  vector<const vector<SCIP_Real>* >::iterator iterP = points_.begin();
  vector<SCIP_SOL*>::iterator iterSol = solutions_.begin();
  vector<const vector<SCIP_Real>* >::iterator iterW = sol_weights_.begin();
  while (iterP != points_.end()) {
    if (isDominated(iterP)) {
      bool negate = (objSense_ == SCIP_OBJSENSE_MAXIMIZE) ? true : false;
      printVec(std::cout, "\n#deleting weakly nondom point: ", *iterP, negate);
      iterP = points_.erase(iterP);
      iterSol = solutions_.erase(iterSol);
      iterW = sol_weights_.erase(iterW);
    }
    else {
      ++iterP;
      ++iterSol;
      ++iterW;
    }
  }
}

/** return true if p is dominated; return false otherwise */
bool WeightedSolver::isDominated(vector<const vector<SCIP_Real>* >::const_iterator iterP) {
  for (vector<const vector< SCIP_Real>* >::iterator iterQ=points_.begin(); 
       iterQ!=points_.end(); ++iterQ) {
    if (iterQ == iterP)
      continue;
    else if (std::equal((*iterQ)->begin(), (*iterQ)->end(), (*iterP)->begin(), 
			std::less_equal<SCIP_Real>()))
      return true;
  }
  return false;
}
