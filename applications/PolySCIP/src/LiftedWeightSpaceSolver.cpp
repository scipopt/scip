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

/**@file   LiftedWeightSpaceSolver.cpp
 * @brief  Main class of the algorithm
 * @author Sebastian Schenker, Timo Strunk
 *
 * Realization of a weighted solver using the lifted weight space polyhedron to calculate weights.
 * It gets the weight from the skeleton, then changes the objective in the SCIP problem to the weighted
 * objective.  Then it solves the problem and uses the skeleton to determine if the new solution is
 * an extremal supported nondominated point.
 */

#include "LiftedWeightSpaceSolver.h"
#include "ProbDataObjectives.h"
#include "ReaderMOP.h"
#include "objscip/objscip.h"
#include "scip/struct_sol.h"
#include "Skeleton.h"

#include <algorithm>   // std::transform, std::fill
#include <cmath>       // std::abs
#include <functional>  // std::plus
#include <iostream>
#include <limits>   // std::numeric_limits<>::max()
#include <map>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>

using std::cout;
using std::fill;
using std::map;
using std::string;
using std::transform;
using std::vector;

/** SCIP style constructor */
LiftedWeightSpaceSolver::LiftedWeightSpaceSolver()
  : WeightedSolver(),
    solving_stage_(MULTIOPT_UNSOLVED),
    timeLimit_(false),
    timeLimitVal_(std::numeric_limits<unsigned long>::max()),
    mip_status_(SCIP_STATUS_UNKNOWN),
    initialWeight_(nullptr),
    nNodesLastRun_(0),
    nIterationsLastRun_(0),
    durationLastRun_(0),
    skeleton_(nullptr)
{
  SCIPcreateClock(scip_, &clock_iteration_);
  SCIPcreateClock(scip_, &clock_total_);
}

/** default destructor */
LiftedWeightSpaceSolver::~LiftedWeightSpaceSolver()
{
  SCIPfreeClock(scip_, &clock_total_);
  SCIPfreeClock(scip_, &clock_iteration_);
  delete initialWeight_;
  delete skeleton_;
}

/** set time limit for entire computation */
void LiftedWeightSpaceSolver::setTimeLimit(unsigned long timeLimit) {
  timeLimit_ = true;
  timeLimitVal_ = timeLimit;
}

/** returns true if there is a weight left to check */
bool LiftedWeightSpaceSolver::hasUncheckedWeight() const {
  return solving_stage_ != MULTIOPT_SOLVED; 
}

/** solve instance with next weight */
SCIP_RETCODE LiftedWeightSpaceSolver::solveUncheckedWeight(bool beVerbose) {
  foundNewOptimum_ = false;
  foundNewUnbounded_ = false;
  nNodesLastRun_ = 0;
  nIterationsLastRun_ = 0;
  durationLastRun_ = 0;
  cost_vector_ = NULL;
  solution_ = NULL;

  SCIP_CALL( SCIPresetClock(scip_, clock_iteration_) );
  SCIP_CALL( SCIPstartClock(scip_, clock_iteration_) );
  SCIP_CALL( loadWeight(false) );
  SCIP_CALL( solveWeighted() ); 
  SCIP_CALL( evalSolution() ); // solving_stage_ is set
  nNodesLastRun_ = SCIPgetNNodes(scip_);
  nIterationsLastRun_ = SCIPgetNLPIterations(scip_);
  SCIP_CALL( SCIPstopClock(scip_, clock_iteration_) );
  durationLastRun_ = SCIPgetClockTime(scip_, clock_iteration_);
  ++nRuns_;
  if (solving_stage_ == MULTIOPT_SOLVED)
    SCIP_CALL( SCIPstopClock(scip_, clock_total_) );
  if (beVerbose)
    printVerboseInfo(true);
  
  return SCIP_OKAY;
}

/** prints weight, objective vector and solution 
 * (and skeleton information) to standard output 
 */
void LiftedWeightSpaceSolver::printVerboseInfo(bool withSkeletonInfo) {
  printVec(std::cout, "\nConsidered weight: ", weight_);
  if (foundNewUnbounded_)
    cout << " Weight yielded new ray.\n";
  else if (foundNewOptimum_)
    cout << " Weight yielded new solution.\n";
  else 
    cout << " Weight yielded already known solution.\n";
  
  cout << "\nNo of SCIP nodes used in this run: " 
       << getNNodesLastRun() << std::endl;
  cout << "No of lp iterations in this run: " 
       << getNLPIterationsLastRun() << std::endl;
  cout << "Computation time in seconds in this run: " 
       << getDurationLastRun() << std::endl;

  if (withSkeletonInfo) {
    cout << "No of new verts in 1-skeleton in this run: " 
	 << getNNewSkeletVerts() << std::endl;
    cout << "No of processed verts in 1-skeleton in this run: " 
	 << getNProcSkeletVerts() << std::endl;
  }
}


/** set weight_ and change objective function accordingly */
SCIP_RETCODE LiftedWeightSpaceSolver::loadWeight(bool initPhase) {
  if (initPhase)
    weight_ = initialWeight_;
  else
    weight_ = skeleton_->nextWeight(); 
  if( SCIPisTransformed(scip_) )
    SCIP_CALL( SCIPfreeTransform(scip_) );
  
  ProbDataObjectives* objs = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
  SCIP_VAR** vars = SCIPgetOrigVars(scip_);
  unsigned nVars = SCIPgetNOrigVars(scip_);
  for(unsigned i=0; i<nVars; ++i) {
    SCIP_Real val = objs->getWeightedObjVal(vars[i], weight_);
    SCIP_CALL( SCIPchgVarObj(scip_, vars[i], val) );
  }
  return SCIP_OKAY;
}


/** compute first non-dominated point and initialise weight space or
 * find out that there is no non-dominated point
 */
SCIP_RETCODE LiftedWeightSpaceSolver::initialise(bool beVerbose) {
  bool checkRay = false;
  unsigned objCounter = 0;
  initialWeight_ = new vector<SCIP_Real>(nObjs_,0.);
  SCIP_CALL( SCIPstartClock(scip_, clock_total_) );
  foundNewOptimum_ = false;

  while (objCounter < nObjs_) {
    foundNewUnbounded_ = false;
    SCIPresetClock(scip_, clock_iteration_);
    SCIP_CALL( SCIPstartClock(scip_, clock_iteration_) );
    (*initialWeight_)[objCounter] = 1.;
    SCIP_CALL( loadWeight(true) );
    SCIP_CALL( solveWeighted() ); //mip_status_ is set
    ++nRuns_;
    nNodesLastRun_ = SCIPgetNNodes(scip_);
    nIterationsLastRun_ = SCIPgetNLPIterations(scip_);
    SCIP_CALL( SCIPstopClock(scip_, clock_iteration_) );
    durationLastRun_ = SCIPgetClockTime(scip_, clock_iteration_);
    
    if (mip_status_ == SCIP_STATUS_INFORUNBD) { //separate status: either infeasible or unbounded
      (*initialWeight_)[objCounter] = 0.0;   // re-compute with zero objective
      SCIP_CALL( loadWeight(true) );
      SCIP_CALL( solveWeighted() ); //mip_status_ is set
      assert (mip_status_ == SCIP_STATUS_INFEASIBLE || mip_status_ == SCIP_STATUS_OPTIMAL);
      if (mip_status_ == SCIP_STATUS_OPTIMAL) { // status with non-zero initialWeight was UNBOUNDED
	SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, TRUE) );
	(*initialWeight_)[objCounter] = 1.;
	SCIP_CALL( loadWeight(true) );
	SCIP_CALL( solveWeighted() ); 
	mip_status_ = SCIP_STATUS_UNBOUNDED; 
	checkRay = true;
      }
    }

    if (mip_status_ == SCIP_STATUS_OPTIMAL) {
      cost_vector_ = setSolAndGetCostVector();
      assert( cost_vector_ != nullptr);
      assert( !hasInfiniteComponent(cost_vector_) );
      addSolAndObjValsAndWeight();
      foundNewOptimum_ = true;
      if (beVerbose)
	printVerboseInfo();
      break; 
    }
    else if (mip_status_ == SCIP_STATUS_UNBOUNDED) {
      if (checkRay) {
	assert (SCIPhasPrimalRay(scip_));
	SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_DEFAULT, TRUE) );
	checkRay = false;
      }
      if (!SCIPhasPrimalRay(scip_)) {
	SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, TRUE) );
	checkRay = true;
	continue; // re-compute for same objective without presolving
      }
      cost_vector_ = setPrimalRayAndGetCostVector();
      assert( cost_vector_ != nullptr);
      cost_rays_.push_back(cost_vector_);
      primal_ray_weights_.push_back(new vector<SCIP_Real>(*initialWeight_));
      foundNewUnbounded_ = true;
      if (beVerbose)
	printVerboseInfo(false);
      (*initialWeight_)[objCounter++] = 0.0;
    }
    else if (mip_status_ == SCIP_STATUS_INFEASIBLE) {
      std::cout << "PROBLEM IS INFEASIBLE.\n";
      solving_stage_ = MULTIOPT_SOLVED;
      break;
    }
    else if (mip_status_ == SCIP_STATUS_TIMELIMIT) {
      std::cout << "TIME LIMIT REACHED\n";
      solving_stage_ = MULTIOPT_SOLVED;
      break;
    }
    else {
      std::cout << "UNEXPECTED SCIP STATUS\n";
      solving_stage_ = MULTIOPT_SOLVED;
      break;
    }
  }

  if (foundNewOptimum_ && nObjs_ > 1) {
    skeleton_ = new Skeleton(scip_, nObjs_);
    skeleton_->init(objCounter, cost_vector_, cost_rays_);
    solving_stage_ = MULTIOPT_SOLVING;
  }
  else {
    solving_stage_ = MULTIOPT_SOLVED;
  }

  if (solving_stage_ == MULTIOPT_SOLVED)
    SCIP_CALL( SCIPstopClock(scip_, clock_total_) );    
  
  return SCIP_OKAY;
}

/** set solution, corresp. objective val and corresp. weight */
void LiftedWeightSpaceSolver::addSolAndObjValsAndWeight() {
  solutions_.push_back(solution_);
  points_.push_back(cost_vector_);
  sol_weights_.push_back(new vector<SCIP_Real>(*weight_));
}


/** returns the number of branching nodes used in the last run */
SCIP_Longint LiftedWeightSpaceSolver::getNNodesLastRun() const {
   return nNodesLastRun_;
}

/** returns the number of LP iterations done in the last run */
SCIP_Longint LiftedWeightSpaceSolver::getNLPIterationsLastRun() const {
   return nIterationsLastRun_;
}

/** returns the compuational time taken in the last run */
SCIP_Real LiftedWeightSpaceSolver::getDurationLastRun() const {
  return durationLastRun_;
}

/** sets found solution and returns corresp. objective vector */
vector<SCIP_Real>* LiftedWeightSpaceSolver::setSolAndGetCostVector() {
  SCIP_SOL* bestSol = SCIPgetBestSol(scip_);
  SCIP_SOL* finite = NULL;
  SCIP_Bool copySuccess = FALSE;
  SCIPcreateFiniteSolCopy(scip_, &finite, bestSol, &copySuccess);
  if (!copySuccess) 
    assert (std::abs(SCIPgetSolOrigObj(scip_,bestSol)-SCIPgetSolOrigObj(scip_,finite)) < 1.0e-5);
  solution_ = finite;
  ProbDataObjectives* objs = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
  unsigned int nVars = SCIPgetNOrigVars(scip_);
  SCIP_VAR** vars = SCIPgetOrigVars(scip_);
  vector<SCIP_Real>* result = new vector<SCIP_Real>(nObjs_, 0.0);
  vector<SCIP_Real> tmp(nObjs_, 0.0);
  SCIP_Real val;
  for (unsigned i=0; i<nVars; ++i) {
    val = SCIPgetSolVal(scip_, solution_, vars[i]);
    for (unsigned j=0; j<nObjs_; ++j) 
      tmp[j] = objs->getObjVal(vars[i], j, val);
    transform(result->begin(), result->end(), tmp.begin(), result->begin(), std::plus<SCIP_Real>());
  }
  return result;
}

/** sets found ray and returns corresp. objective vector */
vector<SCIP_Real>* LiftedWeightSpaceSolver::setPrimalRayAndGetCostVector() {
  ProbDataObjectives* objs = dynamic_cast<ProbDataObjectives*>(SCIPgetObjProbData(scip_));
  unsigned int nVars = SCIPgetNOrigVars(scip_);
  SCIP_VAR** vars = SCIPgetOrigVars(scip_);
  map<string,SCIP_Real>* rayVals = new map<string,SCIP_Real>();
  vector<SCIP_Real>* result = new vector<SCIP_Real>(nObjs_, 0.0);
  vector<SCIP_Real> tmp(nObjs_, 0.0);
  assert (SCIPhasPrimalRay(scip_));
  SCIP_Real val;
  for (unsigned i=0; i<nVars; ++i) {
    val = SCIPgetPrimalRayVal(scip_, vars[i]);
    if (!SCIPisZero(scip_,val)) {
      (*rayVals)[SCIPvarGetName(vars[i])] = val;
      for (unsigned j=0; j<nObjs_; ++j) 
	tmp[j] = objs->getObjVal(vars[i], j, val);
    }
    else 
      fill(tmp.begin(), tmp.end(), 0.0);
    transform(result->begin(), result->end(), tmp.begin(), result->begin(), std::plus<double>());
  }
  primal_rays_.push_back(rayVals);
  return result;
}

/** call the mip solver */
SCIP_RETCODE LiftedWeightSpaceSolver::solveWeighted()
{
  if (timeLimit_) // set SCIP timelimit so that total algorithm timelimit is met
    SCIP_CALL( SCIPsetRealParam(scip_, "limits/time", 
				timeLimitVal_ - SCIPgetClockTime(scip_, clock_total_)) );

  SCIP_CALL( SCIPsolve(scip_) );    // actual SCIP solver call 
  mip_status_ = SCIPgetStatus(scip_);
  return SCIP_OKAY;
}

/** get total time for algorithm */
SCIP_Real LiftedWeightSpaceSolver::getTotalDuration() const {
   return SCIPgetClockTime(scip_, clock_total_);
}

/** return true if the given vector has an entry close to infinity */
bool LiftedWeightSpaceSolver::hasInfiniteComponent(const std::vector<SCIP_Real>* cost_vector) {
   for (std::vector<SCIP_Real>::const_iterator it=cost_vector_->begin();
        it!=cost_vector_->end(); ++it) {
     if (*it >= SCIPinfinity(scip_)/1000.)
       return true;
   }
   return false;
}


/** evaluates solution status */
SCIP_RETCODE LiftedWeightSpaceSolver::evalSolution() {
  if (mip_status_ == SCIP_STATUS_OPTIMAL) {
    cost_vector_ = setSolAndGetCostVector();
    assert( cost_vector_ != NULL);
    assert( !hasInfiniteComponent(cost_vector_) );
    foundNewOptimum_ = skeleton_->isExtremal(cost_vector_);
    if (foundNewOptimum_) 
      addSolAndObjValsAndWeight();
    else {
      delete cost_vector_;
      SCIP_CALL( SCIPfreeSol(scip_, &solution_) );
    }
  }
  else if (mip_status_ == SCIP_STATUS_INFORUNBD || mip_status_ == SCIP_STATUS_UNBOUNDED) {
    if (!SCIPhasPrimalRay(scip_)) {
      SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_OFF, TRUE) );
      SCIP_CALL( SCIPfreeTransform(scip_) );
      SCIP_CALL( solveWeighted() ); // mip_status_ is re-set
      assert (SCIPhasPrimalRay(scip_) && 
	      (mip_status_ == SCIP_STATUS_INFORUNBD || mip_status_ == SCIP_STATUS_UNBOUNDED));
      SCIP_CALL( SCIPsetPresolving(scip_, SCIP_PARAMSETTING_DEFAULT, TRUE) );
    }
    cost_vector_ = setPrimalRayAndGetCostVector();
    primal_ray_weights_.push_back(new vector<SCIP_Real>(*weight_));
    assert( cost_vector_ != NULL);
    cost_rays_.push_back(cost_vector_);
    foundNewUnbounded_ = true;
    skeleton_->addPrimalRay(cost_vector_);
  }
  else {
    solving_stage_ = MULTIOPT_SOLVED;
  }

  if (solving_stage_ == MULTIOPT_SOLVING && !skeleton_->hasNextWeight()) 
    solving_stage_ = MULTIOPT_SOLVED;

  return SCIP_OKAY;
}


/** get number of new vertices in the 1-skeleton added in last run */
int LiftedWeightSpaceSolver::getNNewSkeletVerts() const {
  return skeleton_->getNNewVertices();
}

/** get number of vertices in the 1-skeleton processed in last run */
int LiftedWeightSpaceSolver::getNProcSkeletVerts() const {
  return skeleton_->getNProcessedVertices();
}

