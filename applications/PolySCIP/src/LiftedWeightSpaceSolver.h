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

/**@file   LiftedWeightSpaceSolver.h
 * @brief  Main class of the algorithm
 * @author Sebastian Schenker, Timo Strunk
 *
 * Realization of a weighted solver using the lifted weight space polyhedron to calculate weights.
 */

#ifndef LIFTED_WEIGHT_SPACE_SOLVER
#define LIFTED_WEIGHT_SPACE_SOLVER


#include "lpi/lpi.h"
#include "WeightedSolver.h"

#include <iosfwd>
#include <vector>
#include <string>

class Skeleton;

/** Class implementing the weight space based solving algorithm for multi objective
    optimization */
class LiftedWeightSpaceSolver : public WeightedSolver
{
 public:
  /** SCIP style constructor */
  LiftedWeightSpaceSolver();
  
  /** default destructor */
  virtual ~LiftedWeightSpaceSolver();

  /** set time limit for entire computation */
  void setTimeLimit(unsigned long timeLimit);
  
  /** returns true if there is an unchecked weight */
  bool hasUncheckedWeight() const;
  
  /** solve instance with next weight */
  SCIP_RETCODE solveUncheckedWeight(bool beVerbose);
  
  /** find first non-dominated point and initialise skeleton or find out that there is are
      no non-dominated points (only non-dominated rays) **/
  SCIP_RETCODE initialise(bool beVerbose);
  
  /** returns the number of branching nodes used in the last run */  
  SCIP_Longint getNNodesLastRun() const;

  /** returns the number of LP iterations done in the last run */
  SCIP_Longint getNLPIterationsLastRun() const ;

  /** returns the compuational time taken in the last run */
  SCIP_Real getDurationLastRun() const;
  
  /** get total time for algorithm */
  SCIP_Real getTotalDuration() const;

  /** get number of new vertices in the 1-skeleton added in last step*/
  int getNNewSkeletVerts() const;

  /** get number of vertices in the 1-skeleton processed in last step*/
  int getNProcSkeletVerts() const;

 private:
   /** stage of the multi-criteria solving process */
  enum Multiopt_Stage {
    MULTIOPT_UNSOLVED,           /**< solving process has not yet started */
    MULTIOPT_SOLVING,            /**< finding more solution candidates */
    MULTIOPT_SOLVED,             /**< found all solution candidates */
  } solving_stage_;               /**< stage the algorithm is currently in */
  
  bool timeLimit_;                /**< true if time limit was set by user */
  unsigned long timeLimitVal_;    /**< value of given time limit in seconds */
  
  SCIP_Status mip_status_;        /**< status after SCIP optimisaton  */
  
  std::vector<SCIP_Real>* initialWeight_;   /**< weight used to find first (weakly non-dom) 
					       point and initalize weight space */
  
  SCIP_Longint nNodesLastRun_;        /**< number of branch and bound nodes in last 
					 call of solveNext() */

  SCIP_Longint nIterationsLastRun_;   /**< number of lp iterations in last call of solveNext() */
  
  SCIP_Real durationLastRun_;       /**< duration of last call of solveNext() in seconds */
  
  SCIP_CLOCK* clock_iteration_;      /**< clock measuring the time needed for every iteration */
  SCIP_CLOCK* clock_total_;          /**< clock measuring the time needed for the entire program */
  
  Skeleton* skeleton_;             /**< lifted weight space polyhedron data structure */
  
  /** change objective function for next weighted optimization run */
  SCIP_RETCODE loadWeight(bool initPhase);
  
  /** find the optimal solution for the current weight */
  SCIP_RETCODE solveWeighted();
  
  /** get the MIP solution and check wheather it is a new optimum*/
  SCIP_RETCODE evalSolution();
  
  /** adds solution, correps. weight and objective vector to internal data structures */
  void addSolAndObjValsAndWeight();
  
  /** prints weight, objective vector and solution to standard output */
  void printVerboseInfo(bool withSkeletonInfo=false);
  
  /** return true if the given vector has an entry close to infinity */
  bool hasInfiniteComponent(const std::vector<SCIP_Real>* cost_vector);
  
  /** sets found solution and returns corresp. objective vector */
  std::vector<SCIP_Real>* setSolAndGetCostVector();
  
  /** sets found ray and returns corresp. objective vector */
  std::vector<SCIP_Real>* setPrimalRayAndGetCostVector();
  
};

#endif
