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

/**@file   WeightedSolver.h
 * @brief  Class providing all methods for using the algorithm
 * @author Sebastian Schenker, Timo Strunk
 *
 * Abstract superclass for a multi objective solver using weighted objective functions.
 */

#ifndef WEIGHTED_SOLVER
#define WEIGHTED_SOLVER

#include "scip/scip.h"
#include "scip/message_default.h"

#include <iosfwd>
#include <string>
#include <map>
#include <unordered_set>
#include <vector>

using std::vector;

/** generic weight based solver */
class WeightedSolver
{
 public:
  using unor_Set = std::unordered_set<unsigned>;

  /** SCIP style constructor */
  WeightedSolver();

  /** destructor */
  virtual ~WeightedSolver();
  
   /** reads problem data from file */
  SCIP_RETCODE readProblem(const std::string& file, bool beVerbose=false);
  
  /** reads parameter settings file */
  SCIP_RETCODE readParamSettings(const std::string& file);

  /** returns true if there is a weight left to check */
  virtual bool hasUncheckedWeight() const = 0;
  
  /** solves weighted problem */
  virtual SCIP_RETCODE solveUncheckedWeight(bool beVerbose) = 0;
  
  /** get total time for algorithm */
  virtual SCIP_Real getTotalDuration() const=0;
  
  /** returns true if the last weighted run found a new pareto optimum */
  virtual bool foundNewOptimum() const;

  /** returns true if the last weighted run found an unbounded ray */
  virtual bool foundNewUnbounded() const;

  /** returns the number of objective functions */
  unsigned int getNObjs() const;

  /** returns the number of weighted runs so far */
  unsigned long getNRuns() const;

  /** returns the number of efficient solutions so far */
  virtual unsigned long getNSolutions() const;

  /** prints the unbounded ray information to given output stream */
  virtual void printRayInfo(std::ostream& os, bool onlyLastAdded, bool withWeightsconst,
			    const std::string& prelim = "", unor_Set weaklyEff = unor_Set());

  /** prints the soluton information to given output stream */
  virtual void printSolInfo(std::ostream& os, bool onlyLastAdded, bool withWeightsconst,
			    const std::string& prelim = "", unor_Set weaklyEff = unor_Set());

  /** prints the name of given objective to given output stream */
  virtual void printObjName(std::ostream& os, unsigned objNo);

  /** delete weakly non-dominated points and corresp. solutions and weights */
  void deleteWeaklyNondomPointsAndCorrespSols();

 protected:

  SCIP* scip_;                /**< SCIP solver */

  SCIP_Objsense objSense_;   /**< objective sense of given problem */
  
  unsigned nObjs_;            /**< number of objectives */

  unsigned long nRuns_;       /**< number of solveNext() calls */

  bool foundNewOptimum_;      /**< true if last call of solveNext() found a new optimum */

  bool foundNewUnbounded_;    /**< true if last call of solveNext() found an unbounded ray*/

  SCIP_SOL* solution_;        /**< last solution found by solveNext() */

  const vector<SCIP_Real>* weight_;            /**< weight used in last solveNext() call */

  const vector<SCIP_Real>* cost_vector_;       /**< objective vector of last found solution */
  
  vector< SCIP_SOL* > solutions_;              /**< (weakly) efficient SCIP solutions found so far */

  vector<const vector< SCIP_Real>* > points_;  /**< (weakly) non dominated points found so far */

  vector<const vector<SCIP_Real>* > sol_weights_;  /**< weights corresp. to computed solutions */

  vector< const vector<SCIP_Real>* > cost_rays_;   /**< objective vectors of unbounded primal rays*/

  vector< const std::map<std::string,SCIP_Real>* > primal_rays_; /**< found primal rays */

  vector< const vector<SCIP_Real>* > primal_ray_weights_;  /**< weights corres. to primal rays */

  /** print given vector to given output stream */
  void printVec(std::ostream& os, const std::string& str, const vector<SCIP_Real>* vals, 
		bool negate = false, unsigned precision = 3);
  
 private:
  /** return true if p is dominated, false otherwise */
  bool isDominated(std::vector<const std::vector<SCIP_Real>* >::const_iterator p);

};

#endif
