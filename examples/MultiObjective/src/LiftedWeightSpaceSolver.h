/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   LiftedWeightSpaceSolver.h
 * @brief  Main class of the algorithm
 * @author Timo Strunk
 *
 * Realization of a weighted solver using the lifted weight space polyhedron to calculate weights.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef LIFTED_WEIGHT_SPACE_SOLVER
#define LIFTED_WEIGHT_SPACE_SOLVER

#include <vector>
#include <string>

#include "lpi/lpi.h"

#include "WeightedSolver.h"

class Skeleton;

/** Class implementing the weight space based solving algorithm for multi objective
    optimization */
class LiftedWeightSpaceSolver : public WeightedSolver
{
 public:
   /** SCIP style constructor */
   LiftedWeightSpaceSolver(
      const char*        paramfilename       /**< name of file with SCIP parameters */
      );

   /** default destructor */
   virtual ~LiftedWeightSpaceSolver();

   /** returns true if there is a weight left to check */
   bool hasNext() const;

   /** solve instance with next weight */
   SCIP_RETCODE solveNext();

   /** true if the last solved weighted problem is unbounded */
   bool isWeightedUnbounded() const;

   /** get total time for algorithm */
   SCIP_Real getTotalDuration() const;

   /** get number of new vertices in the 1-skeleton added in last step*/
   int getNNewVertices() const;

   /** get number of vertices in the 1-skeleton processed in last step*/
   int getNProcessedVertices() const;

 private:
   Skeleton*             skeleton_;             /**< lifted weight space polyhedron data structure */
   SCIP_CLOCK*           clock_iteration_;      /**< clock measuring the time needed for every iteration */
   SCIP_CLOCK*           clock_total_;          /**< clock measuring the time needed for the entire program */

   /** phase of the weighted solving process */
   enum Multiopt_Stage {
      MULTIOPT_UNSOLVED,                        /**< solving process has not yet started */
      MULTIOPT_INIT_WEIGHTSPACE,                /**< finding first valid solution candidates */
      MULTIOPT_SOLVING,                         /**< finding more solution candidates */
      MULTIOPT_SOLVED                           /**< found all solution candidates */
   }                     solving_stage_;        /**< phase the algorithm is currently in */

   SCIP_LPI*             feasible_weight_lpi_;  /**< LP instance for finding feasible weight */
   SCIP_Real*            feasible_weight_sol_;  /**< solution of feasible weight LP */

   const std::vector<SCIP_Real>*                feasible_weight_; /**< first weight found with bounded solutions */

   SCIP_Status           mip_status_;           /**< status of weighted problem */

   /** prepare to start solving */
   SCIP_RETCODE init();

   /** calculate weight for next weighted optimization run */
   SCIP_RETCODE loadNextWeight();

   /** find the optimal solution for the current weight */
   SCIP_RETCODE solveWeighted();

   /** call the mip solver */
   SCIP_RETCODE doSCIPrun();

   /** get the MIP solution and check wheather it is a new optimum*/
   SCIP_RETCODE evaluateSolution();

   /** initialize lp for feasible weight generation */
   SCIP_RETCODE createFeasibleWeightLPI();

   /** solve feasible weight lp to get next feasible weight candidate */
   SCIP_RETCODE solveFeasibleWeightLPI();

   /** copy feasible weight lp solution to vector */
   const std::vector<SCIP_Real>* getFeasibleWeight(SCIP_Real* sol);

   /** add new cost ray constraint to feasible weight lp */
   SCIP_RETCODE updateFeasibleWeightLPI(const std::vector<SCIP_Real>* cost_ray);

   /** return true if the given vector has an entry close to infinity */
   bool hasInfiniteComponent(const std::vector<SCIP_Real>* cost_vector);

   /** copies best scip solution into sol */
   SCIP_RETCODE copyBestOriginalSolution();
};

#endif
