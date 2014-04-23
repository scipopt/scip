/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
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
 * @desc   Realization of a weighted solver using the lifted weight space polyhedron to calculate weights.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef LIFTED_WEIGHT_SPACE_SOLVER
#define LIFTED_WEIGHT_SPACE_SOLVER

#include <vector>
#include <string>

#include "WeightedSolver.h"

class Skeleton;

/** Class implementing the weight space based solving algorithm for multi objective
    optimization */
class LiftedWeightSpaceSolver : public WeightedSolver
{
 public:
   /** standard constructor */
   LiftedWeightSpaceSolver(
      bool               verbose,            /**< true if scip output should be displayed */
      SCIP_Real          timelimit,          /**< maximum allowed time in seconds */
      int                solstore            /**< number of solutions stored in SCIP */
      );

   /** default destructor */
   virtual ~LiftedWeightSpaceSolver();

   /** returns true if there is a weight left to check */
   bool hasNext() const;

   /** load next weight into solver */
   SCIP_RETCODE next();

   /** solve instance with next weight */
   SCIP_RETCODE solve();

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

   const std::vector<SCIP_Real>* first_weight_; /**< first weight used before skeleton is properly initialized */

   /** get the MIP solution and check wheather it is a new optimum*/
   void evaluateSolution();
};

#endif
