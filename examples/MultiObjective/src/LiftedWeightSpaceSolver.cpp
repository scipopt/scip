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

/**@file   LiftedWeightSpaceSolver.cpp
 * @brief  Main class of the algorithm
 * @author Timo Strunk
 *
 * @desc   Realization of a weighted solver using the lifted weight space polyhedron to calculate weights.
 *         It gets the weight from the skeleton, then changes the objective in the SCIP problem to the weighted
 *         objective.  Then it solves the problem and uses the skeleton to determine if the new solution is
 *         an extremal supported nondominated point.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <sstream>
#include <time.h>

#include "scip/def.h"
#include "scip/struct_sol.h"

#include "Skeleton.h"
#include "reader_mop.h"
#include "Objectives.h"

#include "LiftedWeightSpaceSolver.h"

/** standard constructor */
LiftedWeightSpaceSolver::LiftedWeightSpaceSolver(
   bool                  verbose,            /**< true if scip output should be displayed */
   SCIP_Real             timelimit,          /**< maximum allowed time in seconds */
   int                   solstore            /**< number of solutions stored in SCIP */
      )
   : WeightedSolver(verbose, timelimit, solstore),
     first_weight_(NULL)
{
   SCIPcreateClock(scip_, &clock_iteration_);
   SCIPcreateClock(scip_, &clock_total_);

   skeleton_ = new Skeleton(scip_);
}

/** default destructor */
LiftedWeightSpaceSolver::~LiftedWeightSpaceSolver()
{
   SCIPfreeClock(scip_, &clock_total_);
   SCIPfreeClock(scip_, &clock_iteration_);

   if( first_weight_ != NULL )
   {
      delete first_weight_;
   }

   delete skeleton_;
}

/** returns true if there is a weight left to check */
bool LiftedWeightSpaceSolver::hasNext() const
{
   /* if the last iteration did not yield a solution, the algorithm shall terminate */
   return  (nondom_points_.empty() || skeleton_->hasNextWeight())
      && (status_ == SCIP_STATUS_OPTIMAL || status_ == SCIP_STATUS_UNKNOWN) ;
}

/** load next weight into solver */
SCIP_RETCODE LiftedWeightSpaceSolver::next()
{
   Objectives*           objectives;
   int                   nobjs;

   objectives = SCIPgetProbData(scip_)->objectives;
   nobjs = objectives->getNObjs();

   /* start the clock
    * we measure the time for the iteration from SCIP solving to evaluating the solution */
   SCIPresetClock(scip_, clock_iteration_);
   SCIPstartClock(scip_, clock_iteration_);

   if( nondom_points_.empty() )
   {
      /* this is the first iteration  */
      SCIPstartClock(scip_, clock_total_);
      first_weight_ = new std::vector<SCIP_Real>(nobjs, 1. / nobjs);
      weight_ = first_weight_;
   }
   else
   {
      /* after the first iteration */
      weight_ = skeleton_->nextWeight();
      if( SCIPisTransformed(scip_) )
      {
        SCIP_CALL( SCIPfreeTransform(scip_) );
      }
   }

   SCIP_CALL( objectives->setWeightedObjective(scip_, weight_) );

   found_new_optimum_ = false;
   status_ = SCIP_STATUS_UNKNOWN;

   return SCIP_OKAY;
}

/** solve instance with next weight */
SCIP_RETCODE LiftedWeightSpaceSolver::solve()
{
   /* set SCIP timelimit so that total algorithm timelimit is met*/
   SCIPsetRealParam(scip_, "limits/time", timelimit_ - SCIPgetClockTime(scip_, clock_total_));

   /* actual SCIP solver call */
   SCIP_CALL( SCIPsolve( scip_ ) );
   status_ = SCIPgetStatus(scip_);

   /* process the result of the SCIP run */
   if( status_ == SCIP_STATUS_OPTIMAL )
   {
      evaluateSolution();
   }

   /* update SCIP run information */
   ++nruns_;
   nnodes_last_run_ = SCIPgetNNodes(scip_);

   if( SCIPgetStage(scip_) != SCIP_STAGE_PRESOLVING )
   {
      niterations_last_run_ = SCIPgetNLPIterations(scip_);
   }
   else
   {
      /* SCIP was interrupted before entering solving stage */
      niterations_last_run_ = 0;
   }
   /* stop the clock */
   SCIPstopClock(scip_, clock_iteration_);
   duration_last_run_ = SCIPgetClockTime(scip_, clock_iteration_);

   return SCIP_OKAY;
}

/** get total time for algorithm */
SCIP_Real LiftedWeightSpaceSolver::getTotalDuration() const
{
   return SCIPgetClockTime(scip_, clock_total_);
}

/** get the MIP solution and check wheather it is a new optimum*/
void LiftedWeightSpaceSolver::evaluateSolution()
{
   solution_ = SCIPgetBestSol(scip_);
   assert(solution_ != NULL);
   assert(solution_->vals != NULL);

   cost_vector_ = SCIPgetProbData(scip_)->objectives->calculateCost(scip_, solution_);
   found_new_optimum_ = skeleton_->checkSolution(cost_vector_);

   if( found_new_optimum_ )
   {
      nondom_points_.push_back(cost_vector_);
   }
   else
   {
      delete cost_vector_;
   }
}
