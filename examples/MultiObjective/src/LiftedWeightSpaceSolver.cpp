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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   LiftedWeightSpaceSolver.cpp
 * @brief  Main class of the algorithm
 * @author Timo Strunk
 *
 * Realization of a weighted solver using the lifted weight space polyhedron to calculate weights.
 * It gets the weight from the skeleton, then changes the objective in the SCIP problem to the weighted
 * objective.  Then it solves the problem and uses the skeleton to determine if the new solution is
 * an extremal supported nondominated point.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <sstream>
#include <time.h>

#include "scip/def.h"
#include "scip/struct_sol.h"

#include "Skeleton.h"
#include "reader_mop.h"
#include "Objectives.h"
#include "main.h"

#include "LiftedWeightSpaceSolver.h"

/** SCIP style constructor */
LiftedWeightSpaceSolver::LiftedWeightSpaceSolver(
   const char*           paramfilename       /**< name of file with SCIP parameters */
      )
   : WeightedSolver(paramfilename),
     solving_stage_(MULTIOPT_UNSOLVED),
     feasible_weight_lpi_(NULL),
     feasible_weight_sol_(NULL),
     mip_status_(SCIP_STATUS_UNKNOWN)
{
   SCIPcreateClock(scip_, &clock_iteration_);
   SCIPcreateClock(scip_, &clock_total_);

   skeleton_ = new Skeleton(scip_);
   feasible_weight_ = new std::vector<SCIP_Real>(0,0.);
}

/** default destructor */
LiftedWeightSpaceSolver::~LiftedWeightSpaceSolver()
{
   SCIPfreeClock(scip_, &clock_total_);
   SCIPfreeClock(scip_, &clock_iteration_);

   if( feasible_weight_lpi_ != NULL )
   {
      SCIPlpiFree(&feasible_weight_lpi_);
   }

   for( std::vector< const std::vector<SCIP_Real>* >::iterator it = cost_rays_.begin();
        it != cost_rays_.end();
        ++it )
   {
      delete *it;
   }

   if( feasible_weight_sol_ != NULL )
   {
      delete[] feasible_weight_sol_;
      feasible_weight_sol_ = NULL;
   }

   delete feasible_weight_;
   delete skeleton_;
}

/** returns true if there is a weight left to check */
bool LiftedWeightSpaceSolver::hasNext() const
{
   return solving_stage_ != MULTIOPT_SOLVED;
}

/** solve instance with next weight */
SCIP_RETCODE LiftedWeightSpaceSolver::solveNext()
{
   SCIPresetClock(scip_, clock_iteration_);
   SCIPstartClock(scip_, clock_iteration_);

   found_new_optimum_ = false;
   nnodes_last_run_ = 0;
   niterations_last_run_ = 0;
   duration_last_run_ = 0;
   cost_vector_ = NULL;
   solution_ = NULL;

   /* if this is the first iteration, run initialization code */
   if( solving_stage_ == MULTIOPT_UNSOLVED )
   {
      SCIP_CALL( init() );
   }

   /* solve with next weight and process the result */
   SCIP_CALL( loadNextWeight() );
   SCIP_CALL( solveWeighted() );
   SCIP_CALL( evaluateSolution() );

   ++nruns_;

   SCIPstopClock(scip_, clock_iteration_);
   duration_last_run_ = SCIPgetClockTime(scip_, clock_iteration_);

   return SCIP_OKAY;
}

/** prepare to start solving */
SCIP_RETCODE LiftedWeightSpaceSolver::init()
{
   SCIPstartClock(scip_, clock_total_);
   solving_stage_ = MULTIOPT_INIT_WEIGHTSPACE;
   SCIP_CALL( createFeasibleWeightLPI() );
   SCIP_CALL( solveFeasibleWeightLPI() );

   return SCIP_OKAY;
}

/** calculate weight for next weighted optimization run */
SCIP_RETCODE LiftedWeightSpaceSolver::loadNextWeight()
{
   /* if in weight finding mode, find weight by weight lp */
   if( solving_stage_ == MULTIOPT_INIT_WEIGHTSPACE )
   {
      delete feasible_weight_;
      feasible_weight_ = getFeasibleWeight(feasible_weight_sol_);
      weight_ = feasible_weight_;
   }
   /* if first weight has been found, get weight from skeleton */
   else if( solving_stage_ == MULTIOPT_SOLVING )
   {
      weight_ = skeleton_->nextWeight();
   }

   return SCIP_OKAY;
}

/** find the optimal solution for the current weight */
SCIP_RETCODE LiftedWeightSpaceSolver::solveWeighted()
{
   Objectives* objectives = SCIPgetProbData(scip_)->objectives;

   /* load weight into solver */
   SCIP_CALL( objectives->setWeightedObjective(scip_, weight_) );

   /* optimize with weight */
   SCIP_CALL( doSCIPrun() );

   return SCIP_OKAY;
}

/** call the mip solver */
SCIP_RETCODE LiftedWeightSpaceSolver::doSCIPrun()
{
   /* as the solution data is overwritten, it must be cleared beforehand */
   assert(cost_vector_ == NULL);
   assert(solution_ == NULL);

   Objectives* objectives = SCIPgetProbData(scip_)->objectives;

   /* set SCIP timelimit so that total algorithm timelimit is met*/
   SCIPsetRealParam(scip_, "limits/time", timelimit_ - SCIPgetClockTime(scip_, clock_total_));

   /* actual SCIP solver call */
   SCIP_CALL( SCIPsolve(scip_) );

   /* update SCIP run information */
   nnodes_last_run_ += SCIPgetNNodes(scip_);

   if( SCIPgetStage(scip_) != SCIP_STAGE_PRESOLVING )
   {
      niterations_last_run_ += SCIPgetNLPIterations(scip_);
   }
   else
   {
      /* SCIP was interrupted before entering solving stage */
      niterations_last_run_ = 0;
   }

   mip_status_ = SCIPgetStatus(scip_);

   if( mip_status_ == SCIP_STATUS_OPTIMAL )
   {
      SCIP_CALL( copyBestOriginalSolution() );

      cost_vector_ = objectives->calculateCost(scip_, solution_);

      assert( !hasInfiniteComponent(cost_vector_) );
   }
   else if( mip_status_ == SCIP_STATUS_UNBOUNDED )
   {
      cost_vector_ = SCIPgetProbData(scip_)->objectives->calculateCostRay(scip_);
      cost_rays_.push_back(cost_vector_);
   }

   return SCIP_OKAY;
}

/** get total time for algorithm */
SCIP_Real LiftedWeightSpaceSolver::getTotalDuration() const
{
   return SCIPgetClockTime(scip_, clock_total_);
}

/** return true if the given vector has an entry close to infinity */
bool LiftedWeightSpaceSolver::hasInfiniteComponent(const std::vector<SCIP_Real>* cost_vector)
{
   bool result = false;

   for( std::vector<SCIP_Real>::const_iterator it = cost_vector_->begin();
        it != cost_vector_->end();
        ++it )
   {
      if( *it >= SCIPinfinity(scip_) / 1000. )
      {
         result = true;
         break;
      }
   }

   return result;
}

/** get the MIP solution and check wheather it is a new optimum*/
SCIP_RETCODE LiftedWeightSpaceSolver::evaluateSolution()
{
   if( mip_status_ == SCIP_STATUS_OPTIMAL )
   {
      if( solving_stage_ == MULTIOPT_INIT_WEIGHTSPACE )
      {
         found_new_optimum_ = true;
         skeleton_->init(cost_vector_, &cost_rays_);
         solving_stage_ = MULTIOPT_SOLVING;
      }
      else
      {
         found_new_optimum_ = skeleton_->isExtremal(cost_vector_);
      }

      if( found_new_optimum_ )
      {
         nondom_points_.push_back(cost_vector_);
         solutions_.push_back(solution_);
         cost_to_sol_[cost_vector_] = solution_;
      }
      else
      {
         delete cost_vector_;
         cost_vector_ = NULL;
         SCIP_CALL( SCIPfreeSol(scip_, &solution_) );
         solution_ = NULL;
      }
   }
   else if( mip_status_ == SCIP_STATUS_UNBOUNDED )
   {
      if( solving_stage_ == MULTIOPT_INIT_WEIGHTSPACE )
      {
         SCIP_CALL( updateFeasibleWeightLPI(cost_vector_) );
         SCIP_CALL( solveFeasibleWeightLPI() );
      }
      else
      {
         skeleton_->addPrimalRay(cost_vector_);
      }
   }

   if( multiopt_status_ == SCIP_STATUS_UNBOUNDED )
   {
      solving_stage_ = MULTIOPT_SOLVED;
   }
   else if( solving_stage_ == MULTIOPT_SOLVING && !skeleton_->hasNextWeight() )
   {
      solving_stage_ = MULTIOPT_SOLVED;
      multiopt_status_ = SCIP_STATUS_OPTIMAL;
   }
   else if(    mip_status_ != SCIP_STATUS_OPTIMAL
            && mip_status_ != SCIP_STATUS_UNBOUNDED )
   {
      solving_stage_ = MULTIOPT_SOLVED;
      multiopt_status_ = mip_status_;
   }

   return SCIP_OKAY;
}

/** get number of new vertices in the 1-skeleton added in last step*/
int LiftedWeightSpaceSolver::getNNewVertices() const
{
   return skeleton_->getNNewVertices();
}

/** get number of vertices in the 1-skeleton processed in last step*/
int LiftedWeightSpaceSolver::getNProcessedVertices() const
{
   return skeleton_->getNProcessedVertices();
}

/** initialize lp for feasible weight generation */
SCIP_RETCODE LiftedWeightSpaceSolver::createFeasibleWeightLPI()
{
   assert( feasible_weight_lpi_ == NULL );

   int nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();

   int ncols = nobjs + 1;
   SCIP_Real* obj = new SCIP_Real[ncols];
   SCIP_Real* lb = new SCIP_Real[ncols];
   SCIP_Real* ub = new SCIP_Real[ncols];

   int nrows = nobjs + 1;
   SCIP_Real* lhs = new SCIP_Real[nrows];
   SCIP_Real* rhs = new SCIP_Real[nrows];
   int nnonz = 3 * nobjs;

   int* beg = new int[ncols];
   int* ind = new int[nnonz];
   SCIP_Real* val = new SCIP_Real[nnonz];

   /** init weight columns */
   for( int j = 0; j < nobjs; ++j )
   {
      obj[j] = 0.;
      lb[j]  = 0.;
      ub[j]  = SCIPinfinity(scip_);
      beg[j] = 2 * j;
      ind[2 * j] = j;
      val[2 * j] = 1.;
      ind[2 * j + 1] = nobjs;
      val[2 * j + 1] = 1.;
   }

   /** init slack column */
   obj[nobjs] = 1.;
   lb[nobjs]  = - SCIPinfinity(scip_);
   ub[nobjs]  = SCIPinfinity(scip_);
   beg[nobjs] = 2 * nobjs;
   for( int i = 0; i < nobjs; ++i )
   {
      ind[2 * nobjs + i] = i;
      val[2 * nobjs + i] = -1.;
   }

   /** init slack rows */
   for( int i = 0; i < nobjs; ++i )
   {
      lhs[i] = 0.;
      rhs[i] = SCIPinfinity(scip_);
   }

   /** init normalized row */
   lhs[nobjs] = 1.;
   rhs[nobjs] = 1.;

   SCIP_CALL( SCIPlpiCreate(
         &feasible_weight_lpi_,
         NULL,
         "feasible weight",
         SCIP_OBJSEN_MAXIMIZE
         ) );

   SCIP_CALL(  SCIPlpiLoadColLP(
         feasible_weight_lpi_,
         SCIP_OBJSEN_MAXIMIZE,
         ncols,
         obj,
         lb,
         ub,
         NULL,
         nrows,
         lhs,
         rhs,
         NULL,
         nnonz,
         beg,
         ind,
         val
         ) );

   delete[] obj;
   delete[] lb;
   delete[] ub;

   delete[] lhs;
   delete[] rhs;

   delete[] beg;
   delete[] ind;
   delete[] val;

   return SCIP_OKAY;
}

/** solve feasible weight lp to get next feasible weight candidate */
SCIP_RETCODE LiftedWeightSpaceSolver::solveFeasibleWeightLPI()
{
   assert( feasible_weight_lpi_ != NULL );

   SCIP_CALL( SCIPlpiSolvePrimal(feasible_weight_lpi_) );

   if( SCIPlpiIsPrimalFeasible(feasible_weight_lpi_) )
   {
      int nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();
      if( feasible_weight_sol_ != NULL )
      {
         delete[] feasible_weight_sol_;
         feasible_weight_sol_ = NULL;
      }

      feasible_weight_sol_ = new SCIP_Real[nobjs + 1];

      SCIP_CALL( SCIPlpiGetSol(
          feasible_weight_lpi_,
          NULL,
          feasible_weight_sol_,
          NULL,
          NULL,
          NULL
          ) );
   }
   else
   {
      multiopt_status_ = SCIP_STATUS_UNBOUNDED;
   }

   return SCIP_OKAY;
}

/** copy feasible weight lp solution to vector */
const std::vector<SCIP_Real>* LiftedWeightSpaceSolver::getFeasibleWeight(SCIP_Real* sol)
{
   int nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();
   std::vector<SCIP_Real>* result = new std::vector<SCIP_Real>(nobjs, 0.);

   for( int i = 0; i < nobjs; ++i )
   {
     (*result)[i] = sol[i];
   }

   return result;
}

/** add new cost ray constraint to feasible weight lp */
SCIP_RETCODE LiftedWeightSpaceSolver::updateFeasibleWeightLPI(const std::vector<SCIP_Real>* cost_ray)
{
   int nobjs = SCIPgetProbData(scip_)->objectives->getNObjs();

   SCIP_Real lhs = 0;
   SCIP_Real rhs = SCIPinfinity(scip_);


   int nnonz = nobjs + 1;
   int beg = 0;
   int* ind = new int[nnonz];
   SCIP_Real* val = new SCIP_Real[nnonz];

   for( int i = 0; i < nobjs; ++i )
   {
      ind[i] = i;
      val[i] = (*cost_ray)[i];
   }

   ind[nobjs] = nobjs;
   val[nobjs] = - 1;

   SCIP_CALL( SCIPlpiAddRows(
      feasible_weight_lpi_,
      1,
      &lhs,
      &rhs,
      NULL,
      nnonz,
      &beg,
      ind,
      val
      ) );

   delete[] ind;
   delete[] val;

   return SCIP_OKAY;
}

/** copies best scip solution into sol */
SCIP_RETCODE LiftedWeightSpaceSolver::copyBestOriginalSolution()
{
   assert( solution_ == NULL );

   SCIP_SOL* sol = SCIPgetBestSol(scip_);

   SCIP_SOL* finite = NULL;
   SCIP_Bool copy_success = FALSE;

   SCIP_CALL( SCIPcreateFiniteSolCopy(scip_, &finite, sol, &copy_success) );

   assert( copy_success );

   int        nvars = SCIPgetNOrigVars(scip_);
   SCIP_VAR** vars  = SCIPgetOrigVars(scip_);
   SCIP_Real* vals  = new SCIP_Real[nvars];

   SCIP_CALL( SCIPgetSolVals(scip_, finite, nvars, vars, vals) );

   SCIP_CALL( SCIPfreeSol(scip_, &finite) );
   SCIP_CALL( SCIPfreeTransform(scip_) );

   SCIP_CALL( SCIPcreateOrigSol(scip_, &solution_, NULL) );

   for( int i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip_, solution_, vars[i], vals[i]) );
   }

   delete[] vals;

   return SCIP_OKAY;
}

/** true if the last solved weighted problem is unbounded */
bool LiftedWeightSpaceSolver::isWeightedUnbounded() const
{
   return mip_status_ == SCIP_STATUS_UNBOUNDED;
}
