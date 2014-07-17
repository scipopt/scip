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

/**@file   WeightedSolver.cpp
 * @brief  Class providing all methods for using the algorithm
 * @author Timo Strunk
 * 
 * @desc   Abstract superclass for a multi objective solver using weighted objective functions.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "cstring"
#include <iostream>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "sys/stat.h"
#include "reader_mop.h"
#include "lpi/lpi.h"

#include "WeightedSolver.h"
#include "Objectives.h"
#include "main.h"

/** SCIP style constructor */
WeightedSolver::WeightedSolver(
   const char*           paramfilename       /**< name of file with SCIP parameters */     
   )
   : scip_(NULL),
     multiopt_status_(SCIP_STATUS_UNKNOWN),
     nruns_(0),
     found_new_optimum_(false),
     solution_(NULL),
     weight_(NULL),
     cost_vector_(NULL),
     n_written_sols_(0)
{
   SCIPcreate(&scip_);
   assert(scip_ != NULL);
   SCIPincludeDefaultPlugins(scip_);
   SCIPincludeReaderMop(scip_);

   if( paramfilename != NULL )
   {
      SCIPreadParams(scip_, paramfilename);
   }

   SCIPgetRealParam(scip_, "limits/time", &timelimit_);
}

/** destructor */
WeightedSolver::~WeightedSolver()
{
   if( SCIPisTransformed(scip_) )
   {
      SCIPfreeTransform(scip_);
   }

   for( std::vector< SCIP_SOL* >::iterator it = solutions_.begin();
        it != solutions_.end();
        ++it )
   {
      SCIPfreeSol( scip_, &*it);
   }

   for( std::vector< const std::vector< SCIP_Real>* >::iterator it = nondom_points_.begin();
        it != nondom_points_.end();
        ++it )
   {
      delete *it;
   }

   delete SCIPgetProbData(scip_)->objectives;

   SCIPfree(&scip_);
}

/** reads problem data from file */
SCIP_RETCODE WeightedSolver::readProblem(
   const char*           filename            /**< name of instance file */
   )
{
   unsigned              startname;

   filename_ = std::string(filename);
   startname = filename_.find_last_of("/\\");

   /* prepare path for output files */
   outfilestump_ = "solutions"+filename_.substr(
      startname,
      filename_.length() - 4 - startname
      );

   SCIP_CALL( SCIPreadProb(scip_, filename, NULL) );

   return SCIP_OKAY;
}

/** returns the SCIP problem status of the multiobjective problem */
SCIP_Status WeightedSolver::getStatus() const
{
   return multiopt_status_;
}

/** returns true if the last run lead to a new pareto optimal solution */
bool WeightedSolver::foundNewOptimum() const
{
   return found_new_optimum_;
}

/** returns the last weight used in a weighted optimization run */
const  std::vector<SCIP_Real>*  WeightedSolver::getWeight() const
{
   return weight_;
}

/** returns the cost vector of the last found pareto optimum */
const  std::vector<SCIP_Real>* WeightedSolver::getCost() const
{
   return cost_vector_;
}

/** returns the last found pareto optimal solution */
SCIP_SOL* WeightedSolver::getSolution() const
{
   return solution_;
}

/** returns the number of branching nodes used in the last run */
SCIP_Longint WeightedSolver::getNNodesLastRun() const
{
   return nnodes_last_run_;
}

/** returns the number of LP iterations used in the last run */
SCIP_Longint WeightedSolver::getNLPIterationsLastRun() const
{
   return niterations_last_run_;
}

/** returns the number of all found pareto optima */
int WeightedSolver::getNSolutions() const
{
   return nondom_points_.size();
}

/** returns the name of the file containing the last written solution */
std::string WeightedSolver::getSolutionFileName() const
{
   return std::string(solution_file_name_, 10);
}

/** returns the time needed for the last iteration in seconds */
SCIP_Real WeightedSolver::getDurationLastRun() const
{
  return duration_last_run_;
}

/** returns the number of objective functions */
int WeightedSolver::getNObjs() const
{
   return SCIPgetProbData(scip_)->objectives->getNObjs();
}

/** returns the number of weighted runs so far */
int WeightedSolver::getNRuns() const
{
   return nruns_;
}

/** writes a solution to a file in folder solutions/ with name <instance name>-<solution number>.sol*/
SCIP_RETCODE WeightedSolver::writeSolution(SCIP_SOL* sol)
{
   std::stringstream     s_outfile;
   FILE*                 fsol;

   /* build and save file name */
   s_outfile << outfilestump_
             << "-"
             << ++n_written_sols_
             << ".sol";
   solution_file_name_ = s_outfile.str();

   /* write file */
   fsol = fopen(solution_file_name_.c_str(), "w");
   SCIP_CALL( SCIPprintSol(scip_, sol, fsol, FALSE) );

   return SCIP_OKAY;
}

/** delete non extremal solutions 
 *  for each nondominated point p we try to find a point 
 *  q <= p which is a convex combination of the other nondominated points */
SCIP_RETCODE WeightedSolver::enforceExtremality()
{
   SCIP_CALL( createExtremalityLP() );

   Objectives* objectives = SCIPgetProbData(scip_)->objectives;
   const std::vector<std::string>* objnames = objectives->getObjNames();

   for( std::vector<std::string>::const_iterator it = objnames->begin();
        it != objnames->end();
        ++it )
   {
      std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
   }

   std::cout << "   solution file" << std::endl;

   std::vector< const std::vector<SCIP_Real>* >::iterator it = nondom_points_.begin();

   int npoints = nondom_points_.size();
   for( int point_index = 0; point_index < npoints; ++point_index )
   {
      cost_vector_ = *it;
      SCIP_CALL( solveExtremalityLP(cost_vector_, point_index) );

      if( candidate_is_extremal_ )
      {
         writeSolution(nondom_point_to_sol_[cost_vector_]);

         for( std::vector<SCIP_Real>::const_iterator it = cost_vector_->begin();
              it != cost_vector_->end();
              ++it )
         {
            std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
         }
      
         std::cout << "   " << getSolutionFileName() << std::endl;
         ++it;
      }
      else
      {
         it = nondom_points_.erase(it);
         /* the iterator now points to the next candidate */
      }
   }

   SCIP_CALL( SCIPlpiFree(&extremality_lpi_) );

   return SCIP_OKAY;
}

/** print every unbounded cost ray */
void WeightedSolver::printUnboundedRays()
{
   if( !cost_rays_.empty() )
   {
      std::cout << "unbounded rays" << std::endl;
      Objectives* objectives = SCIPgetProbData(scip_)->objectives;
      const std::vector<std::string>* objnames = objectives->getObjNames();

      for( std::vector<std::string>::const_iterator it = objnames->begin();
           it != objnames->end();
           ++it )
      {
         std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
      }

      std::cout << std::endl;

      for( std::vector< const std::vector<SCIP_Real>* >::iterator 
              it = cost_rays_.begin();
           it != cost_rays_.end();
           ++it
         )
      {
         cost_vector_ = *it;
         for( std::vector<SCIP_Real>::const_iterator it = cost_vector_->begin();
              it != cost_vector_->end();
              ++it )
         {
            std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
         }
      
         std::cout << std::endl;
      }
   }
}
/** prepare the LP for the extremality check
 *  For a given point p* it has the form 
 *  x_1 * p_1 + ... + x_n * p_n <= p* 
 *  x_1 + ... + x_n = 1
 *  x_i >= 0 for all i = 1,...,n */
SCIP_RETCODE WeightedSolver::createExtremalityLP()
{
   int ncols = nondom_points_.size();
   int nrows = getNObjs() + 1; 
   int nnonz = ncols * nrows;

   SCIP_Real* obj = new SCIP_Real[ncols];
   SCIP_Real* lb  = new SCIP_Real[ncols];
   SCIP_Real* ub  = new SCIP_Real[ncols];

   SCIP_Real* lhs = new SCIP_Real[nrows];
   SCIP_Real* rhs = new SCIP_Real[nrows];

   int*       beg = new int[ncols];
   int*       ind = new int[nnonz];
   SCIP_Real* val = new SCIP_Real[nnonz];

   for( int i = 0; i < nrows - 1; ++i )
   {
      lhs[i] = - SCIPinfinity(scip_);
      rhs[i] = 0.;
   } 

   lhs[nrows - 1] = 1.;
   rhs[nrows - 1] = 1.;

   for( int j = 0; j < ncols; ++j )
   {
      obj[j] = 0.;
      lb[j]  = 0.;
      ub[j]  = 1.;
      beg[j] = nrows * j;

      for( int i = 0; i < nrows - 1; ++i )
      {
         ind[nrows * j + i] = i;
         val[nrows * j + i] = nondom_points_[j]->at(i);
      }

      ind[nrows * (j + 1) - 1] = nrows - 1;
      val[nrows * (j + 1) - 1] = 1.;
   }

   SCIP_CALL( SCIPlpiCreate(&extremality_lpi_, NULL, "calculate convex combination", SCIP_OBJSEN_MINIMIZE) );            

   SCIP_CALL( SCIPlpiLoadColLP(
      extremality_lpi_,
      SCIP_OBJSEN_MINIMIZE,
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

   delete obj;
   delete lb;
   delete ub;
   delete lhs;
   delete rhs;
   delete beg;
   delete ind;
   delete val;

   return SCIP_OKAY;
}

/** solve the extremality lp for one particular nondom point */
SCIP_RETCODE WeightedSolver::solveExtremalityLP(const std::vector<SCIP_Real>* nondom_point, int point_index)
{
   int nrows = getNObjs() + 1; 

   SCIP_Real obj_one = 1.;
   SCIP_Real obj_zero = 0.;
   SCIP_Real objval;

   /* change lp according to nondom point */
   int*       ind = new int[nrows];
   SCIP_Real* lhs = new SCIP_Real[nrows];
   SCIP_Real* rhs = new SCIP_Real[nrows];

   for( int i = 0; i < nrows - 1; ++i )
   {
      lhs[i] = - SCIPinfinity(scip_);
      rhs[i] = nondom_point->at(i);
      ind[i] = i;
   }

   SCIP_CALL( SCIPlpiChgObj(extremality_lpi_, 1, &point_index, &obj_one) );
   SCIP_CALL( SCIPlpiChgSides(extremality_lpi_, nrows - 1, ind, lhs, rhs) );

   /* solve the lp */
   SCIP_CALL( SCIPlpiSolvePrimal(extremality_lpi_) );
   assert( SCIPlpiIsOptimal(extremality_lpi_) );

   /* evaluate the solution */
   SCIP_CALL( SCIPlpiGetObjval(extremality_lpi_, &objval) );
   candidate_is_extremal_ = SCIPisGT(scip_, objval, 0.);

   /* clean up */
   SCIP_CALL( SCIPlpiChgObj(extremality_lpi_, 1, &point_index, &obj_zero) );		
   delete ind;
   delete lhs;
   delete rhs;
   
   return SCIP_OKAY;
}


/** return verblevel parameter set in SCIP */
int WeightedSolver::getVerbosity() const
{
   return verbosity_;
}
