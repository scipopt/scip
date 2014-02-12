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

#include "WeightedSolver.h"
#include "Objectives.h"

/** standard constructor */
WeightedSolver::WeightedSolver(
   bool                  verbose,            /**< true if scip output should be displayed */
   SCIP_Real             timelimit,          /**< maximum allowed time in seconds */
   int                   solstore            /**< number of solutions stored in SCIP */
   )
   : scip_(NULL),
     timelimit_(timelimit),
     found_new_optimum_(false),
     duration_last_run_(0),
     status_(SCIP_STATUS_UNKNOWN),
     nruns_(0),
     weight_(NULL),
     cost_vector_(NULL)
{
   SCIPcreate(&scip_);
   assert(scip_ != NULL);
   SCIPsetMessagehdlrQuiet(scip_, !verbose);
   SCIPincludeDefaultPlugins(scip_);
   SCIPincludeReaderMop(scip_);
   SCIPsetIntParam(scip_, "limits/maxorigsol", solstore);
   SCIPsetIntParam(scip_, "limits/maxsol", solstore);
}

/** destructor */
WeightedSolver::~WeightedSolver()
{
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

/** returns the SCIP problem status */
SCIP_Status WeightedSolver::getStatus() const
{
   return status_;
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

/** writes the last solution to a file in folder solutions/ with name <instance name>-<solution number>.sol*/
SCIP_RETCODE WeightedSolver::writeSolution()
{
   std::stringstream     s_outfile;
   FILE*                 fsol;

   /* build and save file name */
   s_outfile << outfilestump_
             << "-"
             << filename_by_point_.size() + 1
             << ".sol";
   solution_file_name_ = s_outfile.str();
   filename_by_point_[cost_vector_] = solution_file_name_;

   /* write file */
   fsol = fopen(solution_file_name_.c_str(), "w");
   SCIP_CALL( SCIPprintSol(scip_, solution_, fsol, FALSE) );

   return SCIP_OKAY;
}
