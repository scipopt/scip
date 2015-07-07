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

/**@file   main.cpp
 * @brief  main file
 * @author Timo Strunk
 *
 * This is the main class of the program responsible for reading the arguments and running the program.
 * Furthermore it is responsible for writing standard output.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "WeightedSolver.h"
#include "LiftedWeightSpaceSolver.h"

#include "main.h"

/** writes a vector to an output stream */
std::ostream& operator<<(
   std::ostream&                   os,       /**< stream the vector should be written to*/
   const std::vector<SCIP_Real>&   v         /**< vector that should be written */
  )
{
   std::stringstream ss;

   ss.precision(WIDTH_VEC_ENTRY-7);
   /* build string representation: leading "[", fixed spaced value strings and ending "]" */
   ss << "[";
   for( std::vector<SCIP_Real>::const_iterator it = v.begin(); it<v.end(); ++it )
   {
      SCIP_Real x = *it;
      ss <<  std::setw(WIDTH_VEC_ENTRY);
      if( x >= SCIP_DEFAULT_INFINITY )
      {
        ss << "inf";
      }
      else if( x <= - SCIP_DEFAULT_INFINITY )
      {
        ss << "-inf";
      }
      else
      {
        ss << x;
      }
   }
   ss << " ]";

   /* stream representation string to output */
   os << ss.str();

   return os;
}

/** return the scalar product of two vectors */
SCIP_Real scalar_product(
   const std::vector<SCIP_Real>&   u,
   const std::vector<SCIP_Real>&   v
   )
{
   assert( u.size() == v.size() );

   SCIP_Real result = 0.;

   std::vector<SCIP_Real>::const_iterator uit = u.begin();
   std::vector<SCIP_Real>::const_iterator vit = v.begin();

   while( uit != u.end() )
   {
      result += (*uit) * (*vit);
      ++uit;
      ++vit;
   }

   return result;
}

/** default constructor */
Main::Main()
   :     solver_(NULL),
         filename_(NULL),
         nnodes_total_(0),
         niterations_total_(0),
         n_v_new_total_(0),
         n_v_proc_total_(0)
{
}

/** default destructor */
Main::~Main()
{
   delete solver_;
}

/** run the program */
SCIP_RETCODE Main::run(
   int                argc,          /**< number of command line arguments */
   char**             argv           /**< array of command line arguments */
   )
{
   if(argc <= 1)
   {
      std::cout << "usage: multiopt <problem file>.mop [<parameter file>.set]"
                << std::endl;

      return SCIP_READERROR;
   }

   const char* paramfilename = readParamfilename(argc, argv);

   solver_ = new LiftedWeightSpaceSolver(paramfilename);

   SCIP_CALL( readProblem(argv[1]) );

   std::cout << "solving" << std::endl;

   printHeadline();

   /* main loop iterating over all weights */
   while( solver_->hasNext() )
   {
      SCIP_CALL( solver_->solveNext() );

      nnodes_total_ += solver_->getNNodesLastRun();
      niterations_total_ += solver_->getNLPIterationsLastRun();
      n_v_new_total_ += solver_->getNNewVertices();
      n_v_proc_total_ += solver_->getNProcessedVertices();

      printRun();
   }

   evaluateStatus();

   SCIP_CALL( solver_->checkAndWriteSolutions() );
   printSolutions();

   printUnboundedRays();

   printBottomLine();

   return SCIP_OKAY;
}

/** determine name of SCIP parameter file from arguments */
const char* Main::readParamfilename(
   int                argc,          /**< number of command line arguments */
   char**             argv           /**< array of command line arguments */
   )
{
   const char* result = NULL;

   if(argc >= 3)
   {
      result = argv[2];
   }
   else
   {
      result = "scipmip.set";
   }

   if( SCIPfileExists(result) )
   {
      std::cout << "reading parameter file <"
                << result << ">"
                << std::endl;
   }
   else
   {
      std::cout << "parameter file "
                << "<" << result << ">"
                << " not found - using default parameters"
                << std::endl;
      result = NULL;
   }

   return result;
}

/** prints comments about the result of file reading */
SCIP_RETCODE Main::readProblem(
   const char*           filename            /**< filename of file to read */
   )
{
   SCIP_RETCODE          result;

   if( filename == NULL )
   {
     result = SCIP_NOFILE;
   }
   else
   {
      /* try to read input file and evaluate output */
      result = solver_->readProblem(filename);
      if( result == SCIP_READERROR )
      {
        std::cout << "file not in multi objective mps format: "
                  << filename
                  << std::endl;
      }
      else if( result == SCIP_NOFILE )
      {
        std::cout << "file not found "<< filename << std::endl;
      }
      else if( result == SCIP_PLUGINNOTFOUND)
      {
        std::cout << "file extension needs to be .mop" << std::endl;
      }
      else if( result != SCIP_OKAY)
      {
        std::cout << "usage: exa_can <instance>.mop" << std::endl;
      }
      else
      {
         std::cout << "reading problem file "
                   << "<" << filename << ">"
                   << std::endl;
      }
   }

   /* make enough space to write vectors */
   width_vec_ = WIDTH_VEC_ENTRY * solver_->getNObjs() + WIDTH_VEC_PADDING;

   return result;
}

/** prints head row for run statistics table */
void Main::printHeadline()
{
   std::cout << std::setw(width_vec_)        << "Weight"
             << std::setw(width_vec_)        << "Cost"
             << std::setw(WIDTH_DEFAULT)     << "Time/s"
             << std::setw(WIDTH_DEFAULT)     << "Nodes"
             << std::setw(WIDTH_DEFAULT)     << "LP Iter"
             << std::setw(WIDTH_DEFAULT)     << "V_new"
             << std::setw(WIDTH_DEFAULT)     << "V_proc"
             << std::endl;
}

/** prints statistics for last run in table row format */
void Main::printRun()
{
   if( solver_->getVerbosity() > 0 )
   {
      std::cout << "----------------------------------------------------"<<std::endl;
      printHeadline();
   }

   /* print weight */
   std::cout << std::setw(width_vec_) << *(solver_->getWeight());

   /* print new optimum if one was found */
   if( solver_->foundNewOptimum() )
   {
      std::cout << std::setw(width_vec_) << *(solver_->getCost());
   }
   else if( solver_->isWeightedUnbounded() )
   {
      std::cout << std::setw(width_vec_) << "unbounded";
   }
   else
   {
      std::cout << std::setw(width_vec_) << "-";
   }

   /* print solver_ statistics */
   std::cout  << std::setprecision(2) << std::fixed
              << std::setw(WIDTH_DEFAULT) << solver_->getDurationLastRun()
              << std::setw(WIDTH_DEFAULT) << solver_->getNNodesLastRun()
              << std::setw(WIDTH_DEFAULT) << solver_->getNLPIterationsLastRun()
              << std::setw(WIDTH_DEFAULT) << solver_->getNNewVertices()
              << std::setw(WIDTH_DEFAULT) << solver_->getNProcessedVertices()
              << std::endl;

   if( solver_->getVerbosity() > 0 )
   {
      std::cout << std::endl;
   }
}

/** prints message according to SCIP status code */
void Main::evaluateStatus()
{
   SCIP_Status solver_status;

   solver_status = solver_->getStatus();
   if( solver_status == SCIP_STATUS_TIMELIMIT )
   {
      std::cout << " ABORTED: time limit reached" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_INFEASIBLE )
   {
      std::cout << " ABORTED: problem is infeasible" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_INFORUNBD )
   {
      std::cout << " ABORTED: problem is infeasible or unbounded" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_USERINTERRUPT )
   {
      std::cout << " ABORTED: user interrupt" << std::endl;
   }
   else if( solver_status != SCIP_STATUS_OPTIMAL
            && solver_status != SCIP_STATUS_UNBOUNDED )
   {
      std::cout << " ABORTED: SCIP_STATUS = " << (int)solver_status << std::endl;
   }
}

/** print all solutions written so files by solver */
void Main::printSolutions()
{
   const std::vector<std::string>* objnames         = solver_->getObjNames();
   const std::map< const std::vector<SCIP_Real>*,
                   const char* >*  cost_to_filename = solver_->getCostToFilename();
   const std::vector<SCIP_Real>*   cost_vector;

   std::cout << "extremal solutions" << std::endl;

   for( std::vector<std::string>::const_iterator it = objnames->begin();
        it != objnames->end();
        ++it )
   {
      std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
   }

   std::cout << "   solution file" << std::endl;

   for( std::map< const std::vector<SCIP_Real>*,
           const char* >::const_iterator it = cost_to_filename->begin();
        it != cost_to_filename->end();
        ++it )
   {
      cost_vector = it->first;

      for( std::vector<SCIP_Real>::const_iterator jt = cost_vector->begin();
           jt != cost_vector->end();
           ++jt )
      {
         std::cout << std::setw(WIDTH_VEC_ENTRY) << *jt;
      }

      std::cout << "   " << it->second << std::endl;
   }
}

/** print every unbounded cost ray */
void Main::printUnboundedRays()
{
   const std::vector< const std::vector<SCIP_Real>* >*
                                   cost_rays = solver_->getCostRays();
   const std::vector<std::string>* objnames  = solver_->getObjNames();
   const std::vector<SCIP_Real>*   cost_vector;

   if( !cost_rays->empty() )
   {
      std::cout << "unbounded rays" << std::endl;

      for( std::vector<std::string>::const_iterator it = objnames->begin();
           it != objnames->end();
           ++it )
      {
         std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
      }

      std::cout << std::endl;

      for( std::vector< const std::vector<SCIP_Real>* >::const_iterator
              it = cost_rays->begin();
           it != cost_rays->end();
           ++it
         )
      {
         cost_vector = *it;
         for( std::vector<SCIP_Real>::const_iterator it = cost_vector->begin();
              it != cost_vector->end();
              ++it )
         {
            std::cout << std::setw(WIDTH_VEC_ENTRY) << *it;
         }

         std::cout << std::endl;
      }
   }
}

/** prints overall statistics */
void Main::printBottomLine()
{
   std::cout << std::endl
             << std::setw(width_vec_)         << "Runs"
             << std::setw(width_vec_)         << "Solutions"
             << std::setw(WIDTH_DEFAULT)      << "Time/s"
             << std::setw(WIDTH_DEFAULT)      << "Nodes"
             << std::setw(WIDTH_DEFAULT)      << "LP Iter"
             << std::setw(WIDTH_DEFAULT)      << "V_new"
             << std::setw(WIDTH_DEFAULT)      << "V_proc"
             << std::endl                     << std::endl
             << std::setw(width_vec_)         << solver_->getNRuns()
             << std::setw(width_vec_)         << solver_->getNSolutions()
             << std::setw(WIDTH_DEFAULT)      << solver_->getTotalDuration()
             << std::setw(WIDTH_DEFAULT)      << nnodes_total_
             << std::setw(WIDTH_DEFAULT)      << niterations_total_
             << std::setw(WIDTH_DEFAULT)      << n_v_new_total_
             << std::setw(WIDTH_DEFAULT)      << n_v_proc_total_
             << std::endl;
}

/** runs the program: reads specified instance file, calls the solver and prints out
 *  relevant information */
int main(
   int                   argc,          /**< number of arguments */
   char**                argv           /**< array of arguments */
   )
{
   SCIP_RETCODE retcode;
   Main main;

   retcode = main.run(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      return -1;
   }

   return 0;
}
