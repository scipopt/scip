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

/**@file   main.cpp
 * @brief  main file
 * @author Timo Strunk
 *
 * @desc   This is the main class of the program responsible for reading the arguments and running the program.
 *         Furthermore it is responsible for writing standard output.
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

/** constructor from command line arguments */
Main::Main()
   :     filename_(NULL),
         verbose_(DEFAULT_VERBOSE),
         timelimit_(DEFAULT_TIMELIMIT),
         solstore_(DEFAULT_SOLSTORE),
         log_(DEFAULT_LOG),
         solver_(NULL),
         os_(NULL),
         nnodes_total_(0),
         niterations_total_(0)
{
   /*
   readParameters(argc,argv );
   solver_ = new LiftedWeightSpaceSolver(verbose_, timelimit_, solstore_);
   makeOutputStream();
   */
}

/** default destructor */
Main::~Main()
{
   delete solver_;

   if( log_ )
   {
     static_cast<std::ofstream*>(os_)->close();
     delete os_;
   }
}

/** run the program */
SCIP_RETCODE Main::run(   
   int                argc,          /**< number of command line arguments */
   char**             argv           /**< array of command line arguments */
   )
{
   const char* paramfilename;
   if(argc <= 1)
   {
      std::cout << "usage: multiopt <problem file>.mop [<parameter file>.set]" << std::endl;
      return SCIP_READERROR;
   }
   else if(argc >= 3)
   {
      paramfilename = argv[2];
   }
   else
   {
      paramfilename = NULL;
   }
   solver_ = new LiftedWeightSpaceSolver(paramfilename);
   os_  = &std::cout;

   SCIP_CALL( readProblem(argv[1]) );

   /* make enough space to write vectors */
   width_vec_ = WIDTH_VEC_ENTRY * solver_->getNObjs() + WIDTH_VEC_PADDING;

   printHeadline();

   /* main loop iterating over all weights */
   while( solver_->hasNext() )
   {
      SCIP_CALL( solver_->next() );
      SCIP_CALL( solver_->solve() );

      if( solver_->foundNewOptimum() )
      {
         SCIP_CALL( solver_->writeSolution() );
      }

      nnodes_total_ += solver_->getNNodesLastRun();
      niterations_total_ += solver_->getNLPIterationsLastRun();
      n_v_new_total_ += solver_->getNNewVertices();
      n_v_proc_total_ += solver_->getNProcessedVertices();

      printRun();
   }

   evaluateStatus();

   SCIP_CALL( solver_->enforceExtremality() );

   printBottomLine();

   return SCIP_OKAY;
}

/** read parameters from command line input */
void Main::readParameters(
   int                   argc,               /**< number of command line arguments */
   char**                argv                /**< command line arguments */
   )
{
   int i;

   /* loop over arguments */
   i = 1;
   while( i < argc )
   {
      if( !strcmp(argv[i], "-v") )
      {
         verbose_ = TRUE;
      }
      else if( !strcmp(argv[i], "-l") )
      {
         log_ = TRUE;
      }
      else if( !strcmp(argv[i], "-t") )
      {
         ++i;
         timelimit_ = atof(argv[i]);
      }
      else if( !strcmp(argv[i], "-s") )
      {
         ++i;
         solstore_ = atoi(argv[i]);
      }
      else if( filename_ == NULL )
      {
         filename_ = argv[i];
      }
      /* else: ignore */
      ++i;
   }
}

/** generates output stream according to log flag */
void Main::makeOutputStream()
{
   std::string           sfilename;
   int                   startname;
   std::string           logname;

   if( log_ )
   {
      /* generate stream to output file */
      os_       = new std::ofstream();
      sfilename = std::string(filename_);
      startname = sfilename.find_last_of("/\\");

      logname   = "solutions"+sfilename.substr(
         startname,
         sfilename.length() - 4 - startname
         ) + ".log";

      static_cast<std::ofstream*>(os_)->open(logname.c_str());
   }
   else
   {
      os_  = &std::cout;
   }
}

/** prints comments about the result of file reading */
SCIP_RETCODE Main::readProblem(const char* filename)
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
        std::cout << "file extension needs to be .momps" << std::endl;
      }
      else if( result != SCIP_OKAY)
      {
        std::cout << "usage: exa_can <instance>.momps" << std::endl;
      }
   }

   return result;
}

/** prints head row for run statistics table */
void Main::printHeadline()
{
   *os_ << std::setw(width_vec_)        << "Weight"
        << std::setw(width_vec_)        << "Cost"
        << std::setw(WIDTH_DEFAULT)     << "Time/s"
        << std::setw(WIDTH_DEFAULT)     << "Nodes"
        << std::setw(WIDTH_DEFAULT)     << "LP Iter"
        << std::setw(WIDTH_DEFAULT)     << "V_new"
        << std::setw(WIDTH_DEFAULT)     << "V_proc"
        << "   Solution File"
        << std::endl;
}

/** prints statistics for last run in table row format */
void Main::printRun()
{
   if( verbose_ )
   {
      *os_ << "----------------------------------------------------"<<std::endl;
      printHeadline();
   }

   /* print weight */
   *os_ << std::setw(width_vec_) << *(solver_->getWeight());

   /* print new optimum if one was found */
   if( solver_->foundNewOptimum() )
   {
      *os_ << std::setw(width_vec_) << *(solver_->getCost());
   }
   else
   {
      *os_ << std::setw(width_vec_) << "-";
   }

   /* print solver_ statistics */
   *os_  << std::setprecision(2) << std::fixed
         << std::setw(WIDTH_DEFAULT) << solver_->getDurationLastRun()
         << std::setw(WIDTH_DEFAULT) << solver_->getNNodesLastRun()
         << std::setw(WIDTH_DEFAULT) << solver_->getNLPIterationsLastRun()
         << std::setw(WIDTH_DEFAULT) << solver_->getNNewVertices()
         << std::setw(WIDTH_DEFAULT) << solver_->getNProcessedVertices();

   /* print name of solution if optimum was found */
   if( solver_->foundNewOptimum() )
   {
      *os_ << "   " << solver_->getSolutionFileName();
   }

   *os_ << std::endl;

   if(verbose_)
   {
      *os_ << std::endl;
   }
}

/** prints message according to SCIP status code */
void Main::evaluateStatus()
{
   SCIP_Status solver_status;

   solver_status = solver_->getStatus();
   if( solver_status == SCIP_STATUS_TIMELIMIT )
   {
      *os_ << " ABORTED: time limit reached" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_INFEASIBLE )
   {
      *os_ << " ABORTED: problem is infeasible" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_UNBOUNDED )
   {
      *os_ << " ABORTED: problem has at least one infeasible objective" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_INFORUNBD )
   {
      *os_ << " ABORTED: problem is infeasible or unbounded" << std::endl;
   }
   else if( solver_status == SCIP_STATUS_USERINTERRUPT )
   {
      *os_ << " ABORTED: user interrupt" << std::endl;
   }
   else if( solver_status != SCIP_STATUS_OPTIMAL )
   {
      *os_ << " ABORTED: SCIP_STATUS = " << (int)solver_status << std::endl;
   }
}

/** prints overall statistics */
void Main::printBottomLine()
{
   *os_ << std::endl
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
   int                   argc,          /** number of arguments */
   char**                argv           /** array of arguments */
   )
{
   Main* main;

   main = new Main();
   SCIP_CALL(main->run(argc, argv));

   return 0;
}
