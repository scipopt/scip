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

/**@file   main.h
 * @brief  main file
 * @author Timo Strunk
 * @desc   This is the main class of the program responsible for reading the arguments and running the program.
 * Furthermore it is responsible for writing standard output.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iomanip>
#include <sstream>
#include <cstring>
#include <fstream>
#include <iostream>

#include "scip/def.h"

#ifndef CLASS_MAIN
#define CLASS_MAIN

/* default parameters */
const SCIP_Real     DEFAULT_TIMELIMIT   = 3600;
const SCIP_Bool     DEFAULT_VERBOSE     = FALSE;
const int           DEFAULT_SOLSTORE    = 1024;
const SCIP_Bool     DEFAULT_WEIGHTED    = TRUE;
const SCIP_Bool     DEFAULT_LOG         = FALSE;

/* spacing constants for tabular printing */
const int           WIDTH_SOLVER        = 13;
const int           WIDTH_TIME          = 10;
const int           WIDTH_NODES         = 10;
const int           WIDTH_ITERATIONS    = 10;
const int           WIDTH_FILE          = 40;
const int           WIDTH_STORE         = 6;
const int           WIDTH_VEC_ENTRY     = 12;
const int           WIDTH_VEC_PADDING   = 4;

/** writes a vector to an output stream */
std::ostream& operator<<(
   std::ostream&                   os,       /** stream the vector should be written to*/
   const std::vector<SCIP_Real>&   v         /** vector that should be written */
   );

class WeightedSolver;

/** class responsible for executing the program */
class Main
{
 public:
  /** constructor from program arguments */
   Main(
      int                argc,          /**< number of command line arguments */
      char**             argv           /**< array of command line arguments */
      );

   /** default destructor */
   ~Main();

   /** execute the program */
   SCIP_RETCODE main();

 private:
   char*                 filename_;          /**< problem file name */
   SCIP_Bool             verbose_;           /**< whether to print more log text */
   SCIP_Real             timelimit_;         /**< maximal time for entire solve */
   int                   solstore_;          /**< maximal solution storage size */
   SCIP_Bool             log_;               /**< whether to print log to file */

   WeightedSolver*       solver_;            /**< actual solver */

   std::ostream*         os_;                /**< log output stream */
   int                   width_vec_;         /**< space for writing a vector to log */

   int                   nnodes_total_;      /**< sum of branch and bound nodes used in all SCIP runs so far */
   int                   niterations_total_; /**< sum of lp iterations used in all SCIP runs so far */

   /** read parameters from command line input */
   void readParameters(  
      int                 argc,               /**< number of command line arguments */
      char**              argv                /**< command line arguments */
      );

   /** generates output stream according to log flag */
   void makeOutputStream();

   /** prints comments about the result of file reading */
   SCIP_RETCODE readFile();

   /** prints head row for run statistics table */
   void printHeadline();

   /** prints statistics for last run in table row format */
   void printRun();

   /** prints message according to solver status code */
   void evaluateStatus();

   /** prints overall statistics */
   void printBottomLine();
};

#endif
