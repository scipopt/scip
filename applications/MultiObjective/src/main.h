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

/**@file   main.h
 * @brief  main file
 * @author Timo Strunk
 *
 * This is the main class of the program responsible for reading the arguments and running the program.
 * Furthermore it is responsible for writing standard output.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iomanip>
#include <sstream>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

#include "scip/def.h"

#ifndef CLASS_MAIN
#define CLASS_MAIN

/* spacing constants for tabular printing */
const int           WIDTH_DEFAULT       = 10;
const int           WIDTH_FILE          = 40;
const int           WIDTH_VEC_ENTRY     = 12;
const int           WIDTH_VEC_PADDING   = 4;

/** writes a vector to an output stream */
std::ostream& operator<<(
   std::ostream&                   os,       /** stream the vector should be written to*/
   const std::vector<SCIP_Real>&   v         /** vector that should be written */
   );

/** return the scalar product of two vectors */
SCIP_Real scalar_product(
   const std::vector<SCIP_Real>&   u,
   const std::vector<SCIP_Real>&   v
   );

class WeightedSolver;

/** class responsible for executing the program */
class Main
{
 public:
  /** default constructor */
   Main();

   /** default destructor */
   ~Main();

   /** execute the program */
   SCIP_RETCODE run(
      int                argc,               /**< number of command line arguments */
      char**             argv                /**< array of command line arguments */
      );

 private:
   WeightedSolver*       solver_;            /**< actual solver */
   char*                 filename_;          /**< problem file name */
   int                   width_vec_;         /**< space for writing a vector to log */

   int                   nnodes_total_;      /**< sum of branch and bound nodes used in all SCIP runs so far */
   int                   niterations_total_; /**< sum of lp iterations used in all SCIP runs so far */
   int                   n_v_new_total_;     /**< sum of 1-skeleton vertices generated so far */
   int                   n_v_proc_total_;    /**< sum of 1-skeleton vertices processed so far */

   /** determine name of SCIP parameter file from arguments */
   const char* readParamfilename(
      int                argc,          /**< number of command line arguments */
      char**             argv           /**< array of command line arguments */
      );

   /** prints comments about the result of file reading */
   SCIP_RETCODE readProblem(const char* filename);

   /** prints head row for run statistics table */
   void printHeadline();

   /** prints statistics for last run in table row format */
   void printRun();

   /** prints message according to solver status code */
   void evaluateStatus();

   /** print all solutions written so files by solver */
   void printSolutions();

   /** print every unbounded cost ray */
   void printUnboundedRays();

   /** prints overall statistics */
   void printBottomLine();
};

#endif
