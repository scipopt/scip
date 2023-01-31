/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   examples/MIPSolver/src/cppmain.cpp
 * @brief  main file for C++ example project using SCIP as a callable library
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iostream>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"


/** reads parameters */
static
SCIP_RETCODE readParams(
   SCIP*                      scip,               /**< SCIP data structure */
   const char*                filename            /**< parameter file name, or NULL */
   )
{
   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         std::cout << "reading parameter file <" << filename << ">" << std::endl;
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
         std::cout << "parameter file <" << filename << "> not found - using default parameters" << std::endl;
   }
   else if( SCIPfileExists("scipmip.set") )
   {
      std::cout << "reading parameter file <scipmip.set>" << std::endl;
      SCIP_CALL( SCIPreadParams(scip, "scipmip.set") );
   }

   return SCIP_OKAY;
}

/** starts SCIP */
static
SCIP_RETCODE fromCommandLine(
   SCIP*                      scip,               /**< SCIP data structure */
   const char*                filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   std::cout << std::endl << "read problem <" << filename << ">" << std::endl;
   std::cout << "============" << std::endl << std::endl;
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   std::cout << "solve problem" << std::endl;
   std::cout << "=============" << std::endl;
   SCIP_CALL( SCIPsolve(scip) );

   std::cout << std::endl << "primal solution:" << std::endl;
   std::cout << "================" << std::endl << std::endl;
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   std::cout << std::endl << "Statistics" << std::endl;
   std::cout << "==========" << std::endl << std::endl;

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** starts user interactive mode */
static
SCIP_RETCODE interactive(
   SCIP*                      scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( SCIPstartInteraction(scip) );

   return SCIP_OKAY;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runSCIP(
   int                        argc,               /**< number of shell parameters */
   char**                     argv                /**< array with shell parameters */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(scip, NULL);
   std::cout << std::endl;


   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**************
    * Parameters *
    **************/

   if( argc >= 3 )
   {
      SCIP_CALL( readParams(scip, argv[2]) );
   }
   else
   {
      SCIP_CALL( readParams(scip, NULL) );
   }
   /*CHECK_OKAY( SCIPwriteParams(scip, "scipmip.set", TRUE) );*/


   /**************
    * Start SCIP *
    **************/

   if( argc >= 2 )
   {
      SCIP_CALL( fromCommandLine(scip, argv[1]) );
   }
   else
   {
      printf("\n");

      SCIP_CALL( interactive(scip) );
   }

   
   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
