/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cppmain.cpp,v 1.1 2004/10/06 10:47:28 bzfpfend Exp $"

/**@file   cppmain.cpp
 * @brief  main file for C++ example project using SCIP as a callable library
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iostream>

#include "objscip.h"

extern "C"
{
#include "scipdefplugins.h"
}


static
RETCODE readParams(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< parameter file name, or NULL */
   )
{
   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         std::cout << "reading parameter file <" << filename << ">" << std::endl;
         CHECK_OKAY( SCIPreadParams(scip, filename) );
      }
      else
         std::cout << "parameter file <" << filename << "> not found - using default parameters" << std::endl;
   }
   else if( SCIPfileExists("example.set") )
   {
      std::cout << "reading parameter file <example.set>" << std::endl;
      CHECK_OKAY( SCIPreadParams(scip, "example.set") );
   }

   return SCIP_OKAY;
}

static
RETCODE fromCommandLine(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   std::cout << std::endl << "read problem <" << filename << ">" << std::endl;
   std::cout << "============" << std::endl << std::endl;
   CHECK_OKAY( SCIPreadProb(scip, filename) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   std::cout << "solve problem" << std::endl;
   std::cout << "=============" << std::endl;
   CHECK_OKAY( SCIPsolve(scip) );

   std::cout << std::endl << "primal solution:" << std::endl;
   std::cout << "================" << std::endl << std::endl;
   CHECK_OKAY( SCIPprintBestSol(scip, NULL) );


   /**************
    * Statistics *
    **************/

   std::cout << std::endl << "Statistics" << std::endl;
   std::cout << "==========" << std::endl << std::endl;

   CHECK_OKAY( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

static
RETCODE interactive(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   /* start user interactive mode */
   CHECK_OKAY( SCIPstartInteraction(scip) );

   return SCIP_OKAY;
}

static
RETCODE runSCIP(
   int              argc,
   char**           argv
   )
{
   SCIP* scip = NULL;


   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(NULL);
   std::cout << std::endl;


   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   CHECK_OKAY( SCIPcreate(&scip) );


   /* include default SCIP plugins */
   CHECK_OKAY( SCIPincludeDefaultPlugins(scip) );


   /**************
    * Parameters *
    **************/

   if( argc >= 3 )
   {
      CHECK_OKAY( readParams(scip, argv[2]) );
   }
   else
   {
      CHECK_OKAY( readParams(scip, NULL) );
   }
   /*CHECK_OKAY( SCIPwriteParams(scip, "example.set", TRUE) );*/


   /**************
    * Start SCIP *
    **************/

   if( argc >= 2 )
   {
      CHECK_OKAY( fromCommandLine(scip, argv[1]) );
   }
   else
   {
      printf("\n");

      CHECK_OKAY( interactive(scip) );
   }

   
   /********************
    * Deinitialization *
    ********************/

   CHECK_OKAY( SCIPfree(&scip) );


   /*****************************
    * Local Memory Deallocation *
    *****************************/

#ifndef NDEBUG
   memoryCheckEmpty();
#endif

   return SCIP_OKAY;
}

int
main(
   int              argc,
   char**           argv
   )
{
   RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode, stderr);
      return -1;
   }

   return 0;
}
