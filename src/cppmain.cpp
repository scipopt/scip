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
#pragma ident "@(#) $Id: cppmain.cpp,v 1.1 2004/11/01 13:49:56 bzfpfend Exp $"

/**@file   cppmain.cpp
 * @brief  main file for C++ compilation
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>

extern "C"
{
#include "scip.h"
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
         printf("reading parameter file <%s>\n", filename);
         CHECK_OKAY( SCIPreadParams(scip, filename) );
      }
      else
         printf("parameter file <%s> not found - using default parameters\n", filename);
   }
   else if( SCIPfileExists("scip.set") )
   {
      printf("reading parameter file <scip.set>\n");
      CHECK_OKAY( SCIPreadParams(scip, "scip.set") );
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

   printf("\nread problem <%s>\n", filename);
   printf("============\n\n");
   CHECK_OKAY( SCIPreadProb(scip, filename) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n\n");
   CHECK_OKAY( SCIPsolve(scip) );

   printf("\nprimal solution:\n");
   printf("================\n\n");
   CHECK_OKAY( SCIPprintBestSol(scip, NULL) );


   /**************
    * Statistics *
    **************/

   printf("\nStatistics\n");
   printf("==========\n\n");

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
   printf("\n");


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
   /*CHECK_OKAY( SCIPwriteParams(scip, "scip.set", TRUE) );*/


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
