/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cmain.c,v 1.69 2003/11/24 12:12:42 bzfpfend Exp $"

/**@file   cmain.c
 * @brief  main file for C compilation
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>

#include "scip.h"
#include "scipdefplugins.h"



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
      if( SCIPfileExists(argv[2]) )
      {
         printf("reading parameter file <%s>\n", argv[2]);
         CHECK_OKAY( SCIPreadParams(scip, argv[2]) );
      }
      else
         printf("parameter file <%s> not found - using default parameters\n", argv[2]);
   }
   else if( SCIPfileExists("scip.set") )
   {
      printf("reading parameter file <scip.set>\n");
      CHECK_OKAY( SCIPreadParams(scip, "scip.set") );
   }
   /*CHECK_OKAY( SCIPwriteParams(scip, "scip.set", TRUE) );*/


   /********************
    * Interactive mode *
    ********************/

   if( argc < 2 )
   {
      printf("\n");

      /* start user interactive mode */
      CHECK_OKAY( SCIPstartInteraction(scip) );

      return SCIP_OKAY;
   }

   
   /********************
    * Problem Creation *
    ********************/

   printf("\nread problem <%s>\n", argv[1]);
   printf("============\n\n");
   CHECK_OKAY( SCIPreadProb(scip, argv[1]) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n\n");
   CHECK_OKAY( SCIPsolve(scip) );

#if 1
   printf("\nprimal solution:\n");
   printf("================\n\n");
   CHECK_OKAY( SCIPprintBestSol(scip, NULL) );
#endif

#ifndef NDEBUG
   /*SCIPdebugMemory(scip);*/
#endif


   /**************
    * Statistics *
    **************/

   printf("\nStatistics\n");
   printf("==========\n\n");

   CHECK_OKAY( SCIPprintStatistics(scip, NULL) );


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

   /**@todo implement remaining events */
   /**@todo avoid addition of identical rows */
   /**@todo avoid addition of identical constraints */
   /**@todo pricing for pseudo solutions */
   /**@todo integrality check on objective function, abort if gap is below 1.0 */
   /**@todo statistics: count domain reductions and constraint additions of constraint handlers */
   /**@todo it's a bit ugly, that user call backs may be called before the nodequeue was processed */
   /**@todo unboundness detection in presolving -> convert problem into feasibility problem to decide
    *       unboundness/infeasibility */
   /**@todo variable event PSSOLCHANGED, update pseudo activities in constraints to speed up checking of pseudo solutions */
   /**@todo branching rule acting as a filter by temporary changing the branching priority of variables and returning 
    *       SCIP_DIDNOTFIND to let the next branching rule select the branching variable */
   /**@todo in the documentation, make a clear distinction between external "user" methods and internal methods in
    *       the .h files */

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode, stderr);
      return -1;
   }

   return 0;
}
