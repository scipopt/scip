/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for healthcare pricer example
 * @author Tobias Achterberg
 * @author Arne Nielsen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* user defined includes */
#include "pricer_healthcare.h"
#include "probdata_healthcare.h"


static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name, or NULL */
   )
{
   if( filename != NULL )
   {
      if( SCIPfileExists(filename) )
      {
         printf("reading parameter file <%s>\n", filename);
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
         printf("parameter file <%s> not found - using default parameters\n", filename);
   }
   else if( SCIPfileExists("pmedian.set") )
   {
      printf("reading parameter file <pmedian.set>\n");
      SCIP_CALL( SCIPreadParams(scip, "pmedian.set") );
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   printf("\nread problem <%s>\n", filename);
   printf("============\n\n");
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n\n");
   SCIP_CALL( SCIPsolve(scip) );

   printf("\nprimal solution:\n");
   printf("================\n\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   printf("\nStatistics\n");
   printf("==========\n\n");

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE runSCIP(
   int                   argc,
   char**                argv
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
   SCIP_CALL( SCIPcreate(&scip) );

   /* include healthcare pricer */
   SCIP_CALL( HCPincludePricerHealthcare(scip) );
   
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

      /* read problem data */
      SCIP_CALL( HCPcreateProbHealthcare(scip, "healthcare") );

      /* generate the model constraints */
      HCPgenerateModel(scip);

      /* activate the pricer */
      SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "healthcare")) );
   }

   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
   

   /***
       Solve:
   ***/
   SCIP_CALL( SCIPsolve(scip) );

   /***
       Statistics:
    ***/
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPprintTransProblem(scip, NULL, NULL, FALSE) );

   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );

   
   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                   argc,
   char**                argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode, stderr);
      return -1;
   }

   return 0;
}
