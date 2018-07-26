/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   GMI/src/cmain.c
 * @brief  main file for GMI cut example
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include "sepa_gmi.h"

/** reads parameters */
static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name, or NULL */
   )
{
   if ( filename != NULL )
   {
      if ( SCIPfileExists(filename) )
      {
         SCIPinfoMessage(scip, NULL, "reading parameter file <%s> ...\n", filename);
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
         SCIPinfoMessage(scip, NULL, "parameter file <%s> not found - using default parameters.\n", filename);
   }
   else if ( SCIPfileExists("scipgmi.set") )
   {
      SCIPinfoMessage(scip, NULL, "reading parameter file <scipgmi.set> ...\n");
      SCIP_CALL( SCIPreadParams(scip, "scipgmi.set") );
   }

   return SCIP_OKAY;
}

/** starts SCIP */
static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< input file name */
   )
{
   /********************
    * Problem Creation *
    ********************/

   SCIPinfoMessage(scip, NULL, "read problem <%s> ...\n\n", filename);
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   SCIPinfoMessage(scip, NULL, "solve problem ...\n\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIPinfoMessage(scip, NULL, "primal solution:\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   SCIPinfoMessage(scip, NULL, "Statistics:\n");
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** starts user interactive mode */
static
SCIP_RETCODE interactive(
   SCIP*                 scip                /**< SCIP data structure */
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
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array with shell parameters */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /***********************
    * Version information *
    ***********************/

   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");


   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include GMI cut separator */
   SCIP_CALL( SCIPincludeSepaGMI(scip) );

   /**************
    * Parameters *
    **************/

   if ( argc >= 3 )
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

   if ( argc >= 2 )
   {
      SCIP_CALL( fromCommandLine(scip, argv[1]) );
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "\n");

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
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
