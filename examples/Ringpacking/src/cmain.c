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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Ringpacking/src/cmain.c
 * @brief  Main file for ringpacking pricing example
 * @author Benjamin Mueller
 *
 *  This the file contains the \ref main() main function of the projects. This includes all the default plugins of
 *  \SCIP and the once which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */
#include <stdio.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "cons_rpa.h"
#include "reader_rpa.h"
#include "pricer_rpa.h"

/* parameters */
#define DEFAULT_VERIFICATION_NLPTILIM        10.0      /**< time limit for verification NLP */
#define DEFAULT_VERIFICATION_NLPNODELIM      10000L    /**< node limit for verification NLP */
#define DEFAULT_VERIFICATION_HEURTILIM       10.0      /**< time limit for heuristic verification */
#define DEFAULT_VERIFICATION_HEURITERLIM     1000      /**< iteration limit for heuristic verification */
#define DEFAULT_VERIFICATION_NLPTILIMSOFT    1.0       /**< soft time limit for verification NLP */
#define DEFAULT_VERIFICATION_NLPNODELIMSOFT  1000L     /**< soft node limit for verification NLP */
#define DEFAULT_VERIFICATION_HEURTILIMSOFT   1.0       /**< soft time limit for heuristic verification */
#define DEFAULT_VERIFICATION_HEURITERLIMSOFT 100       /**< soft iteration limit for heuristic verification */
#define DEFAULT_PRICING_NLPTILIM             10.0      /**< time limit for each pricing NLP */
#define DEFAULT_PRICING_NLPNODELIM           SCIP_LONGINT_MAX /**< node limit for each pricing NLP */
#define DEFAULT_PRICING_HEURTILIM            60.0     /**< time limit for each heuristic pricing */
#define DEFAULT_PRICING_HEURITERLIM          1000      /**< iteration limit for each heuristic pricing */

#define DEFAULT_TEXFILENAME                  ""        /**< filename for tex output for the best found solution (\"\": disable) */

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
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

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include reader for ringpacking instances */
   SCIP_CALL( SCIPincludeReaderRpa(scip) );

   /* include ringpacking constraint handler */
   SCIP_CALL( SCIPincludeConshdlrRpa(scip) );

   /* include ringpacking pricer  */
   SCIP_CALL( SCIPincludePricerRingpacking(scip) );

   /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts",0) );

   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /*
    * parameters for verification
    */

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/verification/nlptilim",
         "time limit for verification NLP",
         NULL, FALSE, DEFAULT_VERIFICATION_NLPTILIM, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "ringpacking/verification/nlpnodelim",
         "node limit for verification NLP",
         NULL, FALSE, DEFAULT_VERIFICATION_NLPNODELIM, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/verification/heurtilim",
         "time limit for heuristic verification",
         NULL, FALSE, DEFAULT_VERIFICATION_HEURTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "ringpacking/verification/heuriterlim",
         "iteration limit for heuristic verification",
         NULL, FALSE, DEFAULT_VERIFICATION_HEURITERLIM, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/verification/nlptilimsoft",
         "soft time limit for verification NLP",
         NULL, FALSE, DEFAULT_VERIFICATION_NLPTILIMSOFT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "ringpacking/verification/nlpnodelimsoft",
         "soft node limit for verification NLP",
         NULL, FALSE, DEFAULT_VERIFICATION_NLPNODELIMSOFT, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/verification/heurtilimsoft",
         "soft time limit for heuristic verification",
         NULL, FALSE, DEFAULT_VERIFICATION_HEURTILIMSOFT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "ringpacking/verification/heuriterlimsoft",
         "soft iteration limit for heuristic verification",
         NULL, FALSE, DEFAULT_VERIFICATION_HEURITERLIMSOFT, 0, INT_MAX, NULL, NULL) );

   /*
    * parameters for pricing
    */

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/pricing/nlptilim",
         "time limit for each pricing NLP",
         NULL, FALSE, DEFAULT_PRICING_NLPTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "ringpacking/pricing/nlpnodelim",
         "node limit for each pricing NLP",
         NULL, FALSE, DEFAULT_PRICING_NLPNODELIM, 0L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "ringpacking/pricing/heurtilim",
         "time limit for each heuristic pricing",
         NULL, FALSE, DEFAULT_PRICING_HEURTILIM, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "ringpacking/pricing/heuriterlim",
         "iteration limit for each heuristic pricing",
         NULL, FALSE, DEFAULT_PRICING_HEURITERLIM, 0, INT_MAX, NULL, NULL) );

   /*
    * miscellaneous parameters
    */

   SCIP_CALL( SCIPaddStringParam(scip,
         "ringpacking/texfilename",
         "filename for tex output for the best found solution (\"\": disable)",
         NULL, FALSE, DEFAULT_TEXFILENAME, NULL, NULL) );

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runShell(argc, argv, "scip.set");
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
