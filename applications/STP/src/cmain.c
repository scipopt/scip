/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   STP/src/cmain.c
 * @brief  Main file for SCIP-Jack
 * @author Gerald Gamrath
 * @author Daniel Rehfeldt
 *
 *  This the file contains the \ref main() main function of the projects. This includes all the default plugins of
 *  \SCIP and the once which belong to that projects. After that is starts the interactive shell of \SCIP or processes
 *  the shell arguments if given.
 */
#include <stdio.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "reader_stp.h"
#include "reader_gr.h"
#include "cons_stp.h"
#include "heur_tm.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "heur_ascendprune.h"
#include "heur_slackprune.h"
#include "heur_lurkprune.h"
#include "heur_rec.h"
#include "pricer_stp.h"
#include "event_bestsol.h"
#include "probdata_stp.h"
#include "dialog_stp.h"
#include "prop_stp.h"
#include "branch_stp.h"

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array containing shell parameters */
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

   SCIP_CALL( SCIPincludePricerStp(scip) );

   SCIP_CALL( SCIPincludeReaderStp(scip) );
   SCIP_CALL( SCIPincludeReaderGr(scip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPincludeDialogStp(scip) );

   SCIP_CALL( SCIPincludeConshdlrStp(scip) );

   SCIP_CALL( SCIPStpIncludeHeurTM(scip) );

   SCIP_CALL( SCIPStpIncludeHeurLocal(scip) );

   SCIP_CALL( SCIPStpIncludeHeurRec(scip) );

   SCIP_CALL( SCIPStpIncludeHeurPrune(scip) );

   SCIP_CALL( SCIPStpIncludeHeurAscendPrune(scip) );

   SCIP_CALL( SCIPStpIncludeHeurSlackPrune(scip) );

   SCIP_CALL( SCIPStpIncludeHeurLurkPrune(scip) );

   /* include event handler for printing primal solution development */
   SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );

   SCIP_CALL( SCIPincludeBranchruleStp(scip) );

   SCIP_CALL( SCIPincludePropStp(scip) );

   /* set hard-coded default parameters */
   SCIP_CALL( SCIPprobdataSetDefaultParams(scip) );

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
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array containing shell parameters */
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
