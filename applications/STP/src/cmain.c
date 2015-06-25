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

/**@file   cmain.c
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
#include "cons_stp.h"
#include "heur_tm.h"
#include "heur_local.h"
#include "heur_rec.h"
#include "heur_rs.h"
#include "pricer_stp.h"
#include "event_bestsol.h"
#include "probdata_stp.h"
#include "dialog_stp.h"
#include "prop_stp.h"

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

   /* include stp pricer */
   SCIP_CALL( SCIPincludePricerStp(scip) );

   /* include steiner tree reader */
   SCIP_CALL( SCIPincludeReaderStp(scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include STP dialog */
   SCIP_CALL( SCIPincludeDialogStp(scip) );

   /* include steiner tree constraint handler */
   SCIP_CALL( SCIPincludeConshdlrStp(scip) );

   /* include Takahashi Matsuyama heuristic */
   SCIP_CALL( SCIPincludeHeurTM(scip) );

   /* include local heuristics */
   SCIP_CALL( SCIPincludeHeurLocal(scip) );

   /* include recombination heuristic */
   SCIP_CALL( SCIPincludeHeurRec(scip) );

   /* include event handler for printing primal solution development */
   SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropStp(scip) );

   /* set STP-specific default parameters */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxsol", 400) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/rowagelimit", 30) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxstallrounds", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/maxproprounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/coefdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/guideddiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/linesearchdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/nlpdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/objpscostdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/pscostdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/randrounding/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rootsoldiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shiftandpropagate/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/shifting/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/subnlp/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/undercover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/veclendiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/zirounding/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/pseudoobj/timingmask", 5) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/redcost/freq", -1) );
   SCIP_CALL( SCIPsetRealParam(scip, "branching/relpscost/maxreliable", 1.0) );
   SCIP_CALL( SCIPsetRealParam(scip, "separating/minefficacyroot", 0.01) );

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
