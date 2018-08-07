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
#include "reader_gr.h"
#include "cons_stp.h"
#include "heur_tm.h"
#include "heur_local.h"
#include "heur_prune.h"
#include "heur_ascendprune.h"
#include "heur_slackprune.h"
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

   /* include stp pricer */
   SCIP_CALL( SCIPincludePricerStp(scip) );

   /* include Steiner tree reader */
   SCIP_CALL( SCIPincludeReaderStp(scip) );
   SCIP_CALL( SCIPincludeReaderGr(scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include STP dialog */
   SCIP_CALL( SCIPincludeDialogStp(scip) );

   /* include Steiner tree constraint handler */
   SCIP_CALL( SCIPincludeConshdlrStp(scip) );

   /* include shortest path heuristic */
   SCIP_CALL( SCIPStpIncludeHeurTM(scip) );

   /* include local heuristics */
   SCIP_CALL( SCIPStpIncludeHeurLocal(scip) );

   /* include recombination heuristic */
   SCIP_CALL( SCIPStpIncludeHeurRec(scip) );

   /* include pruning heuristic */
   SCIP_CALL( SCIPStpIncludeHeurPrune(scip) );

   /* include ascend-and-prune heuristic */
   SCIP_CALL( SCIPStpIncludeHeurAscendPrune(scip) );

   /* include slack-and-prune heuristic */
   SCIP_CALL( SCIPStpIncludeHeurSlackPrune(scip) );

   /* include event handler for printing primal solution development */
   SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleStp(scip) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropStp(scip) );

   SCIP_CALL( SCIPsetSubscipsOff(scip, FALSE) );

   /* set STP-specific default parameters */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(scip, "limits/maxsol", 400) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/rowagelimit", 30) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxroundsroot", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxstallroundsroot", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxcutsroot", 100000) );
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxcuts", 1000) );   // todo tune

   SCIP_CALL( SCIPsetRealParam(scip, "separating/minefficacyroot", 0.01) ); // todo tune
   SCIP_CALL( SCIPsetRealParam(scip, "separating/minorthoroot", 0.4) ); // todo tune > 0.4
   SCIP_CALL( SCIPsetRealParam(scip, "separating/minortho", 0.4) ); // todo tune > 0.4 best soplex: 0.8
   SCIP_CALL( SCIPsetRealParam(scip, "separating/objparalfac", 0.01) ); // todo tune < 0.1

   SCIP_CALL( SCIPsetRealParam(scip, "separating/intsupportfac", 0.0) );
   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/maxproprounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/alns/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/coefdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/feaspump/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/fracdiving/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/farkasdiving/freq", -1) );
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
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/oneopt/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/rounding/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/locks/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/probing/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/pseudoobj/timingmask", 5) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/redcost/freq", -1) );
   SCIP_CALL( SCIPsetRealParam(scip, "branching/relpscost/maxreliable", 1.0) );

   // todo test properly; normal dfs?
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/restartdfs/stdpriority", 400000) );


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
