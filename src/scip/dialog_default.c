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
#pragma ident "@(#) $Id: dialog_default.c,v 1.22 2004/06/24 15:34:36 bzfpfend Exp $"

/**@file   dialog_default.c
 * @brief  default user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "dialog_default.h"



/** executes a menu dialog */
static
RETCODE dialogExecMenu(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          dialog,             /**< dialog menu */
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG**         nextdialog          /**< pointer to store next dialog to execute */
   )
{
   const char* command;
   Bool again;
   int nfound;

   do
   {
      again = FALSE;

      /* get the next word of the command string */
      command = SCIPdialoghdlrGetWord(dialoghdlr, dialog, NULL);
      
      /* exit to the root dialog, if command is empty */
      if( command[0] == '\0' )
      {
         CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );
         *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
         return SCIP_OKAY;
      }
      else if( strcmp(command, "..") == 0 )
      {
         *nextdialog = SCIPdialogGetParent(dialog);
         if( *nextdialog == NULL )
            *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);
         return SCIP_OKAY;
      }
      
      /* find command in dialog */
      nfound = SCIPdialogFindEntry(dialog, command, nextdialog);
      
      /* check result */
      if( nfound == 0 )
      {
         printf("command <%s> not available\n", command);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
         *nextdialog = dialog;
      }
      else if( nfound >= 2 )
      {
         printf("\npossible completions:\n");
         CHECK_OKAY( SCIPdialogDisplayCompletions(dialog, scip, command) );
         printf("\n");
         SCIPdialoghdlrClearBuffer(dialoghdlr);
         again = TRUE;
      }
   }
   while( again );

   return SCIP_OKAY;
}

/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
DECL_DIALOGEXEC(SCIPdialogExecMenu)
{  /*lint --e{715}*/
   /* if remaining command string is empty, display menu of available options */
   if( SCIPdialoghdlrIsBufferEmpty(dialoghdlr) )
   {
      printf("\n");
      CHECK_OKAY( SCIPdialogDisplayMenu(dialog, scip) );
      printf("\n");
   }

   CHECK_OKAY( dialogExecMenu(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** standard menu dialog execution method, that doesn't display it's help screen */
DECL_DIALOGEXEC(SCIPdialogExecMenuLazy)
{  /*lint --e{715}*/
   CHECK_OKAY( dialogExecMenu(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** dialog execution method for the display branching command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayBranching)
{  /*lint --e{715}*/
   BRANCHRULE** branchrules;
   BRANCHRULE** sorted;
   int nbranchrules;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   branchrules = SCIPgetBranchrules(scip);
   nbranchrules = SCIPgetNBranchrules(scip);

   /* copy branchrules array into temporary memory for sorting */
   CHECK_OKAY( SCIPduplicateBufferArray(scip, &sorted, branchrules, nbranchrules) );

   /* sort the branching rules */
   SCIPbsortPtr((void**)sorted, nbranchrules, SCIPbranchruleComp);

   /* display sorted list of branching rules */
   printf("\n");
   printf(" branching rule       priority maxdepth  description\n");
   printf(" --------------       -------- --------  -----------\n");
   for( i = 0; i < nbranchrules; ++i )
   {
      printf(" %-20s ", SCIPbranchruleGetName(sorted[i]));
      if( strlen(SCIPbranchruleGetName(sorted[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%8d %8d  ", SCIPbranchruleGetPriority(sorted[i]), SCIPbranchruleGetMaxdepth(sorted[i]));
      printf(SCIPbranchruleGetDesc(sorted[i]));
      printf("\n");
   }
   printf("\n");

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &sorted) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display conflict command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayConflict)
{  /*lint --e{715}*/
   CONFLICTHDLR** conflicthdlrs;
   CONFLICTHDLR** sorted;
   int nconflicthdlrs;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   conflicthdlrs = SCIPgetConflicthdlrs(scip);
   nconflicthdlrs = SCIPgetNConflicthdlrs(scip);

   /* copy conflicthdlrs array into temporary memory for sorting */
   CHECK_OKAY( SCIPduplicateBufferArray(scip, &sorted, conflicthdlrs, nconflicthdlrs) );

   /* sort the conflict handlers */
   SCIPbsortPtr((void**)sorted, nconflicthdlrs, SCIPconflicthdlrComp);

   /* display sorted list of conflict handlers */
   printf("\n");
   printf(" conflict handler     priority  description\n");
   printf(" ----------------     --------  -----------\n");
   for( i = 0; i < nconflicthdlrs; ++i )
   {
      printf(" %-20s ", SCIPconflicthdlrGetName(sorted[i]));
      if( strlen(SCIPconflicthdlrGetName(sorted[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%8d  ", SCIPconflicthdlrGetPriority(sorted[i]));
      printf(SCIPconflicthdlrGetDesc(sorted[i]));
      printf("\n");
   }
   printf("\n");

   /* free temporary memory */
   CHECK_OKAY( SCIPfreeBufferArray(scip, &sorted) );

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display conshdlrs command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayConshdlrs)
{  /*lint --e{715}*/
   CONSHDLR** conshdlrs;
   int nconshdlrs;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);

   /* display list of constraint handlers */
   printf("\n");
   printf(" constraint handler   chckprio enfoprio sepaprio sepaf propf eager  description\n");
   printf(" ------------------   -------- -------- -------- ----- ----- -----  -----------\n");
   for( i = 0; i < nconshdlrs; ++i )
   {
      printf(" %-20s ", SCIPconshdlrGetName(conshdlrs[i]));
      if( strlen(SCIPconshdlrGetName(conshdlrs[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%8d %8d %8d %5d %5d %5d  ",
         SCIPconshdlrGetCheckPriority(conshdlrs[i]),
         SCIPconshdlrGetEnfoPriority(conshdlrs[i]),
         SCIPconshdlrGetSepaPriority(conshdlrs[i]),
         SCIPconshdlrGetSepaFreq(conshdlrs[i]),
         SCIPconshdlrGetPropFreq(conshdlrs[i]),
         SCIPconshdlrGetEagerFreq(conshdlrs[i]));
      printf(SCIPconshdlrGetDesc(conshdlrs[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display displaycols command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayDisplaycols)
{  /*lint --e{715}*/
   DISP** disps;
   int ndisps;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   disps = SCIPgetDisps(scip);
   ndisps = SCIPgetNDisps(scip);

   /* display list of display columns */
   printf("\n");
   printf(" display column       header           position width priority status  description\n");
   printf(" --------------       ------           -------- ----- -------- ------  -----------\n");
   for( i = 0; i < ndisps; ++i )
   {
      printf(" %-20s ", SCIPdispGetName(disps[i]));
      if( strlen(SCIPdispGetName(disps[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%-16s ", SCIPdispGetHeader(disps[i]));
      if( strlen(SCIPdispGetHeader(disps[i])) > 16 )
         printf("\n %20s %16s ", "", "-->");
      printf("%8d ", SCIPdispGetPosition(disps[i]));
      printf("%5d ", SCIPdispGetWidth(disps[i]));
      printf("%8d ", SCIPdispGetPriority(disps[i]));
      switch( SCIPdispGetStatus(disps[i]) )
      {
      case SCIP_DISPSTATUS_OFF:
         printf("%6s  ", "off");
         break;
      case SCIP_DISPSTATUS_AUTO:
         printf("%6s  ", "auto");
         break;
      case SCIP_DISPSTATUS_ON:
         printf("%6s  ", "on");
         break;
      default:
         printf("%6s  ", "???");
         break;
      }
      printf(SCIPdispGetDesc(disps[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display heuristics command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayHeuristics)
{  /*lint --e{715}*/
   HEUR** heurs;
   int nheurs;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   heurs = SCIPgetHeurs(scip);
   nheurs = SCIPgetNHeurs(scip);

   /* display list of primal heuristics */
   printf("\n");
   printf(" primal heuristic     c priority freq ofs  description\n");
   printf(" ----------------     - -------- ---- ---  -----------\n");
   for( i = 0; i < nheurs; ++i )
   {
      printf(" %-20s ", SCIPheurGetName(heurs[i]));
      if( strlen(SCIPheurGetName(heurs[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%c ", SCIPheurGetDispchar(heurs[i]));
      printf("%8d ", SCIPheurGetPriority(heurs[i]));
      printf("%4d ", SCIPheurGetFreq(heurs[i]));
      printf("%3d  ", SCIPheurGetFreqofs(heurs[i]));
      printf(SCIPheurGetDesc(heurs[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display nodeselectors command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayNodeselectors)
{  /*lint --e{715}*/
   NODESEL** nodesels;
   int nnodesels;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   nodesels = SCIPgetNodesels(scip);
   nnodesels = SCIPgetNNodesels(scip);

   /* display list of node selectors */
   printf("\n");
   printf(" node selector        std priority memsave prio  description\n");
   printf(" -------------        ------------ ------------  -----------\n");
   for( i = 0; i < nnodesels; ++i )
   {
      printf(" %-20s ", SCIPnodeselGetName(nodesels[i]));
      if( strlen(SCIPnodeselGetName(nodesels[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%12d ", SCIPnodeselGetStdPriority(nodesels[i]));
      printf("%12d  ", SCIPnodeselGetMemsavePriority(nodesels[i]));
      printf(SCIPnodeselGetDesc(nodesels[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display presolvers command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayPresolvers)
{  /*lint --e{715}*/
   PRESOL** presols;
   int npresols;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   presols = SCIPgetPresols(scip);
   npresols = SCIPgetNPresols(scip);

   /* display list of presolvers */
   printf("\n");
   printf(" presolver            priority  description\n");
   printf(" ---------            --------  -----------\n");
   for( i = 0; i < npresols; ++i )
   {
      printf(" %-20s ", SCIPpresolGetName(presols[i]));
      if( strlen(SCIPpresolGetName(presols[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%8d  ", SCIPpresolGetPriority(presols[i]));
      printf(SCIPpresolGetDesc(presols[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display problem command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayProblem)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintOrigProblem(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display readers command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayReaders)
{  /*lint --e{715}*/
   READER** readers;
   int nreaders;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   readers = SCIPgetReaders(scip);
   nreaders = SCIPgetNReaders(scip);

   /* display list of readers */
   printf("\n");
   printf(" file reader          extension  description\n");
   printf(" -----------          ---------  -----------\n");
   for( i = 0; i < nreaders; ++i )
   {
      printf(" %-20s ", SCIPreaderGetName(readers[i]));
      if( strlen(SCIPreaderGetName(readers[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%9s  ", SCIPreaderGetExtension(readers[i]));
      printf(SCIPreaderGetDesc(readers[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display separators command */
DECL_DIALOGEXEC(SCIPdialogExecDisplaySeparators)
{  /*lint --e{715}*/
   SEPA** sepas;
   int nsepas;
   int i;

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   sepas = SCIPgetSepas(scip);
   nsepas = SCIPgetNSepas(scip);

   /* display list of separators */
   printf("\n");
   printf(" separator            priority freq  description\n");
   printf(" ---------            -------- ----  -----------\n");
   for( i = 0; i < nsepas; ++i )
   {
      printf(" %-20s ", SCIPsepaGetName(sepas[i]));
      if( strlen(SCIPsepaGetName(sepas[i])) > 20 )
         printf("\n %20s ", "-->");
      printf("%8d ", SCIPsepaGetPriority(sepas[i]));
      printf("%4d  ", SCIPsepaGetFreq(sepas[i]));
      printf(SCIPsepaGetDesc(sepas[i]));
      printf("\n");
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display solution command */
DECL_DIALOGEXEC(SCIPdialogExecDisplaySolution)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintBestSol(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display statistics command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintStatistics(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display transproblem command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayTransproblem)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintTransProblem(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the display varbranchstatistics command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayVarbranchstatistics)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintBranchingStatistics(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the help command */
DECL_DIALOGEXEC(SCIPdialogExecHelp)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPdialogDisplayMenu(SCIPdialogGetParent(dialog), scip) );
   printf("\n");

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the display transsolution command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayTranssolution)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintBestTransSol(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the free command */
DECL_DIALOGEXEC(SCIPdialogExecFree)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   CHECK_OKAY( SCIPfreeProb(scip) );

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the newstart command */
DECL_DIALOGEXEC(SCIPdialogExecNewstart)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   CHECK_OKAY( SCIPfreeSolve(scip) );

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the optimize command */
DECL_DIALOGEXEC(SCIPdialogExecOptimize)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   switch( SCIPstage(scip) )
   {
   case SCIP_STAGE_INIT:
      printf("no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPsolve(scip) );
      break;

   case SCIP_STAGE_SOLVED:
      printf("problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the presolve command */
DECL_DIALOGEXEC(SCIPdialogExecPresolve)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   switch( SCIPstage(scip) )
   {
   case SCIP_STAGE_INIT:
      printf("no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
      CHECK_OKAY( SCIPpresolve(scip) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      printf("problem is already presolved\n");
      break;

   case SCIP_STAGE_SOLVED:
      printf("problem is already solved\n");
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      errorMessage("invalid SCIP stage\n");
      return SCIP_INVALIDCALL;
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the quit command */
DECL_DIALOGEXEC(SCIPdialogExecQuit)
{  /*lint --e{715}*/
   printf("\n");

   *nextdialog = NULL;

   return SCIP_OKAY;
}

/** dialog execution method for the read command */
DECL_DIALOGEXEC(SCIPdialogExecRead)
{  /*lint --e{715}*/
   RETCODE retcode;
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      if( SCIPfileExists(filename) )
      {
         CHECK_OKAY( SCIPfreeProb(scip) );
         retcode = SCIPreadProb(scip, filename);
         if( retcode == SCIP_READERROR || retcode == SCIP_NOFILE || retcode == SCIP_PARSEERROR )
         {
            printf("error reading file <%s>\n", filename);
            CHECK_OKAY( SCIPfreeProb(scip) );
         }
         else
         {
            CHECK_OKAY( retcode );
         }
      }
      else
      {
         printf("file <%s> not found\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set load command */
DECL_DIALOGEXEC(SCIPdialogExecSetLoad)
{  /*lint --e{715}*/
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      if( SCIPfileExists(filename) )
      {
         CHECK_OKAY( SCIPreadParams(scip, filename) );
         printf("loaded parameter file <%s>\n", filename);
      }
      else
      {
         printf("file <%s> not found\n", filename);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
      }
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set save command */
DECL_DIALOGEXEC(SCIPdialogExecSetSave)
{  /*lint --e{715}*/
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      CHECK_OKAY( SCIPwriteParams(scip, filename, TRUE, FALSE) );
      printf("saved parameter file <%s>\n", filename);
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set diffsave command */
DECL_DIALOGEXEC(SCIPdialogExecSetDiffsave)
{  /*lint --e{715}*/
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      CHECK_OKAY( SCIPwriteParams(scip, filename, TRUE, TRUE) );
      printf("saved non-default parameter settings to file <%s>\n", filename);
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set parameter command */
DECL_DIALOGEXEC(SCIPdialogExecSetParam)
{  /*lint --e{715}*/
   RETCODE retcode;
   PARAM* param;
   char prompt[MAXSTRLEN];
   const char* valuestr;
   int intval;
   Longint longintval;
   Real realval;
   char charval;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* get the parameter to set */
   param = (PARAM*)SCIPdialogGetData(dialog);
   
   /* depending on the parameter type, request a user input */
   switch( SCIPparamGetType(param) )
   {
   case SCIP_PARAMTYPE_BOOL:
      snprintf(prompt, MAXSTRLEN, "current value: %s, new value (TRUE/FALSE): ",
         SCIPparamGetBool(param) ? "TRUE" : "FALSE");
      valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
      CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      switch( valuestr[0] )
      {
      case 'f':
      case 'F':
      case '0':
      case 'n':
      case 'N':
         CHECK_OKAY( SCIPparamSetBool(param, scip, FALSE) );
         printf("parameter <%s> set to FALSE\n", SCIPparamGetName(param));
         break;
      case 't':
      case 'T':
      case '1':
      case 'y':
      case 'Y':
         CHECK_OKAY( SCIPparamSetBool(param, scip, TRUE) );
         printf("parameter <%s> set to TRUE\n", SCIPparamGetName(param));
         break;
      default:
         printf("\ninvalid parameter value <%s>\n\n", valuestr);
         break;
      }
      break;

   case SCIP_PARAMTYPE_INT:
      snprintf(prompt, MAXSTRLEN, "current value: %d, new value [%d,%d]: ",
         SCIPparamGetInt(param), SCIPparamGetIntMin(param), SCIPparamGetIntMax(param));
      valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
      CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      if( sscanf(valuestr, "%d", &intval) != 1 )
      {
         printf("\ninvalid input <%s>\n\n", valuestr);
         return SCIP_OKAY;
      }
      retcode = SCIPparamSetInt(param, scip, intval);
      if( retcode != SCIP_PARAMETERWRONGVAL )
      {
         CHECK_OKAY( retcode );
      }
      printf("parameter <%s> set to %d\n", SCIPparamGetName(param), SCIPparamGetInt(param));
      break;

   case SCIP_PARAMTYPE_LONGINT:
      snprintf(prompt, MAXSTRLEN, "current value: %lld, new value [%lld,%lld]: ",
         SCIPparamGetLongint(param), SCIPparamGetLongintMin(param), SCIPparamGetLongintMax(param));
      valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
      CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      if( sscanf(valuestr, LONGINT_FORMAT, &longintval) != 1 )
      {
         printf("\ninvalid input <%s>\n\n", valuestr);
         return SCIP_OKAY;
      }
      retcode = SCIPparamSetLongint(param, scip, longintval);
      if( retcode != SCIP_PARAMETERWRONGVAL )
      {
         CHECK_OKAY( retcode );
      }
      printf("parameter <%s> set to %lld\n", SCIPparamGetName(param), SCIPparamGetLongint(param));
      break;

   case SCIP_PARAMTYPE_REAL:
      snprintf(prompt, MAXSTRLEN, "current value: %g, new value [%g,%g]: ",
         SCIPparamGetReal(param), SCIPparamGetRealMin(param), SCIPparamGetRealMax(param));
      valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
      CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      if( sscanf(valuestr, REAL_FORMAT, &realval) != 1 )
      {
         printf("\ninvalid input <%s>\n\n", valuestr);
         return SCIP_OKAY;
      }
      retcode = SCIPparamSetReal(param, scip, realval);
      if( retcode != SCIP_PARAMETERWRONGVAL )
      {
         CHECK_OKAY( retcode );
      }
      printf("parameter <%s> set to %g\n", SCIPparamGetName(param), SCIPparamGetReal(param));
      break;

   case SCIP_PARAMTYPE_CHAR:
      snprintf(prompt, MAXSTRLEN, "current value: <%c>, new value: ", SCIPparamGetChar(param));
      valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
      CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      if( sscanf(valuestr, "%c", &charval) != 1 )
      {
         printf("\ninvalid input <%s>\n\n", valuestr);
         return SCIP_OKAY;
      }
      retcode = SCIPparamSetChar(param, scip, charval);
      if( retcode != SCIP_PARAMETERWRONGVAL )
      {
         CHECK_OKAY( retcode );
      }
      printf("parameter <%s> set to <%c>\n", SCIPparamGetName(param), SCIPparamGetChar(param));
      break;

   case SCIP_PARAMTYPE_STRING:
      snprintf(prompt, MAXSTRLEN, "current value: <%s>, new value: ", SCIPparamGetString(param));
      valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
      CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
      if( valuestr[0] == '\0' )
         return SCIP_OKAY;

      retcode = SCIPparamSetString(param, scip, valuestr);
      if( retcode != SCIP_PARAMETERWRONGVAL )
      {
         CHECK_OKAY( retcode );
      }
      printf("parameter <%s> set to <%s>\n", SCIPparamGetName(param), SCIPparamGetString(param));
      break;

   default:
      errorMessage("invalid parameter type\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** dialog description method for the set parameter command */
DECL_DIALOGDESC(SCIPdialogDescSetParam)
{  /*lint --e{715}*/
   PARAM* param;
   char valuestr[MAXSTRLEN];

   /* get the parameter to set */
   param = (PARAM*)SCIPdialogGetData(dialog);
   
   /* retrieve parameter's current value */
   switch( SCIPparamGetType(param) )
   {
   case SCIP_PARAMTYPE_BOOL:
      if( SCIPparamGetBool(param) )
         sprintf(valuestr, "TRUE");
      else
         sprintf(valuestr, "FALSE");
      break;

   case SCIP_PARAMTYPE_INT:
      sprintf(valuestr, "%d", SCIPparamGetInt(param));
      break;

   case SCIP_PARAMTYPE_LONGINT:
      sprintf(valuestr, LONGINT_FORMAT, SCIPparamGetLongint(param));
      break;

   case SCIP_PARAMTYPE_REAL:
      sprintf(valuestr, "%g", SCIPparamGetReal(param));
      if( strchr(valuestr, '.') == NULL && strchr(valuestr, 'e') == NULL )
         sprintf(valuestr, "%.1f", SCIPparamGetReal(param));
      break;

   case SCIP_PARAMTYPE_CHAR:
      sprintf(valuestr, "%c", SCIPparamGetChar(param));
      break;

   case SCIP_PARAMTYPE_STRING:
      snprintf(valuestr, MAXSTRLEN, "%s", SCIPparamGetString(param));
      break;

   default:
      errorMessage("invalid parameter type\n");
      return SCIP_INVALIDDATA;
   }
   valuestr[MAXSTRLEN-1] = '\0';

   /* display parameter's description */
   printf(SCIPparamGetDesc(param));

   /* display parameter's current value */
   printf(" [%s]", valuestr);
   
   return SCIP_OKAY;
}

/** dialog execution method for the set limits objective command */
DECL_DIALOGEXEC(SCIPdialogExecSetLimitsObjective)
{  /*lint --e{715}*/
   char prompt[MAXSTRLEN];
   const char* valuestr;
   Real objlim;

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   /* objective limit cannot be set, if no problem was created */
   if( SCIPstage(scip) == SCIP_STAGE_INIT )
   {
      printf("cannot set objective limit before problem was created\n");
      return SCIP_OKAY;
   }

   /* get new objective limit from user */
   snprintf(prompt, MAXSTRLEN, "current value: %g, new value: ", SCIPgetObjlimit(scip));
   valuestr = SCIPdialoghdlrGetWord(dialoghdlr, dialog, prompt);
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, valuestr) );
   if( valuestr[0] == '\0' )
      return SCIP_OKAY;

   if( sscanf(valuestr, REAL_FORMAT, &objlim) != 1 )
   {
      printf("\ninvalid input <%s>\n\n", valuestr);
      return SCIP_OKAY;
   }
   
   /* check, if new objective limit is valid */
   if( SCIPstage(scip) > SCIP_STAGE_PROBLEM
      && SCIPtransformObj(scip, objlim) > SCIPtransformObj(scip, SCIPgetObjlimit(scip)) )
   {
      printf("\ncannot relax objective limit from %g to %g after problem was transformed\n\n",
         SCIPgetObjlimit(scip), objlim);
      return SCIP_OKAY;
   }

   /* set new objective limit */
   CHECK_OKAY( SCIPsetObjlimit(scip, objlim) );
   printf("objective value limit set to %g\n", SCIPgetObjlimit(scip));

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** dialog execution method for the debug memory command */
DECL_DIALOGEXEC(SCIPdialogExecDebugMemory)
{  /*lint --e{715}*/
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   SCIPdebugMemory(scip);
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}
#endif

/** includes or updates the default dialog menus in SCIP */
RETCODE SCIPincludeDialogDefault(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   DIALOG* root;
   DIALOG* submenu;
   DIALOG* dialog;

   /* root menu */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &root, SCIPdialogExecMenuLazy, NULL,
                     "SCIP", "SCIP's main menu", TRUE, NULL) );
      CHECK_OKAY( SCIPsetRootDialog(scip, root) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &root) );
      root = SCIPgetRootDialog(scip);
   }

   /* display */
   if( !SCIPdialogHasEntry(root, "display") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "display", "display information", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "display", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;
   
   /* display branching */
   if( !SCIPdialogHasEntry(submenu, "branching") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayBranching, NULL,
                     "branching", "display branching rule priorities", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display conflict */
   if( !SCIPdialogHasEntry(submenu, "conflict") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayConflict, NULL,
                     "conflict", "display conflict handler priorities", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display conshdlrs */
   if( !SCIPdialogHasEntry(submenu, "conshdlrs") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayConshdlrs, NULL,
                     "conshdlrs", "display constraint handler settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display displaycols */
   if( !SCIPdialogHasEntry(submenu, "displaycols") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayDisplaycols, NULL,
                     "displaycols", "display display column settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display heuristics */
   if( !SCIPdialogHasEntry(submenu, "heuristics") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayHeuristics, NULL,
                     "heuristics", "display primal heuristics settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display nodeselectors */
   if( !SCIPdialogHasEntry(submenu, "nodeselectors") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayNodeselectors, NULL,
                     "nodeselectors", "display node selectors settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display presolvers */
   if( !SCIPdialogHasEntry(submenu, "presolvers") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayPresolvers, NULL,
                     "presolvers", "display presolvers settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display problem */
   if( !SCIPdialogHasEntry(submenu, "problem") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayProblem, NULL,
                     "problem", "display original problem", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display readers */
   if( !SCIPdialogHasEntry(submenu, "readers") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayReaders, NULL,
                     "readers", "display file readers settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display separators */
   if( !SCIPdialogHasEntry(submenu, "separators") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplaySeparators, NULL,
                     "separators", "display cut separators settings", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display solution */
   if( !SCIPdialogHasEntry(submenu, "solution") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplaySolution, NULL,
                     "solution", "display best primal solution", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display statistics */
   if( !SCIPdialogHasEntry(submenu, "statistics") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayStatistics, NULL,
                     "statistics", "display problem and optimization statistics", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display transproblem */
   if( !SCIPdialogHasEntry(submenu, "transproblem") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayTransproblem, NULL,
                     "transproblem", "display transformed/preprocessed problem", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display varbranchstatistics */
   if( !SCIPdialogHasEntry(submenu, "varbranchstatistics") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayVarbranchstatistics, NULL,
                     "varbranchstatistics", "display statistics for branching on variables", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* display transsolution */
   if( !SCIPdialogHasEntry(submenu, "transsolution") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayTranssolution, NULL,
                     "transsolution", "display best primal solution in transformed variables", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   
   /* free */
   if( !SCIPdialogHasEntry(root, "free") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecFree, NULL,
                     "free", "free current problem from memory", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* help */
   if( !SCIPdialogHasEntry(root, "help") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecHelp, NULL,
                     "help", "display this help", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* newstart */
   if( !SCIPdialogHasEntry(root, "newstart") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecNewstart, NULL,
                     "newstart", "reset branch and bound tree to start again from root", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* optimize */
   if( !SCIPdialogHasEntry(root, "optimize") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecOptimize, NULL,
                     "optimize", "solve the problem", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* presolve */
   if( !SCIPdialogHasEntry(root, "presolve") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecPresolve, NULL,
                     "presolve", "solve the problem, but stop after presolving stage", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* quit */
   if( !SCIPdialogHasEntry(root, "quit") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecQuit, NULL,
                     "quit", "leave SCIP", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* read */
   if( !SCIPdialogHasEntry(root, "read") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecRead, NULL,
                     "read", "read a problem", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set */
   CHECK_OKAY( SCIPincludeDialogDefaultSet(scip) );

#ifndef NDEBUG
   /* debug */
   if( !SCIPdialogHasEntry(root, "debug") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "debug", "debugging information", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(root, "debug", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;
   
   /* debug memory */
   if( !SCIPdialogHasEntry(submenu, "memory") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDebugMemory, NULL,
                     "memory", "display memory diagnostics", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
#endif

   return SCIP_OKAY;
}

/** if a '/' occurs in the parameter's name, adds a sub menu dialog to the given menu and inserts the parameter dialog
 *  recursively in the sub menu; if no '/' occurs in the name, adds a parameter change dialog into the given dialog menu
 */
static
RETCODE addParamDialog(
   SCIP*            scip,               /**< SCIP data structure */
   DIALOG*          menu,               /**< dialog menu to insert the parameter into */
   PARAM*           param,              /**< parameter to add a dialog for */
   char*            paramname           /**< parameter name to parse */
   )
{
   char* slash;
   char* dirname;

   assert(paramname != NULL);

   /* check for a '/' */
   slash = strchr(paramname, '/');

   if( slash == NULL )
   {
      /* check, if the corresponding dialog already exists */
      if( !SCIPdialogHasEntry(menu, paramname) )
      {
         DIALOG* paramdialog;

         /* create a parameter change dialog */
         CHECK_OKAY( SCIPcreateDialog(scip, &paramdialog, SCIPdialogExecSetParam, SCIPdialogDescSetParam, 
                        paramname, SCIPparamGetDesc(param), FALSE, (DIALOGDATA*)param) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, menu, paramdialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &paramdialog) );
      }
   }
   else
   {
      DIALOG* submenu;

      /* split the parameter name into dirname and parameter name */
      dirname = paramname;
      paramname = slash+1;
      *slash = '\0';

      /* if not yet existing, create a corresponding sub menu */
      if( !SCIPdialogHasEntry(menu, dirname) )
      {
         char desc[MAXSTRLEN];

         sprintf(desc, "parameters for <%s>", dirname);
         CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL, dirname, desc, TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
      }

      /* find the corresponding sub menu */
      (void)SCIPdialogFindEntry(menu, dirname, &submenu);
      if( submenu == NULL )
      {
         errorMessage("dialog sub menu not found\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* recursively call add parameter method */
      CHECK_OKAY( addParamDialog(scip, submenu, param, paramname) );
   }

   return SCIP_OKAY;
}

/** includes or updates the "set" menu for each available parameter setting */
RETCODE SCIPincludeDialogDefaultSet(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   DIALOG* root;
   DIALOG* setmenu;
   DIALOG* submenu;
   DIALOG* dialog;
   PARAM** params;
   const char* pname;
   char* paramname;
   int nparams;
   int i;

   /* get root dialog */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      errorMessage("root dialog not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* find (or create) the "set" menu of the root dialog */
   if( !SCIPdialogHasEntry(root, "set") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &setmenu, SCIPdialogExecMenu, NULL,
                     "set", "load/save/change parameters", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, root, setmenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &setmenu) );
   }
   if( SCIPdialogFindEntry(root, "set", &setmenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   /* set load */
   if( !SCIPdialogHasEntry(setmenu, "load") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecSetLoad, NULL,
                     "load", "load parameter settings from a file", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set save */
   if( !SCIPdialogHasEntry(setmenu, "save") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecSetSave, NULL,
                     "save", "save parameter settings to a file", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set diffsave */
   if( !SCIPdialogHasEntry(setmenu, "diffsave") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecSetDiffsave, NULL,
                     "diffsave", "save non-default parameter settings to a file", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }

   /* set branching */
   if( !SCIPdialogHasEntry(setmenu, "branching") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "branching", "change parameters for branching rules", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "branching", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNBranchrules(scip); ++i )
   {
      BRANCHRULE* branchrule = SCIPgetBranchrules(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPbranchruleGetName(branchrule)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPbranchruleGetName(branchrule), SCIPbranchruleGetDesc(branchrule), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set conflict */
   if( !SCIPdialogHasEntry(setmenu, "conflict") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "conflict", "change parameters for conflict handlers", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "conflict", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNConflicthdlrs(scip); ++i )
   {
      CONFLICTHDLR* conflicthdlr = SCIPgetConflicthdlrs(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPconflicthdlrGetName(conflicthdlr)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPconflicthdlrGetName(conflicthdlr), SCIPconflicthdlrGetDesc(conflicthdlr), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set constraints */
   if( !SCIPdialogHasEntry(setmenu, "constraints") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "constraints", "change parameters for constraint handlers", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "constraints", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNConshdlrs(scip); ++i )
   {
      CONSHDLR* conshdlr = SCIPgetConshdlrs(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPconshdlrGetName(conshdlr)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPconshdlrGetName(conshdlr), SCIPconshdlrGetDesc(conshdlr), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set display */
   if( !SCIPdialogHasEntry(setmenu, "display") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "display", "change parameters for display columns", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "display", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNDisps(scip); ++i )
   {
      DISP* disp = SCIPgetDisps(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPdispGetName(disp)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPdispGetName(disp), SCIPdispGetDesc(disp), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set heuristics */
   if( !SCIPdialogHasEntry(setmenu, "heuristics") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "heuristics", "change parameters for primal heuristics", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "heuristics", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNHeurs(scip); ++i )
   {
      HEUR* heur = SCIPgetHeurs(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPheurGetName(heur)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPheurGetName(heur), SCIPheurGetDesc(heur), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set limits */
   if( !SCIPdialogHasEntry(setmenu, "limits") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "limits", "change parameters for time, memory, objective value, and other limits", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecSetLimitsObjective, NULL,
                     "objective", "set limit on objective value", FALSE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set lp */
   if( !SCIPdialogHasEntry(setmenu, "lp") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "lp", "change parameters for linear programming relaxations", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set memory */
   if( !SCIPdialogHasEntry(setmenu, "memory") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "memory", "change parameters for memory management", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set nodeselection */
   if( !SCIPdialogHasEntry(setmenu, "nodeselection") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "nodeselection", "change parameters for node selectors", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "nodeselection", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNNodesels(scip); ++i )
   {
      NODESEL* nodesel = SCIPgetNodesels(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPnodeselGetName(nodesel)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPnodeselGetName(nodesel), SCIPnodeselGetDesc(nodesel), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set numerics */
   if( !SCIPdialogHasEntry(setmenu, "numerics") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "numerics", "change parameters for numerical values", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set presolving */
   if( !SCIPdialogHasEntry(setmenu, "presolving") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "presolving", "change parameters for presolving", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "presolving", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNPresols(scip); ++i )
   {
      PRESOL* presol = SCIPgetPresols(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPpresolGetName(presol)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPpresolGetName(presol), SCIPpresolGetDesc(presol), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set pricing */
   if( !SCIPdialogHasEntry(setmenu, "pricing") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "pricing", "change parameters for pricing variables", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "pricing", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNPricers(scip); ++i )
   {
      PRICER* pricer = SCIPgetPricers(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPpricerGetName(pricer)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPpricerGetName(pricer), SCIPpricerGetDesc(pricer), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set reading */
   if( !SCIPdialogHasEntry(setmenu, "reading") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "reading", "change parameters for problem file readers", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "reading", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNReaders(scip); ++i )
   {
      READER* reader = SCIPgetReaders(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPreaderGetName(reader)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPreaderGetName(reader), SCIPreaderGetDesc(reader), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set separating */
   if( !SCIPdialogHasEntry(setmenu, "separating") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "separating", "change parameters for cut separators", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }
   if( SCIPdialogFindEntry(setmenu, "separating", &submenu) != 1 )
      return SCIP_PLUGINNOTFOUND;

   for( i = 0; i < SCIPgetNSepas(scip); ++i )
   {
      SEPA* sepa = SCIPgetSepas(scip)[i];

      if( !SCIPdialogHasEntry(submenu, SCIPsepaGetName(sepa)) )
      {
         CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                        SCIPsepaGetName(sepa), SCIPsepaGetDesc(sepa), TRUE, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
      }
   }

   /* set timing */
   if( !SCIPdialogHasEntry(setmenu, "timing") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "timing", "change parameters for timing issues", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }

   /* set misc */
   if( !SCIPdialogHasEntry(setmenu, "misc") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                     "misc", "change parameters for miscellaneous stuff", TRUE, NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, setmenu, submenu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
   }

   /* get SCIP's parameters */
   params = SCIPgetParams(scip);
   nparams = SCIPgetNParams(scip);

   /* insert each parameter into the set menu */
   for( i = 0; i < nparams; ++i )
   {
      pname = SCIPparamGetName(params[i]);
      ALLOC_OKAY( duplicateMemoryArray(&paramname, pname, strlen(pname)+1) );
      CHECK_OKAY( addParamDialog(scip, setmenu, params[i], paramname) );
      freeMemoryArray(&paramname);
   }

   return SCIP_OKAY;
}

