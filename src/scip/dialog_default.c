/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
{
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
{
   CHECK_OKAY( dialogExecMenu(scip, dialog, dialoghdlr, nextdialog) );

   return SCIP_OKAY;
}

/** dialog execution method for the display statistics command */
DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics)
{
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPprintStatistics(scip, NULL) );
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the help command */
DECL_DIALOGEXEC(SCIPdialogExecHelp)
{
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   CHECK_OKAY( SCIPdialogDisplayMenu(SCIPdialogGetParent(dialog), scip) );
   printf("\n");

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the free command */
DECL_DIALOGEXEC(SCIPdialogExecFree)
{
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   CHECK_OKAY( SCIPfreeProb(scip) );

   *nextdialog = SCIPdialogGetParent(dialog);

   return SCIP_OKAY;
}

/** dialog execution method for the optimize command */
DECL_DIALOGEXEC(SCIPdialogExecOptimize)
{
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   switch( SCIPstage(scip) )
   {
   case SCIP_STAGE_INIT:
      printf("no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      CHECK_OKAY( SCIPsolve(scip) );
      break;

   case SCIP_STAGE_SOLVED:
      printf("problem is already solved\n");
      break;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_FREESOLVE:
   default:
      errorMessage("invalid SCIP stage");
      return SCIP_INVALIDCALL;
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the presolve command */
DECL_DIALOGEXEC(SCIPdialogExecPresolve)
{
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   switch( SCIPstage(scip) )
   {
   case SCIP_STAGE_INIT:
      printf("no problem exists\n");
      break;

   case SCIP_STAGE_PROBLEM:
      CHECK_OKAY( SCIPpresolve(scip) );
      break;

   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      printf("problem is already presolved\n");
      break;

   case SCIP_STAGE_SOLVED:
      printf("problem is already solved\n");
      break;

   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_FREESOLVE:
   default:
      errorMessage("invalid SCIP stage");
      return SCIP_INVALIDCALL;
   }
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the quit command */
DECL_DIALOGEXEC(SCIPdialogExecQuit)
{
   *nextdialog = NULL;

   return SCIP_OKAY;
}

/** dialog execution method for the read command */
DECL_DIALOGEXEC(SCIPdialogExecRead)
{
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      if( SCIPfileExists(filename) )
      {
         CHECK_OKAY( SCIPfreeProb(scip) );
         CHECK_OKAY( SCIPreadProb(scip, filename) );
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

/** dialog execution method for the set read command */
DECL_DIALOGEXEC(SCIPdialogExecSetRead)
{
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      if( SCIPfileExists(filename) )
      {
         CHECK_OKAY( SCIPreadParams(scip, filename) );
         printf("read parameter file <%s>\n", filename);
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

/** dialog execution method for the set write command */
DECL_DIALOGEXEC(SCIPdialogExecSetWrite)
{
   const char* filename;

   filename = SCIPdialoghdlrGetWord(dialoghdlr, dialog, "enter filename: ");

   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, filename) );

   if( filename[0] != '\0' )
   {
      CHECK_OKAY( SCIPwriteParams(scip, filename, TRUE) );
      printf("written parameter file <%s>\n", filename);
   }

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}

/** dialog execution method for the set <parameter> command */
DECL_DIALOGEXEC(SCIPdialogExecSetParam)
{
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
      break;

   default:
      errorMessage("invalid parameter type");
      return SCIP_INVALIDDATA;
   }


   return SCIP_OKAY;
}

/** dialog description method for the set <parameter> command */
DECL_DIALOGDESC(SCIPdialogDescSetParam)
{
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
      errorMessage("invalid parameter type");
      return SCIP_INVALIDDATA;
   }
   valuestr[MAXSTRLEN-1] = '\0';

   /* display parameter's description */
   printf(SCIPparamGetDesc(param));

   /* display parameter's current value */
   printf(" [%s]", valuestr);
   
   return SCIP_OKAY;
}

#ifndef NDEBUG
/** dialog execution method for the debug memory command */
DECL_DIALOGEXEC(SCIPdialogExecDebugMemory)
{
   CHECK_OKAY( SCIPdialoghdlrAddHistory(dialoghdlr, dialog, NULL) );

   printf("\n");
   SCIPdebugMemory(scip);
   printf("\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}
#endif

/** includes the default dialog menus in SCIP */
RETCODE SCIPincludeDialogDefault(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   DIALOG* root;
   DIALOG* menu;
   DIALOG* submenu;
   DIALOG* dialog;
   int i;

   /* root menu */
   CHECK_OKAY( SCIPcreateDialog(scip, &root, SCIPdialogExecMenuLazy, NULL,
                  "SCIP", "SCIP's main menu", NULL) );
   CHECK_OKAY( SCIPsetRootDialog(scip, root) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &root) );
   
   /* display */
   CHECK_OKAY( SCIPcreateDialog(scip, &menu, SCIPdialogExecMenu, NULL,
                  "display", "display information", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, menu) );

   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDisplayStatistics, NULL,
                  "statistics", "display solution statistics", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );

   CHECK_OKAY( SCIPreleaseDialog(scip, &menu) );
   
   /* free */
   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecFree, NULL,
                  "free", "free current problem from memory", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   
   /* help */
   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecHelp, NULL,
                  "help", "display this help", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   
   /* optimize */
   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecOptimize, NULL,
                  "optimize", "solve the problem", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   
   /* presolve */
   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecPresolve, NULL,
                  "presolve", "solve the problem, but stop after presolving stage", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   
   /* quit */
   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecQuit, NULL,
                  "quit", "leave SCIP", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   
   /* read */
   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecRead, NULL,
                  "read", "read a problem", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );

   /* set */
   CHECK_OKAY( SCIPcreateDialog(scip, &menu, SCIPdialogExecMenu, NULL,
                  "set", "read/write/change parameters", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, menu) );

   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecSetRead, NULL,
                  "read", "read parameter settings from a file", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );

   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecSetWrite, NULL,
                  "write", "save parameter settings to a file", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );

   /* set branching */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "branching", "change parameters for branching rules", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNBranchrules(scip); ++i )
   {
      BRANCHRULE* branchrule = SCIPgetBranchrules(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPbranchruleGetName(branchrule), SCIPbranchruleGetDesc(branchrule), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set constraints */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "constraints", "change parameters for constraint handlers", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNConshdlrs(scip); ++i )
   {
      CONSHDLR* conshdlr = SCIPgetConshdlrs(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPconshdlrGetName(conshdlr), SCIPconshdlrGetDesc(conshdlr), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set display */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "display", "change parameters for display columns", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNDisps(scip); ++i )
   {
      DISP* disp = SCIPgetDisps(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPdispGetName(disp), SCIPdispGetDesc(disp), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set heuristics */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "heuristics", "change parameters for primal heuristics", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNHeurs(scip); ++i )
   {
      HEUR* heur = SCIPgetHeurs(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPheurGetName(heur), SCIPheurGetDesc(heur), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set limits */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "limits", "change parameters for time, memory, and other limits", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set lp */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "lp", "change parameters for linear programming relaxations", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set memory */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "memory", "change parameters for memory management", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set nodeselection */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "nodeselection", "change parameters for node selectors", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNNodesels(scip); ++i )
   {
      NODESEL* nodesel = SCIPgetNodesels(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPnodeselGetName(nodesel), SCIPnodeselGetDesc(nodesel), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set numerics */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "numerics", "change parameters for numerical values", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set presolving */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "presolving", "change parameters for presolving", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNPresols(scip); ++i )
   {
      PRESOL* presol = SCIPgetPresols(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPpresolGetName(presol), SCIPpresolGetDesc(presol), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set pricing */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "pricing", "change parameters for pricing variables", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   todoMessage("include pricers in standard set dialog");
#if 0
   for( i = 0; i < SCIPgetNPricers(scip); ++i )
   {
      PRICER* pricer = SCIPgetPricers(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPpricerGetName(pricer), SCIPpricerGetDesc(pricer), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
#endif
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set reading */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "reading", "change parameters for problem file readers", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNReaders(scip); ++i )
   {
      READER* reader = SCIPgetReaders(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPreaderGetName(reader), SCIPreaderGetDesc(reader), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set separating */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "separating", "change parameters for cut separators", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   for( i = 0; i < SCIPgetNSepas(scip); ++i )
   {
      SEPA* sepa = SCIPgetSepas(scip)[i];

      CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecMenu, NULL,
                     SCIPsepaGetName(sepa), SCIPsepaGetDesc(sepa), NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, submenu, dialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );
   }
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set timing */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "timing", "change parameters for timing issues", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   /* set misc */
   CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL,
                  "misc", "change parameters for miscellaneous stuff", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );

   CHECK_OKAY( SCIPreleaseDialog(scip, &menu) );

#ifndef NDEBUG
   /* debug */
   CHECK_OKAY( SCIPcreateDialog(scip, &menu, SCIPdialogExecMenu, NULL,
                  "debug", "debugging information", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, menu) );

   CHECK_OKAY( SCIPcreateDialog(scip, &dialog, SCIPdialogExecDebugMemory, NULL,
                  "memory", "display memory diagnostics", NULL) );
   CHECK_OKAY( SCIPaddDialogEntry(scip, menu, dialog) );
   CHECK_OKAY( SCIPreleaseDialog(scip, &dialog) );

   CHECK_OKAY( SCIPreleaseDialog(scip, &menu) );
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
      DIALOG* paramdialog;

      /* create a parameter change dialog */
      CHECK_OKAY( SCIPcreateDialog(scip, &paramdialog, SCIPdialogExecSetParam, SCIPdialogDescSetParam, 
                     paramname, SCIPparamGetDesc(param), (DIALOGDATA*)param) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, menu, paramdialog) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &paramdialog) );
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
         CHECK_OKAY( SCIPcreateDialog(scip, &submenu, SCIPdialogExecMenu, NULL, dirname, desc, NULL) );
         CHECK_OKAY( SCIPaddDialogEntry(scip, menu, submenu) );
         CHECK_OKAY( SCIPreleaseDialog(scip, &submenu) );
      }

      /* find the corresponding sub menu */
      (void)SCIPdialogFindEntry(menu, dirname, &submenu);
      if( submenu == NULL )
      {
         errorMessage("dialog sub menu not found");
         return SCIP_PLUGINNOTFOUND;
      }

      /* recursively call add parameter method */
      CHECK_OKAY( addParamDialog(scip, submenu, param, paramname) );
   }

   return SCIP_OKAY;
}

/** includes a dialog in the "set" menu for each available parameter setting */
RETCODE SCIPincludeDialogParams(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   DIALOG* rootmenu;
   DIALOG* setmenu;
   DIALOG* dialog;
   DIALOG* menu;
   PARAM** params;
   const char* pname;
   char* paramname;
   int nparams;
   int i;

   /* find (or create) the "set" menu of the root dialog */
   rootmenu = SCIPgetRootDialog(scip);
   if( rootmenu == NULL )
   {
      errorMessage("root dialog not found");
      return SCIP_PLUGINNOTFOUND;
   }

   if( !SCIPdialogHasEntry(rootmenu, "set") )
   {
      CHECK_OKAY( SCIPcreateDialog(scip, &setmenu, SCIPdialogExecMenu, NULL,
                     "set", "read/write/change parameters", NULL) );
      CHECK_OKAY( SCIPaddDialogEntry(scip, NULL, menu) );
      CHECK_OKAY( SCIPreleaseDialog(scip, &setmenu) );
   }

   (void)SCIPdialogFindEntry(rootmenu, "set", &setmenu);
   if( setmenu == NULL )
   {
      errorMessage("set menu not found");
      return SCIP_PLUGINNOTFOUND;
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

