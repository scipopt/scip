/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog_stp.c
 * @brief  stp user interface dialog
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_stp.h"

/** executes a menu dialog */
static
SCIP_RETCODE dialogExecMenu(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog menu */
   SCIP_DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SCIP_DIALOG**         nextdialog          /**< pointer to store next dialog to execute */
   )
{
   char* command;
   SCIP_Bool again;
   SCIP_Bool endoffile;
   int nfound;

   do
   {
      again = FALSE;

      /* get the next word of the command string */
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, NULL, &command, &endoffile) );
      if( endoffile )
      {
         *nextdialog = NULL;
         return SCIP_OKAY;
      }

      /* exit to the root dialog, if command is empty */
      if( command[0] == '\0' )
      {
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
         SCIPdialogMessage(scip, NULL, "command <%s> not available\n", command);
         SCIPdialoghdlrClearBuffer(dialoghdlr);
         *nextdialog = dialog;
      }
      else if( nfound >= 2 )
      {
         SCIPdialogMessage(scip, NULL, "\npossible completions:\n");
         SCIP_CALL( SCIPdialogDisplayCompletions(dialog, scip, command) );
         SCIPdialogMessage(scip, NULL, "\n");
         SCIPdialoghdlrClearBuffer(dialoghdlr);
         again = TRUE;
      }
   }
   while( again );

   return SCIP_OKAY;
}


/* parse the given string to detect a Boolean value and returns it */
static
SCIP_Bool parseBoolValue(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           valuestr,           /**< string to parse */
   SCIP_Bool*            error               /**< pointer to store the error result */
   )
{
   assert( scip  != NULL );
   assert( valuestr != NULL );
   assert( error != NULL );

   *error = FALSE;

   switch( valuestr[0] )
   {
   case 'f':
   case 'F':
   case '0':
   case 'n':
   case 'N':
      return FALSE;
   case 't':
   case 'T':
   case '1':
   case 'y':
   case 'Y':
      return TRUE;
   stp:
      SCIPdialogMessage(scip, NULL, "\ninvalid parameter value <%s>\n\n", valuestr);
      *error = TRUE;
      break;
   }

   return FALSE;
}


/** dialog execution method for the write LP command */
static
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteStpsol)
{  /*lint --e{715}*/
   SCIPdialogMessage(scip, NULL, "\n");

   SCIP_CALL( SCIPprobdataWriteLogfileEnd(scip) );

   SCIPdialogMessage(scip, NULL, "written STP solution\n");

   *nextdialog = SCIPdialoghdlrGetRoot(dialoghdlr);

   return SCIP_OKAY;
}



/** includes or updates the stp dialog menus in SCIP */
SCIP_RETCODE SCIPincludeDialogStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DIALOG* root;
   SCIP_DIALOG* submenu;
   SCIP_DIALOG* dialog;

   /* root menu */
   root = SCIPgetRootDialog(scip);
   if( root == NULL )
   {
      SCIP_CALL( SCIPcreateRootDialog(scip, &root) );
   }

   /* write */
   assert(SCIPdialogHasEntry(root, "write"));
   if( SCIPdialogFindEntry(root, "write", &submenu) != 1 )
   {
      SCIPerrorMessage("write sub menu not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* write STP solution */
   if( !SCIPdialogHasEntry(submenu, "stpsolution") )
   {
      SCIP_CALL( SCIPincludeDialog(scip, &dialog,
            NULL,
            SCIPdialogExecWriteStpsol, NULL, NULL,
            "stpsolution", "write solution to STP logfile", FALSE, NULL) );
      SCIP_CALL( SCIPaddDialogEntry(scip, submenu, dialog) );
      SCIP_CALL( SCIPreleaseDialog(scip, &dialog) );
   }

   return SCIP_OKAY;
}
