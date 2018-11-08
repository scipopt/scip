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

/**@file   dialog_stp.c
 * @brief  stp user interface dialog
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "dialog_stp.h"
#include "probdata_stp.h"


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
   assert(root != NULL);

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
