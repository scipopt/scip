/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_dialog.h,v 1.4 2005/01/21 09:17:02 bzfpfend Exp $"

/**@file   pub_dialog.h
 * @brief  public methods for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_DIALOG_H__
#define __PUB_DIALOG_H__


#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"
#include "type_dialog.h"



/*
 * dialog handler
 */

/** returns the root dialog of the dialog handler */
extern
DIALOG* SCIPdialoghdlrGetRoot(
   DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** clears the input command buffer of the dialog handler */
extern
void SCIPdialoghdlrClearBuffer(
   DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** returns TRUE iff input command buffer is empty */
extern
Bool SCIPdialoghdlrIsBufferEmpty(
   DIALOGHDLR*      dialoghdlr          /**< dialog handler */
   );

/** returns the next word in the handler's command buffer; if the buffer is empty, displays the given prompt or the 
 *  current dialog's path and asks the user for further input; the user must not free or modify the returned string
 */
extern
const char* SCIPdialoghdlrGetWord(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG*          dialog,             /**< current dialog */
   const char*      prompt              /**< prompt to display, or NULL to display the current dialog's path */
   );

/** adds a command to the command history of the dialog handler; if a dialog is given, the command is preceeded
 *  by the dialog's command path; if no command is given, only the path to the dialog is added to the command history
 */
extern
RETCODE SCIPdialoghdlrAddHistory(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG*          dialog,             /**< current dialog, or NULL */
   const char*      command             /**< command string to add to the command history, or NULL */
   );




/*
 * dialog
 */

/** returns TRUE iff a dialog entry matching exactly the given name is existing in the given dialog */
extern
Bool SCIPdialogHasEntry(
   DIALOG*          dialog,             /**< dialog */
   const char*      entryname           /**< name of the dialog entry to find */
   );

/** searches the dialog for entries corresponding to the given name;
 *  If a complete match is found, the entry is returned as "subdialog" and
 *  the return value is 1.
 *  If no dialog entry completely matches the given "entryname", the number
 *  of entries with names beginning with "entryname" is returned. If this
 *  number is 1, the single match is returned as "subdialog". Otherwise,
 *  "subdialog" is set to NULL.
 */
extern
int SCIPdialogFindEntry(
   DIALOG*          dialog,             /**< dialog */
   const char*      entryname,          /**< name of the dialog entry to find */
   DIALOG**         subdialog           /**< pointer to store the found dialog entry */
   );

/** displays the dialog's menu */
extern
RETCODE SCIPdialogDisplayMenu(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** displays the entry for the dialog in it's parent's menu */
extern
RETCODE SCIPdialogDisplayMenuEntry(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** displays all dialog entries with names starting with the given "entryname" */
extern
RETCODE SCIPdialogDisplayCompletions(
   DIALOG*          dialog,             /**< dialog */
   SCIP*            scip,               /**< SCIP data structure */   
   const char*      entryname           /**< name of the dialog entry to find */
   );

/** gets the name of the current path in the dialog tree, separated by the given character */
extern
void SCIPdialogGetPath(
   DIALOG*          dialog,             /**< dialog */
   const char       sepchar,            /**< separation character to insert in path */
   char*            path                /**< string buffer to store the path */
   );

/** gets the command name of the dialog */
extern
const char* SCIPdialogGetName(
   DIALOG*          dialog              /**< dialog */
   );

/** gets the description of the dialog */
extern
const char* SCIPdialogGetDesc(
   DIALOG*          dialog              /**< dialog */
   );

/** returns whether the dialog is a sub menu */
extern
Bool SCIPdialogIsSubmenu(
   DIALOG*          dialog              /**< dialog */
   );

/** gets the parent dialog of the given dialog */
extern
DIALOG* SCIPdialogGetParent(
   DIALOG*          dialog              /**< dialog */
   );

/** gets the array of subdialogs associated with the given dialog */
extern
DIALOG** SCIPdialogGetSubdialogs(
   DIALOG*          dialog              /**< dialog */
   );

/** gets the number of subdialogs associated with the given dialog */
extern
int SCIPdialogGetNSubdialogs(
   DIALOG*          dialog              /**< dialog */
   );

/** gets the user defined data associated with the given dialog */
extern
DIALOGDATA* SCIPdialogGetData(
   DIALOG*          dialog              /**< dialog */
   );


#endif
