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
#pragma ident "@(#) $Id: dialog.h,v 1.11 2005/02/08 16:13:23 bzfpfend Exp $"

/**@file   dialog.h
 * @brief  internal methods for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DIALOG_H__
#define __DIALOG_H__


#include "def.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_dialog.h"
#include "pub_dialog.h"




/*
 * dialog handler
 */

/** creates a dialog handler */
extern
RETCODE SCIPdialoghdlrCreate(
   DIALOGHDLR**     dialoghdlr          /**< pointer to store dialog handler */
   );

/** frees a dialog handler and it's dialog tree */
extern
RETCODE SCIPdialoghdlrFree(
   DIALOGHDLR**     dialoghdlr          /**< pointer to dialog handler */
   );

/** executes the root dialog of the dialog handler */
extern
RETCODE SCIPdialoghdlrExec(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   SET*             set                 /**< global SCIP settings */
   );

/** makes given dialog the root dialog of dialog handler; captures dialog and releases former root dialog */
extern
RETCODE SCIPdialoghdlrSetRoot(
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG*          dialog              /**< dialog to be the root */
   );




/*
 * dialog
 */

/** creates and captures a user interface dialog */
extern
RETCODE SCIPdialogCreate(
   DIALOG**         dialog,             /**< pointer to store the dialog */
   DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   const char*      name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*      desc,               /**< description of dialog used if description output method is NULL */
   Bool             issubmenu,          /**< is the dialog a submenu? */
   DIALOGDATA*      dialogdata          /**< user defined dialog data */
   );

/** captures a dialog */
extern
void SCIPdialogCapture(
   DIALOG*          dialog              /**< dialog */
   );

/** releases a dialog */
extern
RETCODE SCIPdialogRelease(
   DIALOG**         dialog              /**< pointer to dialog */
   );

/** executes dialog */
extern
RETCODE SCIPdialogExec(
   DIALOG*          dialog,             /**< dialog */
   SET*             set,                /**< global SCIP settings */
   DIALOGHDLR*      dialoghdlr,         /**< dialog handler */
   DIALOG**         nextdialog          /**< pointer to store the next dialog to process */
   );

/** adds a sub dialog to the given dialog as menu entry and captures the sub dialog */
extern
RETCODE SCIPdialogAddEntry(
   DIALOG*          dialog,             /**< dialog */
   SET*             set,                /**< global SCIP settings */
   DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   );


#endif
