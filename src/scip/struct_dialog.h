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
#pragma ident "@(#) $Id: struct_dialog.h,v 1.6 2005/01/25 12:46:22 bzfpfend Exp $"

/**@file   struct_dialog.h
 * @brief  datastructures for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_DIALOG_H__
#define __STRUCT_DIALOG_H__


#include "def.h"
#include "type_dialog.h"



/** user interface dialog */
struct Dialog
{
   DECL_DIALOGEXEC  ((*dialogexec));    /**< execution method of dialog */
   DECL_DIALOGDESC  ((*dialogdesc));    /**< description output method of dialog, or NULL */
   char*            name;               /**< name of dialog: command name appearing in parent's dialog menu */
   char*            desc;               /**< description of dialog used if description output method is NULL */
   DIALOG*          parent;             /**< parent dialog of dialog */
   DIALOG**         subdialogs;         /**< sub dialogs of dialog */
   DIALOGDATA*      dialogdata;         /**< user defined dialog data */
   int              nsubdialogs;        /**< number of sub dialogs */
   int              subdialogssize;     /**< size of subdialogs array */
   int              nuses;              /**< number of times, the dialog is used */
   Bool             issubmenu;          /**< is the dialog a submenu? */
};

/** dialog handler */
struct Dialoghdlr
{
   DIALOG*          rootdialog;         /**< main (root) dialog */
   char*            buffer;             /**< command buffer */
   int              buffersize;         /**< size of command buffer */
   int              bufferpos;          /**< position of first unprocessed character in buffer */
   int              nprotectedhistelems;/**< number of history entries protected from cleaning up */
};


#endif
