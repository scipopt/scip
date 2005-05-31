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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_dialog.h,v 1.6 2005/05/31 17:20:24 bzfpfend Exp $"

/**@file   type_dialog.h
 * @brief  type definitions for user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_DIALOG_H__
#define __TYPE_DIALOG_H__


typedef struct Dialog DIALOG;           /**< user interface dialog */
typedef struct DialogData DIALOGDATA;   /**< user defined dialog data */
typedef struct Dialoghdlr DIALOGHDLR;   /**< dialog handler */


/** execution method of dialog
 *
 *  This method is invoked, if the user selected the dialog's command name in the parent's dialog menu.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - dialoghdlr      : dialog handler to call for user interaction
 *  - dialog          : the dialog itself
 *
 *  output:
 *  - *nextdialog     : next dialog to process (or NULL to quit dialog processing)
 */
#define DECL_DIALOGEXEC(x) RETCODE x (SCIP* scip, DIALOG* dialog, DIALOGHDLR* dialoghdlr, DIALOG** nextdialog)

/** description output method of dialog
 *
 *  This method should output (usually a single line of) information describing the meaning of the dialog.
 *  The method is called, when the help menu of the parent's dialog is displayed.
 *  If no description output method is given, the description string of the dialog is displayed instead.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - *dialog         : the dialog itself
 */
#define DECL_DIALOGDESC(x) RETCODE x (SCIP* scip, DIALOG* dialog)



#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"


#endif
