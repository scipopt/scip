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

/**@file   dialog_default.h
 * @brief  default user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DIALOG_DEFAULT_H__
#define __DIALOG_DEFAULT_H__


#include "scip.h"



/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
extern
DECL_DIALOGEXEC(SCIPdialogExecMenu);

/** standard menu dialog execution method, that doesn't display it's help screen */
extern
DECL_DIALOGEXEC(SCIPdialogExecMenuLazy);

/** dialog execution method for the help command */
extern
DECL_DIALOGEXEC(SCIPdialogExecHelp);

/** dialog execution method for the quit command */
extern
DECL_DIALOGEXEC(SCIPdialogExecQuit);

/** includes the default dialog menus in SCIP */
extern
RETCODE SCIPincludeDialogDefault(
   SCIP*            scip                /**< SCIP data structure */
   );

/** includes a dialog in the "set" menu for each available parameter setting */
extern
RETCODE SCIPincludeDialogParams(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
