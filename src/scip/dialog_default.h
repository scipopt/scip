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

/** dialog execution method for the display statistics command */
extern
DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics);

/** dialog execution method for the help command */
extern
DECL_DIALOGEXEC(SCIPdialogExecHelp);

/** dialog execution method for the free command */
extern
DECL_DIALOGEXEC(SCIPdialogExecFree);

/** dialog execution method for the optimize command */
extern
DECL_DIALOGEXEC(SCIPdialogExecOptimize);

/** dialog execution method for the presolve command */
extern
DECL_DIALOGEXEC(SCIPdialogExecPresolve);

/** dialog execution method for the quit command */
extern
DECL_DIALOGEXEC(SCIPdialogExecQuit);

/** dialog execution method for the read command */
extern
DECL_DIALOGEXEC(SCIPdialogExecRead);

/** dialog execution method for the set read command */
extern
DECL_DIALOGEXEC(SCIPdialogExecSetRead);

/** dialog execution method for the set write command */
extern
DECL_DIALOGEXEC(SCIPdialogExecSetWrite);

/** dialog execution method for the set parameter command */
extern
DECL_DIALOGEXEC(SCIPdialogExecSetParam);

/** dialog description method for the set parameter command */
extern
DECL_DIALOGDESC(SCIPdialogDescSetParam);

#ifndef NDEBUG
/** dialog execution method for the debug memory command */
extern
DECL_DIALOGEXEC(SCIPdialogExecDebugMemory);
#endif

/** includes or updates the default dialog menus in SCIP */
extern
RETCODE SCIPincludeDialogDefault(
   SCIP*            scip                /**< SCIP data structure */
   );

/** includes or updates the "set" menu for each available parameter setting */
extern
RETCODE SCIPincludeDialogDefaultSet(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
