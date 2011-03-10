/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dialog_default.h
 * @brief  default user interface dialog
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DIALOG_DEFAULT_H__
#define __SCIP_DIALOG_DEFAULT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** standard menu dialog execution method, that displays it's help screen if the remaining command line is empty */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenu);

/** standard menu dialog execution method, that doesn't display it's help screen */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecMenuLazy);

/** dialog execution method for the change add constraint */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeAddCons);

/** dialog execution method for the change bounds command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChangeBounds);

/** dialog execution method for the checksol command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecChecksol);

/** dialog execution method for the conflictgraph command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecConflictgraph);

/** dialog execution method for the display branching command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayBranching);

/** dialog execution method for the display conflict command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConflict);

/** dialog execution method for the display conshdlrs command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayConshdlrs);

/** dialog execution method for the display displaycols command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayDisplaycols);

/** dialog execution method for the display heuristics command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayHeuristics);

/** dialog execution method for the display memory command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayMemory);

/** dialog execution method for the display nodeselectors command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNodeselectors);

/** dialog execution method for the display nlpi command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayNlpi);

/** dialog execution method for the display parameters command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayParameters);

/** dialog execution method for the display presolvers command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPresolvers);

/** dialog execution method for the display pricer command */
extern  
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPricers);

/** dialog execution method for the display problem command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayProblem);

/** dialog execution method for the display propagators command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayPropagators);

/** dialog execution method for the display readers command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayReaders);

/** dialog execution method for the display separators command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySeparators);

/** dialog execution method for the display solution command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplaySolution);

/** dialog execution method for the display statistics command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayStatistics);

/** dialog execution method for the display transproblem command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTransproblem);

/** dialog execution method for the display value command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayValue);

/** dialog execution method for the display varbranchstatistics command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayVarbranchstatistics);

/** dialog execution method for the display LP solution quality command */
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayLPSolutionQuality);

/** dialog execution method for the display transsolution command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecDisplayTranssolution);

/** dialog execution method for the help command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecHelp);

/** dialog execution method for the free command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecFree);

/** dialog execution method for the newstart command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecNewstart);

/** dialog execution method for the optimize command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecOptimize);

/** dialog execution method for the presolve command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecPresolve);

/** dialog execution method for the quit command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecQuit);

/** dialog execution method for the read command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecRead);

/** dialog execution method for the set default command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDefault);

/** dialog execution method for the set load command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLoad);

/** dialog execution method for the set save command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSave);

/** dialog execution method for the set diffsave command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetDiffsave);

/** dialog execution method for the set parameter command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetParam);

/** dialog description method for the set parameter command */
extern
SCIP_DECL_DIALOGDESC(SCIPdialogDescSetParam);

/** dialog execution method for the set branching direction command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingDirection);

/** dialog execution method for the set branching priority command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetBranchingPriority);

/** dialog execution method for the set heuristics aggressive command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsEmphasisAggressive);

/** dialog execution method for the set heuristics fast command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsEmphasisFast);

/** dialog execution method for the set heuristics off command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetHeuristicsEmphasisOff);

/** dialog execution method for the set presolving aggressive command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingEmphasisAggressive);

/** dialog execution method for the set presolving fast command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingEmphasisFast);

/** dialog execution method for the set presolving off command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetPresolvingEmphasisOff);

/** dialog execution method for the set separating aggressive command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingEmphasisAggressive);

/** dialog execution method for the set separating fast command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingEmphasisFast);

/** dialog execution method for the set separating off command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetSeparatingEmphasisOff);

/** dialog execution method for the set emphasis counter command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCounter);

/** dialog execution method for the set emphasis cpsolver command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisCpsolver);

/** dialog execution method for the set emphasis easy CIP command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisEasycip);

/** dialog execution method for the set emphasis feasibility command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisFeasibility);

/** dialog execution method for the set emphasis hard LP command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisHardlp);

/** dialog execution method for the set emphasis optimality command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetEmphasisOptimality);

/** dialog execution method for the set limits objective command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecSetLimitsObjective);

/** creates a root dialog */
extern
SCIP_RETCODE SCIPcreateRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         root                /**< pointer to store the root dialog */
   );

/** includes or updates the default dialog menus in SCIP */
extern
SCIP_RETCODE SCIPincludeDialogDefault(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes or updates the "set" menu for each available parameter setting */
extern
SCIP_RETCODE SCIPincludeDialogDefaultSet(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
