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

/**@file   scip_dialog.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for dialog handler plugins
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_DIALOG_H__
#define __SCIP_SCIP_DIALOG_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "scip/type_timing.h"
#include "scip/type_paramset.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_nlp.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_scip.h"

#include "scip/type_bandit.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_cons.h"
#include "scip/type_dialog.h"
#include "scip/type_disp.h"
#include "scip/type_heur.h"
#include "scip/type_compr.h"
#include "scip/type_history.h"
#include "scip/type_nodesel.h"
#include "scip/type_presol.h"
#include "scip/type_pricer.h"
#include "scip/type_reader.h"
#include "scip/type_relax.h"
#include "scip/type_sepa.h"
#include "scip/type_table.h"
#include "scip/type_prop.h"
#include "nlpi/type_nlpi.h"
#include "scip/type_concsolver.h"
#include "scip/type_syncstore.h"
#include "scip/type_benders.h"
#include "scip/type_benderscut.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicDialogMethods
 *
 * @{
 */

/** creates and includes dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPincludeDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog,             /**< pointer to store the dialog */
   SCIP_DECL_DIALOGCOPY  ((*dialogcopy)),    /**< copy method of dialog or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   SCIP_DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   SCIP_DECL_DIALOGFREE  ((*dialogfree)),    /**< destructor of dialog to free user data, or NULL */
   const char*           name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*           desc,               /**< description of dialog used if description output method is NULL */
   SCIP_Bool             issubmenu,          /**< is the dialog a submenu? */
   SCIP_DIALOGDATA*      dialogdata          /**< user defined dialog data */
   );

/** returns if the dialog already exists
 *
 *  @return TRUE is returned if the dialog exists, otherwise FALSE.
 */
EXTERN
SCIP_Bool SCIPexistsDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** captures a dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPcaptureDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   );

/** releases a dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPreleaseDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog              /**< pointer to the dialog */
   );

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPsetRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog to be the root */
   );

/** returns the root dialog of SCIP's interactive user shell
 *
 *  @return the root dialog of SCIP's interactive user shell is returned.
 */
EXTERN
SCIP_DIALOG* SCIPgetRootDialog(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds a sub dialog to the given dialog as menu entry and captures it
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPaddDialogEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   SCIP_DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   );

/** adds a single line of input which is treated as if the user entered the command line
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPaddDialogInputLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   );

/** adds a single line of input to the command history which can be accessed with the cursor keys
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 */
EXTERN
SCIP_RETCODE SCIPaddDialogHistoryLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   );

/** starts interactive mode of SCIP by executing the root dialog
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. otherwise a suitable error code is passed. see \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_FREE
 *
 *  @post After calling this method \SCIP reaches one of the following stages depending on if and when the
 *        interactive shell was closed:
 *        - \ref SCIP_STAGE_PROBLEM if the interactive shell was closed after the problem was created
 *        - \ref SCIP_STAGE_TRANSFORMED if the interactive shell was closed after the problem was transformed
 *        - \ref SCIP_STAGE_PRESOLVING if the interactive shell was closed  during presolving
 *        - \ref SCIP_STAGE_PRESOLVED if the interactive shell was closed after presolve
 *        - \ref SCIP_STAGE_SOLVING if the interactive shell was closed during the tree search
 *        - \ref SCIP_STAGE_SOLVED if the interactive shell was closed after the problem was solved
 *        - \ref SCIP_STAGE_FREE if the interactive shell was closed after the problem was freed
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
EXTERN
SCIP_RETCODE SCIPstartInteraction(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
