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

/**@file   scip_bandit.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for bandit algorithms
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

#ifndef __SCIP_SCIP_BANDIT_H__
#define __SCIP_SCIP_BANDIT_H__


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

/**@addtogroup PublicBanditMethods
 *
 * @{
 */

/** includes a bandit algorithm virtual function table  */
EXTERN
SCIP_RETCODE SCIPincludeBanditvtable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDITVTABLE**   banditvtable,       /**< bandit algorithm virtual function table */
   const char*           name,               /**< a name for the algorithm represented by this vtable */
   SCIP_DECL_BANDITFREE  ((*banditfree)),    /**< callback to free bandit specific data structures */
   SCIP_DECL_BANDITSELECT((*banditselect)),  /**< selection callback for bandit selector */
   SCIP_DECL_BANDITUPDATE((*banditupdate)),  /**< update callback for bandit algorithms */
   SCIP_DECL_BANDITRESET ((*banditreset))    /**< update callback for bandit algorithms */
   );

/** returns the bandit virtual function table of the given name, or NULL if not existing */
EXTERN
SCIP_BANDITVTABLE* SCIPfindBanditvtable(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of bandit algorithm virtual function table */
   );

/** calls destructor and frees memory of bandit algorithm */
EXTERN
SCIP_RETCODE SCIPfreeBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT**         bandit              /**< pointer to bandit algorithm data structure */
   );

/** reset the bandit algorithm */
EXTERN
SCIP_RETCODE SCIPresetBandit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BANDIT*          bandit,             /**< pointer to bandit algorithm data structure */
   SCIP_Real*            priorities,         /**< priorities for every action, or NULL if not needed */
   unsigned int          seed                /**< initial random seed for bandit selection */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
