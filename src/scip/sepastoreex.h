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

/**@file   sepastoreex.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing separated exact cuts
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPASTOREEX_H__
#define __SCIP_SEPASTOREEX_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_implics.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_lpex.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_sepastore.h"
#include "scip/type_branch.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates separation storage */
extern
SCIP_RETCODE SCIPsepastoreexCreate(
   SCIP_SEPASTOREEX**    sepastoreex,          /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees separation storage */
extern
SCIP_RETCODE SCIPsepastoreexFree(
   SCIP_SEPASTOREEX**    sepastoreex,          /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** informs separation storage that the setup of the initial LP starts now */
extern
void SCIPsepastoreexStartInitialLP(
   SCIP_SEPASTOREEX*     sepastoreex           /**< separation storage */
   );

/** informs separation storage that the setup of the initial LP is now finished */
extern
void SCIPsepastoreexEndInitialLP(
   SCIP_SEPASTOREEX*     sepastoreex           /**< separation storage */
   );

/** adds cut to separation storage and captures it */
extern
SCIP_RETCODE SCIPsepastoreexAddCut(
   SCIP_SEPASTOREEX*     sepastoreex,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEX*            lp,                 /**< LP data */
   SCIP_ROWEX*           cut,                /**< separated cut */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut is infeasible */
   );

/** adds cuts to the LP and clears separation storage */
extern
SCIP_RETCODE SCIPsepastoreexApplyCuts(
   SCIP_SEPASTOREEX*     sepastoreex,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LPEX*            lp,                 /**< LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_Bool             root,               /**< are we at the root node? */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   );

/** clears the separation storage without adding the cuts to the LP */
extern
SCIP_RETCODE SCIPsepastoreexClearCuts(
   SCIP_SEPASTOREEX*     sepastoreex,          /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEX*            lp                  /**< LP data */
   );

/** indicates whether a cut is applicable
 *
 *  A cut is applicable if it is modifiable, not a bound change, or a bound change that changes bounds by at least epsilon.
 */
extern
SCIP_Bool SCIPsepastoreexIsCutApplicable(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             cut                 /**< cut to check */
   );

/** get cuts in the separation storage */
extern
SCIP_ROWEX** SCIPsepastoreexGetCuts(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   );

/** get number of cuts in the separation storage */
extern
int SCIPsepastoreexGetNCuts(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   );

/** get total number of cuts found so far */
extern
int SCIPsepastoreexGetNCutsFound(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   );

/** get number of cuts found so far in current separation round */
extern
int SCIPsepastoreexGetNCutsFoundRound(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   );

/** get total number of cuts applied to the LPs */
extern
int SCIPsepastoreexGetNCutsApplied(
   SCIP_SEPASTOREEX*       sepastoreex           /**< separation storage */
   );

#ifdef __cplusplus
}
#endif

#endif
