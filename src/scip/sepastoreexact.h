/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepastoreexact.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for storing separated exact cuts
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPASTOREEXACT_H__
#define __SCIP_SEPASTOREEXACT_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_implics.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_lpexact.h"
#include "scip/type_prob.h"
#include "scip/type_tree.h"
#include "scip/type_reopt.h"
#include "scip/type_sepastore.h"
#include "scip/type_branch.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates separation storage */
SCIP_RETCODE SCIPsepastoreExactCreate(
   SCIP_SEPASTOREEXACT** sepastoreexact,     /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees separation storage */
SCIP_RETCODE SCIPsepastoreExactFree(
   SCIP_SEPASTOREEXACT** sepastoreexact,     /**< pointer to store separation storage */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );


/** adds cut to separation storage and captures it */
SCIP_RETCODE SCIPsepastoreexAddCut(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEXACT*         lp,                 /**< LP data */
   SCIP_ROWEXACT*        cut,                /**< separated cut */
   SCIP_Bool*            infeasible          /**< pointer to store whether the cut is infeasible */
   );

/** adds cuts to the LP and clears separation storage */
SCIP_RETCODE SCIPsepastoreExactSyncLPs(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_LPEXACT*         lp,                 /**< LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter         /**< global event filter */
   );

/** clears the separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreExactClearCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global events */
   SCIP_LPEXACT*         lp                  /**< LP data */
   );

/** get cuts in the separation storage */
SCIP_ROWEXACT** SCIPsepastoreExactGetCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   );

/** get number of cuts in the separation storage */
int SCIPsepastoreExactGetNCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   );

/** get total number of cuts found so far */
int SCIPsepastoreExactGetNCutsFound(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   );

/** get number of cuts found so far in current separation round */
int SCIPsepastoreExactGetNCutsFoundRound(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   );

/** get total number of cuts applied to the LPs */
int SCIPsepastoreExactGetNCutsApplied(
   SCIP_SEPASTOREEXACT*  sepastoreexact      /**< separation storage */
   );

#ifdef __cplusplus
}
#endif

#endif
