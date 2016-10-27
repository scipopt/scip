/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   conflictstore.h
 * @brief  internal methods for storing conflicts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICTSTORE_H__
#define __SCIP_CONFLICTSTORE_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_conflictstore.h"
#include "scip/type_retcode.h"
#include "scip/type_cons.h"
#include "scip/type_event.h"
#include "scip/type_conflict.h"
#include "scip/type_prob.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates separation storage */
extern
SCIP_RETCODE SCIPconflictstoreCreate(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict storage */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees separation storage */
extern
SCIP_RETCODE SCIPconflictstoreFree(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** adds a constraint to the pool of dual rays */
extern
SCIP_RETCODE SCIPconflictstoreAddDualray(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_CONS*            dualray,            /**< dual ray to add */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob           /**< transformed problem */
   );

/** adds a conflict to the conflict storage */
extern
SCIP_RETCODE SCIPconflictstoreAddConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_EVENTFILTER*     eventfilter,        /**< eventfiler */
   SCIP_CONS*            cons,               /**< constraint representing the conflict */
   SCIP_NODE*            node,               /**< node to add conflict (or NULL if global) */
   SCIP_NODE*            validnode,          /**< node at whichaddConf the constraint is valid (or NULL) */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             cutoffinvolved,     /**< is a cutoff bound invaled in this conflict */
   SCIP_Real             primalbound         /**< primal bound the conflict depend on (or -SCIPinfinity) */
   );

/** delete all conflicts depending a cutoff bound larger than the given bound */
extern
SCIP_RETCODE SCIPconflictstoreCleanNewIncumbant(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem*/
   SCIP_Real             cutoffbound         /**< current cutoff bound */
   );

/** return the maximal size of the conflict pool */
extern
int SCIPconflictstoreGetMaxPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   );

/** return the initial size of the conflict pool */
extern
int SCIPconflictstoreGetInitPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   );

/** returns the number of stored conflicts on the conflict pool
 *
 *  note: the number of active conflicts can be less
 */
extern
int SCIPconflictstoreGetNConflictsInStore(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   );

/** returns all active conflicts stored in the conflict store */
extern
SCIP_RETCODE SCIPconflictstoreGetConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_CONS**           conflicts,          /**< array to store conflicts */
   int                   conflictsize,       /**< site of the conflict array */
   int*                  nconflicts          /**< pointer to store the number of conflicts */
   );

#ifdef __cplusplus
}
#endif

#endif
