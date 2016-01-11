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
   SCIP_CONFLICTSTORE**  conflictstore       /**< pointer to store conflict storage */
   );

/** frees separation storage */
extern
SCIP_RETCODE SCIPconflictstoreFree(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict storage */
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** adds a conflict to the conflict storage */
extern
SCIP_RETCODE SCIPconflictstoreAddConflict(
   SCIP_CONFLICTSTORE*   conflictstore,
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,
   SCIP_CONS*            cons,
   SCIP_NODE*            node,
   SCIP_NODE*            validnode,
   SCIP_Bool             global,
   SCIP_CONFTYPE         conftype,
   SCIP_Bool             cutoffinvolved,     /**< is a cutoff bound invaled in this conflict */
   SCIP_Real             primalbound
   );

#ifdef __cplusplus
}
#endif

#endif
