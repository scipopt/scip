/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
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
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** frees separation storage */
SCIP_RETCODE SCIPsepastoreExactFree(
   SCIP_SEPASTOREEXACT** sepastoreexact      /**< pointer to store separation storage */
   );


/** adds cut to separation storage and captures it */
SCIP_RETCODE SCIPsepastoreExactAddCut(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_ROWEXACT*        cut                 /**< separated cut */
   );

/** clears the separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPsepastoreExactClearCuts(
   SCIP_SEPASTOREEXACT*  sepastoreexact,     /**< separation storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
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
