/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   conflict_resolution.h
 * @ingroup OTHER_CFILES
 * @brief   Methods for generalized resolution conflict analysis.
 * @author  Gioni Mexi
 *
 * This file implements a conflict analysis method based on generalized resolution,
 * as detailed in the paper:
 *
 * Gioni Mexi, et al. "Cut-based Conflict Analysis in Mixed Integer Programming."
 * arXiv preprint arXiv:2410.15110 (2024).
 *
 * The generalized resolution conflict analysis procedure starts with an initial
 * conflict row and it iteratively aggregates this row with a reason row—the row
 * that propagated the bound change causing the conflict. The aggregation
 * cancels the coefficient of the resolving variable. This process continues
 * until a first unique implication point (FUIP) is reached. If the aggregation
 * does not yield a valid (infeasible) row, the algorithm attempts to reduce the
 * reason row (e.g., using MIR reduction) and retries the aggregation. Once a
 * valid conflict row with negative slack is generated, a conflict constraint is
 * constructed and added to the problem.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONFLICT_RESOLUTION_H__
#define __SCIP_CONFLICT_RESOLUTION_H__

#include "scip/type_conflict.h"
#include "scip/type_reopt.h"
#include "scip/type_implics.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"
#include "scip/type_branch.h"
#include "scip/type_var.h"
#include "scip/type_prob.h"
#include "scip/type_event.h"

#ifdef __cplusplus
extern "C" {
#endif

/** return TRUE if generalized resolution conflict analysis is applicable */
SCIP_Bool SCIPconflictResolutionApplicable(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** gets number of conflict constraints found during propagation with the generalized resolution conflict analysis */
SCIP_Longint SCIPconflictGetNResConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to generalized resolution conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNResSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to generalized resolution conflict analysis terminating because of large coefficients */
SCIP_Longint SCIPconflictGetNResLargeCoefs(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to generalized resolution conflict analysis terminating because of long conflicts */
SCIP_Longint SCIPconflictGetNResLongConflicts(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** gets number of calls to generalized resolution conflict analysis */
SCIP_Longint SCIPconflictGetNResCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   );

/** create constraints out of the conflict row and add them to the problem */
SCIP_RETCODE SCIPconflictAddConflictCon(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LP*              lp,                 /**< current LP data */
   SCIP_BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< global event filter */
   SCIP_CLIQUETABLE*     cliquetable,        /**< clique table data structure */
   SCIP_CONFLICTROW*     conflictrow,        /**< conflict row to add to the tree */
   SCIP_Bool*            success             /**< true if the conflict is added to the problem */
   );

/** creates and clears the conflict rows */
SCIP_RETCODE SCIPconflictInitRows(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   );

/** frees a conflict row */
void SCIPconflictRowFree(
   SCIP_CONFLICTROW**    row,                /**< conflict row */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

#ifdef __cplusplus
}
#endif


#endif
