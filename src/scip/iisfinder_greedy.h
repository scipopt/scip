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

/**@file   iisfinder_greedy.h
 * @ingroup IISFINDERS
 * @brief  greedy deletion and addition filter heuristic to compute IISs
 * @author Paul Meinhold
 * @author Marc Pfetsch
 * @author Mark Turner
 *
 * An irreducible infeasible subsystem (IIS) is a subset of the constraints and bounds from a problem
 * that is infeasible.
 * The greedy filter heuristic has two different approaches for producing an (I)IS:
 * - Remove constraints such that the remaining problem is still infeasible.
 * - Add constraints from an empty problem until the problem becomes infeasible.
 * Both these approaches can be augmented to include bounds. That is, existing bounds can be deleted
 * while the IS remains infeasible. A common approach is to also apply the deletion based
 * filter after applying the additive based filter.
 * This greedy algorithm is based on
 *
 * O. Guieu and J. Chinneck, Analyzing infeasible mixed-integer and integer linear programs,@p
 * INFORMS J. Comput. 11, no. 1 (1999), pp. 63–77.
 *
 * If the appropriate parameters are set then we can guarantee that the result is minimal, i.e.,
 * an irreducible infeasible subsystem (IIS). Otherwise we may only obtain an infeasible subsystem (IS).
 * For no settings can we guarantee the smallest possible infeasible subsystem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __IISFINDER_GREEDY_H__
#define __IISFINDER_GREEDY_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** creates the greedy IIS finder rule and includes it in SCIP
 *
 * @ingroup IISfinderIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeIISfinderGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup IISFINDERS
 *
 * @{
 */

/** perform the greedy deletion algorithm with singleton batches to obtain an irreducible infeasible subsystem (IIS) */
SCIP_EXPORT
SCIP_RETCODE SCIPiisGreedyMinimize(
   SCIP_IIS*             iis                 /**< IIS data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
