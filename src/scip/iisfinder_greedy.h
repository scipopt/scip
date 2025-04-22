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

/** perform a greedy addition or deletion algorithm to obtain an infeasible subsystem (IS).
 *
 *  This is the generation method for the greedy IIS finder rule.
 *  Depending on the parameter choices, constraints are either greedily added from an empty problem,
 *  or deleted from a valid problem state. In the case of constraints being added, this is done until the problem
 *  becomes infeasible, after which one can then begin deleting constraints. In the case of deleting constraints,
 *  this is done until no more constraints (or batches of constraints) can be deleted without making
 *  the problem feasible.
 *  The algorithm also extends to variable bounds.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexecIISfinderGreedy(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_Real             timelim,            /**< The global time limit on the IIS call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS call */
   SCIP_Bool             removebounds,       /**< Whether the algorithm should remove bounds as well as constraints */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */

   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             additive,           /**< whether an additive approach instead of deletion based approach should be used */
   SCIP_Bool             conservative,       /**< should a hit limit (e.g. node / time) solve be counted as feasible when deleting constraints */
   SCIP_Bool             delafteradd,        /**< should the deletion routine be performed after the addition routine (in the case of additive) */
   SCIP_Bool             dynamicreordering,  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */

   int                   initbatchsize,      /**< the initial batchsize for the first iteration */
   SCIP_Real             initrelbatchsize,   /**< the initial batchsize relative to the original problem for the first iteration */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   SCIP_Real             maxrelbatchsize,    /**< the maximum batchsize relative to the original problem per iteration */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied each iteration */
   int                   batchingoffset,     /**< the offset with which the batchsize is summed each iteration */

   SCIP_RESULT*          result              /**< pointer to store the result of the IIS run */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
