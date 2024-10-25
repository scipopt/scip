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

/**@file   iis_greedy.h
 * @brief  greedy deletion and addition filter heuristic to compute IISs
 * @author Marc Pfetsch
 * @author Mark Turner
 *
 * An irreducibly infeasible subsystem (IIS) is a subset of the constraints that is infeasible and (set-wise) minimial
 * in this respect.
 * The deletion filter heuristic greedily removes constraints and checks whether the remaining problem
 * is still infeasible.
 * The addition filter heuristic greedily adds constraints from an empty problem until the problem becomes
 * infeasible. The deletion filter can then be applied as a post-processing step.
 * The method is based on
 *
 * O. Guieu and J. Chinneck, Analyzing infeasible mixed-integer and integer linear programs,@p
 * INFORMS J. Comput. 11, no. 1 (1999), pp. 63–77.
 *
 * If the appropriate parameters are set then we can guarantee that the result is minimal, i.e.,
 * an irreducible infeasible subsystem (IIS). Otherwise we may only obtain an infeasible subsystem (IS).
 * For no settings can guarantee the smallest possible infeasible subsystem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __IIS_GREEDY_H__
#define __IIS_GREEDY_H__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** creates the greedy IIS rule and includes it in SCIP
 *
 * @ingroup IISIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeIISGreedy(
   SCIP*                 scip                /**< SCIP data structure */
);

/**@addtogroup IIS
 *
 * @{
 */

/** perform a greedy addition or deletion algorithm to obtain an infeasible subsystem (IS).
 *
 *  This is the generation method for the greedy IIS rule.
 *  Depending on the parameter choices, constraints are either greedily added from an empty problem,
 *  or deleted from a complete problem. In the case of constraints being added, this is done until the problem
 *  becomes infeasible, after which one can then begin deleting constraints. In the case of deleting constraints,
 *  this is done until no more constraints (or batches of constraints) can be deleted withough making
 *  the problem feasible.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgenerateIISGreedy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Bool             minify,             /**< whether the computed IS should undergo a final deletion round to ensure an IIS */
   SCIP_Bool             additive,           /**< whether an additive approach instead of deletion based approach should be used */
   SCIP_Bool             conservative,       /**< should a node or time limit solve be counted as feasible when deleting constraints */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */
   SCIP_Longint          maxnnodesperiter,   /**< maximum number of nodes per individual solve call */
   int                   batchsize           /**< the number of constraints to delete or add per iteration */
);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
