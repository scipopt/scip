/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   heur_localsearch.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Local search heuristic
 * @author Dominik Kamp
 * @author Gennesaret Tjusila
 *
 * Local search heuristic based on weighted constraint satisfaction with tabu control.
 *
 * The algorithm operates in two modes:
 * - Feasibility mode (cutoff bound infinite): starts from variable bounds, searches for a feasible solution using
 *   constraint-guided moves.
 * - Optimality mode (cutoff bound finite): starts from the best known solution, searches for an improving solution with
 *   objective breakpoint target set to cutoff bound minus cutoff bound delta.
 *
 * The search uses five move operators:
 * - UnsatTightMove: samples violated constraints and moves toward satisfaction.
 * - SatTightMove: samples satisfied constraints and moves toward constraint boundaries (optimality mode only).
 * - FlipMove: flips binary variables with bounds checking.
 * - RandomTightMove: escapes local minima using soft aspiration tabu and accepting negative-score moves.
 * - LiftMove: improves objective variables while maintaining feasibility (optimality mode only).
 *
 * Constraints are normalized to sum(coeff * var) <= rhs form. Move scores are integer sums of constraint weights based
 * on feasibility transitions. Constraint weights grow for persistent violations and decay for satisfied constraints.
 * Objective and constraint comparisons use infinity guards to handle values beyond the infinity range. Terminates on
 * infinite objective value.
 *
 * Based on:
 *  Peng Lin, Shaowei Cai, Mengchuan Zou, Jinkun Lin.
 *  "Local-MIP: Efficient Local Search for Mixed Integer Programming"
 *  Artificial Intelligence, Volume 348, 2025, 104405.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCALSEARCH_H__
#define __SCIP_HEUR_LOCALSEARCH_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the localsearch primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurLocalsearch(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
