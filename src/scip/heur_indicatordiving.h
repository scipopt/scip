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

/**@file   heur_indicatordiving.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP diving heuristic that fixes indicator variables controlling semicontinuous variables
 * @author Katrin Halbig
 * @author Alexander Hoen
 *
 * A diving heuristic iteratively rounds some fractional variables or variables determined by constraint handlers,
 * and resolves the LP relaxation. Thereby simulating a depth-first-search in the tree.
 *
 * Indicatordiving focuses on indicator variables, which control semicontinuous variables.
 * If the semicontinuous variable is unbounded, the indicator constraint is not part of the LP and,
 * therefore, the indicator variable is not set to an useful value in the LP solution.
 *
 * For these indicator variables the score depends on the LP value and the bounds of the corresponding semicontinuous variable.
 * If parameter usevarbounds=TRUE, also varbound constraints modeling semicontinuous variables are considered.
 * For all other variables the Farkas score (scaled) is returned.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_INDICATORDIVING_H__
#define __SCIP_HEUR_INDICATORDIVING_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the indicatordiving heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurIndicatordiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
