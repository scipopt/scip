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

/**@file   branch_gomory.h
 * @ingroup BRANCHINGRULES
 * @brief  Gomory cut branching rule
 * @author Mark Turner
 *
 * The approach is based on the following papers.
 *
 * M. Turner, T. Berthold, M. Besancon, T. Koch@n
 * Branching via Cutting Plane Selection: Improving Hybrid Branching,@n
 * arXiv preprint arXiv:2306.06050
 *
 * The Gomory cut branching rule selects a candidate integer variable $j$ with a fractional solution value.
 * Each candidate variable must be a basic variable in the LP Tableau (if not then it would have to be at its bound
 * that is integer-valued)
 * This branching rule calculates the GMI cut for the aggregated row of the LP tableau associated with the
 * candidate variable.
 * The generated cut is then scored using a weighted sum rule.
 * The branching candidate whose cut is highest scoring is then selected.
 * For more details on the method, see:
 *
 * @par
 * Mark Turner, Timo Berthold, Mathieu Besançon, Thorsten Koch@n
 * Branching via Cutting Plane Selection: Improving Hybrid Branching@n
 * 2023@n
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_GOMORY_H__
#define __SCIP_BRANCH_GOMORY_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the Gomory cut branching rule and includes it in SCIP
 *
 *  @ingroup BranchingRuleIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleGomory(
   SCIP*                 scip                /**< SCIP data structure */
);

#ifdef __cplusplus
}
#endif

#endif
