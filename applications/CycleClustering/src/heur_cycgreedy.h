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

/**@file   heur_cycgreedy.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Greedy primal heuristic. States are assigned to clusters iteratively. At each iteration all possible
 * assignments are computed and the one with the best change in objective value is selected.
 * @author Leon Eifler
 *
 */

#ifndef __SCIP_HEUR_CYCGREEDY_H__
#define __SCIP_HEUR_CYCGREEDY_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the CycGreedy primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurCycGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
