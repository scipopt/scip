/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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
EXTERN
SCIP_RETCODE SCIPincludeHeurCycGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
