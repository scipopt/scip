/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_dps.h
 * @ingroup PRIMALHEURISTICS
 * @brief  dynamic partition search
 * @author Katrin Halbig
 *
 * The dynamic partition search (DPS) is a construction heuristic which additionally needs a
 * user decomposition with linking constraints only.
 *
 * This heuristic splits the problem into several sub-SCIPs according to the given decomposition. Thereby the linking constraints
 * with their right-hand and left-hand sides are also split. DPS searches for a partition of the sides on the blocks
 * so that a feasible solution is obtained.
 * For each block the parts of the original linking constraints are extended by slack variables. Moreover, the objective function
 * is replaced by the sum of these additional variables weighted by penalty parameters lambda. If all blocks have an optimal solution
 * of zero, the algorithm terminates with a feasible solution for the main problem. Otherwise, the partition and the penalty parameters
 * are updated, and the sub-SCIPs are solved again.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_DPS_H__
#define __SCIP_HEUR_DPS_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dps primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurDps(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
