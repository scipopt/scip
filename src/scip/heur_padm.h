/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_padm.h
 * @ingroup PRIMALHEURISTICS
 * @brief  PADM primal heuristic based on ideas published in the paper
 *         "A Decomposition Heuristic for Mixed-Integer Supply Chain Problems"
 *         by Martin Schmidt, Lars Schewe, and Dieter Weninger
 * @author Dieter Weninger
 * @author Katrin Halbig
 * 
 * The penalty alternating direction method (PADM) heuristic is a construction heuristic which additionally needs a
 * user decomposition with linking variables only.
 * 
 * PADM splits the problem into several sub-SCIPs according to the decomposition, whereby the linking variables get
 * copied and the difference is penalized. Then the sub-SCIPs are solved on an alternating basis until they arrive at
 * the same values of the linking variables (ADM-loop). If they don't reconcile after a couple of iterations,
 * the penalty parameters are increased (penalty-loop) and the sub-SCIPs are solved again on an alternating basis.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_PADM_H__
#define __SCIP_HEUR_PADM_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the PADM primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurPADM(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
