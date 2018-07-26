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

/**@file   heur_cyckerlin.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Improvement heuristic that trades bin-variables between clusters
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_CYCKERLIN_H__
#define __SCIP_HEUR_CYCKERLIN_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the oneopt primal heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeHeurCycKerlin(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** External method that adds a solution to the list of candidate-solutions that should be improved */
EXTERN
SCIP_RETCODE addCandSolCyckerlin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< The given solution */
   );

#ifdef __cplusplus
}
#endif

#endif
