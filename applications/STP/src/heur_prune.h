/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_prune.h
 * @ingroup PRIMALHEURISTICS
 * @brief  reduction-based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reducion based heuristic for Steiner problems. It is based on an approach
 * described in T. Polzin's "Algorithms for the Steiner problem in networks".
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_PRUNE_H__
#define __SCIP_HEUR_PRUNE_H__


#include "scip/scip.h"
#include "grph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the prune primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPStpIncludeHeurPrune(
   SCIP*                 scip                /**< SCIP data structure */
   );

extern
/** execute prune heuristic on given graph */
SCIP_RETCODE SCIPStpHeurPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   GRAPH*                g,                  /**< the graph */
   int*                  soledge,            /**< array to store primal solution (if no solution is provided,
                                                solgiven must be set to FALSE) */
   SCIP_Bool*            success,            /**< feasible solution found? */
   const SCIP_Bool       solgiven,           /**< solution given? */
   const SCIP_Bool       reducegraph         /**< try to reduce graph initially? */
	     );

#ifdef __cplusplus
}
#endif

#endif
