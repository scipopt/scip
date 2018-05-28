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

/** updates solutions for pruned graph */
extern
SCIP_RETCODE SCIPStpHeurPruneUpdateSols(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                g,                  /**< graph data structure */
   GRAPH*                prunegraph,         /**< pruned graph data structure */
   PATH*                 path,               /**< shortest path struct */
   int*                  nodearrint,         /**< array */
   int*                  edgearrint,         /**< array */
   int*                  solnode,            /**< array for best solution nodes wrt prunegraph */
   int*                  soledge,            /**< array for best solution edges wrt prunegraph */
   int*                  globalsoledge,      /**< array storing best solution wrt g */
   STP_Bool*             nodearrchar,        /**< array */
   SCIP_Real*            globalobj,          /**< pointer to objective value of best solution wrt g */
   SCIP_Bool             incumbentgiven,     /**< incumbent solution for pruned graph given? */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
   );

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
   const SCIP_Bool       withinitialsol,     /**< solution given? */
   const SCIP_Bool       reducegraph         /**< try to reduce graph initially? */
	     );

#ifdef __cplusplus
}
#endif

#endif
