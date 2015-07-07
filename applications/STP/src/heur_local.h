/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_local.h
 * @brief  Improvement heuristic for STP
 * @author Daniel Rehfeldt
 *
 * This file implements three local heuristics, namely vertex insertion, key-path exchange and key-vertex elimination,
 * see "Fast Local Search for Steiner Trees in Graphs" by Uchoa and Werneck.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCAL_H__
#define __SCIP_HEUR_LOCAL_H__


#include "scip/scip.h"
#include "grph.h"
#ifdef __cplusplus
extern "C" {
#endif

/** creates the local primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurLocal(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** perform local heuristics on a given Steiner tree */
extern
SCIP_RETCODE SCIPheurImproveSteinerTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   const SCIP_Real*      cost,               /**< arc cost array */
   const SCIP_Real*      costrev,            /**< reversed arc cost array */
   int*                  best_result         /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   );

#ifdef __cplusplus
}
#endif

#endif
