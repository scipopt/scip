/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_lurkprune.h
 * @ingroup PRIMALHEURISTICS
 * @brief  reduction based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reduction based heuristic for Steiner problems that makes use of lurking bounds.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LURKPRUNE_H__
#define __SCIP_HEUR_LURKPRUNE_H__


#include "scip/scip.h"
#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the lurk prune primal heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPStpIncludeHeurLurkPrune(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** execute lurk-and-prune heuristic on given graph */
EXTERN
SCIP_RETCODE SCIPStpHeurLurkPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables or NULL */
   const SCIP_Real*      lurkingbounds,      /**< lurking edge bounds */
   GRAPH*                g,                  /**< graph data structure */
   SCIP_Bool             initialreduce,      /**< try to reduce graph initially? */
   SCIP_Bool             ascendprune,        /**< use ascend-prune? */
   int*                  soledge,            /**< array to 1. provide and 2. return primal solution */
   SCIP_Bool*            solimproved         /**< could a better solution be found? */
   );

#ifdef __cplusplus
}
#endif

#endif
