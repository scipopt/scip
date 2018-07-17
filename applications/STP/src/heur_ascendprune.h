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

/**@file   heur_ascendprune.h
 * @ingroup PRIMALHEURISTICS
 * @brief  reduction and dual-cost based primal heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a reducion and dual-cost based heuristic for Steiner problems. It is based on an approach
 * described in T. Polzin's "Algorithms for the Steiner problem in networks".
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_ASCENTPRUNE_H__
#define __SCIP_HEUR_ASCENTPRUNE_H__


#include "scip/scip.h"
#include "grph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the prune primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPStpIncludeHeurAscendPrune(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** ascent and prune */
extern
SCIP_RETCODE SCIPStpHeurAscendPruneRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure or NULL */
   const GRAPH*          g,                  /**< the graph */
   const SCIP_Real*      redcosts,           /**< the reduced costs */
   int*                  edgearrint,         /**< int edges array to store solution */
   int*                  nodearrint,         /**< int vertices array for internal computations */
   int                   root,               /**< the root (used for dual ascent) */
   STP_Bool*             nodearrchar,        /**< char vertices array for internal computations */
   SCIP_Bool*            solfound,           /**< has a solution been found? */
   SCIP_Bool             addsol              /**< should the solution be added to SCIP by this method? */
   );


#ifdef __cplusplus
}
#endif

#endif
