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

/**@file   heur_rec.h
 * @brief  Primal recombination heuristic for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file implements a recombination heuristic for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_REC_H__
#define __SCIP_HEUR_REC_H__

#include "scip/scip.h"
#include "grph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the rec primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurRec(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** heuristic to exclude vertices or edges from a given solution (and inserting other edges) to improve objective */
extern
SCIP_RETCODE SCIPheurExclusion(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   const SCIP_Real*      costrev,            /**< reversed edge costs */
   const int*            result,             /**< edge solution array (UNKNOWN/CONNECT) */
   int*                  newresult,          /**< new edge solution array (UNKNOWN/CONNECT) */
   int*                  dnodemap,           /**< node array for internal use */
   int*                  nodearrint,         /**< node array for internal use */
   STP_Bool*             stvertex,           /**< node array for internally marking solution vertices */
   SCIP_Bool*            success             /**< solution improved? */
   );

#ifdef __cplusplus
}
#endif

#endif
