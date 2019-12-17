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
#include "graph.h"
#include "solpool.h"

#ifdef __cplusplus
extern "C" {
#endif


/** run REC heuristic */
SCIP_RETCODE SCIPStpHeurRecRun(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL*           pool,               /**< solution pool or NULL */
   SCIP_HEUR*            heur,               /**< heuristic or NULL */
   SCIP_HEURDATA*        heurdata,           /**< heuristic data or NULL */
   const GRAPH*          graph,              /**< graph data */
   SCIP_VAR**            vars,               /**< variables or NULL */
   int*                  newsolindex,        /**< index of new solution */
   int                   runs,               /**< number of runs */
   int                   nsols,              /**< number of solutions */
   SCIP_Bool             restrictheur,       /**< use restricted version of heur? */
   SCIP_Bool*            solfound            /**< new solution found? */
);


/** creates the rec primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPStpIncludeHeurRec(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** heuristic to exclude vertices or edges from a given solution (and inserting other edges) to improve objective */
SCIP_RETCODE SCIPStpHeurRecExclude(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph structure */
   const int*            result,             /**< edge solution array (UNKNOWN/CONNECT) */
   int*                  newresult,          /**< new edge solution array (UNKNOWN/CONNECT) */
   int*                  dnodemap,           /**< node array for internal use */
   STP_Bool*             stvertex,           /**< node array for internally marking solution vertices */
   SCIP_Bool*            success             /**< solution improved? */
   );

#ifdef __cplusplus
}
#endif

#endif
