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

/**@file   heur_tm.h
 * @brief  shortest paths based primal heuristics for Steiner problems
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file implements several shortest paths based primal heuristics for Steiner problems, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_TM_H__
#define __SCIP_HEUR_TM_H__

#include "scip/scip.h"
#include "grph.h"

#define DEFAULT_HOPFACTOR 0.33

#ifdef __cplusplus
extern "C" {
#endif

/** compute starting points among marked (w.r.t. g->mark) vertices for constructive heuristics */
EXTERN
void SCIPStpHeurTMCompStarts(
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< starting points array */
   int*                  runs                /**< pointer to number of runs */
   );

/** creates the TM primal heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPStpIncludeHeurTM(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** execute shortest paths heuristic to obtain a Steiner tree */
EXTERN
SCIP_RETCODE SCIPStpHeurTMRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata,           /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   int*                  starts,             /**< array containing start vertices (NULL to not provide any) */
   int*                  bestnewstart,       /**< pointer to the start vertex resulting in the best soluton */
   int*                  best_result,        /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int                   runs,               /**< number of runs */
   int                   bestincstart,       /**< best incumbent start vertex */
   SCIP_Real*            cost,               /**< arc costs */
   SCIP_Real*            costrev,            /**< reversed arc costs */
   SCIP_Real*            hopfactor,          /**< edge cost multiplicator for HC problems */
   SCIP_Real*            nodepriority,       /**< vertex priorities for vertices to be starting points (NULL for no priorities) */
   SCIP_Real             maxcost,            /**< maximal edge cost (only for HC) */
   SCIP_Bool*            success,            /**< pointer to store whether a solution could be found */
   SCIP_Bool             pcmwfull            /**< use full computation of tree (i.e. connect all terminals and prune), only for prize-collecting variants */
   );

/** run shortest path heuristic, but bias edge costs towards best current LP solution */
EXTERN
SCIP_RETCODE SCIPStpHeurTMRunLP(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< graph data structure */
   SCIP_HEUR*            heur,               /**< heuristic or NULL */
   int*                  result,             /**< array indicating whether an arc is part of the solution (CONNECTED/UNKNOWN) */
   int                   runs,               /**< number of runs */
   SCIP_Real*            cost,               /**< arc costs (uninitialized) */
   SCIP_Real*            costrev,            /**< reversed arc costs (uninitialized) */
   SCIP_Bool*            success             /**< pointer to store whether a solution could be found */
   );

/** prune a Steiner tree in such a way that all leaves are terminals */
EXTERN
SCIP_RETCODE SCIPStpHeurTMPrune(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int                   layer,              /**< layer, @note: should be set to 0 */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   );

/** prune the (rooted) prize collecting Steiner tree in such a way that all leaves are terminals */
EXTERN
SCIP_RETCODE SCIPStpHeurTMPrunePc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   const SCIP_Real*      cost,               /**< edge costs */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   );

/** build (rooted) prize collecting Steiner tree in such a way that all leaves are positive-weight vertices */
EXTERN
SCIP_RETCODE SCIPStpHeurTMBuildTreePcMw(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH*                 mst,                /**< path data structure array */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            objresult,          /**< pointer to store objective value of result */
   int*                  connected           /**< CONNECT/UNKNOWN */
   );

/** build Steiner tree in such a way that all leaves are terminals */
EXTERN
SCIP_RETCODE SCIPStpHeurTMBuildTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   PATH*                 mst,                /**< path data structure array */
   const SCIP_Real*      cost,               /**< edge costs */
   SCIP_Real*            objresult,          /**< pointer to store objective value of result */
   int*                  connected           /**< CONNECT/UNKNOWN */
   );

/** prune a degree constrained Steiner tree in such a way that all leaves are terminals */
EXTERN
SCIP_RETCODE SCIPStpHeurTMBuildTreeDc(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph structure */
   int*                  result,             /**< ST edges */
   STP_Bool*             connected           /**< ST nodes */
   );

#ifdef __cplusplus
}
#endif

#endif
