/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_subnlp.h
 * @ingroup PRIMALHEURISTICS
 * @brief  NLP local search primal heuristic using sub-SCIPs
 * @author Stefan Vigerske
 *
 * This heuristic applies a NLP local search to a nonlinear CIP after fixing all discrete variables.
 * That is, the CIP is copied, all discrete variables are fixed, presolving is applied,
 * and if the resulting CIP has a nonlinear relaxation, then it is tried to solve this relaxation
 * by an NLP solver.
 * The heuristic only runs if continuous nonlinearities are present (@ref SCIPhasNLPContinuousNonlinearity()).
 *
 * Fixing values for discrete values are either taken from a solution of the LP relaxation which
 * satisfies all integrality constraints, or are provided by SCIPupdateStartpointHeurSubNlp().
 *
 * This heuristic is orthogonal to the undercover heuristic (@ref heur_undercover.h), which fixes
 * variables in a nonlinear CIP in a way that a (possibly mixed-integer) linear subproblem is obtained.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef HEUR_SUBNLP_H_
#define HEUR_SUBNLP_H_

#include "scip/def.h"
#include "scip/type_heur.h"
#include "scip/type_result.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the NLP local search primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurSubNlp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup PRIMALHEURISTICS
  *
  * @{
  */

/** main procedure of the subNLP heuristic */
SCIP_EXPORT
SCIP_RETCODE SCIPapplyHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< pointer to store result of: solution found, no solution found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_SOL*             resultsol           /**< a solution where to store found solution values, if any, or NULL if to try adding to SCIP */
   );

/** updates the starting point for the NLP heuristic
 *
 * Is called, for example, by a constraint handler that handles nonlinear constraints when a check on feasibility of a solution fails.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPupdateStartpointHeurSubNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< subNLP heuristic */
   SCIP_SOL*             solcand,            /**< solution candidate */
   SCIP_Real             violation           /**< constraint violation of solution candidate */
   );

/** gets startpoint candidate to be used in next call to NLP heuristic, or NULL if none */
SCIP_EXPORT
SCIP_SOL* SCIPgetStartCandidateHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /*HEUR_SUBNLP_H_*/
