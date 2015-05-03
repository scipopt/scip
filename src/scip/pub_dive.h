
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

/**@file   scip.c
 * @brief  library methods for diving heuristics
 * @author Gregor Hendel
 *
 * @todo check all checkStage() calls, use bit flags instead of the SCIP_Bool parameters
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

#ifndef PUB_DIVE_H_
#define PUB_DIVE_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** performs a diving within the limits of the diveset parameters
 *
 *  This method performs a diving according to the settings defined by the diving settings @p diveset; Contrary to the
 *  name, SCIP enters probing mode (not diving mode) and dives along a path into the tree. Domain propagation
 *  is applied at every node in the tree, whereas probing LPs might be solved less frequently.
 *
 *  Starting from the current LP candidates, the algorithm determines a fraction of the candidates that should be
 *  branched on; if a single candidate should be fixed, the algorithm selects a candidate which minimizes the score
 *  defined by the @p diveset.  If more than one candidate should be selected, the candidates are sorted in
 *  non-decreasing order of their score.
 *
 *  The algorithm iteratively selects the the next (unfixed) candidate in the list, until the
 *  targeted depth is reached, or the last node is proven to be infeasible. It optionally backtracks and tries the
 *  other branching direction.
 *
 *  After the set of remaining candidates is empty or the targeted depth is reached, the node LP is
 *  solved, and the old candidates are replaced by the new LP candidates.
 *
 *  @see heur_guideddiving.c for an example implementation of a dive set controlling the diving algorithm.
 *
 *  @see the parameter @p heuristics/startdivefrac to determine the fraction of candidates that should be dived on at the
 *       beginning. Setting this parameter to 0.0 will result in an LP solved after every candidate selection.
 *
 *  @note the fraction of candidate variables is subject to change during solving. It is decreased by a factor of
 *        2 every time the algorithm could not dive half as deep as desired. However, if it succeeded, the fraction
 *        is multiplied by a factor of 1.1.
 *
 *  @note the node from where the algorithm is called is checked for a basic LP solution. If the solution
 *        is non-basic, e.g., when barrier without crossover is used, the method returns without performing a dive.
 *
 *  @note currently, when multiple diving heuristics call this method and solve an LP at the same node, only the first
 *        call will be executed, @see SCIPgetLastDiveNode().
 */
EXTERN
SCIP_RETCODE SCIPperformGenericDivingAlgorithm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIVESET*         diveset,            /**< settings for diving */
   SCIP_SOL*             worksol,            /**< non-NULL working solution */
   SCIP_HEUR*            heur,               /**< the calling primal heuristic */
   SCIP_RESULT*          result,             /**< SCIP result pointer */
   SCIP_Bool             nodeinfeasible      /**< is the current node known to be infeasible? */
   );

#ifdef __cplusplus
}
#endif

#endif /* PUB_DIVE_H_ */
