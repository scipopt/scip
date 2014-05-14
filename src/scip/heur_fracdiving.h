/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_fracdiving.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP diving heuristic that chooses fixings w.r.t. the fractionalities
 * @author Tobias Achterberg
 *
 * Diving heuristic: Iteratively fixes some fractional variable and resolves the LP-relaxation, thereby simulating a
 * depth-first-search in the tree. Fractional Diving chooses the variable with the highest fractionality and rounds it to the
 * nearest integer. One-level backtracking is applied: If the LP gets infeasible, the last fixing is undone, and the
 * opposite fixing is tried. If this is infeasible, too, the procedure aborts.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_FRACDIVING_H__
#define __SCIP_HEUR_FRACDIVING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the fracdiving heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeHeurFracdiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** score candidate variable
 *
 *  if candidate cannot be trivially rounded,
 *  score is the difference between frac and [frac] (rounding to the nearest integer)
 *
 *  if candidate can be trivially rounded in at least one direction, the objective gain is used to score the variables
 *
 *  candidate which cannot be rounded trivially always have a lower score than roundable candidates
 */
EXTERN
SCIP_Real SCIPgetFracdivingVarScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             cand,               /**< diving candidate for score */
   SCIP_Real             frac,               /**< fractionality of candidate in (last) LP solution */
   SCIP_Bool             roundup             /**< should the variable be rounded up? */
   );

#ifdef __cplusplus
}
#endif

#endif
