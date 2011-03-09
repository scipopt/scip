/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_trysol.h
 * @brief  primal heuristic that tries a given solution
 * @author Marc Pfetsch
 *
 * This heuristic is just a stand in for other components of SCIP that should not or cannot submit
 * solutions. For example, this heuristic might take a solution from a constraint handler that knows
 * how to produce feasible solutions. The heuristic then tries to submit the solution to SCIP.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_TRYSOL_H__
#define __SCIP_HEUR_TRYSOL_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the trysol primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurTrySol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** pass solution to trysol heuristic */
extern
SCIP_RETCODE SCIPheurPassSolTrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< trysol heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   );

#ifdef __cplusplus
}
#endif

#endif
