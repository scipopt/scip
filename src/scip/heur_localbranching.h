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

/**@file   heur_localbranching.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Local branching heuristic according to Fischetti and Lodi
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCALBRANCHING_H__
#define __SCIP_HEUR_LOCALBRANCHING_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates local branching primal heuristic and includes it in SCIP
 *
 *  @ingroup PrimalHeuristicIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeHeurLocalbranching(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the working limits for the auxiliary problem */
EXTERN
SCIP_RETCODE SCIPlocalbranchingSetWorkingLimits(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< the subproblem created by localbranching */
   SCIP_Longint          nsubnodes,          /**< nodelimit for subscip */
   int                   bestsollimit,       /**< the limit on the number of best solutions found */
   SCIP_Bool             useuct              /**< should the uct node selector be used? */
   );

/** checks the solutions from the subscip and adds them to the master SCIP is feasible */
EXTERN
SCIP_RETCODE SCIPlocalbranchingCheckAndAddSolution(
   SCIP*                 scip,               /**< the SCIP data structure */
   SCIP*                 subscip,            /**< the subproblem created by localbranching */
   SCIP_HEUR*            heur,               /**< localbranching heuristic */
   SCIP_VAR**            subvars,            /**< the variables from the subproblem */
   SCIP_RESULT*          result              /**< result pointer */
   );

#ifdef __cplusplus
}
#endif

#endif
