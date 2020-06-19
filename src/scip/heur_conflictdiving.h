/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_conflictdiving.h
 * @ingroup PRIMALHEURISTICS
 * @brief  LP diving heuristic that chooses fixings w.r.t. conflict locks
 * @author Jakob Witzig
 *
 * Diving heuristic: Iteratively fixes some fractional variable and resolves the LP-relaxation, thereby simulating a
 * depth-first-search in the tree. Conflict Diving chooses the variable with the fewest conflict locking number in any
 * direction and rounds it into this direction.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_CONFLICTDIVING_H__
#define __SCIP_HEUR_CONFLICTDIVING_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the conflictdiving heuristic and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurConflictdiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
