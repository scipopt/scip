/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_hailmary.h
 * @ingroup PRIMALHEURISTICS
 * @brief  hailmary primal heuristic
 * @author Timo Berthold
 *
 * template file for primal heuristic plugins
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_HAILMARY_H__
#define __SCIP_HEUR_HAILMARY_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** main procedure of the hailmary heuristic, creates and solves a sub-SCIP */
extern
SCIP_RETCODE SCIPapplyHailmary(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minimprove,         /**< factor by which hailmary should at least improve the incumbent          */
   SCIP_Longint          nnodes              /**< node limit for the subproblem                                       */
   );

/** creates the hailmary primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurHailmary(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
