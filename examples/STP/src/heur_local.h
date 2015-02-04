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

/**@file   heur_local.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Improvement heuristic for STP
 * @author Daniel Rehfeldt
 *
 * This is an improvement heuristic.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCAL_H__
#define __SCIP_HEUR_LOCAL_H__


#include "scip/scip.h"
#include "grph.h"
#ifdef __cplusplus
extern "C" {
#endif

   /** creates the local primal heuristic and includes it in SCIP */
   extern
   SCIP_RETCODE SCIPincludeHeurLocal(
      SCIP*                 scip                /**< SCIP data structure */
      );

   extern
   SCIP_RETCODE do_local(
      SCIP*                 scip,               /**< SCIP data structure */
      const GRAPH*  graph,
      const SCIP_Real* cost,
      const SCIP_Real* costrev,
      int*          best_result
      );

#ifdef __cplusplus
}
#endif

#endif
