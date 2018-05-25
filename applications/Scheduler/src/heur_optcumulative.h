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

/**@file   heur_optcumulative.h
 * @brief  heuristic for cumulative scheduling with optional activities
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_OPTCUMULATIVE_H__
#define __SCIP_HEUR_OPTCUMULATIVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the clique primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurOptcumulative(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** initialize the heuristics data structure */
extern
SCIP_RETCODE SCIPinitHeurOptcumulative(
   SCIP*                 scip,               /**< original SCIP data structure */
   int                   nmachines,          /**< number of machines */
   int                   njobs,              /**< number of njobs */
   int*                  machines,           /**< number of jobs for each machines */
   SCIP_VAR***           binvars,            /**< machnine job matrix (choice variables) */
   SCIP_VAR***           vars,               /**< machnine job matrix (start time variables) */
   int**                 durations,          /**< machnine job duration matrix */
   int**                 demands,            /**< machnine job demands matrix */
   int*                  capacities          /**< machine capacities */
   );

#ifdef __cplusplus
}
#endif

#endif
