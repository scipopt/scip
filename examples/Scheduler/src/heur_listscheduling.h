/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_listscheduling.h
 * @brief  scheduling specific primal heuristic which is based on bidirectional serial generation scheme.
 * @author Eamonn Coughlan
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LISTSCHEDULING_H__
#define __SCIP_HEUR_LISTSCHEDULING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the list scheduling primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurListScheduling(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** initialize heuristic */
extern
SCIP_RETCODE SCIPinitializeHeurListScheduling(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIGRAPH*         precedencegraph,    /**< precedence graph */
   SCIP_VAR**            vars,               /**< start time variables */
   int*                  durations,          /**< duration of the jobs independent of the resources */
   int**                 resourcedemands,    /**< resource demand matrix */
   int*                  capacities,         /**< resource capacities */
   int                   njobs,              /**< number if jobs */
   int                   nresources          /**< number of resources */
   );

#ifdef __cplusplus
}
#endif

#endif
