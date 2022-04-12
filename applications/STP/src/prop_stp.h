/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_stp.h
 * @brief  propagator for Steiner tree problems, using the LP reduced costs
 * @author Daniel Rehfeldt
 *
 * This propagator makes use of the reduced cost of an optimally solved LP relaxation to propagate the variables, see
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by
 * Gamrath, Koch, Maher, Rehfeldt and Shinano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_STP_H__
#define __SCIP_PROP_STP_H__

#include <stdio.h>
#include <stdlib.h>
#include "scip/scip.h"
#include "graph.h"
#include "probdata_stp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the stp propagator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropStp(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** fix a variable (corresponding to an edge) to 0 */
SCIP_EXPORT
SCIP_RETCODE SCIPStpFixEdgeVarTo0(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   SCIP_Bool*            success             /**< could variable be fixed? */
   );


/** fix a variable (corresponding to an edge) to 1 */
SCIP_EXPORT
SCIP_RETCODE SCIPStpFixEdgeVarTo1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   SCIP_Bool*            success             /**< could variable be fixed? */
   );


/** return total number of arcs fixed by 'fixedgevar' method of this propagator */
SCIP_EXPORT
int SCIPStpNfixedEdges(
   SCIP*                 scip                /**< SCIP data structure */
   );


/** checks whether problem has become infeasible at current node */
SCIP_EXPORT
SCIP_RETCODE SCIPStpPropCheckForInfeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            probisinfeas        /**< is infeasible? */
   );

/** gets propagator graph  */
SCIP_EXPORT
SCIP_RETCODE SCIPStpPropGetGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data */
   SCIP_Longint*         graphnodenumber,    /**< point to b&b node for which graph is valid */
   SCIP_Bool*            probisinfeas,       /**< infeasible problem? */
   SCIP_Real*            offset              /**< needed for PC/MW */
   );

/** gives array indicating which nodes are degree-2 bounded */
SCIP_EXPORT
const SCIP_Bool* SCIPStpPropGet2BoundedArr(
   SCIP*                 scip                /**< SCIP data structure */
);

#ifdef __cplusplus
}
#endif

#endif
