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
#include "grph.h"
#include "probdata_stp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the stp propagator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludePropStp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** fix a variable (corresponding to an edge) to zero */
EXTERN
SCIP_RETCODE fixedgevar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             edgevar,            /**< the variable to be fixed */
   int*                  nfixed              /**< counter that is incriminated if variable could be fixed */
   );

/** return total number of arcs fixed by 'fixedgevar' method of this propagator */
EXTERN
int SCIPStpNfixedEdges(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets propagator graph  */
EXTERN
void SCIPStpPropGetGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH**               graph,              /**< graph data */
   SCIP_Longint*         graphnodenumber     /**< pointer to b&b node for which graph is valid */
   );

#ifdef __cplusplus
}
#endif

#endif
