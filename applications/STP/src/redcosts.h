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

/**@file   redcosts.h
 * @brief  Reduced cost based routines for Steiner problems
 * @author Daniel Rehfeldt
 * *
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_REDCOSTS_H_
#define APPLICATIONS_STP_SRC_REDCOSTS_H_

#include "scip/scip.h"
#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** reduced costs available? */
EXTERN
SCIP_Bool redcosts_forLPareAvailable(
   SCIP*                 scip                /**< SCIP structure */
   );

/** initialize reduced costs */
EXTERN
void redcosts_forLPget(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables */
   const GRAPH*          graph,              /**< graph data */
   SCIP_Real*            redcosts            /**< reduced costs */
   );


/** reduced cost result data */
typedef struct reduce_costs_data
{
   SCIP_Real*            redEdgeCost;        /**< reduced costs */
   SCIP_Real*            rootToNodeDist;     /**< shortest path distances from root  */
   PATH*                 nodeTo3TermsPaths;  /**< paths to three nearest terminals */
   int*                  nodeTo3TermsBases;  /**< three nearest terminals */
   SCIP_Real             cutoff;             /**< reduced cost cutoff value or -1.0 if not used */
   SCIP_Real             dualBound;          /**< dual bound or -1.0 if not used */
   int                   redCostRoot;        /**< graph root for reduced cost calculation */
#ifndef NDEBUG
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
#endif
} REDCOST;


#ifdef __cplusplus
}
#endif

#endif /* APPLICATIONS_STP_SRC_REDCOSTS_H_ */
