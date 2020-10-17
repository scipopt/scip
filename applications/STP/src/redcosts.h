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
 *
 * Reduced costs are stored level-wise
 *
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef APPLICATIONS_STP_SRC_REDCOSTS_H_
#define APPLICATIONS_STP_SRC_REDCOSTS_H_

#include "scip/scip.h"
#include "graph.h"

#ifdef __cplusplus
extern "C" {
#endif

/** reduced cost data */
typedef struct reduce_costs_data REDCOST;


/** returns number of nodes for which reduced costs are stored */
EXTERN
int redcosts_getNnodes(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns number of edges for which reduced costs are stored */
EXTERN
int redcosts_getNedges(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns top level reduced costs */
EXTERN
SCIP_Real* redcosts_getEdgeCostsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns root to node distances */
EXTERN
SCIP_Real* redcosts_getRootToNodeDistTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns paths from nodes to closes terms */
EXTERN
PATH* redcosts_getNodeToTermsPathsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns closest terms to nodes */
EXTERN
int* redcosts_getNodeToTermsBasesTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns cutoff */
EXTERN
SCIP_Real redcosts_getCutoffTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns dual-bound */
EXTERN
SCIP_Real redcosts_getDualBoundTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns root used for reduced cost computation */
EXTERN
int redcosts_getRootTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** sets cutoff */
EXTERN
void redcosts_setCutoffTop(
   SCIP_Real           cutoff,             /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   );


/** sets dual-bound */
EXTERN
void redcosts_setDualBoundTop(
   SCIP_Real           dualbound,          /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   );


/** sets root used for reduced cost computation */
EXTERN
void redcosts_setRootTop(
   int                   root,               /**< the root */
   REDCOST*              redcostdata         /**< reduced costs data */
   );


/** initializes reduced costs data structure */
EXTERN
SCIP_RETCODE redcosts_init(
   SCIP*                 scip,               /**< SCIP */
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   SCIP_Real             cutoff,             /**< reduced cost cutoff value or -1.0 if not used */
   int                   redCostRoot,        /**< graph root for reduced cost calculation */
   REDCOST**             redcostdata         /**< reduced costs data */
   );


/** frees */
EXTERN
void redcosts_free(
   SCIP*                 scip,               /**< SCIP */
   REDCOST**             redcostdata         /**< data to free */
);


/** sets cutoff */
EXTERN
void redcosts_setAndReturnCutoffFromBound(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata,        /**< reduced cost data */
  SCIP_Real*            cutoffbound         /**< cutoff */
);


/** sets cutoff */
void redcosts_setCutoffFromBound(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata         /**< reduced cost data */
);

/** increases reduced cost for deleted arcs */
EXTERN
void redcosts_increaseOnDeletedArcs(
   const GRAPH*          graph,              /**< graph */
   const STP_Bool*       arcsdeleted,        /**< array to mark deleted arcs */
   REDCOST*              redcostdata         /**< reduced cost data */
);


/* initialize distances from reduced costs */
EXTERN
SCIP_RETCODE redcosts_initializeDistances(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata         /**< reduced cost data */
   );


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



#ifdef __cplusplus
}
#endif

#endif /* APPLICATIONS_STP_SRC_REDCOSTS_H_ */
