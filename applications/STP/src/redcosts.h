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


/** parameters */
typedef struct reduced_costs_parameters
{
   SCIP_Real             cutoff;             /**< FOR FIRST LEVEL: reduced cost cutoff value or -1.0 if not used */
   int                   nLevels;            /**< number of of levels */
   int                   nCloseTerms;        /**< number of close terminals: 1, 2, or 3 */
   int                   nnodes;             /**< number of nodes */
   int                   nedges;             /**< number of edges */
   int                   redCostRoot;        /**< FOR FIRST LEVEL: graph root for reduced cost calculation */
} RCPARAMS;



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


/** returns reduced costs */
EXTERN
SCIP_Real* redcosts_getEdgeCosts(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get reduced costs for */
   );

/** returns top level reduced costs */
EXTERN
SCIP_Real* redcosts_getEdgeCostsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns root to node distances */
EXTERN
SCIP_Real* redcosts_getRootToNodeDist(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get distances for */
   );

/** returns root to node distances for top level */
EXTERN
SCIP_Real* redcosts_getRootToNodeDistTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns paths from nodes to closes terms */
EXTERN
PATH* redcosts_getNodeToTermsPaths(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get reduced costs for */
   );


/** returns paths from nodes to closes terms for top level */
EXTERN
PATH* redcosts_getNodeToTermsPathsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns closest terminals to nodes */
int* redcosts_getNodeToTermsBases(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get terminals for */
   );


/** returns closest terms to nodes for top level */
EXTERN
int* redcosts_getNodeToTermsBasesTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns cutoff */
EXTERN
SCIP_Real redcosts_getCutoff(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get cutoff for */
   );


/** returns cutoff for top level */
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
int redcosts_getRoot(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get root for */
   );


/** returns root used for reduced cost computation for top level */
EXTERN
int redcosts_getRootTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns current (top) level; 0-indexed */
EXTERN
int redcosts_getLevel(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns current number of levels */
int redcosts_getNlevels(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );

/** sets cutoff */
EXTERN
void redcosts_setCutoff(
   int                 level,               /**< level to set cutoff for */
   SCIP_Real           cutoff,             /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   );


/** sets cutoff for top level */
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


/** adds a new level */
EXTERN
void redcosts_addLevel(
   REDCOST*              redcostdata         /**< reduced costs data */
   );


/** initializes reduced costs data structure from given parameter struct */
EXTERN
SCIP_RETCODE redcosts_initFromParams(
   SCIP*                 scip,               /**< SCIP */
   const RCPARAMS*       parameters,         /**< parameters for initialization */
   REDCOST**             redcostdata         /**< data to initialize */
);

/** initializes reduced costs data structure.
 *  DEPRECATED! Use redcosts_initFromParams */
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
SCIP_RETCODE redcosts_initializeDistancesTop(
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
