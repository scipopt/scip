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
#include "graphdefs.h"

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
SCIP_EXPORT
int redcosts_getNnodes(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns number of edges for which reduced costs are stored */
SCIP_EXPORT
int redcosts_getNedges(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns reduced costs */
SCIP_EXPORT
SCIP_Real* redcosts_getEdgeCosts(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get reduced costs for */
   );

/** returns top level reduced costs */
SCIP_EXPORT
SCIP_Real* redcosts_getEdgeCostsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns root to node distances */
SCIP_EXPORT
SCIP_Real* redcosts_getRootToNodeDist(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get distances for */
   );

/** returns root to node distances for top level */
SCIP_EXPORT
SCIP_Real* redcosts_getRootToNodeDistTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns paths from nodes to closes terms */
SCIP_EXPORT
PATH* redcosts_getNodeToTermsPaths(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get reduced costs for */
   );


/** returns paths from nodes to closes terms for top level */
SCIP_EXPORT
PATH* redcosts_getNodeToTermsPathsTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns closest terminals to nodes */
int* redcosts_getNodeToTermsBases(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get terminals for */
   );


/** returns closest terms to nodes for top level */
SCIP_EXPORT
int* redcosts_getNodeToTermsBasesTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns cutoff */
SCIP_EXPORT
SCIP_Real redcosts_getCutoff(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get cutoff for */
   );


/** returns cutoff for top level */
SCIP_EXPORT
SCIP_Real redcosts_getCutoffTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns dual-bound */
SCIP_EXPORT
SCIP_Real redcosts_getDualBound(
   int                   level,              /**< level */
   const REDCOST*        redcostdata         /**< reduced costs data */
   );

/** returns dual-bound for top level*/
SCIP_EXPORT
SCIP_Real redcosts_getDualBoundTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns root used for reduced cost computation */
int redcosts_getRoot(
   const REDCOST*        redcostdata,        /**< reduced costs data */
   int                   level               /**< level to get root for */
   );


/** returns root used for reduced cost computation for top level */
SCIP_EXPORT
int redcosts_getRootTop(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns current (top) level; 0-indexed */
SCIP_EXPORT
int redcosts_getLevel(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );


/** returns current number of levels */
int redcosts_getNlevels(
   const REDCOST*        redcostdata         /**< reduced costs data */
   );

/** sets cutoff */
SCIP_EXPORT
void redcosts_setCutoff(
   SCIP_Real           cutoff,             /**< the value */
   int                 level,               /**< level to set cutoff for */
   REDCOST*            redcostdata         /**< reduced costs data */
   );


/** sets cutoff for top level */
SCIP_EXPORT
void redcosts_setCutoffTop(
   SCIP_Real           cutoff,             /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   );



/** sets dual-bound */
SCIP_EXPORT
void redcosts_setDualBound(
   SCIP_Real           dualbound,          /**< the value */
   int                 level,              /**< level to set dual bound for */
   REDCOST*            redcostdata         /**< reduced costs data */
   );


/** sets dual-bound for top level */
SCIP_EXPORT
void redcosts_setDualBoundTop(
   SCIP_Real           dualbound,          /**< the value */
   REDCOST*            redcostdata         /**< reduced costs data */
   );


/** sets root used for reduced cost computation */
SCIP_EXPORT
void redcosts_setRoot(
   int                   root,               /**< the root */
   int                   level,              /**< level to set dual bound for */
   REDCOST*              redcostdata         /**< reduced costs data */
   );


/** sets root used for reduced cost computation for top level */
SCIP_EXPORT
void redcosts_setRootTop(
   int                   root,               /**< the root */
   REDCOST*              redcostdata         /**< reduced costs data */
   );


/** adds a new level */
SCIP_EXPORT
void redcosts_addLevel(
   REDCOST*              redcostdata         /**< reduced costs data */
   );


/** initializes reduced costs data structure from given parameter struct */
SCIP_EXPORT
SCIP_RETCODE redcosts_initFromParams(
   SCIP*                 scip,               /**< SCIP */
   const RCPARAMS*       parameters,         /**< parameters for initialization */
   REDCOST**             redcostdata         /**< data to initialize */
);

/** initializes reduced costs data structure.
 *  DEPRECATED! Use redcosts_initFromParams */
SCIP_EXPORT
SCIP_RETCODE redcosts_init(
   SCIP*                 scip,               /**< SCIP */
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   SCIP_Real             cutoff,             /**< reduced cost cutoff value or -1.0 if not used */
   int                   redCostRoot,        /**< graph root for reduced cost calculation */
   REDCOST**             redcostdata         /**< reduced costs data */
   );


/** frees */
SCIP_EXPORT
void redcosts_free(
   SCIP*                 scip,               /**< SCIP */
   REDCOST**             redcostdata         /**< data to free */
);


/** sets cutoff for top level */
SCIP_EXPORT
void redcosts_setAndReturnCutoffFromBoundTop(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata,        /**< reduced cost data */
  SCIP_Real*            cutoffbound         /**< cutoff */
);


/** sets cutoff */
SCIP_EXPORT
void redcosts_setCutoffFromBound(
  SCIP_Real             upperbound,         /**< bound */
  int                   level,              /**< level */
  REDCOST*              redcostdata         /**< reduced cost data */
);


/** sets cutoff for top level */
SCIP_EXPORT
void redcosts_setCutoffFromBoundTop(
  SCIP_Real             upperbound,         /**< bound */
  REDCOST*              redcostdata         /**< reduced cost data */
);



/** increases reduced cost for deleted arcs */
SCIP_EXPORT
void redcosts_increaseOnDeletedArcs(
   const GRAPH*          graph,              /**< graph */
   const STP_Bool*       arcsdeleted,        /**< array to mark deleted arcs */
   int                   level,              /**< the level */
   REDCOST*              redcostdata         /**< reduced cost data */
);


/** unifies costs */
SCIP_EXPORT
void redcosts_unifyBlockedEdgeCosts(
   const GRAPH*          graph,              /**< graph */
   int                   level,              /**< the level */
   REDCOST*              redcostdata         /**< reduced cost data */
);


/** increases reduced cost for deleted arcs for top level */
SCIP_EXPORT
void redcosts_increaseOnDeletedArcsTop(
   const GRAPH*          graph,              /**< graph */
   const STP_Bool*       arcsdeleted,        /**< array to mark deleted arcs */
   REDCOST*              redcostdata         /**< reduced cost data */
);


/* initialize distances from reduced costs */
SCIP_EXPORT
SCIP_RETCODE redcosts_initializeDistances(
   SCIP*                 scip,               /**< SCIP */
   int                   level,              /**< level to inizialize for*/
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata         /**< reduced cost data */
   );


/* initialize top distances from reduced costs */
SCIP_EXPORT
SCIP_RETCODE redcosts_initializeDistancesTop(
   SCIP*                 scip,               /**< SCIP */
   GRAPH*                g,                  /**< graph data structure */
   REDCOST*              redcostdata         /**< reduced cost data */
   );


/** reduced costs available? */
SCIP_EXPORT
SCIP_Bool redcosts_forLPareAvailable(
   SCIP*                 scip                /**< SCIP structure */
   );


/** are reduced costs reliable? */
SCIP_EXPORT
SCIP_Bool redcosts_forLPareReliable(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_VAR**            vars,               /**< variables (in) */
   const GRAPH*          graph               /**< graph data */
   );

/** initialize reduced costs */
SCIP_EXPORT
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
