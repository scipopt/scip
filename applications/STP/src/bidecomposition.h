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

/**@file   bidecomposition.h
 * @brief  several decomposition methods for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_BIDECOMPOSITION_H_
#define APPLICATIONS_STP_SRC_BIDECOMPOSITION_H_

#include "scip/scip.h"
#include "graph.h"
#include "stpvector.h"


#ifdef __cplusplus
extern "C" {
#endif


//#define CUTTREE_PRINT_STATISTICS


/** decompose */
typedef struct biconnected_component_decomposition
{
   SUBINOUT*             subinout;           /**< helper */
   int*                  nodes;              /**< nodes */
   int*                  starts;             /**< starts array per component */
   int                   nbicomps;           /**< number of components */
} BIDECOMP;



/** for internal computation */
typedef struct biconnected_stack_node STACK_NODE;


/** cut nodes/ articulation points todo: hide */
typedef struct cut_nodes
{
   SCIP*                 scip;               /**< SCIP data structure */
#ifdef CUTTREE_PRINT_STATISTICS
   int*                  childcount_nodes;   /**< number of nodes below each node */
   int*                  childcount_terms;   /**< number of terminals below each node */
#endif
   STACK_NODE*           stack_nodes;        /**< data for iterative computation */
   STP_Vectype(int)      biconn_stack;       /**< stack for marking bi-connected component */
   int*                  biconn_nodesmark;   /**< marks in which component each node is 0, 1,.., biconn_ncomps - 1 */
   int*                  biconn_comproots;   /**< root of each component with index 0,1,...,biconn_ncomps - 1 */
   STP_Vectype(int)      artpoints;          /**< cut nodes */
   int*                  nodes_hittime;      /**< hit time 0,1,... */
   int                   stack_size;         /**< size of stack */
   int                   biconn_ncomps;      /**< number of components */
   int                   dfsroot;            /**< root */
   int                   nrootcomps;         /**< number of root components */
   int                   curr_lowpoint;      /**< current low-point */
   int                   curr_hittime;       /**< current hit time */
} CUTNODES;


extern SCIP_RETCODE bidecomposition_cutnodesInit(SCIP*, const GRAPH*, CUTNODES**);
extern void         bidecomposition_cutnodesFree(SCIP*, CUTNODES**);
extern void         bidecomposition_cutnodesCompute(const GRAPH*, CUTNODES*);

extern SCIP_RETCODE bidecomposition_init(SCIP*, const CUTNODES*, const GRAPH*, BIDECOMP**);
extern SCIP_RETCODE bidecomposition_initSubInOut(SCIP*, const GRAPH*, BIDECOMP*);
extern void         bidecomposition_free(SCIP*, BIDECOMP**);
extern void         bidecomposition_markSub(const BIDECOMP*, int, GRAPH*);
extern SCIP_Bool    bidecomposition_componentIsTrivial(const BIDECOMP*, int);
extern SCIP_Bool    bidecomposition_isPossible(const GRAPH*);
extern SCIP_RETCODE bidecomposition_getMarkedSubRoot(SCIP*, const BIDECOMP*, const GRAPH*, const GRAPH*, int*);
extern SCIP_Real    bidecomposition_getCompNodeRatio(const BIDECOMP*, int);
extern SCIP_Real    bidecomposition_getMaxcompNodeRatio(const BIDECOMP*);



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_BIDECOMPOSITION_H_ */
