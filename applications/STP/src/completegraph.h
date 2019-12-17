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

/**@file   completegraph.c
 * @brief  includes complete graph methods, in particular for MST calculation
 * @author Daniel Rehfeldt
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_COMPLETEGRAPH_H_
#define APPLICATIONS_STP_SRC_COMPLETEGRAPH_H_

#include "scip/scip.h"
#include "graph.h"

/** complete, undirected graph storage */
typedef struct complete_graph
{
   SCIP_Real*            edgecosts;          /**< edge cost array; of size nnodes_max * (nnodes_max - 1) */
   int*                  nodeids;            /**< node ids; of size nnodes_max */
   int                   nnodes_max;         /**< maximum number of edges */
   int                   nnodes_curr;        /**< current number of nodes */
} CGRAPH;


/** complete, undirected graph MST storage */
typedef struct complete_mst
{
   DHEAP*                heap;               /**< heap needed for MST computation */
   int*                  predecessors;       /**< predecessor array of size nnodes_max */
   int                   nnodes_max;         /**< maximum number of edges */
   int                   nnodes_curr;        /**< current number of nodes */
} CMST;



#endif /* APPLICATIONS_STP_SRC_COMPLETEGRAPH_H_ */
