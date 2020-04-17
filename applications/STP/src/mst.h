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

/**@file   mst.h
 *
 *
 * @brief  minimum spanning tree based algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various minimum spanning tree based algorithms.
 * Note: This file is supposed to (partly) replace graph_path.c in the long run, as it includes the faster implementations.
 *
 */
#ifndef APPLICATIONS_STP_SRC_MST_H_
#define APPLICATIONS_STP_SRC_MST_H_

#include "graph.h"
#include "completegraph.h"


/** information for (sparse) MST computations */
typedef
struct minimum_spanning_tree
{
   const CSR*            csr;                /**< CSR */
   DHEAP*                dheap;              /**< Dijkstra heap */
   SCIP_Real* RESTRICT   nodes_dist;         /**< distance array (on vertices) */
   int* RESTRICT         nodes_predEdge;     /**< predecessor edge array (on vertices);
                                                NOTE: with respect to original graph edge IDs
                                                NOTE: might contain uninitialized values in opt mode! */
} MST;


void mst_computeOnMarked(const GRAPH*, int, MST*);


#endif /* APPLICATIONS_STP_SRC_MST_H_ */
