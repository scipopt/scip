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

/**@file   shortestpaths.c
 * @brief  Shortest path based algorithms for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file encompasses various shortest path based algorithms.
 * Note: This file is supposed to replace graph_path.c in the long run, as it includes the faster implementations.
 *
 */


#ifndef APPLICATIONS_STP_SRC_SHORTESTPATH_H_
#define APPLICATIONS_STP_SRC_SHORTESTPATH_H_

#include "graph.h"

/** information needed for prize-collecting problems */
typedef
struct stortest_path_prizecollecting
{
   SCIP_Real*            orderedprizes;      /**< ordered prizes for (pseudo) terminals */
   int*                  orderedprizes_id;   /**< ordered prizes IDs */
   SCIP_Real             maxoutprize;        /**< maximum prize of not yet connected vertex */
   int                   maxoutprize_idx;    /**< index */
} SPATHSPC;


/** information for shortest paths */
typedef
struct stortest_paths
{
   const CSR*            csr;                /**< CSR */
   DHEAP*                dheap;              /**< Dijkstra heap */
   SCIP_Real* RESTRICT   nodes_dist;         /**< distance array (on vertices) */
   int* RESTRICT         nodes_pred;         /**< predecessor node array (on vertices)
                                                  NOTE: might contain uninitialized values in opt mode! */
   STP_Bool* RESTRICT    nodes_isConnected;  /**< array to mark whether a vertex is part of computed Steiner tree */
} SPATHS;

/* todo should be made static, currently kept global for legacy reasons */
void shortestpath_pcStart(SPATHSPC*);
void shortestpath_pcConnectNode(const GRAPH*, const STP_Bool*, int, SPATHSPC*);

void shortestpath_computeSteinerTree(const GRAPH*, int, SPATHS*);
void shortestpath_computeSteinerTreePcMw(const GRAPH*, int, const SCIP_Real*, SCIP_Bool, SPATHSPC*, SPATHS*);
void shortestpath_computeSteinerTreePcMwFull(const GRAPH*, int, SPATHS*);


#endif /* APPLICATIONS_STP_SRC_SHORTESTPATH_H_ */
