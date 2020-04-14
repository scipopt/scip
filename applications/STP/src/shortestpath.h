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


void shortestpath_computeSteinerTree(const GRAPH*, int, SCIP_Real* RESTRICT, int* RESTRICT, DHEAP*, STP_Bool* RESTRICT);


#endif /* APPLICATIONS_STP_SRC_SHORTESTPATH_H_ */
