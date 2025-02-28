/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_strongcoloring.h
 * @brief  branching rule performing strong branching for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements an additional branching rule for the coloring algorithm.
 *
 * We are looking for two nodes v and w, which are not adjacent in the current graph, and consider
 * the following two constraints: SAME(v,w) and DIFFER(v,w). More information about the meaning of
 * these constraints can be found in the documentation of the branching rule in branch_coloring.c.
 *
 * This branching rule puts some more effort into the choice of the two nodes and performs a
 * strongbranching. This means that for every possible choice of two nodes, it solves the LPs of the
 * created children and computes a score with respect to the increase of the lower bound in both
 * nodes. After that, it takes the combination of nodes yielding the best score. The interesting
 * point is that the strongbranching is not performed for each variable, as it is done in some
 * default branching rules of SCIP and supported by the LP-solver, but is done for a constraint,
 * since we are branching on constraints. Look at executeStrongBranching() to see how it is
 * done. There are also some improvements, since testing all possible combination of nodes is very
 * expensive.  The first possibility to avoid this is to stop the computation of scores once a
 * possible branching is found that has only one feasible child. This results in more restrictions
 * in this child without increasing the number of unprocessed nodes.
 *
 * The second improvement is to compute a priority for all possible combinations, w.r.t. the
 * fractional values of the variables. Then, only the first best k combinations are investigated by
 * strongbranching.
 *
 * This code is not optimized and in most cases inferior to the standard branching rule. It is only
 * a demonstration of how to perform strongbranching on constraints!
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_STRONGCOLORING_H__
#define __SCIP_BRANCH_STRONGCOLORING_H__


#include "scip/scip.h"
#include "probdata_coloring.h"
#include "cons_storeGraph.h"
#include "scip/cons_linear.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the coloring branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleStrongcoloring(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
