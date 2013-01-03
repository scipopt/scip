/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_coloring.h
 * @brief  default branching rule for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements the standard branching rule for the coloring algorithm.
 *
 * As we use column generation, we may not branch on the variables themselves,
 * but on some sort of constraints that we introduce in the pricing problem.
 *
 * In our case, we choose two nodes v and w, which are not adjacent in the current graph, and
 * consider the following two constraints: SAME(v,w) and DIFFER(v,w).  SAME(v,w) requires that both
 * nodes v and w get the same color, whereas DIFFER(v,w) forbids this. For each pair of nodes, each
 * feasible solution fulfills exactly one of these constraints. Hence, splitting the solution space
 * into two parts, one fulfilling SAME(v,w) and the other DIFFER(v,w), does not cut off any feasible
 * solution and can therefore be used as the branching rule.
 *
 * The branching is done as follows: Given the optimal (fractional) solution of the current
 * branch-and-bound node, choose the most fractional variable and the corresponding stable set
 * s1. Now choose two nodes v, w and another stable set s2, such that v is part of both stable sets,
 * whereas w is part of exactly one of the stable sets.  Create two children of the current node,
 * one with the restriction SAME(v,w), the other one with restriction DIFFER(v,w). Therefore, each
 * node gets a constraint of type @c cons_storeGraph, which enforces the branching decision and
 * assures that each coloring of the nodes in the respective subgraph assigns to both nodes the same
 * color/different colors by fixing stable sets to 0 that violate this constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_COLORING_H__
#define __SCIP_BRANCH_COLORING_H__


#include "scip/scip.h"
#include "probdata_coloring.h"
#include "cons_storeGraph.h"
#include "scip/cons_setppc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the coloring branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleColoring(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
