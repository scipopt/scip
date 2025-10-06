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

/**@file   heur_init.h
 * @brief  initial primal heuristic for the vertex coloring problem
 * @author Gerald Gamrath
 *
 * This file implements a heuristic which computes a starting solution for the coloring problem. It
 * therefore computes maximal stable sets and creates one variable for each set, which is added to the
 * LP.
 *
 * The heuristic is called only one time: before solving the root node.
 *
 * It checks whether a solution-file was read in and there already is a starting solution.  If this
 * is not the case, an initial possible coloring is computed by a greedy method.  After that, a
 * tabu-search is called, which tries to reduce the number of colors needed. The tabu-search algorithm
 * follows the description in
 *
 * "A Survey of Local Search Methods for Graph Coloring"@n
 * by P. Galinier and A. Hertz@n
 * Computers & Operations Research, 33 (2006)
 *
 * The tabu-search works as follows: given the graph and a number of colors it tries to color the
 * nodes of the graph with at most the given number of colors.  It starts with a random coloring. In
 * each iteration, it counts the number of violated edges, that is, edges for which both incident
 * nodes have the same color. It now switches one node to another color in each iteration, taking
 * the node and color, that cause the greatest reduction of the number of violated edges, or if no
 * such combination exists, the node and color that cause the smallest increase of that number.  The
 * former color of the node is forbidden for a couple of iterations in order to give the possibility
 * to leave a local minimum.
 *
 * As long as the tabu-search finds a solution with the given number of colors, this number is reduced
 * by 1 and the tabu-search is called another time. If no coloring was found after a given number
 * of iterations, the tabu-search is stopped and variables for all sets of the last feasible coloring
 * are created and added to the LP (after possible extension to maximal stable sets).
 *
 * The variables of these sets result in a feasible starting solution of the coloring problem.
 *
 * The tabu-search can be deactivated by setting the parameter <heuristics/initcol/usetabu> to
 * FALSE.  The number of iterations after which the tabu-search stops if no solution was yet found
 * can be changed by the param <heuristics/initcol/maxiter>. A great effect is also obtained by
 * changing the parameters <heuristics/initcol/tabubase> and <heuristics/initcol/tabugamma>, which
 * distinguish the number of iterations for which the former color of a node is forbidden; more
 * precisely, this number is \<tabubase\> + ncritical * \<tabugamma\>, where ncritical is the number
 * of nodes, which are incident to violated edges.  Finally, the level of output and the frequency of
 * status lines can be changed by <heuristics/initcol/output> and <heuristics/initcol/dispfreq>.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_INIT_H__
#define __SCIP_HEUR_INIT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the initial primal heuristic for coloring and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurInit(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
