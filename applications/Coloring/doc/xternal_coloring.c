/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   xternal_coloring.c
 * @brief  main document page
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page COLORING_MAIN Coloring
 * @version  0.1
 * @author   Gerald Gamrath
 *
 * This branch-and-price graph coloring code gives an example for
 * a pricer and associated modules.
 *
 * It implements the approach described in
 *
 * "A column generation approach for graph coloring"@n
 * by Anuj Mehrotra and Micheal A. Trick,@n
 * INFORMS J. Comput. 8, no. 4, 1995, pp. 344-354.
 *
 * The input format for the graph files is the DIMACS standard format; the name of the file must end with ".col".
 *
 * The graph coloring problem is the following:
 *
 * Given a graph \f$G = (V,E)\f$, the goal is to assign a color to each vertex, such that no
 * adjacent vertices have the same color; the number of colors needed should be minimized.
 *
 * We use the following integer programming model: We have binary
 * variables \f$ x_{s}, s \in \mathcal{S}\f$ where \f$\mathcal{S}\f$
 * is the set of all stable sets in the graph \f$G\f$.
 *
 * The basic model is then:
 * \f[
 *  \begin{array}[t]{rl}
 *    \min & \displaystyle \sum_{s \in \mathcal{S}} x_{s} \\
 *         & \\
 *    s.t. & \displaystyle \sum_{s \in \mathcal{S}, v \in s} x_{s} \ge 1 \quad \forall v \in V \\
 *  \end{array}
 * \f]
 *
 * Since the number of stable sets can be exponential in the size of the graph, the algorithm starts
 * with some subset \f$ \bar{\mathcal{S}} \subseteq \mathcal{S}\f$ of the stable sets and adds
 * further stable sets during the solution process. This way it tries to improve the current LP
 * solution.
 *
 * Further information about particular modules like the pricing routine and the
 * branching rule can be found in the documentation of the corresponding files.
 *
 * The pricer pricer_coloring.c shows how to perform column generation in SCIP.  The constraint
 * handler cons_storeGraph.c demonstrates how to store branching decisions at nodes und enforce them
 * by propagation.  The default branching rule branch_coloring.c describes how these constraints are
 * added to the branch-and-bound nodes.  Some more sophisticated approaches for the branching,
 * especially a strongbranching on these constraints, can be found in the second branching rule
 * branch_strongcoloring.c. An initial solution is computed by a start heuristic which is
 * described in heur_init.c.  The organization of the data for the problem is described in the
 * problem data file (probdata_coloring.c). The file readers are described in reader_col.c and
 * reader_csol.c.
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 *
 *
 */
