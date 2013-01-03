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

/**@file   xternal.c
 * @brief  main document page
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
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
 */
