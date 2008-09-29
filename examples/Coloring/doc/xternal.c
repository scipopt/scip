/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2008 Konrad-Zuse-Zentrum                            */
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

/**@mainpage Coloring Example
 * @version  0.1
 * @author   Gerald Gamrath
 *
 * This branch-and-price coloring code gives an example for 
 * a pricer and associated modules.
 *
 * It implements the approach described in
 *
 * "A column generation approach for graph coloring"@n
 * by Anuj Mehrotra and Micheal A. Trick,@n
 * INFORMS J. Comput. 8, no. 4, 1995, pp. 344-354.
 *
 * The input format for the graph files is the DIMACS standard format, the file must end with ".col".
 *
 * The coloring problem is the following:
 *
 * Given a graph \f$G = (V,E)\f$, the goal is to assign a color to each vertex,
 * such that no adjacent vertices have the same color. At the same time
 * minimize the number of colors needed.
 *
 * We use the following integer programming model: We have binary
 * variables \f$ x_{s}, s \in \mathcal{S}\f$ where \f$\mathcal{S}\f$ 
 * is the set of all (inclusion-) maximal stable sets in the graph \f$G\f$.
 * 
 * 
 * The basic model is then:
 * \f[
 *  \begin{array}[t]{rl}
 *    \min & \sum_{s \in \mathcal{S}} x_{s} \\
 *         & \\
 *    s.t. & \sum_{s \in \mathcal{S}, v \in s} x_{s} \ge 1 \quad \forall v \in V \\
 *  \end{array}
 * \f]
 *
 * Since the number of stable sets is exponential in most cases, the algorithm starts
 * with some subset \f$ \bar{\mathcal{S}} \subseteq \mathcal{S}\f$ of the stable sets and adds further
 * stable sets during the solution process if they can improve the current LP solution.
 *
 * Further information about particular modules like the pricing routine and the 
 * branching rule can be found in the corresponding (pricer_coloring.c and branch_coloring.c).
 *
 */
