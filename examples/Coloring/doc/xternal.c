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
 * It implements the approach described in "A column generation approach
 * for graph coloring" by Anuj Mehrotra and Micheal A. Trick, April 11, 1995.
 * The input format for the graph files is the DIMACS standard format, the file must end with .col.
 *
 * The coloring problem is the following:
 *
 * Given a graph G = (V,E), the goal is to assign a color to each vertex,
 * such that no adjacent vertices have the same color. At the same time
 * minimize the number of colors needed.
 *
 * We use the following integer programming model: We have binary
 * variables \f$ x_{s}, s \in S\f$ where \f$S \subseteq \mathcal{P}(V)\f$ 
 * is the set of all (inclusion-) maximal stable sets in the graph G.
 * 
 * 
 * The basic model is then:
 * \f[
 *  \begin{array}[t]{rl}
 *    \min & \sum_{s \in S} x_{s} \\
 *         & \\
 *    s.t. & \sum_{s \in S, v \in s} x_{s} \ge 1 \quad \forall v \in V \\
 *  \end{array}
 * \f]
 *
 * Since the number of stable sets is exponential in most cases, the algorithm starts
 * with some subset \f$ \bar{S} \subseteq S\f$ of the stable sets and adds further
 * stable sets during the solution process if they can improve the current LP solution.
 *
 * Further information about particular modules like the pricing routine and the 
 * branching rule can be found in the corresponding .c-files.
 *
 */
