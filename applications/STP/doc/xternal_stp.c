/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal_stp.c
 * @brief  main document page
 * @author Daniel Rehfeldt
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Yuji Shinano
 * @author Michael Winkler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page STP_MAIN STP
 * @author Daniel Rehfeldt
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Yuji Shinano
 * @author Michael Winkler
 *
 * This application contains the branch-and-cut based solver SCIP-Jack for Steiner tree problems, realized within the framework
 * \SCIP, see: "A generic approach to solving the Steiner tree problem and variants" by D. Rehfeldt. The following plugins are implemented:
 *
 * - a problem reader, which parses the problem out of a .stp file
 *   (reader_stp.c)
 * - a (global) problem data structure, containing all necessary information including the graph, which creates the model within \SCIP (probdata_stp.c)
 * - shortest path based construction heuristics (heur_tm.c)
 * - reduction based construction heuristics (heur_prune.c, heur_ascendprune, heur_slackprune)
 * - improvment heuristics (heur_local.c)
 * - a recombination heuristic (heur_rec.c)
 * - basic reduction techniques (reduce_simple.c)
 * - alternative-based reduction techniques (reduce_alt.c)
 * - bound-based reduction techniques (reduce_bnd.c)
 * - a constraint handler, which checks solutions for feasibility and separates any violated model constraints (cons_stp.c)
 * - a propagator, which attempts to fix (edge) variables to zero utilizing their reduced costs (prop_stp.c)
 * - an event handler, which simply writes each incumbent solution to a file -- if activated (event_bestsol.c)
 *
 *
 *
 * In the following the problem is introduced and the solving process is delineated. Furthermore, the two main plugins are
 * sketched.
 *
 * -# \subpage STP_PROBLEM "Problem description and solving approach"
 * -# \subpage STP_PROBLEMDATA "Main problem data, creating the problem"
 * -# \subpage STP_CONS "Separating violated constraints"
 * -# \subpage STP_MINCUT "Graph minimum cut routine"
 * -# \subpage STP_MISC "Miscellaneous methods used for Steiner tree problems"
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 */

/**@page STP_PROBLEM Problem description and solving approach
 *
 * The Steiner tree problem in graphs (SPG) can be described as follows: Given an undirected connected graph
 * \f$ G=(V,E)\f$, costs \f[ c: E \rightarrow  \mathcal{Q}^+ \f] and a set \f$ T \subset V \f$ of terminals,
 * the problem is to find a minimum-weight tree \f$ S\subseteq G \f$ that spans \f$ T \f$. Each tree \f$ S \f$ that spans \f$ T \f$, called Steiner tree, is
 * a feasible solution to the problem.
 * The following picture shows an SPG instance with the terminals given as squares:
 *
 * <CENTER>
 * \image html stp.png
 * </CENTER>
 *
 *
 * The solving approach of SCIP-Jack can be dissected into three major components:
 *
 * First, problem specific preprocessing is extremely important. Apart from some pathological instances
 * specifically constructed to defy presolving techniques, preprocessing is often able to significantly
 * reduce instances or even solve them.
 *
 * Second, heuristics are needed, especially for hard instances, to find good or even optimal solutions.
 *
 * Finally, the core of our approach is constituted by a branch-and-cut procedure used to compute lower bounds and prove optimality.
 *
 * The problem can be formulated using the directed equivalent of the SPG, the Steiner arborescence problem (SAP):
 * Given a directed graph \f$ D=(V,A) \f$, a root \f$ r \in V \f$, costs \f$ c: A \rightarrow \mathcal{Q}^+ \f$
 * and a set \f$ T \subset V \f$ of terminals, a subgraph \f$ S\subseteq D \f$ such that
 * for all \f$ t \in T \f$, \f$ S \f$ contains exactly one directed path from \f$ r \f$ to \f$ t \f$ is called Steiner arborescence.
 * Thereupon, a Steiner arborescence \f$ S = (V_S, A_S) \f$ is required that minimizes \f$ \sum_{a \in A_S} c_a \f$.  An SPG can be
 * transformed to an SAP by replacing each edge with two anti-parallel arcs of the same cost and distinguishing an arbitrary
 * terminal as the root. This transformation results in a one-to-one correspondence between the respective solution sets.
 * Introducing variables \f$y_a\f$ for \f$a\in A\f$ with the interpretation that \f$y_a:=1\f$ if and only if \f$a\f$ is in the
 * Steiner arborecence, and \f$y_a:=0\f$ otherwise, we obtain the integer program:
 *
 * \f[
 *  \begin{array}[t]{rll}
 *   \min {c}^T y\\
 *  \\
 *   y(\delta^+(W))&\geq& 1, ~~~ \forall  , W\subset V, r\in W, (V\setminus W)\cap T\neq \emptyset\\
 *   y(\delta^-(v))&
 *   \left\{{\begin{array}{l} = \\
 *       = \\
 *       \leq
 *     \end{array}}\right.
 *   &
 *   {\begin{array}{l}
 *     0, \mbox{if } v=r;\\
 *     1, \mbox{if } v\in T\setminus{r};\\
 *     1, \mbox{if } v\in N; \end{array}} \hspace{2.9mm}\forall  v \in V
 *   \\
 *   y(\delta^-(v))&\leq& y(\delta^+(v)), \hspace{10.5mm}\forall  v\in N;\\
 *   y(\delta^-(v))&\geq& y_a, \hspace{20.2mm}\forall  a\in\delta^+(v), v\in N;\\
 *  0\leq y_a&\leq& 1, \hspace{22mm}\forall  a\in A;\\
 *   y_a&\in& \{0,1\}, \hspace{15.1mm}\forall  a\in A,
 *  \end{array}
 * \f]
 * where \f$N=V\setminus T\f$, \f$\delta^+(X):=\{(u,v)\in A| u\in X, v\in V\setminus X\}\f$, \f$\delta^-(X):= \delta^+(V \setminus X)\f$ for
 * \f$X\subset V\f$ i.e., \f$\delta^+(X)\f$ is the set of all arcs going out of and \f$\delta^-(X)\f$ the set of all arcs going into \f$X\f$.
 *
 * Since the model potentially contains an exponential number of constraints, a separation approach is employed.
 * Violated constraints are separated during the execution of the branch-and-cut algorithm.
 *
 * In addition to Steiner problems in graphs there exist several variations. The following Steiner problem variants can be solved by SCIP-JACK,
 * by transforming them to a Steiner arborescence problem, and in some cases introducing additional constraints:
 *
 * -Steiner arborescence problems,
 *
 * -rectilinear Steiner minimum tree problems,
 *
 * -node-weighted Steiner tree problems,
 *
 * -prize-collecting Steiner tree problems,
 *
 * -rooted prize-collecting Steiner tree problems,
 *
 * -maximum-weight connected subgraph problems,
 *
 * -degree-constrained Steiner tree problems,
 *
 * -group Steiner tree problems, and
 *
 * -hop-honstrained directed Steiner tree problems.
 *
 * A more detailed description of SCIP-Jack and its various components can be found in
 * "SCIP-Jack – A solver for STP and variants with parallelization extensions" by Gerald Gamrath et. al.
 *
 *
 */
