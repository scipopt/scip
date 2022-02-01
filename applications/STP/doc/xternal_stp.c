/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
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
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 * @author Yuji Shinano
 * @author Michael Winkler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page STP_MAIN STP
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Daniel Rehfeldt
 * @author Yuji Shinano
 * @author Michael Winkler
 *
 * This application contains the branch-and-cut based solver SCIP-Jack for Steiner tree and related problems, realized within the framework
 * \SCIP. See the doctoral thesis: "Faster algorithms for Steiner tree and related problems: From theory to practice"
 * by Daniel Rehfeldt for more details. See also https://scipjack.zib.de. This application includes (among others):
 *
 * - problem readers, which parse the problem out of .stp or .gr files
 *   (reader_stp.c, reader_gr.c)
 * - a (global) problem data structure, containing all necessary information including the graph, which creates the model within \SCIP (probdata_stp.c)
 * - various algorithms and data structures for handling Steiner tree graphs (graph_base.c, graph_edge.c, graph_history.c, graph_mcut.c, graph_path.c, graph_save.c, graph_stats.c, graph_tpath.c, graph_util.c
 *   graph_delpseudo.c, graph_grid.c, graph_load.c, graph_node.c, graph_pcbase.c, graph_sdpath.c, graph_sub.c, graph_trans.c, graph_vnoi.c)
 * - decomposition methods (bidecomposition.c, cons_stpcomponents.c)
 * - algorithms to solve fixed-parameter tractable (sub-)problems (dpborder_base.c, dpborder_core.c, dpborder_util.c, dpterms_base.c, dpterms_core.c,
 *   dpterms_util.c, enumeration.c)
 * - various primal heuristics (heur_tm.c, heur_prune.c, heur_ascendprune.c, heur_slackprune.c, heur_local.c, heur_rec.c)
 * - a dual heuristic (dualascent.c)
 * - various reduction techniques (reduce_alt.c, reduce_bnd.c, reduce_ext.c,  reduce_pcsimple.c, reduce_sdcomp.c, reduce_sdutil.c, reduce_simple.c,
 *   reduce_termsepa.c,reduce_termsepafull.c, reduce_base.c, reduce_da.c, reduce_path.c, reduce_sd.c, reduce_sdgraph.c, reduce_sepa.c,
 *   reduce_sol.c, reduce_termsepada.c, reduce_util.c)
 * - extended reduction techniques (extreduce_base.c, extreduce_contract.c, extreduce_data.c, extreduce_dist.c, extreduce_extmst.c, extreduce_mldists.c,
 *   extreduce_util.c, extreduce_bottleneck.c, extreduce_core.c, extreduce_dbg.c, extreduce_extmstbiased.c, extreduce_extspg.c, extreduce_redcosts.c)
 * - a constraint handler, which checks solutions for feasibility and separates any violated model constraints (cons_stp.c)
 * - a propagator, which attempts to fix (edge) variables, or detects infeasiblity (prop_stp.c)
 * - relaxation handlers, for example to detect and efficiently solve fixed-parameter tractable instances (relax_stp.c, relax_stpdp.c, relax_stpenum.c)
 * - a branching rule (branch_stp.c)
 * - generic data-structures and algorithms (stpvector.h, stpbitset.h, stpprioqueue.c, misc_stp.c, mst.c)
 *
 *
 * In the following, the problem is introduced and the solving process is delineated. Next, two main plugins are
 * sketched. Finally, information on reading a problem and on writing a solution found by SCIP-Jack is given.
 *
 * -# \subpage STP_PROBLEM "Problem description and solving approach"
 * -# \subpage STP_PROBLEMDATA "Main problem data, creating the problem"
 * -# \subpage STP_CONS "Separating violated constraints"
 * -# \subpage STP_READWRITE "Reading and writing"
 * -# \subpage STP_MISC "Miscellaneous methods used for Steiner tree problems"
 * -# \subpage STP_MINCUT "Graph minimum cut routine"
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file", application "stp" (for Make) or "scipstp" (for CMake)
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
 * The solving approach of SCIP-Jack is based on three major components:
 *
 * First, (problem-specific) reduction techniqes are extremely important. Apart from some pathological instances
 * specifically constructed to defy presolving techniques, preprocessing is often able to significantly
 * reduce instances or even solve them. Additionally, reduction techniques are used for domain propagation.
 *
 * Second, primal heuristics are needed, especially for hard instances, to find good or even optimal solutions.
 * In particular, primal heuristics are important for bound-based reduction methods.
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
 * -rooted maximum-weight connected subgraph problems,
 *
 * -maximum-weight connected subgraph problems with budget constraints,
 *
 * -partial-terminal node-weighted Steiner tree problems,
 *
 * -degree-constrained Steiner tree problems,
 *
 * -group Steiner tree problems, and
 *
 * -hop-honstrained directed Steiner tree problems.
 *
 * A more detailed description of SCIP-Jack and its various components can be found in
 * "Faster algorithms for Steiner tree and related problems: From theory to practice" by Daniel Rehfeldt.
 *
 *
 */

/**@page STP_READWRITE Reading and writing
 *
 * Problem instances can be read in both in the .stp format described here: http://steinlib.zib.de/format.php,
 * and the (somewhat simpler) .gr format described here: https://pacechallenge.org/2018/steiner-tree/. To run SCIP-Jack from the command line use
 *
 * bin/stp -f filename.stp
 *
 * Note that the solution shown by SCIP does not correspond to the actual solution computed by SCIP-Jack,
 * because SCIP-Jack performs preprocessing prior to building the SCIP problem (for performance and memory reasons).
 * The solution value shown by SCIP is correct, however.
 * To get the solution computed by SCIP-Jack, create a file called (for example) "write.set" in folder ./settings with content
 *
 * stp/logfile = "use_probname"
 *
 * Then call
 *
 * bin/stp -f filename.stp -s settings/write.set
 *
 * The solution will be written to the file "filename.stplog". Alternatively, you can set the logfile parameter via the interactive shell.
 *
 */
