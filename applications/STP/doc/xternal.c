/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
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
 * @author Daniel Rehfeldt
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Stephen Maher
 * @author Yuji Shinano
 * @author Michael Winkler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 * @author Daniel Rehfeldt
 *
 * This application contains a (by default) branch-and-cut based solver for Steiner problems, realized within the framework
 * \SCIP, see: "SCIP-Jack - A solver for STP and variants with parallelization extensions" by G. Gamrath et al. The following plugins are implemented:
 *
 * - a problem reader which parses the problem out of an .stp file
 *   (reader_stp.c)
 * - a (global) problem data structure, containing all necessary information including the graph, which creates the model within \SCIP (probdata_stp.c)
 * - a construction heuristic (heur_tm.c)
 * - an improvment heuristic (heur_local.c)
 * - a recombination heuristic (heur_rec.c)
 * - a, by default not used, pricer which generates new variables/columns during the search (pricer_stp.c)
 * - a constraint handler which checks solutions for feasibility and separates any violated model constraints (cons_stp.c)
 * - a propagator which attempts to fix (edge) variables to zero utilizing their reduced costs (prop_stp.c)
 * - an event handler which simple writes each incumbent solution to a file, if activated (event_bestsol.c)
 *
 * In the following the problem is introduced and the solving process delineated. Afterwards, two plugins are
 * sketched.
 *
 * -# \ref PROBLEM "Problem description and solving approach"
 * -# \ref PROBLEMDATA "Main problem data, creating the problem"
 * -# \ref CONS "Separating violated constraints"
 * -# \ref MAKEFILE "The Makefile"
 *
 */

/**@page PROBLEM Problem description and solving approach
 *
 * A more intricate account of the following can be found in
 * "SCIP-Jack - A solver for STP and variants with parallelization extensions" by G. Gamrath et al.
 * The Steiner tree problem in graphs (SPG) can be described as follows: Given an undirected connected graph
 * \f$ G=(V,E)\f$, costs \f[ c: E \rightarrow  \mathcal{Q}^+ \f] and a set \f$ T \subset V \f$ of \f$ \textit{terminals} \f$,
 * the problem is to find a minimum weight tree \f$ S\subseteq G \f$ which spans \f$ T \f$.
 * The following picture shows a SPG instance with the terminals given as squares:
 *
 * <CENTER>
 * \image html stp.png
 * </CENTER>
 *
 *
 * Our branch-and-cut based Steiner tree solver can be dissected into three major components.
 *
 * First, problem specific preprocessing is extremely important. Apart from some pathological instances
 * specifically constructed to defy presolving techniques, preprocessing is often able to significantly
 * reduce instances.
 *
 * Second, heuristics are needed to find good or even optimal solutions.
 *
 * Finally, at the core is the branch-and-cut procedure used to compute a lower bound and prove optimality:
 *
 * The problem can be formulated using the directed equivalent of the STP, the Steiner arborescence problem (SAP):
 * Given a directed graph \f$ D=(V,A) \f$, a root \f$ r \in V \f$, costs \f$ c: A \rightarrow \mathcal{Q}^+ \f$
 * and a set \f$ T \subset V \f$ of terminals, a directed tree \f$ S\subseteq D \f$ is required such that
 * for all \f$ t \in T \f$, \f$ (V_S,A_S) \f$ contains exactly one directed path from \f$ r \f$ to \f$ t \f$.  Each STP can be
 * transformed to an SAP replacing each edge by two anti-parallel arcs of the same cost and distinguishing an arbirtrary
 * terminal as the root. This results in a one-to-one correspondence between the respective solution sets
 * Introducing variables \f$y_a\f$ for \f$a\in A\f$ with the interpretation \f$y_a:=1\f$, if \f$a\f$ is in the
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
 *    \label{flowcons1}
 *     0, \mbox{if } v=r;\\
 *     1, \mbox{if } v\in T\setminus{r};\\
 *     1, \mbox{if } v\in N; \end{array}} \hspace{2.9mm}\forall  v \in V
 *   \\
 *   \label{flowcons2}
 *   y(\delta^-(v))&\leq& y(\delta^+(v)), \hspace{10.5mm}\forall  v\in N;\\
 *   \label{flowcons3}
 *   y(\delta^-(v))&\geq& y_a, \hspace{20.2mm}\forall  a\in\delta^+(v), v\in N;\\
 *  0\leq y_a&\leq& 1, \hspace{22mm}\forall  a\in A;\\
 *   y_a&\in& \{0,1\}, \hspace{15.1mm}\forall  a\in A,
 *  \end{array}
 * \f]
 * where \f$N=V\setminus T\f$, \f$\delta^+(X):=\{(u,v)\in A| u\in X, v\in V\setminus X\}\f$, \f$\delta^-(X):= \delta^+(V \setminus X)\f$ for
 * \f$X\subset V\f$ i.e., \f$\delta^+(X)\f$ is the set of all arcs going out of and \f$\delta^-(X)\f$ the set of all arcs going into \f$X\f$.
 *
 * Since the model potentially contains an exponential number of constraints, a separation approach is employed.
 * Violated constraints are separated during the execution of the branch-and-cut algorithm
 *
 * In addition to Steiner Problems in Graphs there exist several variants of which the following can be solved by SCIP-JACK,
 * transforming them to a Steiner Arborescence Problem, and in some cases introducing additional constraints:
 *
 * -Steiner Arborescence Problems,
 *
 * -Rectilinear Steiner Minimum Tree Problems,
 *
 * -Node-Weighted Steiner Tree Problems,
 *
 * -Prize-Collecting Steiner Tree Problems,
 *
 * -Rooted Prize-Collecting Steiner Tree Problems,
 *
 * -Maximum-Weight Connected Subgraph Problems,
 *
 * -Degree-Constrained Steiner Tree Problems,
 *
 * -Group Steiner Tree Problems, and
 *
 * -Hop-Constrained Directed Steiner Tree Problems.
 */



/**@page MAKEFILE The Makefile
 *
 * The Makefile is based on the main \SCIP Makefile. This means that all compiling options which are
 * available for \SCIP are also available for the stp project. Below, you find a list
 * of the most important compiling flags, the values they can take, and a short description. The
 * values in bold face are the default values.
 *
 * - <code>LPS={clp | cpx | none | <b>spx</b>}</code>
 *   <br>
 *   Defines the linear program solver to use:
 *   - <code>clp</code> use COIN-OR CLP
 *   - <code>cpx</code> use IBM CPLEX
 *   - <code>none</code> no LP solver
 *   - <code><b>spx</b></code> use SoPlex
 *
 * - <code>OPT={dbg | <b>opt</b>}</code>
 *   <br>
 *   Defines if the projects gets compiled in debug (<code>dbg</code>) mode or
 *   optimized (<code><b>opt</b></code>) mode. In the debug mode all assertions are checked.
 *
 * - <code>ZIMPL={false | <b>true</b>}</code>
 *   <br>
 *   Defines if the modeling language ZIMPL should be linked to binary or not.
 *
 * In the following we explain the all <b>Makefile targets</b>.
 *
 * - <b>lint</b>
 *   <br>
 *   Statically checks the code for uninitialized variables and many other possible problems,
 *   which even do not lead to compiler errors. For this, the
 *   the external tool flexelint is needed. The call produces the file <code>lint.out</code>
 *   which contains all the detected warnings. From the experience when developing \SCIP, we strongly
 *   recommend to use such a code checker. It is always a surprising the stuff such tools detect.
 *   <br>
 *
 * - <b>doc</b>
 *   <br>
 *   Generates a html documentation. For this call the external tool
 *   <a href="http://doxygen.org">doyxgen</a> is needed.
 *   After generating the documentation, you can use your favorite browser to open the main page
 *   <code>doc/html/index.shtml</code>.
 *   <br>
 *
 * - <b>clean</b>
 *   <br>
 *   Remove all objective files, libraries, and binaries.
 *   <br>
 *
 * - <b>test</b>
 *   <br>
 *   Starts an automated test run based on the SCIP test runs (see <a
 *   href="http://scip.zib.de/doc/html/TEST.php">How to run automated tests with SCIP</a>).
 *   <br>
 *
 * - <b>tags</b>
 *   <br>
 *   Generates tags which can be used in the editor <b>emacs</b> and <b>xemacs</b>.
 *   <br>
 *
 * - <b>depend</b>
 *   <br>
 *   Generates the dependencies for the compiling process.
 */
