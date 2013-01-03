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
 * @author Timo Berthold
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This example contains a branch-and-price approach for the binpacking problem which is realized with the framework 
 * \SCIP. Therefore, the following plugins are implemented:
 *
 * - a problem reader which parses the problem out of file and creates the corresponding problem within \SCIP
 *   (reader_bpa.c)
 * - a (global) problem data structure which contains all necessary information (probdata_binpacking.c)
 * - a pricer which generates new variables/columns during the search (pricer_binpacking.c)
 * - the Ryan/Foster branching rule (branch_ryanfoster.c)
 * - a constraint handler which handles the branching decisions of the Ryan/Foster branching (cons_samediff.c)
 * - a variable data structure which stores information for each variable and is needed to perform the Ryan/Foster
 *   branching (vardata_binpacking.c)
 *
 * In the following we introduce the problem, explain the use of the reader plugin and pricer plugin. Finally, we
 * introduce the Ryan/Foster branching rule and briefly discuss how that specific branching rule can be realized within
 * the framework \SCIP.
 *
 * -# \ref PROBLEM "Problem description"
 * -# \ref READER "Parsing the input format and creating the problem"
 * -# \ref PROBLEMDATA "Main problem data"
 * -# \ref PRICER "Pricing new variables"
 * -# \ref BRANCHING "Ryan/Foster branching"
 * -# \ref MAKEFILE "The Makefile"
 *
 */

/**@page PROBLEM Problem description
 *
 * The binpacking problem consists of the task to distribute a given set of items \f$ [n] := \{1, \dots, n\}\f$ with
 * nonnegative size \f$s_i\f$ to a minimal number of bins, all of the same capacity \f$\kappa\f$.  As example consider 9
 * items with sizes: 2, 1, 2, 1, 1, 2, 3, 2, and 1 and a bin capacity of \f$\kappa\f$ of 4. The following pictures show
 * a feasible solution which needs 5 bins. The minimum number of bins needed for that example is 3.
 *
 * <CENTER>
 * \image html binpacking.png
 * </CENTER>
 *
 * This problem can be formulated as a set covering problem. Therefore, we introduce a binary variable \f$x_{S}\f$ for
 * each feasible packing \f$S\f$. A <b>packing</b> \f$S\f$ is an assignment vector \f$ \lambda_{S}\in\{0,1\}^n \f$ which
 * states the items belonging to that packing. It is <b>feasible</b>, if and only if the total size of the items
  contained in this assignment is not greater than the given capacity \f$\kappa\f$. Let \f$\mathcal{S}\f$ be the set of
 * all feasible packing, this measns:
 *
 * \f[
 *    \mathcal{S} := \{S\subseteq [n] \mid \sum_{i:i\in S} s_{i} \leq \kappa \}
 * \f]
 *
 * An integer program can be formulated as follows:
 * 
 * \f[
 *  \begin{array}[t]{rll}
 *    \min & \displaystyle \sum_{S \in \mathcal{S}} x_{S} \\
 *         & \\
 *    subject \ to & \displaystyle \sum_{S \in \mathcal{S}} (\lambda_{S})_{i}x_{S} \ge 1 & \quad \forall i \in \{1,\dots,n\} \\
 *         & \\
 *         & x_{S} \in \{0,1\} & \quad \forall S \in \mathcal{S} \\
 *  \end{array}
 * \f]
 *
 * This means we are searching for a set of packings such that each item is contained in at least one of the selected
 * packings. Since the objective is to minimize the number of used packings, each optimal solution to the above problem
 * can be transformed into a solution where each item is packed exactly once with the same cost.
 *
 *
 * Since \f$\mathcal{S}\f$ can be of exponential size, we are using a column generation approach to solve this
 * problem. We initialize the (master) problem with a set of \f$ n \f$ variables representing packings of a single item
 * per bin.  Now, we have to iteratively search for variables representing "better" packings, i.e., a packing pattern
 * which reduces the overall cost. For a given solution \f$y^\star\f$ of the (restricted) dual linear program, we have
 * to ﬁnd a variable/packing \f$ \lambda_{S} \f$ for which the reduced costs is negative. This means:
 *
 * \f[
 *     c_{S} - \sum_{i=0}^n (\lambda_S)_i y_i^\star < 0.
 * \f]
 *
 * Since all variables \f$ \lambda_{S} \f$ have an objective coefficient \f$ c_{S} = 1 \f$ the above condition is
 * equivalent to
 *
 * \f[
 *     \sum_{i=0}^n  (\lambda_S)_i y_i^\star > 1.
 * \f]
 *
 *
 * To find such a variable/packing we solve the following integer program:
 *
 *  \f[
 *  \begin{array}[t]{rll}
 *       \max & \displaystyle \sum_{i=1}^n (\lambda_S)_i y^\star_i\\
 *        & \\
 *        subject \ to & \displaystyle \sum_{i=0}^n (\lambda_S)_i s_i \leq \kappa \\
 *        & \\
 *        & (\lambda_S)_i \in \{0,1\} & \quad \forall i \in \{ 1, \dots , n \} \\
 *  \end{array}
 * \f]
 *
 * where \f$ (\lambda_S)_i \f$ for \f$i\in\{1,\dots,n\}\f$ are binary variables and \f$y^\star_i\f$ given by the dual
 * solution of the restricted master problem.
 *
 * The above problem is a knapsack problem which can be solved via dynamic programming or by solving the above integer
 * program. In this example we implemented a pricer which solves the integer program.
 *
 */



/**@page MAKEFILE The Makefile
 *
 * The Makefile is based on the main \SCIP Makefile. This means that all compiling options which are
 * available for \SCIP are also available for the binpacking project. Below, you find a list
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
 *   <code>doc/html/index.html</code>.
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
 *   href="http://scip.zib.de/doc/html/TEST.html">How to run automated tests with SCIP</a>).  
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
 
