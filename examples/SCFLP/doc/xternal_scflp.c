/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal_scflp.c
 * @brief  main document page
 * @author Stephen J. Maher
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page SCFLP_MAIN Stochastic capacitated facility location problem example
 * @version  0.9
 * @author   Stephen J. Maher

 * This is an example of using the Benders' decomposition framework of SCIP to solve the stochastic capacitated facility
 * location problem (abbreviated to SCFLP). The instances used for this problem are taken from the OR-Library CAP
 * instances.  These instances describe the deterministic capacitated facility location problem. The customer demands of
 * the deterministic problem are used as the mean of the normal distribution in the stochastic program.
 *
 * To use the Benders' decomposition framework to solve the SCFLP instances requires the implementation of two plugins:
 *
 * - a \ref reader_scflp.c "problem reader" which parses the data from the CAP instance files and provides it to the
 *   probdata plugin in a convenient format to build the problem within \SCIP.
 * - a \ref probdata_scflp.c "problem data structure" which builds the problem and stores the global information. The
 *   storage of global information is not absolutely necessary in this example, but it can be useful in post processing
 *   of the solutions and checking their correctness.
 *
 * The SCFLP example formulates the problem as the determinstic equivalent, which can be solved directly by SCIP and by
 * Benders' decomposition. Initially, we will describe how to build the deterministic equivalent problem. Second, we
 * will describe how to build the problem so that the Benders' decomposition framework can be used.
 *
 * -# @subpage SCFLP_PROBLEM "Problem description"
 * -# @subpage SCFLP_READER "Parsing the input format"
 * -# @subpage SCFLP_SOLVEPROB "Solving the deterministic equivalent using SCIP
 *    - @subpage SCFLP_DETEQUIV "Directly as a monolithic MIP"
 *    - @subpage SCFLP_BENDERS "Applying Benders' decomposition"
 */

/**@page SCFLP_PROBLEM Problem description
 *
 * In the following we describe the CIP model that we use: both the monolithic mixed integer program and the decomposed
 * problem (using Benders' decomposition).
 *
 * Given: a set of facilities \f$ I \f$ and a set of customers \f$ J \f$. The set of scenarios is given by \f$ S \f$,
 * which are defined by a set of different customer demands.
 *
 * Task: Find the minimum cost facilities to open such that the customer demand can be satisfied in all scenarios.
 *
 * Variables:
 *  - \f$ x_i \in \{0,1\} \quad \forall i \in I \f$:
 *    - \f$ x_i = 1 \f$ if facility \f$ i \f$ is opened.
 *  - \f$ y^{s}_{ij} \ge 0 \quad \forall i \in I, \forall j \in J, \forall s \in S \f$
 *    - \f$ y^{s}_{ij} \f$ is the level of demand for customer \f$ j \f$ satisfied by facility \f$ i \f$ in scenario
 *    \f$ s \f$.
 *
 * Parameters:
 *  - \f$ f_i \f$ the fixed cost for opening facility \f$ i \f$,
 *  - \f$ q_{ij} \f$ the cost of servicing customer \f$ j \f$ from facility \f$ i \f$,
 *  - \f$ \lambda^{s}_{j} \f$ the demand of customer \f$ j \f$ in scenario \f$ s \f$,
 *  - \f$ k_i \f$ the capacity of facility \f$ i \f$.
 *
 * @section SCFLP_DETEQUIVMODEL The deterministic equivalent
 *
 * The deterministic equivalent can be formulated as:
 *
 * \f[
 *  \begin{array}[t]{rll}
 *    \min & \displaystyle \sum_{i \in I} f_{i} x_{i} + \frac{1}{|S|}\sum_{s \in S}\sum_{i \in I}\sum_{j \in J}q_{ij}y^{s}_{ij} \\
 *         & \\
 *    subject \ to & \displaystyle \sum_{i \in I} y^{s}_{ij} \ge \lambda^{s}_{j} & \quad \forall j \in J, \forall s \in  S \\
 *         & \\
 *         & \displaystyle \sum_{j \in J} y^{s}_{ij} \le k_{i}x_{i} & \quad \forall i \in I, \forall s \in  S \\
 *         & \\
 *         & \displaystyle \sum_{i \in I} k_{i}x_{i} \le \max_{s \in S}\sum_{j \in J}\lambda^{s}_{j} & \\
 *         & \\
 *         & \displaystyle x_{i} \in \{0, 1\} & \quad \forall i \in I \\
 *         & \\
 *         & \displaystyle y^{s}_{ij} \ge 0 & \quad \forall i \in I, \forall j \in J, \forall s \in S \\
 *  \end{array}
 * \f]
 *
 * It can be seen that as the number of scenarios increases, the size of the deterministic equivalent increases
 * significantly. For large numbers of scenarios, the resulting deterministic equivalent can in intractable. This
 * limitation can be addressed by applying decomposition techniques. In this example, Benders' decomposition is applied
 * to solve the stochastic capacitated facility location problem.
 *
 * @section SCFLP_BENDERSMODEL Applying Benders' decomposition to the SCFLP
 *
 * The application of Benders' decomposition forms a master problem, consisting of only the facility location variables,
 * and a subproblem for each scenario, consisting of the customer servicing variables for the given secnario.
 *
 * The master problem is given by:
 *
 * \f[
 *  \begin{array}[t]{rll}
 *    \min & \displaystyle \sum_{i \in I} f_{i} x_{i} + \frac{1}{|S|}\sum_{s \in S}\varphi^{s} \\
 *         & \\
 *    subject \ to & \displaystyle \sum_{i \in I} k_{i}x_{i} \le \max_{s \in S}\sum_{j \in J}\lambda^{s}_{j} & \\
 *         & \\
 *         & \displaystyle \varphi^{s} \geq \sum_{j \in J}\lambda^{s}_{j}u^{p}_{j} + \sum_{i \in I}k_{i}x_{i}v^{p}_{i} & \quad \forall s \in S, \forall p \in P^{s} \\
 *         & \\
 *         & \displaystyle 0 \geq \sum_{j \in J}\lambda^{s}_{j}u^{r}_{j} + \sum_{i \in I}k_{i}x_{i}v^{r}_{i} & \quad \forall s \in S, \forall r \in R^{s} \\
 *         & \\
 *         & \displaystyle x_{i} \in \{0, 1\} & \quad \forall i \in I \\
 *         & \\
 *         & \displaystyle \varphi^{s} \geq 0 & \quad \forall s \in S \\
 *  \end{array}
 * \f]
 *
 * where \f$ \varphi^{s} \f$ is the auxiliary variable for each scenario \f$ s \f$ that is an underestimator of the
 * optimal subproblem objective function value. The second and third constraint of the master problem are the Benders'
* optimality and feasibility cuts. Given a solution to the master problem, an optimality cut for scenario \f$s\f$ is
* generated from the optimal dual solution to the corresponding subproblem. Similarly, if the solution to the master
* problem induces an infeasibl instance of subproblem \f$s\f$, then the resulting dual ray is used to generate a
* feasibility cut.
 *
 * The subproblem for scenario \f$ s \f$ that are solved to generate optimality and feasibility cuts are:
 *
 * \f[
 *  \begin{array}[t]{rll}
 *    z^{s}(\bar{x}) = \min & \displaystyle \sum_{i \in I}\sum_{j \in J}q_{ij}y^{s}_{ij} \\
 *         & \\
 *    subject \ to & \displaystyle \sum_{i \in I} y^{s}_{ij} \ge \lambda^{s}_{j} & \quad \forall j \in J \\
 *         & \\
 *         & \displaystyle \sum_{j \in J} y^{s}_{ij} \le k_{i}\bar{x}_{i} & \quad \forall i \in I \\
 *         & \\
 *         & \displaystyle y^{s}_{ij} \ge 0 & \quad \forall i \in I, \forall j \in J \\
 *  \end{array}
 * \f]
 *
 * The solution \f$\bar{x}\f$ is the candidate solution that is verified by solving the subproblem. As explained above,
 * if the subproblem is infeasible, then the corresponding dual ray is used to generate a Benders' feasibility cut.  If
 * the subproblem is optimal and \f$ z^{s}(\bar{x}) > \varphi^{s} \f$, then an optimality cut is generated from the
 * corresponding dual solution. If \f$ z^{s}(\bar{x}) \le \varphi^{s} \f$, then the subproblem is optimal for the given
 * solution \f$ \bar{x} \f$. If \f$ z^{s}(\bar{x}) > \varphi^{s} \f$ for all scenario subproblems, then \f$ \bar{x} \f$
 * is the optimal solution to the original problem.
 */
