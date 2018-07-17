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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 * @version  0.1
 * @author   Marc Pfetsch
 *
 * This code uses Benders decomposition to solve the Minimum IIS-Cover Problem (MinIISC). Here, one is given an
 * infeasible linear system and wants to compute a subsystem of smallest size whose removal leaves a feasible
 * subsystem. This corresponds to removing at least one constraint from each irreducible infeasible subsystem (IIS).
 *
 * The approach is described in:
 *
 * "Finding the Minimum Weight IIS Cover of an Infeasible System of Linear Inequalities"@n
 * by Mark Parker and Jennifer Ryan,@n
 * Ann. Math. Artif. Intell. 17, no. 1-2, 1996, pp. 107-126.
 *
 * @section Appr Solution Approach
 *
 * This approach works as follows. MinIISC can be formulated as:
 * \f{eqnarray*}{
 *    \min && \sum_{i=1}^m y_i \\
 *    s.t. && \sum_{i \in I} y_i \geq 1 \quad \forall \mbox{ IISs }I \\
 *         && y_i \in \{0,1\} \quad \forall i.
 * \f}
 *
 * We begin with a subset of IISs and solve the above set covering problem (using SCIP as a MIP solver). We then check
 * whether the resulting vector \f$y^*\f$ corresponds to an IIS-cover. If this is the case, we are done. Otherwise, we
 * check for an uncovered IIS and add its inequality to the set covering problem. We then repeat the process.
 *
 * Checking for an uncovered IIS can be done using a so-called alternative polyhedron. We explain the approach for the case
 * in which the linear system is \f$Dx \leq d\f$ where \f$D\f$ is an \f$m \times n\f$ matrix. The alternative polyhedron
 * is then
 * \f[
 * \{z : D^T z = 0,\; d^T z \leq -1,\; z \geq 0\}.
 * \f]
 * Gleeson and Ryan [1990] proved that the vertices of this polyhedron are in 1-to-1 correspondence with the IISs of the
 * original system. If we are then given the solution \f$y^* \in \{0,1\}^m\f$ from the set covering problem, we can
 * consider
 * \f[
 * \{z : D^T z = 0,\; d^T z \leq -1,\; z_i = 0 \mbox{ for all i with }y^*_i = 1,\; z \geq 0\}.
 * \f]
 * A vertex of this polyhedron corresponds to an IIS is the system remaining from \f$D x \leq d\f$ when the inequalities
 * given by \f$y^* = 1\f$ are deleted.
 *
 * @section Impl Implementation
 *
 * The implementation uses several tricks to speed up the solution process:
 * - Several IISs are generated in one round using the technique described in
 *      Branch-And-Cut for the Maximum Feasible Subsystem Problem,@n
 *      Marc Pfetsch, SIAM Journal on Optimization 19, No.1, 21-38 (2008)
 * - The master problem can be solved approximately (using a gap limit) or using a stall limit (the final MIP has to be
 *   solved exactly).
 * - Moreover, the master problem can be tackled using reoptimization.
 *
 * The input to the code should be an infeasible linear program (the objective is ignored) in any format that SCIP can
 * handle. The basic benders algorithm is implemented in the file benders.c using a call back for the cut generation.
 *
 * Compiling the MinIISC application
 * ----------------------------------
 *
 * See the @ref INSTALL "Install file"
 */
