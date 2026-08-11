/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   xternal_lj.c
 * @brief  main documentation page of the Lennard-Jones Cluster example
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page LJ_MAIN Nonlinear Handler for Lennard-Jones Cluster problem
 * @author   Stefan Vigerske
 *
 * To find globally optimal solutions for (mixed-integer) nonlinear optimization problems, SCIP relies on a linear relaxations to estimate the gap between the objective function value of a best known feasible solution and the optimal value of the problem.
 * To construct such a relaxation, linear under- or overestimators of nonlinear functions are often used.
 * In the following, we illustrate how to extend the techniques available in SCIP by using the \ref NLHDLR "nonlinear handler plugin".
 * This plugin type allows close interaction with the code that constructs the linear relaxation and performs variable bound tightening for nonlinear constraints.
 *
 * # Example Application
 *
 * As an exemplary application from computational chemistry, consider the solution of the [Lennard-Jones Cluster problem](https://en.wikipedia.org/wiki/Lennard-Jones_potential).
 * Given a fixed number of particles, the task is to position them in 3-dimensional space such that the sum of the Lennard-Jones potentials between pairs of particles is minimized.
 * The Lennard-Jones potential combines a repulsive and an attractive potential that depend on the Euclidean distance of two particles.
 * A widely used formula is \f[4\epsilon\left[\left(\frac{\sigma}{d}\right)^{12} - \left(\frac{\sigma}{d}\right)^{6}\right],\f]
 * where \f$d\f$ is the distance between the particles, and \f$\epsilon\f$, \f$\sigma>0\f$ are parameters.
 * The minimization of the Lennard-Jones potential for \f$\epsilon=\sigma=1\f$ is a popular [test problem for global optimization](https://www-wales.ch.cam.ac.uk/~wales/CCD/jon/structures/LJ.html).
 *
 * The problem for \f$N\f$ particles can be formulated in SCIP as follows:
 * \f{equation}{
 * \begin{aligned}
 * \min\; & 4\sum_{i=1}^N \sum_{j=i+1}^N p_{ij} \\
 * \text{s.t.}\; & p_{ij} \geq (\Vert x^i - x^j\Vert_2^2)^{-6} - (\Vert x^i - x^j\Vert_2^2)^{-3} && \forall i,j\in\{1,\ldots,N\}, i < j, \\
 * & x^i \in [-9,9]^3 && i\in\{1,\ldots,N\}.
 * \end{aligned}
 * \f}
 * Bounding the particles' positions by the box \f$[-9,9]^3\f$ is an artificial constraint that is required to ensure convergence of SCIP's algorithm.
 *
 * # Nonlinear Handler Implementation
 *
 * To solve this problem, SCIP needs to compute linear underestimators for the function \f$f(r) := r^{-6} - r^{-3}\f$, where \f$r\f$ is a variable that SCIP introduces to stand for the distance of two particles (\f$\Vert x^i - x^j\Vert_2^2\f$).
 * In the standard approach, \f$r^{-6}\f$ is recognized as convex and \f$-r^{-3}\f$ as concave.
 * A convex underestimator for \f$f(r)\f$ is then derived by replacing \f$-r^{-3}\f$ by the secant on the function between current lower and upper bounds on \f$r\f$.
 * While this gives a valid underestimator, it is far from being tight, which results in a weak dual bound.
 *
 * A tighter relaxation could be obtained if \f$r^{-6}\f$ and \f$-r^{-3}\f$ were not underestimated separately, but the convex envelope of \f$f(r)\f$ were used instead.
 * This is achieved by means of the nonlinear handler in this example.
 * It implements the following central callbacks:
 *
 * ## Detect Structure
 *
 * Given an algebraic expression, check whether it is of the form \f$r^{-6} - r^{-3} - p\f$.
 * Here, \f$r\f$ and \f$p\f$ could be any expressions, but will correspond to a quadratic expression \f$\Vert x^i - x^j\Vert_2^2\f$ and a variable expression \f$p_{ij}\f$, respectively, for some particles \f$i,j\f$ in this example.
 * If an expression of this structure is found, the callback informs the handler for nonlinear constraints that it will provide underestimators and participate in bound tightening.
 * To produce underestimators that can be added to the LP relaxation, the nonlinear handler requires variables that stand for \f$r\f$ and \f$p\f$.
 * SCIP ensures that these are made available.
 *
 * In the actual implementation, the handler checks for \f$a r^{-6} - a r^{-3} + b p\f$ for some coefficients \f$a\in\{-1,1\}\f$ and \f$b\neq 0\f$, because SCIP's handler for nonlinear constraints may multiply \f$p_{ij} \geq r_{ij}^{-6} - r_{ij}^{-3}\f$ by \f$-1\f$.
 * If \f$a<0\f$, the handler will provide overestimators instead of underestimators.
 * For the remainder of this documentation, we will assume \f$a=1\f$ and \f$b=-1\f$, though.
 *
 * ## Linear Estimators
 *
 * At a node of the branch-and-bound tree, let \f$(r',p')\f$ be the value of the expressions \f$r\f$ and \f$p\f$ at the current node and \f$\ell\f$ and \f$u\f$ be the bounds on the expression \f$r\f$ (SCIP computes these via interval arithmetics).
 * The nonlinear handler needs to compute a linear underestimator for \f$f(r)-p\f$ that is as tight as possible at \f$(r',p')\f$.
 * Note that \f$f(r)\f$ attains its minimum at \f$r_{\min}=\sqrt[3]{2}\approx 1.2599\f$ and has an inflection point at \f$r_{\inf}=\sqrt[3]{\frac{7}{2}}\approx 1.51829\f$.
 * The function is convex for \f$0<r\leq r_{\inf}\f$ and concave for \f$r\geq r_{\inf}\f$.
 *
 * - If \f$r'\leq r_{\min}\f$ or \f$u \leq r_{\inf}\f$, a supporting hyperplane at \f$r'\f$ is used to underestimate \f$f(r)\f$: \f$f(r') + f'(r') (r-r')\f$.
 * - If \f$\ell \geq r_{\inf}\f$, then \f$f(r)\f$ is concave on \f$[\ell,u]\f$ and a linear underestimator is given by the secant between \f$\ell\f$ and \f$u\f$: \f$f(\ell) + \frac{f(u)-f(\ell)}{u-l}(r-\ell)\f$.
 * - In the remaining case, \f$\ell < r_{\inf} < u\f$ with \f$r' > r_{\min}\f$, \f$f(r)\f$ changes curvature within \f$[\ell,u]\f$.
 *   Similar to the technique for [convex envelopes of monomials of odd degree](https://doi.org/10.1023/A:1021924706467), a valid linear underestimator of \f$f(r)\f$ is given by the secant between a point in the convex region, \f$\tilde r \in [r_{\min},r_{\inf}]\f$, and \f$u\f$.
 *   \f$\tilde r\f$ needs to be chosen such that the slope of the secant equals the slope of the tangent at \f$\tilde r\f$, i.e.,
 *   \f[ f'(\tilde r) = \frac{f(u) - f(\tilde r)}{u - \tilde r}. \f]
 *   This equation is solved by using Newton's Method to find a root of \f$f(u) - f(\tilde r) - f'(\tilde r) ( u - \tilde r)\f$, using \f$\frac{1}{2}(r_{\min}+r_{\inf})\f$ as starting pointing.
 *   If \f$\tilde r < \ell\f$, then the secant between \f$\ell\f$ and \f$u\f$ can be used as linear underestimator.
 *
 * ## Interval Evaluation
 *
 * Given bounds on \f$r\f$ and \f$p\f$, compute the range of \f$f(r)-p\f$.
 * The minimum of \f$f(r)\f$ is attained either at \f$\sqrt[3]{2}\f$, the lower, or the upper bound on \f$r\f$.
 * The maximum of \f$f(r)\f$ is attained at the lower or upper bound on \f$r\f$, if \f$r>0\f$.
 *
 * ## Domain Propagation
 *
 * Given an interval \f$[\underline{g},\overline{g}]\f$ for \f$f(r)-p\f$ and bounds \f$[\underline{p},\overline{p}]\f$ for \f$p\f$, compute bounds on \f$r\f$.
 * The code first considers the equation \f$z^2 - z\in [\underline{g},\overline{g}] + [\underline{p},\overline{p}]\f$ and finds an interval \f$[\underline{z},\overline{z}]\f$ such that \f$[\underline{z},\overline{z}]^2-[\underline{z},\overline{z}]\supseteq [\underline{g},\overline{g}] + [\underline{p},\overline{p}]\f$.
 * Next, bounds on \f$r\f$ are given by \f$[\underline{z},\overline{z}]^{-\frac{1}{3}}\f$.
 *
 * # Installation
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 */
