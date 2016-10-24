/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_nlobbt.h
 * @ingroup PROPAGATORS
 * @brief  nonlinear OBBT propagator
 * @author Benjamin Mueller
 *
 * In Nonlinear Optimization-Based Bound Tightening (NLOBBT), we solve auxiliary NLPs of the form
 * \f[
 *      \min / \max \, x_i \\
 * \f]
 * \f[
 *      s.t. \; g_j(x) \le 0 \, \forall j=1,\ldots,m \\
 * \f]
 * \f[
 *      c'x \le \mathcal{U}
 * \f]
 * \f[
 *      x \in [\ell,u]
 * \f]
 *
 * where each \f$ g_j \f$ is a convex function and \f$ \mathcal{U} \f$ the solution value of the current
 * incumbent. Clearly, the optimal objective value of this nonlinear program provides a valid lower/upper bound on
 * variable \f$ x_i \f$.
 *
 * The propagator sorts all variables w.r.t. their occurrences in convex nonlinear constraints and solves sequentially
 * all convex NLPs. Variables which could be successfully tightened by the propagator will be prioritized in the next
 * call of a new node in the branch-and-bound tree.
 *
 * By default, NLOBBT is only applied for non-binary variables. Variables which do not appear non-linearly in the
 * nonlinear constraints will not be considered even though they might lead to additional tightenings.
 *
 * After solving an NLP we try to exploit the dual information to generate a globally valid inequality, called
 * Generalized Variable Bound (see @ref prop_genvbounds.h). Let \f$ \lambda_j \f$, \f$ \mu \f$, \f$ \alpha \f$, and
 * \f$ \beta \f$ be the dual multipliers for the constraints of the NLP where \f$ \alpha \f$ and \f$ \beta \f$
 * correspond to the variable bound constraints. Because of the convexity of \f$ g_j \f$ we know that
 *
 * \f[
 *      g_j(x) \ge g_j(x^*) + \nabla g_j(x^*)(x-x^*)
 * \f]
 *
 * holds for every \f$ x^* \in [\ell,u] \f$. In an optimal solution \f$ x^* \f$ of the NLP we know that the KKT
 * conditions
 *
 * \f[
 *      +/- e_i + \nabla g(x^*) + c + \alpha - \beta = 0
 * \f]
 * \f[
 *      \lambda_j g_j(x) = 0
 * \f]
 *
 * hold. Applying these identities to a redundant but valid inequality
 *
 * \f[
 *      x_i \ge x_i + \sum_{j} g_j(x^*_i)
 * \f]
 *
 * leads to a globally valid inequality
 *
 * \f[
 *      x_i \ge (\beta - \alpha)'x + (e_i + \alpha - \beta) x^* + \mu \mathcal{U} .
 * \f]
 *
 * which is passed to the genvbounds propagator.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_NLOBBT_H__
#define __SCIP_PROP_NLOBBT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the nlobbt propagator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludePropNlobbt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
