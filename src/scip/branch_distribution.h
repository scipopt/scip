/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_distribution.h
 * @ingroup BRANCHINGRULES
 * @brief  probability based branching rule based on an article by J. Pryor and J.W. Chinneck
 * @author Gregor Hendel
 *
 * The distribution branching rule selects a variable based on its impact on row activity distributions. More formally,
 * let \f$ a(x) = a_1 x_1 + \dots + a_n x_n \leq b \f$ be a row of the LP. Let further \f$ l_i, u_i \in R\f$ denote the
 * (finite) lower and upper bound, respectively, of the \f$ i \f$-th variable \f$x_i\f$.
 * Viewing every variable \f$x_i \f$ as (continuously) uniformly distributed within its bounds, we can approximately
 * understand the row activity \f$a(x)\f$ as a gaussian random variate with mean value \f$ \mu = E[a(x)] = \sum_i a_i\frac{l_i + u_i}{2}\f$
 * and variance \f$ \sigma^2 = \sum_i a_i^2 \sigma_i^2 \f$, with \f$ \sigma_i^2 = \frac{(u_i - l_i)^2}{12}\f$ for
 * continuous and \f$ \sigma_i^2 = \frac{(u_i - l_i + 1)^2 - 1}{12}\f$ for discrete variables.
 * With these two parameters, we can calculate the probability to satisfy the row in terms of the cumulative distribution
 * of the standard normal distribution: \f$ P(a(x) \leq b) = \Phi(\frac{b - \mu}{\sigma})\f$.
 *
 * The impact of a variable on the probability to satisfy a constraint after branching can be estimated by altering
 * the variable contribution to the sums described above. In order to keep the tree size small,
 * variables are preferred which decrease the probability
 * to satisfy a row because it is more likely to create infeasible subproblems.
 *
 * The selection of the branching variable is subject to the parameter @p scoreparam. For both branching directions,
 * an individual score is calculated. Available options for scoring methods are:
 * - @b d: select a branching variable with largest difference in satisfaction probability after branching
 * - @b l: lowest cumulative probability amongst all variables on all rows (after branching).
 * - @b h: highest cumulative probability amongst all variables on all rows (after branching).
 * - @b v: highest number of votes for lowering row probability for all rows a variable appears in.
 * - @b w: highest number of votes for increasing row probability for all rows a variable appears in.
 *
 * If the parameter @p usescipscore is set to @a TRUE, a single branching score is calculated from the respective
 * up and down scores as defined by the SCIP branching score function (see advanced branching parameter @p scorefunc),
 * otherwise, the variable with the single highest score is selected, and the maximizing direction is assigned
 * higher branching priority.
 *
 * The original idea of probability based branching schemes appeared in:
 *
 * J. Pryor and J.W. Chinneck:@n
 * Faster Integer-Feasibility in Mixed-Integer Linear Programs by Branching to Force Change@n
 * Computers and Operations Research, vol. 38, 2011, p. 1143–1152@n
 * (http://www.sce.carleton.ca/faculty/chinneck/docs/PryorChinneck.pdf)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_DISTRIBUTION_H__
#define __SCIP_BRANCH_DISTRIBUTION_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the distribution branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleDistribution(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
