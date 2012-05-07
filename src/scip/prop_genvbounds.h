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

/**@file    prop_genvbounds.h
 * @ingroup PROPAGATORS
 * @brief   generalized variable bounds propagator
 * @author  Stefan Weltge
 *
 *  A generalized variable bound is a linear inequality of the form
 *
 *     c * x_i >= sum(a_j * x_j) + d * primal_bound + const,
 *
 *  where x_i's coefficient c is either 1 or -1 and "primal_bound" is an upper bound on the optimal objective value,
 *  which may improve during the solving process. In SCIP, generalized variable bounds are used for providing bounds on
 *  the LHS's variable x_i. If the above inequality is valid, the following bounds, depending on x_i's coefficient, are
 *  also valid:
 *
 *     c = 1   =>   x_i >=   minactivity(sum(a_j * x_j)) + d * primal_bound + const
 *
 *     c = -1  =>   x_i <= - minactivity(sum(a_j * x_j)) - d * primal_bound - const.
 *
 *  Note that for feasible problems, d <= 0 must hold. If d < 0 a decrease of the primal bound causes an improvement of
 *  the provided bound. Similarly, if a_j > 0 (< 0), a tightened lower (upper) bound of a variable x_j also yields a
 *  better bound for x_i.
 *
 *  The genvbounds propagator sorts its stored generalized variable bounds topologically in the following order: A
 *  generalized variable bound A (c * x_i >= ...) preceeds a generalized variable bound B if the left-hand side variable
 *  of A appears in the right-hand side of B with sign of its coefficient equal to c; i.e., if A is propagated and
 *  tightens the corresponding bound of x_i, then the minactivity on the right-hand side of B increases. We assume that
 *  this order is acyclic for the generalized variable bounds added. Under this condition, propagating the generalized
 *  variable bounds in a topological order ensures that all propagations are found in one round.
 *
 *  Both global and local propagation is applied: If the primal bound improves, generalized variable bounds with a
 *  nonzero coefficient d are enforced in order to tighten global bounds using the global variable bounds for computing
 *  the minactivity. Independently, the genvbounds propagator catches events SCIP_EVENTTYPE_LBTIGHTENED and
 *  SCIP_EVENTTYPE_UBTIGHTENED, i.e., locally tightened bounds of variables that occur in the right-hand sides of
 *  generalized variable bounds, in order to perform an efficient local propagation when called.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_GENVBOUNDS_H__
#define __SCIP_PROP_GENVBOUNDS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** adds a generalized variable bound to the genvbounds propagator; if there is already a genvbound for the bound
 *  "boundtype" of variable "var", it will be replaced
 */
extern
SCIP_RETCODE SCIPgenVBoundAdd(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_PROP*            genvboundprop,       /**< genvbound propagator */
   SCIP_VAR**            vars,                /**< array of RHSs variables */
   SCIP_VAR*             var,                 /**< LHSs variable */
   SCIP_Real*            coefs,               /**< array of coefficients for the RHSs variables */
   int                   ncoefs,              /**< size of coefs array */
   SCIP_Real             coefprimalbound,     /**< nonpositive value of the primal bounds multiplier */
   SCIP_Real             constant,            /**< constant term */
   SCIP_BOUNDTYPE        boundtype            /**< type of bound provided by the genvbound */
   );

/** creates the genvbounds propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePropGenvbounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
