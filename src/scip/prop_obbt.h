/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_obbt.h
 * @brief  optimality-based bound tightening propagator
 * @author Stefan Weltge
 *
 * In Optimality-Based Bound Tightening (OBBT), we solve auxiliary LPs of the form
 *
 *    min/max { x_i : x in P' },
 *
 * where P' is the current LP relaxation restricted by the primal cutoff constraint c*x <= z, z the current cutoff
 * bound. Trivially, the optimal objective value of this LP provides a valid lower/upper bound on variable x_i.
 *
 * Since solving LPs may be expensive, the propagator inspects solutions x in P' and does not run for variable bounds
 * which are tight at x: First, we check SCIP's last LP relaxation solution. Second, we solve a sequence of filtering
 * LP's min/max{ sum( w_i * x_i ) : x in P' } in order to push several variables towards one of their bounds in one LP
 * solve. Third, we inspect all solutions of the auxiliary LPs solved along the way.
 *
 * By default, OBBT is only applied for nonbinary variables that occur in nonlinear constraints.
 *
 * Additionally, the propagator uses the dual solution of the auxiliary LPs to construct globally valid generalized
 * variable bounds which may be propagated during the branch-and-bound search.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_OBBT_H__
#define __SCIP_PROP_OBBT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the obbt propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePropObbt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
