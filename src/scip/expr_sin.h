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

/**@file   expr_sin.h
 * @brief  handler for sin expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_SIN_H__
#define __SCIP_EXPR_SIN_H__

#include "scip/scip.h"
#include "scip/type_expr.h"

#define NEWTON_NITERATIONS    100
#define NEWTON_PRECISION      1e-12

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sin expressions and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprHdlrSin(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a sin expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** helper function to compute the new interval for child in reverse propagation */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeRevPropIntervalSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         parentbounds,       /**< bounds for sine expression */
   SCIP_INTERVAL         childbounds,        /**< bounds for child expression */
   SCIP_INTERVAL*        newbounds           /**< buffer to store new child bounds */
);

/** helper function to compute coefficients and constant term of a linear estimator at a given point
 *
 *  The function will try to compute the following estimators in that order:
 *  - soltangent: tangent at specified refpoint
 *  - secant: secant between the points (lb,sin(lb)) and (ub,sin(ub))
 *  - lmidtangent: tangent at some other point that goes through (lb,sin(lb))
 *  - rmidtangent: tangent at some other point that goes through (ub,sin(ub))
 *
 *  They are ordered such that a successful computation for one of them cannot be improved by following ones in terms
 *  of value at the reference point
 */
SCIP_EXPORT
SCIP_Bool SCIPcomputeEstimatorsTrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< sin or cos expression */
   SCIP_Real*            lincoef,            /**< buffer to store the linear coefficient */
   SCIP_Real*            linconst,           /**< buffer to store the constant term */
   SCIP_Real             refpoint,           /**< point at which to underestimate (can be SCIP_INVALID) */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate       /**< whether the estimator should be underestimating */
);

/** helper function to create initial cuts for sine and cosine separation
 *
 *  The following 5 cuts can be generated:
 *  - secant: secant between the points (lb,sin(lb)) and (ub,sin(ub))
 *  - ltangent/rtangent: tangents at the points (lb,sin(lb)) or (ub,sin(ub))
 *  - lmidtangent/rmidtangent: tangent at some other point that goes through (lb,sin(lb)) or (ub,sin(ub))
 *
 *  If one of the computations fails or turns out to be irrelevant, the respective argument pointer is set to NULL.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeInitialCutsTrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< sin or cos expression */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate,      /**< whether the cuts should be underestimating */
   SCIP_Real**           coefs,              /**< buffer to store coefficients of computed estimators */
   SCIP_Real*            constant,           /**< buffer to store constant of computed estimators */
   int*                  nreturned           /**< buffer to store number of estimators that have been computed */
   );

/* helper function that computs the curvature of a sine expression for given bounds and curvature of child */
SCIP_EXPORT
SCIP_EXPRCURV SCIPcomputeCurvatureSin(
   SCIP_EXPRCURV         childcurvature,     /**< curvature of child */
   SCIP_Real             lb,                 /**< lower bound of child */
   SCIP_Real             ub                  /**< lower bound of child */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_SIN_H__ */
