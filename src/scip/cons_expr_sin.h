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

/**@file   cons_expr_sin.h
 * @brief  handler for sin expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_SIN_H__
#define __SCIP_CONS_EXPR_SIN_H__

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sin.h"

#define SCIP_PRIVATE_ROWPREP
#include "scip/cons_quadratic.h"

#define NEWTON_NITERATIONS    100
#define NEWTON_PRECISION      1e-12

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sin expressions and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConsExprExprHdlrSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates a sin expression */
EXTERN
SCIP_RETCODE SCIPcreateConsExprExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   );

/** helper function to compute the new interval for child in reverse propagation */
EXTERN
SCIP_RETCODE SCIPcomputeRevPropIntervalSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTERVAL         parentbounds,       /**< bounds for sine expression */
   SCIP_INTERVAL         childbounds,        /**< bounds for child expression */
   SCIP_INTERVAL*        newbounds           /**< buffer to store new child bounds */
);

/** helper function to create cuts for point- or initial separation for sine and cosine expressions
 *
 *  A total of 6 different cuts can be generated. All except soltangent are independent of a specific solution and
 *  use only the bounds of the child variable. If their pointers are passed with NULL, the respective computation
 *  is not performed at all. If one of the computations fails or turns out to be irrelevant, the respective argument
 *  pointer is set to NULL
 */
EXTERN
SCIP_RETCODE SCIPcomputeCutsSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_ROWPREP**        secant,             /**< pointer to store the secant */
   SCIP_ROWPREP**        ltangent,           /**< pointer to store the left tangent */
   SCIP_ROWPREP**        rtangent,           /**< pointer to store the right tangent */
   SCIP_ROWPREP**        lmidtangent,        /**< pointer to store the left middle tangent */
   SCIP_ROWPREP**        rmidtangent,        /**< pointer to store the right middle tangent */
   SCIP_ROWPREP**        soltangent,         /**< pointer to store the solution tangent */
   SCIP_Real             refpoint,           /**< point that is to be seperated (can be SCIP_INVALID) */
   SCIP_Real             childlb,            /**< lower bound of child variable */
   SCIP_Real             childub,            /**< upper bound of child variable */
   SCIP_Bool             underestimate       /**< whether the cuts should be underestimating */
   );

/* helper function that computs the curvature of a sine expression for given bounds and curvature of child */
EXTERN
SCIP_EXPRCURV SCIPcomputeCurvatureSin(
   SCIP_EXPRCURV         childcurvature,     /**< curvature of child */
   SCIP_Real             lb,                 /**< lower bound of child */
   SCIP_Real             ub                  /**< lower bound of child */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_SIN_H__ */
