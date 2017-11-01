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

/** computes the secant if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
EXTERN
SCIP_Bool SCIPcomputeSecantSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   );

/** computes the tangent at lower bound if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
EXTERN
SCIP_Bool SCIPcomputeLeftTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb                  /**< lower bound of argument variable */
   );

/** computes the tangent at upper bound if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
EXTERN
SCIP_Bool SCIPcomputeRightTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   );

/** computes the tangent at solution point if it is a feasible underestimating cut
 *  returns true if the cut was computed successfully
 */
EXTERN
SCIP_Bool SCIPcomputeSolTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub,                 /**< upper bound of argument variable */
   SCIP_Real             solpoint            /**< solution point to be separated */
   );

/** computes the tangent at some other point that goes through (lb,sin(lb)) and is underestimating
 *  returns true if the cut was computed successfully
 */
EXTERN
SCIP_Bool SCIPcomputeLeftMidTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Bool*            issecant,           /**< buffer to store whether cut is actually a secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   );

/** computes the tangent at some other point that goes through (ub,sin(ub)) and is underestimating
 *  returns true if the cut was computed successfully
 */
EXTERN
SCIP_Bool SCIPcomputeRightMidTangentSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            lincoef,            /**< buffer to store linear coefficient of secant */
   SCIP_Real*            linconst,           /**< buffer to store linear constant of secant */
   SCIP_Bool*            issecant,           /**< buffer to stroe whether cut is actually a secant */
   SCIP_Real             lb,                 /**< lower bound of argument variable */
   SCIP_Real             ub                  /**< upper bound of argument variable */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_SIN_H__ */
