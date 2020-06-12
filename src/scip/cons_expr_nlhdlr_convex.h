/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_convex.h
 * @brief  nonlinear handlers for convex and concave expressions, respectively
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_NLHDLR_CONVEX_H__
#define __SCIP_CONS_EXPR_NLHDLR_CONVEX_H__

#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes convex nonlinear handler to consexpr */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** includes concave nonlinear handler to consexpr */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrConcave(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** checks whether a given expression is convex or concave w.r.t. the original variables
 *
 * This function uses the methods that are used in the detection algorithm of the convex nonlinear handler.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPhasConsExprExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_EXPRCURV         curv,               /**< curvature to check for */
   SCIP_Bool*            success,            /**< buffer to store whether expression has curvature curv (w.r.t. original variables) */
   SCIP_HASHMAP*         assumevarfixed      /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_NLHDLR_CONVEX_H__ */
