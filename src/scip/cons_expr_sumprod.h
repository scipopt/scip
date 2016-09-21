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

/**@file   cons_expr_sumprod.h
 * @brief  sum and product expression handlers
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_SUMPROD_H__
#define __SCIP_CONS_EXPR_SUMPROD_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sum expressions and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConsExprExprHdlrSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates a sum expression */
EXTERN
SCIP_RETCODE SCIPcreateConsExprExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real*            coefficients,       /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real             constant            /**< constant term of sum */
   );

/** gets the coefficients of a summation expression */
EXTERN
SCIP_Real* SCIPgetConsExprExprSumCoefs(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   );

/** gets the constant of a summation expression */
EXTERN
SCIP_Real SCIPgetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   );

EXTERN
void SCIPsetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant */
   );

/** creates the handler for product expressions and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConsExprExprHdlrProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates a product expression */
EXTERN
SCIP_RETCODE SCIPcreateConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real             constant            /**< constant coefficient of product */
   );

/** gets the constant coefficient of a product expression */
EXTERN
SCIP_Real SCIPgetConsExprExprProductCoef(
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   );

/** appends an expression to a sum expression */
EXTERN
SCIP_RETCODE SCIPappendConsExprExprSumExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_CONSEXPR_EXPR*   child,              /**< expression to be appended */
   SCIP_Real             childcoef           /**< child's coefficient */
   );

/** appends an expression to a product expression */
EXTERN
SCIP_RETCODE SCIPappendConsExprExprProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< product expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
   );

/** multiplies given sum expr by a constant */
void SCIPmultiplyConsExprExprSumByConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant that multiplies sum expression */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_SUMPROD_H__ */
