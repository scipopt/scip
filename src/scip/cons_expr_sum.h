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

/**@file   cons_expr_sum.h
 * @brief  sum expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_SUM_H__
#define __SCIP_CONS_EXPR_SUM_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sum expressions and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprExprHdlrSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates a sum expression */
SCIP_EXPORT
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
SCIP_EXPORT
SCIP_Real* SCIPgetConsExprExprSumCoefs(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   );

/** gets the constant of a summation expression */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   );

/** sets the constant of a summation expression */
SCIP_EXPORT
void SCIPsetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant */
   );

/** appends an expression to a sum expression */
SCIP_EXPORT
SCIP_RETCODE SCIPappendConsExprExprSumExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_CONSEXPR_EXPR*   child,              /**< expression to be appended */
   SCIP_Real             childcoef           /**< child's coefficient */
   );

/** multiplies given sum expression by a constant */
void SCIPmultiplyConsExprExprSumByConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant that multiplies sum expression */
   );

/** reverse propagate a weighted sum of expressions in the given interval */
SCIP_EXPORT
SCIP_RETCODE
SCIPreverseConsExprExprPropagateWeightedSum(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nexprs,             /**< number of expressions to propagate */
   SCIP_CONSEXPR_EXPR**  exprs,              /**< expressions to propagate */
   SCIP_Real*            weights,            /**< weights of expressions in sum */
   SCIP_Real             constant,           /**< constant in sum */
   SCIP_INTERVAL         interval,           /**< constant + sum weight_i expr_i \in interval */
   SCIP_QUEUE*           reversepropqueue,   /**< queue used in reverse prop, pass to SCIPtightenConsExprExprInterval */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions,        /**< buffer to store the number of interval reductions */
   SCIP_Bool             force               /**< to force tightening */
   );
#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_SUM_H__ */
