/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   expr_sum.h
 * @ingroup EXPRHDLRS
 * @brief  sum expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_SUM_H__
#define __SCIP_EXPR_SUM_H__


#include "scip/scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sum expressions and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrSum(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup EXPRHDLRS
 *
 * @{
 *
 * @name Sum expression
 *
 * This expression handler provides the sum function, that is,
 * \f[
 *   x \mapsto c + \sum_{i=1}^n a_i x_i
 * \f]
 * for some constant c and constant coefficients \f$a_i\f$.
 *
 * @{
 */

/** creates a sum expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children */
   SCIP_Real*            coefficients,       /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real             constant,           /**< constant term of sum */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** sets the constant of a summation expression */
SCIP_EXPORT
void SCIPsetConstantExprSum(
   SCIP_EXPR*            expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant */
   );

/** appends an expression to a sum expression */
SCIP_EXPORT
SCIP_RETCODE SCIPappendExprSumExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< sum expression */
   SCIP_EXPR*            child,              /**< expression to be appended */
   SCIP_Real             childcoef           /**< child's coefficient */
   );

/** multiplies given sum expression by a constant */
SCIP_EXPORT
void SCIPmultiplyByConstantExprSum(
   SCIP_EXPR*            expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant that multiplies sum expression */
   );

/** constructs the expanded product of two sum expressions */
SCIP_EXPORT
SCIP_RETCODE SCIPmultiplyBySumExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           product,            /**< buffer where to store multiplied sums (expanded as sum) */
   SCIP_EXPR*            factor1,            /**< first sum */
   SCIP_EXPR*            factor2,            /**< second sum */
   SCIP_Bool             simplify,           /**< whether to simplify created terms and sum */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** constructs the expanded power of a sum expression
 *
 * @attention The number of terms in the expansion grows exponential with the exponent. Be aware of what you wish for.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPpowerExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           result,             /**< buffer where to store expanded power of sum */
   SCIP_EXPR*            base,               /**< sum */
   int                   exponent,           /**< exponent > 1 */
   SCIP_Bool             simplify,           /**< whether to simplify created terms and sum */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_SUM_H__ */
