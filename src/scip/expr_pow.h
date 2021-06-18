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

/**@file   expr_pow.h
 * @ingroup EXPRHDLRS
 * @brief  power and signed power expression handlers
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_POW_H__
#define __SCIP_EXPR_POW_H__

#include "scip/scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup EXPRHDLRS
 *
 * @{
 *
 * @name Power and signed power expression.
 *
 * These expression handler provide the power function, that is,
 * \f[
 *   x \mapsto \begin{case}
 *     x^e & \textrm{if} x \geq 0 or e integral, \\
 *     \textrm{undefined}, & otherwise.
 *     \end{cases}
 * \f]
 * and the signed power function, that is,
 * \f[
 *   x \mapsto \textrm{sign}(e) |x|^e
 * \f]
 * for some exponent e.
 *
 * @{
 */

/** creates a power expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** creates a signpower expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprSignpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_Real             exponent,           /**< exponent of the power expression */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** indicates whether expression is of signpower-type */
SCIP_EXPORT
SCIP_Bool SCIPisExprSignpower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** @}
  * @}
  */

/** creates the handler for power expression and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrPow(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the handler for signed power expression and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrSignpower(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_POW_H__ */
