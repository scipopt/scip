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

/**@file   expr_trig.h
 * @ingroup EXPRHDLRS
 * @brief  handler for sin expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_TRIG_H__
#define __SCIP_EXPR_TRIG_H__

#include "scip/scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sin expressions and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrSin(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the handler for cos expressions and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrCos(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup EXPRHDLRS
 *
 * @{
 *
 * @name Sine and Cosine expression
 *
 * These expression handler provide the sine and cosine functions, that is,
 * \f[
 *   x \mapsto \sin(x)
 * \f]
 * and
 * \f[
 *   x \mapsto \cos(x).
 * \f]
 *
 * @{
 */

/** creates a sin expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprSin(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** creates a cos expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_TRIG_H__ */
