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

/**@file   expr_exp.h
 * @ingroup EXPRHDLRS
 * @brief  exponential expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_EXP_H__
#define __SCIP_EXPR_EXP_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup EXPRHDLRS
 *
 * @{
 *
 * @name Exponential value expression
 *
 * This expression handler provides the exponential function, that is,
 * \f[
 *   x \mapsto \exp(x).
 * \f]
 *
 * @{
 */

/** creates an exponential expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPR*            child,              /**< single child */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** indicates whether expression is of exp-type */
SCIP_EXPORT
SCIP_Bool SCIPisExprExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** @}
  * @}
  */

/** creates the handler for exponential expressions and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrExp(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_EXP_H__ */
