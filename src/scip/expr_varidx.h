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

/**@file   expr_varidx.h
 * @ingroup EXPRHDLRS
 * @brief  handler for variable index expressions
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EXPR_VARIDX_H__
#define __SCIP_EXPR_VARIDX_H__


#include "scip/type_scip.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for variable index expressions and includes it into SCIP
 *
 * @ingroup ExprhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprhdlrVaridx(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup EXPRHDLRS
 *
 * @{
 *
 * @name Index variable expression
 *
 * This expression handler handles a variable that is given by a variable index. It cannot have children.
 *
 * This expression handler is used for expressions that are passed to a NLP solver via the NLPI.
 * @{
 */

/** creates a variable index expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateExprVaridx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   varidx,             /**< variable index to represent */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   );

/** indicates whether expression is varidx expression */
SCIP_EXPORT
SCIP_Bool SCIPisExprVaridx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   );

/** gives the index stored in a varidx expression */
SCIP_EXPORT
int SCIPgetIndexExprVaridx(
   SCIP_EXPR*            expr                /**< varindex expression */
   );

/** sets the index stored in a varidx expression */
SCIP_EXPORT
void SCIPsetIndexExprVaridx(
   SCIP_EXPR*            expr,               /**< varindex expression */
   int                   newindex            /**< new index */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_EXPR_VARIDX_H__ */
