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

/**@file   cons_expr_nlhdlr_bilinear.h
 * @brief  bilinear nonlinear handler
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_NLHDLR_BILINEAR_H__
#define __SCIP_CONS_EXPR_NLHDLR_BILINEAR_H__

#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns an array of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR** SCIPgetConsExprNlhdlrBilinearExprs(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   );

/** returns an array of nonlinear handler expressions data of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPORT
SCIP_CONSEXPR_NLHDLREXPRDATA** SCIPgetConsExprNlhdlrBilinearExprsdata(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   );

/** returns the total number of expressions that have been detected by the bilinear nonlinear handler */
SCIP_EXPORT
int SCIPgetConsExprNlhdlrBilinearNExprs(
   SCIP_CONSEXPR_NLHDLR* nlhdlr              /**< nonlinear handler */
   );

/** adds a globally valid inequality of the form xcoef x <= ycoef y + constant to a product expression of the form x*y */
SCIP_EXPORT
SCIP_RETCODE SCIPaddConsExprNlhdlrBilinearIneq(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLR* nlhdlr,             /**< nonlinear handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< product expression */
   SCIP_Real             xcoef,              /**< x coefficient */
   SCIP_Real             ycoef,              /**< y coefficient */
   SCIP_Real             constant,           /**< constant part */
   SCIP_Bool*            success             /**< buffer to store whether inequality has been accepted */
   );

/** includes bilinear nonlinear handler to consexpr */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrBilinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_NLHDLR_BILINEAR_H__ */
