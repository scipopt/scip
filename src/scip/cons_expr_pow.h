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

/**@file   cons_expr_pow.h
 * @brief  power expression handler
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_POW_H__
#define __SCIP_CONS_EXPR_POW_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a power expression */
EXTERN
SCIP_RETCODE SCIPcreateConsExprExprPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child,              /**< single child */
   SCIP_Real             exponent            /**< exponent of the power expression */
   );

/** creates the handler for power expression and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConsExprExprHdlrPow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** gets the exponent of a power expression */
EXTERN
SCIP_Real SCIPgetConsExprExprPowExponent(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_POW_H__ */
