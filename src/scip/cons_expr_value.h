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

/**@file   cons_expr_value.h
 * @brief  constant value expression handler
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_VALUE_H__
#define __SCIP_CONS_EXPR_VALUE_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for constant value expression and includes it into the expression constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprExprHdlrValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates constant value expression */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsExprExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_Real             value               /**< value to be stored */
   );

/** gets the value of a constant value expression */
SCIP_EXPORT
SCIP_Real SCIPgetConsExprExprValueValue(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   );


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_VALUE_H__ */
