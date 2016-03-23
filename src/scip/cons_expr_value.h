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

/**@file   cons_expr_constant.h
 * @brief  constant value operand handler
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_CONSTANT_H__
#define __SCIP_CONS_EXPR_CONSTANT_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for constant value operands and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeExprOperandValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates the data of a constant value operand */
EXTERN
SCIP_RETCODE SCIPcreateExprOperandValue(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< constant value operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   SCIP_Real                   value               /**< value to be stored */
   );

/** gets the value of a constant value operand */
EXTERN
SCIP_Real SCIPgetValueOperandValue(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< operand data */
   );


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_CONSTANT_H__ */
