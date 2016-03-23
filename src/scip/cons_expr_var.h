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

/**@file   cons_expr_var.h
 * @brief  variable operand handler
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_VAR_H__
#define __SCIP_CONS_EXPR_VAR_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for variable operands and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeExprOperandVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates the data of a variable operand */
EXTERN
SCIP_RETCODE SCIPcreateExprOperandVar(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   SCIP_VAR*                   var                 /**< variable to be stored */
   );

/** gets the variable of a variable operand */
EXTERN
SCIP_VAR* SCIPgetVarOperandVar(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< operand data */
   );


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_VAR_H__ */
