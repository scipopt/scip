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

/**@file   cons_expr_sumprod.h
 * @brief  sum and product operand handlers
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_SUMPROD_H__
#define __SCIP_CONS_EXPR_SUMPROD_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for sum operands and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeExprOperandSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates the data of a summation operand */
EXTERN
SCIP_RETCODE SCIPcreateExprOperandSum(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   int                         ncoefficients,      /**< number of coefficients (i.e., number of children) */
   SCIP_Real*                  coefficients,       /**< array with coefficients for all operands (or NULL if all 1.0) */
   SCIP_Real                   constant            /**< constant term of sum */
   );

/** gets the coefficients of a summation operand */
EXTERN
SCIP_Real* SCIPgetCoefsOperandSum(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< sum operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   );

/** gets the constant of a summation operand */
EXTERN
SCIP_Real SCIPgetConstantOperandSum(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< sum operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   );


/** creates the handler for product operands and includes it into the expression constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeExprOperandProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** creates the data of a product operand */
EXTERN
SCIP_RETCODE SCIPcreateExprOperandProduct(
   SCIP*                       scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*              consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_OPERANDHDLR*  operandhdlr,        /**< variable operand handler */
   SCIP_CONSEXPR_OPERANDDATA** operanddata,        /**< pointer where to store data of operand */
   int                         nexponents,         /**< number of exponents (i.e., number of children) */
   SCIP_Real*                  exponents,          /**< array with exponents for all operands (or NULL if all 1.0) */
   SCIP_Real                   constant            /**< constant coefficient of product */
   );

/** gets the exponents of a product operand */
EXTERN
SCIP_Real* SCIPgetCoefsOperandProduct(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< product operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   );

/** gets the constant coefficient of a product operand */
EXTERN
SCIP_Real SCIPgetConstantOperandProduct(
   SCIP_CONSEXPR_OPERANDHDLR* operandhdlr,   /**< product operand handler */
   SCIP_CONSEXPR_OPERANDDATA* operanddata    /**< data of operand */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_SUMPROD_H__ */
