/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_expression.h,v 1.3 2010/05/05 16:20:13 bzfviger Exp $"

/**@file   type_expression.h
 * @brief  type definitions for expressions and expression trees
 * @ingroup TYPEDEFINITIONS
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_EXPRESSION_H__
#define __SCIP_TYPE_EXPRESSION_H__

#ifdef __cplusplus
extern "C" {
#endif

/** Operators of expressions.
 */
enum SCIP_ExprOp {
   /**@name Terminals (Leaves) */
   /**@{ */
   SCIP_EXPR_VARIDX    =  1,  /**< variable given by index (stored in data.idx) */
   SCIP_EXPR_CONST     =  2,  /**< constant (value stored in data.dbl) */
   SCIP_EXPR_PARAM     =  3,  /**< parameter = a constant that can be modified (should not be simplified away) */
   /**@} */

   /**@name Simple Operands */
   /**@{ */
   SCIP_EXPR_PLUS      =  8,  /**< addition (2 operands) */
   SCIP_EXPR_MINUS     =  9,  /**< substraction (2 operands) */
   SCIP_EXPR_MUL       = 10,  /**< multiplication (2 operands) */
   SCIP_EXPR_DIV       = 11,  /**< division (2 operands) */
   SCIP_EXPR_SQUARE    = 12,  /**< square (1 operand) */
   SCIP_EXPR_SQRT      = 13,  /**< square root (1 operand) */
   SCIP_EXPR_POWER     = 14,  /**< power (x^y, 2 operands) */
   SCIP_EXPR_EXP       = 15,  /**< exponential (e^x, 1 operand) */
   SCIP_EXPR_LOG       = 16,  /**< natural logarithm (ln(x), 1 operand) */
   SCIP_EXPR_SIN       = 17,  /**< sinus (1 operand) */
   SCIP_EXPR_COS       = 18,  /**< cosinus (1 operand) */
   SCIP_EXPR_TAN       = 19,  /**< tangent (1 operand) */
   SCIP_EXPR_ERF       = 20,  /**< gaussian error function (1 operand) */
   SCIP_EXPR_ERFI      = 21,  /**< imaginary part of gaussian error function (1 operand) */
   SCIP_EXPR_MIN       = 22,  /**< minimum (2 operands) */
   SCIP_EXPR_MAX       = 23,  /**< maximum (2 operands) */
   SCIP_EXPR_ABS       = 24,  /**< absolute value (1 operand) */
   SCIP_EXPR_SIGN      = 25,  /**< sign of value (1 operand) */
   SCIP_EXPR_SIGNPOWER = 26,  /**< signed power (sign(x)|x|^y, 2 operands) */
   SCIP_EXPR_INTPOWER  = 27,  /**< power with integer exponent (1 operand!, exponent stored in expression data) */
   /**@} */

   /**@name Complex Operands
    * @TODO these expressions are not implemented yet!
    */
   /**@{ */
   SCIP_EXPR_SUM       = 64,  /**< summation sum_{i=1}^n op_i (n operands) */
   SCIP_EXPR_PRODUCT   = 65,  /**< product prod_{i=1}^n op_i (n operands) */
   SCIP_EXPR_LINEAR    = 66,  /**< linear term sum_{i=1}^n a_i op_i (n operands) */
   SCIP_EXPR_SIGNOMIAL = 67,  /**< signomial term prod_{i=1}^n op_i^c_i (n operands, for much later...) */
   SCIP_EXPR_QUADRATIC = 68,  /**< quadratic term sum_{i,j=1}^n a_{i,j} op_i op_j (n operands) */
   SCIP_EXPR_POLYNOM   = 69,  /**< polynomial term sum_{I} a_{I}ops^I (I a multiindex, n operands) */
   /**@} */

   SCIP_EXPR_LAST      = 99   /**< no expression, used for counting reasons */
};

typedef enum   SCIP_ExprOp     SCIP_EXPROP;     /**< expression operand */
typedef union  SCIP_ExprOpData SCIP_EXPROPDATA; /**< expression operand data */
typedef struct SCIP_Expr       SCIP_EXPR;       /**< expression */
typedef struct SCIP_ExprTree   SCIP_EXPRTREE;   /**< expression tree */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_EXPRESSION_H__ */
