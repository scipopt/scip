/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_expr.h
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
   SCIP_EXPR_REALPOWER = 14,  /**< power with real exponent (1 operand!, assumed to be nonnegative, exponent is stored in expression data) */
   SCIP_EXPR_INTPOWER  = 15,  /**< power with integer exponent (1 operand!, exponent stored in expression data) */
   SCIP_EXPR_SIGNPOWER = 16,  /**< signed power (sign(x)|x|^a, 1 operand!, exponent is stored in expression data) */
   SCIP_EXPR_EXP       = 17,  /**< exponential (e^x, 1 operand) */
   SCIP_EXPR_LOG       = 18,  /**< natural logarithm (ln(x), 1 operand) */
   SCIP_EXPR_SIN       = 19,  /**< sinus (1 operand) */
   SCIP_EXPR_COS       = 20,  /**< cosinus (1 operand) */
   SCIP_EXPR_TAN       = 21,  /**< tangent (1 operand) */
   /* SCIP_EXPR_ERF       = 22, */  /**< gaussian error function (1 operand) */
   /* SCIP_EXPR_ERFI      = 23, */  /**< imaginary part of gaussian error function (1 operand) */
   SCIP_EXPR_MIN       = 24,  /**< minimum (2 operands) */
   SCIP_EXPR_MAX       = 25,  /**< maximum (2 operands) */
   SCIP_EXPR_ABS       = 26,  /**< absolute value (1 operand) */
   SCIP_EXPR_SIGN      = 27,  /**< sign of value (1 operand) */
   /**@} */

   /**@name Complex Operands
    */
   /**@{ */
   SCIP_EXPR_SUM       = 64,  /**< summation sum_{i=1}^n op_i (n operands) */
   SCIP_EXPR_PRODUCT   = 65,  /**< product prod_{i=1}^n op_i (n operands) */
   SCIP_EXPR_LINEAR    = 66,  /**< linear term sum_{i=1}^n a_i op_i (n operands) */
   SCIP_EXPR_QUADRATIC = 67,  /**< quadratic term sum_{i,j=1}^n a_{i,j} op_i op_j (n operands) */
   SCIP_EXPR_POLYNOMIAL= 68,  /**< polynomial term sum_{I} a_{I}ops^I (I a multiindex, n operands) */
   /**@} */

   SCIP_EXPR_LAST      = 69   /**< no expression, used for counting reasons */
};

/** Curvature types */
enum SCIP_ExprCurv
{
   SCIP_EXPRCURV_UNKNOWN    = 0,             /**< unknown curvature (or indefinite) */
   SCIP_EXPRCURV_CONVEX     = 1,             /**< convex */
   SCIP_EXPRCURV_CONCAVE    = 2,             /**< concave */
   SCIP_EXPRCURV_LINEAR     = SCIP_EXPRCURV_CONVEX | SCIP_EXPRCURV_CONCAVE/**< linear = convex and concave */
};

typedef enum   SCIP_ExprOp      SCIP_EXPROP;     /**< expression operand */
typedef union  SCIP_ExprOpData  SCIP_EXPROPDATA; /**< expression operand data */
typedef struct SCIP_Expr        SCIP_EXPR;       /**< expression */
typedef struct SCIP_ExprTree    SCIP_EXPRTREE;   /**< expression tree */
typedef enum   SCIP_ExprCurv    SCIP_EXPRCURV;   /**< curvature types */

/** An element of a quadratic term: two variable indices and a coefficient.
 * The convention is to have idx1 <= idx2.
 */
struct SCIP_QuadElement
{
   int                   idx1;             /**< index of first variable */
   int                   idx2;             /**< index of second variable */
   SCIP_Real             coef;             /**< value of coefficient at position (idx1, idx2) */
};
/* We have defined struct SCIP_QuadElement here (instead of type_expression.h) to allow fast access, allocation, and copying. (similar to SCIP_INTERVAL) */

typedef struct SCIP_QuadElement        SCIP_QUADELEM;           /**< element of a quadratic term */
typedef struct SCIP_ExprData_Quadratic SCIP_EXPRDATA_QUADRATIC; /**< the data of a quadratic expression (SCIP_EXPR_QUADRATIC) */

typedef struct SCIP_ExprData_Monomial   SCIP_EXPRDATA_MONOMIAL;   /**< a monomial as part of the data in a polynomial expression */
typedef struct SCIP_ExprData_Polynomial SCIP_EXPRDATA_POLYNOMIAL; /**< the data of a polynomial expression (SCIP_EXPR_POLYNOMIAL) */

#define SCIP_EXPR_DEGREEINFINITY 65535       /**< value that stands for an infinite degree of an expression (see SCIPexprGetMaxDegree) */

/** signature of an expression (pointwise) evaluation function
 * The function should return nan, inf, or -inf in result if the function is undefined for the given arguments.
 *
 * - opdata    operand data
 * - nargs     number of arguments
 * - argvals   values of arguments
 * - varvals   values for variables
 * - paramvals values for parameters
 * - result    buffer where to store result of evaluation
 */
#define SCIP_DECL_EXPREVAL(x) SCIP_RETCODE x (SCIP_EXPROPDATA opdata, int nargs, SCIP_Real* argvals, SCIP_Real* varvals, SCIP_Real* paramvals, SCIP_Real* result)

/** signature of an expression (interval) evaluation function
 * The function should return an empty interval if the function is undefined for the given arguments.
 *
 * - infinity  value for infinity
 * - opdata    operand data
 * - nargs     number of arguments
 * - argvals   interval values of arguments
 * - varvals   interval values for variables
 * - paramvals values for parameters
 * - result    buffer where to store result of evaluation
 */
#define SCIP_DECL_EXPRINTEVAL(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_EXPROPDATA opdata, int nargs, SCIP_INTERVAL* argvals, SCIP_INTERVAL* varvals, SCIP_Real* paramvals, SCIP_INTERVAL* result)

/** signature of a simple expression curvature check function
 *
 * - infinity  value for infinity
 * - opdata    operand data
 * - nargs     number of arguments
 * - argbounds bounds on value of arguments
 * - argcurv   curvature of arguments
 * - paramvals values for parameters
 * - result    buffer where to store result of curvature check
 */
#define SCIP_DECL_EXPRCURV(x) SCIP_RETCODE x (SCIP_Real infinity, SCIP_EXPROPDATA opdata, int nargs, SCIP_INTERVAL* argbounds, SCIP_EXPRCURV* argcurv, SCIP_EXPRCURV* result)

/** signature of a expression data copy function
 *
 * - blkmem       block memory
 * - nchildren    number of children in expression
 * - opdatasource source expression data
 * - opdatatarget pointer to target expression data
 */
#define SCIP_DECL_EXPRCOPYDATA(x) SCIP_RETCODE x (BMS_BLKMEM* blkmem, int nchildren, SCIP_EXPROPDATA opdatasource, SCIP_EXPROPDATA* opdatatarget)

/** signature of a expression data free function
 *
 * - blkmem       block memory
 * - nchildren    number of children in expression
 * - opdata       expression data to free
 */
#define SCIP_DECL_EXPRFREEDATA(x) void x (BMS_BLKMEM* blkmem, int nchildren, SCIP_EXPROPDATA opdata)

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_EXPRESSION_H__ */
