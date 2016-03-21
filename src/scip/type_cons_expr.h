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

/**@file   type_cons_expr.h
 * @brief  (public) types of expression constraint
 * @author Stefan Vigerske
 *
 * These are in particular types that define the expressions in cons_expr
 * and that need to be accessed by the linear estimation plugins of cons_expr.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CONS_EXPR_H__
#define __SCIP_TYPE_CONS_EXPR_H__

#ifdef __cplusplus
extern "C" {
#endif


/** Operators of expressions. */
enum SCIP_ConsExpr_Operand
{
   /**@name Terminals (Leaves) */
   /**@{ */
   SCIP_CONSEXPR_VAR       =  1,  /**< variable given by SCIP_VAR* (stored in data.ptr) */
   SCIP_CONSEXPR_CONST     =  2,  /**< floating-point constant (value stored in data.dbl) */
   /* SCIP_CONSEXPR_PARAM     =  3, */  /**< parameter = a constant that can be modified (should not be simplified away) */
   /**@} */

   /**@name Import and frequently used operands */
   /**@{ */
   SCIP_CONSEXPR_SUM       =  4,  /**< summation with coefficients and a constant */
   SCIP_CONSEXPR_PRODUCT   =  5,  /**< a signomial, i.e., product a with a coefficient and rational-valued exponents */
   SCIP_CONSEXPR_INBOUNDS  =  6,  /**< indicator (0 or 1) whether child is within some constant (possibly infinite) bounds */
   /**@} */

   /**@name Further operands
    */
   /**@{ */
   SCIP_CONSEXPR_EXP       = 8,  /**< exponential (e^x, 1 operand) */
   SCIP_CONSEXPR_LOG       = 9,  /**< natural logarithm (ln(x), 1 operand) */
   /**@} */

   SCIP_CONSEXPR_LAST      = 10   /**< no expression, used for counting reasons */
};

/** variability of expression operands */
enum SCIP_ConsExpr_Variability
{
   SCIP_CONSEXPR_INVARIATE    = 0,   /**< no children */
   SCIP_CONSEXPR_UNIVARIATE   = 1,   /**< exactly one child */
   SCIP_CONSEXPR_BIVARIATE    = 2,   /**< exactly two children */
   SCIP_CONSEXPR_MULTIVARIATE = 3    /**< arbitrary number of children */
};

typedef enum  SCIP_ConsExpr_Variability  SCIP_CONSEXPR_VARIABILITY; /**< variability of expression operands */
typedef union SCIP_ConsExpr_Children     SCIP_CONSEXPR_CHILDREN;    /**< storage type for children of an expression */
typedef enum  SCIP_ConsExpr_Operand      SCIP_CONSEXPR_OPERAND;     /**< expression operand */
typedef union SCIP_ConsExpr_OperatorData SCIP_CONSEXPR_OPERANDDATA; /**< expression operand data */

typedef struct SCIP_ConsExpr_Expr        SCIP_CONSEXPR_EXPR;        /**< expression */


#ifdef __cplusplus
}
#endif

#endif /* __SCIP_TYPE_CONS_EXPR_H__ */
