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

/**@file   struct_expr.h
 * @brief  data definitions for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_EXPRESSION_H__
#define __SCIP_STRUCT_EXPRESSION_H__

#include "scip/def.h"
#include "nlpi/type_expr.h"
#include "nlpi/type_exprinterpret.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** operator data of an expression */
union SCIP_ExprOpData 
{
   int                   intval;             /**< index of a variable or parameter or a constant integer value */
   SCIP_Real             dbl;                /**< a constant double value */
   void*                 data;               /**< pointer to some data structure */
};

/** arithmetic expression node */
struct SCIP_Expr 
{
   SCIP_EXPROP           op;                 /**< operator of the node */
   int                   nchildren;          /**< number of children */
   SCIP_EXPR**           children;           /**< children nodes */
   SCIP_EXPROPDATA       data;               /**< operator data */
};

/** expression tree */
struct SCIP_ExprTree 
{
   BMS_BLKMEM*           blkmem;             /**< block memory data structure */
   SCIP_EXPR*            root;               /**< root node expression of expression tree */
   int                   nvars;              /**< number of variables */
   void**                vars;               /**< mapping of variable indices to user variables, may be NULL */
   int                   nparams;            /**< number of parameters (modifiable constants) in expression */
   SCIP_Real*            params;             /**< current values for parameters, or NULL if no parameters */
   SCIP_EXPRINTDATA*     interpreterdata;    /**< data of expression interpreter (evaluator) */
};

/** data of quadratic expression: sum_i coef_i x_i y_i */
struct SCIP_ExprData_Quadratic 
{
   SCIP_Real             constant;           /**< constant term */
   SCIP_Real*            lincoefs;           /**< linear coefficients of children */
   SCIP_QUADELEM*        quadelems;          /**< quadratic elements */
   int                   nquadelems;         /**< number of quadratic elements */
   SCIP_Bool             sorted;             /**< are the quadratic elements sorted? */
};

/** data of polynomial expression: constant + sum_i monom_i */
struct SCIP_ExprData_Polynomial 
{
   SCIP_Real             constant;           /**< constant term of polynomial */
   SCIP_EXPRDATA_MONOMIAL** monomials;       /**< monomials that constitute the polynomial */
   int                   monomialssize;      /**< size of monomials array */
   int                   nmonomials;         /**< number of monomials */
   SCIP_Bool             sorted;             /**< are the monomials sorted? */
};

/** data of monomial in polynomial expression: coef * prod_i child_i^exponent_i
 * we allow for real values exponents here 
 */
struct SCIP_ExprData_Monomial 
{
   SCIP_Real             coef;               /**< coefficient of monomial */
   int                   factorssize;        /**< size of factors arrays */
   int                   nfactors;           /**< number of factors */
   int*                  childidxs;          /**< children corresponding to factors */
   SCIP_Real*            exponents;          /**< value of exponent for each factor */
   SCIP_Bool             sorted;             /**< are the factors sorted (by childidx)? */
};

/* element in table of expression operands */
struct SCIPexprOpTableElement
{
  const char*           name;               /**< name of operand (used for printing) */
  int                   nargs;              /**< number of arguments (negative if not fixed) */
  SCIP_DECL_EXPREVAL    ((*eval));          /**< evaluation function */
  SCIP_DECL_EXPRINTEVAL ((*inteval));       /**< interval evaluation function */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_EXPRESSION_H__ */
