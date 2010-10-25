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
#pragma ident "@(#) $Id: struct_expr.h,v 1.1 2010/10/25 04:27:34 bzfviger Exp $"

/**@file   struct_expression.h
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

#ifndef __SCIP_TYPE_VAR_H__
#define SCIP_VAR void
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** operator data of an expression */
union SCIP_ExprOpData {
   int                   intval;             /**< index of a variable or parameter or a constant integer value */
   SCIP_Real             dbl;                /**< a constant double value */
   void*                 data;               /**< pointer to some data structure */
};

/** arithmetic expression node */
struct SCIP_Expr {
   SCIP_EXPROP           op;         /**< operator of the node */
   int                   nchildren;  /**< number of children */
   SCIP_EXPR**           children;   /**< children nodes */
   SCIP_EXPROPDATA       data;       /**< operator data */
};

/** expression tree */
struct SCIP_ExprTree {
   BMS_BLKMEM*           blkmem;           /**< block memory data structure */
   SCIP_EXPR*            root;             /**< root node expression of expression tree */
   int                   nvars;            /**< number of variables */
   SCIP_VAR**            vars;             /**< mapping of variable indices to SCIP variables, may be NULL, e.g., if not used in context of SCIP */
   int                   nparams;          /**< number of parameters (modifiable constants) in expression */
   SCIP_Real*            params;           /**< current values for parameters, or NULL if no parameters */
   SCIP_EXPRINTDATA*     interpreterdata;  /**< data of expression interpreter (evaluator) */
};

/** data of quadratic expression: sum_i coef_i x_i y_i */
struct SCIP_ExprData_Quadratic {
   SCIP_QUADELEM*        quadelems;        /**< quadratic elements */
   int                   nquadelems;       /**< number of quadratic elements */
};

/** data of polynomial expression: constant + sum_i monom_i */
struct SCIP_ExprData_Polynom {
   SCIP_Real             constant;         /**< constant term of polynom */
   SCIP_EXPRDATA_MONOM** monoms;           /**< monoms that constitute the polynom */
   int                   nmonoms;          /**< number of monoms */
};

/** data of monom in polynomial expression: coef * prod_i child_i^exponent_i
 * we allow for real values exponents here */
struct SCIP_ExprData_Monom {
   SCIP_Real             coef;             /**< coefficient of monom */
   int                   nfactors;         /**< number of factors */
   int*                  childidxs;        /**< children corresponding to factors */
   SCIP_Real*            exponents;        /**< value of exponent for each factor */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_EXPRESSION_H__ */
