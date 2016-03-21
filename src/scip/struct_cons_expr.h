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

/**@file   struct_cons_expr.h
 * @brief  (public) data structures of expression constraints
 * @author Stefan Vigerske
 *
 * These are in particular data structures to manage the expressions in cons_expr
 * and that need to be accessed by the linear estimation plugins of cons_expr.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CONS_EXPR_H__
#define __SCIP_STRUCT_CONS_EXPR_H__

#include "scip/type_cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** operator data of an expression */
union SCIP_ConsExpr_OperatorData
{
   int                   intval;             /**< index of a parameter or a constant integer value */
   SCIP_Real             dbl;                /**< a constant double value */
   void*                 data;               /**< pointer to some data structure */
};

/** union for storing one, two, or many children */
union SCIP_ConsExpr_Children
{
   SCIP_CONSEXPR_EXPR*   child;              /**< child expression of a univariate expression */
   struct
   {
      SCIP_CONSEXPR_EXPR* child1;            /**< first child of a bivariate expression */
      SCIP_CONSEXPR_EXPR* child2;            /**< second child of a bivariate expression */
   } pair;
   struct
   {
      int                  nchildren;        /**< number of children of a multivariate expression */
      int                  childrensize;     /**< length of children array */
      SCIP_CONSEXPR_EXPR** children;         /**< children expressions of multivariate expression */
   } array;
};


/** a node in the expression graph that is handled by the expression constraint handler */
struct SCIP_ConsExpr_Expr
{
   SCIP_CONSEXPR_OPERAND     op;             /**< operand of expression */
   SCIP_CONSEXPR_OPERANDDATA opdata;         /**< operand data */

   SCIP_CONSEXPR_VARIABILITY variability;    /**< variability of constraint (in-, uni-, bi-, multivariate) */
   SCIP_CONSEXPR_CHILDREN    children;       /**< children of expression, interpretation of union depends on variability */

   int                       nuses;          /**< reference counter */
};




#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_CONS_EXPR_H__ */
