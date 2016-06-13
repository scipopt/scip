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
 * @author Benjamin Mueller
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

/** generic data and callback methods of an expression handler */
struct SCIP_ConsExpr_ExprHdlr
{
   char*                         name;       /**< expression handler name */
   char*                         desc;       /**< expression handler description (can be NULL) */
   SCIP_CONSEXPR_EXPRHDLRDATA*   data;       /**< data of handler */
   unsigned int                  precedence; /**< precedence of expression operation relative to other expression (used for printing) */

   SCIP_DECL_CONSEXPR_EXPRCOPYHDLR((*copyhdlr));  /**< handler copy callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRFREEHDLR((*freehdlr));  /**< handler free callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA((*copydata));  /**< data copy callback, or NULL for expressions that have no data */
   SCIP_DECL_CONSEXPR_EXPRFREEDATA((*freedata));  /**< data free callback, or NULL for expressions that have no data or which data does not need to be freed */
   SCIP_DECL_CONSEXPR_EXPRCMP((*compare));        /**< compare callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print));        /**< print callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRPARSE((*parse));        /**< parse callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPREVAL((*eval));          /**< point evaluation callback (can never be NULL) */
   SCIP_DECL_CONSEXPR_EXPRINTEVAL((*inteval));    /**< interval evaluation callback (can be NULL) */
};

/** a node in the expression graph that is handled by the expression constraint handler */
struct SCIP_ConsExpr_Expr
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;      /**< expression type (as pointer to its handler) */
   SCIP_CONSEXPR_EXPRDATA* exprdata;      /**< expression data */

   int                     nchildren;     /**< number of children */
   int                     childrensize;  /**< length of children array */
   SCIP_CONSEXPR_EXPR**    children;      /**< children expressions */

   int                     nuses;         /**< reference counter */

   /* point-evaluation */
   unsigned int            evaltag;       /**< tag of point for which the expression has been evaluated last, or 0 */
   SCIP_Real               evalvalue;     /**< value of expression from last evaluation (corresponding to evaltag) */

   /* interval-evaluation */
   unsigned int            intevaltag;    /**< tag of domains for which tag for which the expression has been evaluated last, or 0 */
   SCIP_INTERVAL           interval;      /**< interval from the last interval evaluation */

   /* expression walker data */
   SCIP_CONSEXPR_EXPR*     walkparent;    /**< parent expression in expression walk */
   int                     walkcurrentchild; /**< child that is currently visited (or will be visited next) by expression walk */
   SCIP_CONSEXPREXPRWALK_IO walkio;       /**< space for walker callback to store some (temporary) data, e.g., to simulate input or output values of a recursive call */
};




#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_CONS_EXPR_H__ */
