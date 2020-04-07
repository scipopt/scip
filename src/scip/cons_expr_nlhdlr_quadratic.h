/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_quadratic.h
 * @brief  nonlinear handler to handle quadratic expressions
 * @author Felipe Serrano
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_NLHDLR_QUADRATIC_H__
#define __SCIP_CONS_EXPR_NLHDLR_QUADRATIC_H__


#include "scip/scip.h"
#include "scip/cons_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes quadratic nonlinear handler to consexpr */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConsExprNlhdlrQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   );

/** data structure to store a single term associated to a quadratic variable
 */
struct SCIP_QuadExprTerm
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< quadratic expression */
   SCIP_Real             lincoef;            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef;            /**< square coefficient of variable */

   int                   nadjbilin;          /**< number of bilinear terms this variable is involved in */
   int                   adjbilinsize;       /**< size of adjacent bilinear terms array */
   int*                  adjbilin;           /**< indices of associated bilinear terms */
};
typedef struct SCIP_QuadExprTerm SCIP_QUADEXPRTERM;

/** data structure to store a single bilinear term (similar to SCIP_QUADELEM)
 * except for temporary reasons, we assume that the index of var1 is smaller than the index of var2
 */
struct SCIP_BilinExprTerm
{
   SCIP_CONSEXPR_EXPR*   expr1;              /**< first factor of bilinear term expr1 * expr2 */
   SCIP_CONSEXPR_EXPR*   expr2;              /**< second factor of bilinear term expr1 * expr2 */
   SCIP_Real             coef;               /**< coef of bilinear term expr1 * expr2 */
   int                   pos2;               /**< position of expr2's quadexprterm in quadexprterms */
};
typedef struct SCIP_BilinExprTerm SCIP_BILINEXPRTERM;

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_CONS_EXPR_NLHDLR_QUADRATIC_H__ */
