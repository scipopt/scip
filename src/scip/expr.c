/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   expr.c
 * @brief  functions for algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include "scip/scip_expr.h"
#include "scip/expr.h"

typedef struct SCIP_QuadExpr_QuadTerm  SCIP_QUADEXPR_QUADTERM;  /**< a single term associated to a quadratic variable */
typedef struct SCIP_QuadExpr_BilinTerm SCIP_QUADEXPR_BILINTERM; /**< a single bilinear term */

/** data for representation of an expression as quadratic */
struct SCIP_QuadExpr
{
   SCIP_Real                constant;        /**< a constant term */

   int                      nlinexprs;       /**< number of expressions that appear linearly */
   SCIP_EXPR**              linexprs;        /**< expressions that appear linearly */
   SCIP_Real*               lincoefs;        /**< coefficients of expressions that appear linearly */

   int                      nquadexprs;      /**< number of expressions in quadratic terms */
   SCIP_QUADEXPR_QUADTERM*  quadexprterms;   /**< array with quadratic expression terms */

   int                      nbilinexprterms; /**< number of bilinear expressions terms */
   SCIP_QUADEXPR_BILINTERM* bilinexprterms;  /**< bilinear expression terms array */

   SCIP_Bool                allexprsarevars; /**< whether all arguments (linexprs, quadexprterms[.].expr) are variable expressions */

   SCIP_EXPRCURV            curvature;       /**< curvature of the quadratic representation of the expression */
   SCIP_Bool                curvaturechecked;/**< whether curvature has been checked */
   SCIP_Bool                eigeninfostored; /**< whether the eigen information is stored */

   /* eigen decomposition information */
   SCIP_Real*               eigenvalues;     /**< eigenvalues of the Q matrix: size of nquadexprs */
   SCIP_Real*               eigenvectors;    /**< eigenvalues of the Q matrix; size of nquadexprs^2 */
};

/** data structure to store a single term associated to a quadratic variable */
struct SCIP_QuadExpr_QuadTerm
{
   SCIP_EXPR*            expr;               /**< quadratic expression */
   SCIP_Real             lincoef;            /**< linear coefficient of variable */
   SCIP_Real             sqrcoef;            /**< square coefficient of variable */

   int                   nadjbilin;          /**< number of bilinear terms this variable is involved in */
   int                   adjbilinsize;       /**< size of adjacent bilinear terms array */
   int*                  adjbilin;           /**< indices of associated bilinear terms */

   SCIP_EXPR*            sqrexpr;            /**< expression that was found to be the square of expr, or NULL if no square term (sqrcoef==0) */
};

/** data structure to store a single bilinear term coef * expr1 * expr2  (similar to SCIP_QUADELEM)
 * except for temporary reasons, we assume that the index of var1 is smaller than the index of var2
 */
struct SCIP_QuadExpr_BilinTerm
{
   SCIP_EXPR*            expr1;              /**< first factor of bilinear term */
   SCIP_EXPR*            expr2;              /**< second factor of bilinear term */
   SCIP_Real             coef;               /**< coefficient of bilinear term */
   int                   pos2;               /**< position of expr2's quadexprterm in quadexprterms */

   SCIP_EXPR*            prodexpr;           /**< expression that was found to be the product of expr1 and expr2 */
};
