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

/**@file   cons_expr_nlhdlr_quadratic.c
 * @brief  nonlinear handler to handle quadratic expressions
 * @author Felipe Serrano
 *
 * Some definitions:
 * - a BILINEXPRTERM is a product of two expressions
 * - a SCIP_QUADEXPRTERM stores an expression expr that is known to appear in a nonlinear, quadratic term, that is
 *   expr^2 or expr * other_expr. It stores its sqrcoef (that can be 0), its linear coef and all the bilinear expression
 *   terms in which expr appears.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_nlhdlr_quadratic.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_rowprep.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "quadratic"
#define NLHDLR_DESC               "handler for quadratic expressions"
#define NLHDLR_DETECTPRIORITY     100
#define NLHDLR_ENFOPRIORITY       100

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_CONSEXPR_QUADEXPR* quaddata;         /**< data of quadratic form as stored in expr (stored here again for convenient access) */

   SCIP_EXPRCURV         curvature;          /**< curvature of the quadratic representation of the expression */

   SCIP_INTERVAL         linactivity;        /**< activity of linear part */

   /* activities of quadratic parts as defined in nlhdlrIntevalQuadratic */
   SCIP_Real             minquadfiniteact;   /**< minimum activity of quadratic part where only terms with finite min
                                               activity contribute */
   SCIP_Real             maxquadfiniteact;   /**< maximum activity of quadratic part where only terms with finite max
                                               activity contribute */
   int                   nneginfinityquadact;/**< number of quadratic terms contributing -infinity to activity */
   int                   nposinfinityquadact;/**< number of quadratic terms contributing +infinity to activity */
   SCIP_INTERVAL*        quadactivities;     /**< activity of each quadratic term as defined in nlhdlrIntevalQuadratic */
   SCIP_INTERVAL         quadactivity;       /**< activity of quadratic part (sum of quadactivities) */
   unsigned int          activitiestag;      /**< value of activities tag when activities were computed */
};

/*
 * static methods
 */


/** returns whether a quadratic form is "propagable"
 *
 * It is propagable, if a variable (aka child expr) appears at least twice, which is the case if at least two of the following hold:
 * - it appears as a linear term (coef*expr)
 * - it appears as a square term (coef*expr^2)
 * - it appears in a bilinear term
 * - it appears in another bilinear term
 */
static
SCIP_Bool isPropagable(
   SCIP_CONSEXPR_QUADEXPR* quaddata          /**< quadratic representation data */
   )
{
   int nquadexprs;
   int i;

   SCIPgetConsExprQuadraticData(quaddata, NULL, NULL, NULL, NULL, &nquadexprs, NULL);

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      int nadjbilin;

      SCIPgetConsExprQuadraticQuadTermData(quaddata, i, NULL, &lincoef, &sqrcoef, &nadjbilin, NULL, NULL);

      if( (lincoef != 0.0) + (sqrcoef != 0.0) + nadjbilin >= 2 )  /*lint !e514*/ /* actually MIN(2, nadjbilin), but we check >= 2 */
         return TRUE;
   }

   return FALSE;
}

/** returns whether a quadratic term is "propagable"
 *
 * A term is propagable, if its variable (aka child expr) appears at least twice, which is the case if at least two of the following hold:
 * - it appears as a linear term (coef*expr)
 * - it appears as a square term (coef*expr^2)
 * - it appears in a bilinear term
 * - it appears in another bilinear term
 */
static
SCIP_Bool isPropagableTerm(
   SCIP_CONSEXPR_QUADEXPR* quaddata,         /**< quadratic representation data */
   int                     idx               /**< index of quadratic term to consider */
   )
{
   SCIP_Real lincoef;
   SCIP_Real sqrcoef;
   int nadjbilin;

   SCIPgetConsExprQuadraticQuadTermData(quaddata, idx, NULL, &lincoef, &sqrcoef, &nadjbilin, NULL,  NULL);

   return (lincoef != 0.0) + (sqrcoef != 0.0) + nadjbilin >= 2;  /*lint !e514*/ /* actually MIN(2, nadjbilin), but we check >= 2 */
}

/** solves a quadratic equation \f$ a expr^2 + b expr \in rhs \f$ (with b an interval) and reduces bounds on expr or
 * deduces infeasibility if possible; expr is quadexpr.expr
 */
static
SCIP_RETCODE propagateBoundsQuadExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression for which to solve */
   SCIP_Real             sqrcoef,            /**< square coefficient */
   SCIP_INTERVAL         b,                  /**< interval acting as linear coefficient */
   SCIP_INTERVAL         rhs,                /**< interval acting as rhs */
   SCIP_QUEUE*           reversepropqueue,   /**< queue used in reverse prop, pass to SCIPtightenConsExprExprInterval */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions,        /**< buffer to store the number of interval reductions */
   SCIP_Bool             force               /**< to force tightening */
   )
{
   SCIP_INTERVAL a;
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Propagating <expr> by solving a <expr>^2 + b <expr> in rhs, where <expr> is: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "expr in [%g, %g], a = %g, b = [%g, %g] and rhs = [%g, %g]\n",
         SCIPintervalGetInf(SCIPgetConsExprExprActivity(scip, expr)),
         SCIPintervalGetSup(SCIPgetConsExprExprActivity(scip, expr)), sqrcoef, b.inf, b.sup,
         rhs.inf, rhs.sup);
#endif

   /* compute solution of a*x^2 + b*x in rhs */
   SCIPintervalSet(&a, sqrcoef);
   SCIPintervalSolveUnivariateQuadExpression(SCIP_INTERVAL_INFINITY, &newrange, a, b, rhs, SCIPgetConsExprExprActivity(scip, expr));

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Solution [%g, %g]\n", newrange.inf, newrange.sup);
#endif

   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, conshdlr, expr, newrange, force, reversepropqueue, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** solves a linear equation \f$ b expr \in rhs \f$ (with b a scalar) and reduces bounds on expr or deduces infeasibility if possible */
static
SCIP_RETCODE propagateBoundsLinExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression for which to solve */
   SCIP_Real             b,                  /**< linear coefficient */
   SCIP_INTERVAL         rhs,                /**< interval acting as rhs */
   SCIP_QUEUE*           reversepropqueue,   /**< queue used in reverse prop, pass to SCIPtightenConsExprExprInterval */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions,        /**< buffer to store the number of interval reductions */
   SCIP_Bool             force               /**< to force tightening */
   )
{
   SCIP_INTERVAL newrange;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Propagating <expr> by solving %g <expr> in [%g, %g], where <expr> is: ", b, rhs.inf, rhs.sup);
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   /* compute solution of b*x in rhs */
   SCIPintervalDivScalar(SCIP_INTERVAL_INFINITY, &newrange, rhs, b);

#ifdef DEBUG_PROP
   SCIPinfoMessage(scip, NULL, "Solution [%g, %g]\n", newrange.inf, newrange.sup);
#endif

   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, conshdlr, expr, newrange, force, reversepropqueue, infeasible, nreductions) );

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataQuadratic)
{  /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   if( (*nlhdlrexprdata)->quadactivities != NULL )
   {
      int nquadexprs;
      SCIPgetConsExprQuadraticData((*nlhdlrexprdata)->quaddata, NULL, NULL, NULL, NULL, &nquadexprs, NULL);
      SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->quadactivities, nquadexprs);
   }

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree
 *
 * A term is quadratic if:
 * - It is a product expression of two expressions
 * - It is power expression of an expression with exponent 2.0
 *
 * We define a propagable quadratic expression as a quadratic expression whose termwise propagation does not yield the
 * best propagation. In other words, is a quadratic expression that suffers from the dependency problem.
 *
 * Specifically, a propagable quadratic expression is a sum expression such that there is at least one expr that appears
 * at least twice (because of simplification, this means it appears in a quadratic terms and somewhere else). For
 * example: x^2 + y^2 is not a propagable quadratic expression; x^2 + x is a propagable quadratic expression; x^2 + x *
 * y is also a propagable quadratic expression
 *
 * Furthermore, we distinguish between propagable and non-propagable terms. A term is propagable if any the expressions
 * involved in it appear somewhere else. For example, x*y + z^2 + z is a propagable quadratic, the term 'x*y' is
 * non-propagable, and 'z^2' is propagable. For propagation, non-propagable terms are handled as if they were linear
 * terms, that is, we do not use the activity of 'x' and 'y' to compute the activity of 'x*y' but rather we use directly
 * the activity of 'x*y'. Similarly, we do not backward propagate to 'x' and 'y' (the product expr handler will do
 * this), but we backward propagate to 'x*y'. More technically, we register 'x*y' for its activity usage, rather than
 * 'x' and 'y'.
 *
 * For propagation, we store the quadratic in our data structure in the following way: We count how often a variable
 * appears. Then, a bilinear product expr_i * expr_j is stored as expr_i * expr_j if # expr_i appears > # expr_j
 * appears. When # expr_i appears == # expr_j appears, it then it will be stored as expr_i * expr_j if and only if
 * expr_i < expr_j, where '<' is the expression order (see Ordering Rules in cons_expr.c documentation).
 * Heuristically, this should be useful for propagation. The intuition is that by factoring out the variable that
 * appears most often we should be able to take care of the dependency problem better.
 *
 * Simple convex quadratics like x^2 + y^2 are ignored since the default nlhdlr will take care of them.
 * More complicated convex quadratics are handled here. (TODO: extended formulation using eigen-decomposition?)
 * Note that simple convex quadratics are exactly the non-propagable quadratics.
 *
 * TODO: when nonconvex quadratics are handled, detection should be revisited.
 *
 * @note:
 * - the expression needs to be simplified (in particular, it is assumed to be sorted)
 * - common subexpressions are also assumed to have been identified, the hashing will fail otherwise!
 *
 * Sorted implies that:
 *  - expr < expr^2: bases are the same, but exponent 1 < 2
 *  - expr < expr * other_expr: u*v < w holds if and only if v < w (OR8), but here w = u < v, since expr comes before
 *  other_expr in the product
 *  - expr < other_expr * expr: u*v < w holds if and only if v < w (OR8), but here v = w
 *
 *  Thus, if we see somebody twice, it is a propagable quadratic.
 *
 * It also implies that
 *  - expr^2 < expr * other_expr
 *  - other_expr * expr < expr^2
 *
 * It also implies that x^-2 < x^-1, but since, so far, we do not interpret x^-2 as (x^-1)^2, it is not a problem.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectQuadratic)
{  /*lint --e{715,774}*/
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlexprdata;
   SCIP_CONSEXPR_QUADEXPR* quaddata;
   SCIP_Bool propagable;
   SCIP_Bool nlhdlrconvexdoesquadratic;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);
   assert(nlhdlrexprdata != NULL);

   /* don't check if all enforcement methods are already ensured */
   if( (*enforcing & SCIP_CONSEXPR_EXPRENFO_ALL) == SCIP_CONSEXPR_EXPRENFO_ALL )
      return SCIP_OKAY;

   /* if it is not a sum of at least two terms, it is not interesting */
   /* TODO: constraints of the form l<= x*y <= r ? */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(expr) < 2 )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "Nlhdlr quadratic detecting expr %p aka ", (void*)expr);
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "Have to enforce %d\n", *enforcing);
#endif

   /* check whether expression is quadratic (a sum with at least one square or bilinear term) */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, conshdlr, expr, &quaddata) );

   /* not quadratic -> nothing for us */
   if( quaddata == NULL )
   {
      SCIPdebugMsg(scip, "expr %p is not quadratic -> abort detect\n", (void*)expr);
      return SCIP_OKAY;
   }

   propagable = isPropagable(quaddata);

   /* if we are not propagable and are in presolving, return */
   if( !propagable && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIPdebugMsg(scip, "expr %p is not propagable and in presolving -> abort detect\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* XXX: right now we only generate cut for convex quadratics; however, if a convex quadratic is not propagable, we do
    * not want to handle it; in other words, if it is not propagable and convex -> we don't want to do anything, and if
    * it is not propaable and nonconvex -> we can't do anything; thus, there is no point of storing the quadratic and we
    * just return here. However, if the handler is extended later to do something with nonconvex quadratics one should
    * delete this
    */
   if( !propagable )
   {
      SCIPdebugMsg(scip, "expr %p is not propagable -> abort detect\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* store quadratic in nlhdlrexprdata */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   nlexprdata = *nlhdlrexprdata;
   nlexprdata->quaddata = quaddata;

#ifdef DEBUG_DETECT
   SCIPinfoMessage(scip, NULL, "Nlhdlr quadratic detected:\n");
   SCIP_CALL( SCIPprintConsExprQuadratic(scip, conshdlr, quaddata) );
#endif

   /* every propagable quadratic expression will be handled since we can propagate */
   if( propagable )
   {
      SCIP_CONSEXPR_EXPR** linexprs;
      int nlinexprs;
      int nquadexprs;
      int nbilin;
      int i;

      *participating |= SCIP_CONSEXPR_EXPRENFO_ACTIVITY;
      *enforcing |= SCIP_CONSEXPR_EXPRENFO_ACTIVITY;

      SCIPgetConsExprQuadraticData(quaddata, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, &nbilin);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlexprdata->quadactivities, nquadexprs) );

      /* notify children of quadratic that we will need their activity for propagation */
      for( i = 0; i < nlinexprs; ++i )
         SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, linexprs[i], FALSE, TRUE, FALSE, FALSE) );

      for( i = 0; i < nquadexprs; ++i )
      {
         SCIP_CONSEXPR_EXPR* argexpr;
         if( isPropagableTerm(quaddata, i) )
         {
            SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &argexpr, NULL, NULL, &nbilin, NULL, NULL);
            SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, argexpr, FALSE, TRUE, FALSE, FALSE) );

#ifdef DEBUG_DETECT
            SCIPinfoMessage(scip, NULL, "quadterm %d propagable, using %p, unbounded=%d\n", i, (void*)argexpr, nbilin > 0 && SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, SCIPgetConsExprExprActivity(scip, argexpr)));
#endif
         }
         else
         {
            /* non-propagable quadratic is either a single square term or a single bilinear term
             * we should make use nlhdlrs in pow or product for this term, so we register usage of the square or product expr instead of argexpr
             */
            SCIP_CONSEXPR_EXPR* sqrexpr;
            int* adjbilin;

            SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &argexpr, NULL, NULL, &nbilin, &adjbilin, &sqrexpr);

            if( sqrexpr != NULL )
            {
               SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, sqrexpr, FALSE, TRUE, FALSE, FALSE) );
               assert(nbilin == 0);

#ifdef DEBUG_DETECT
               SCIPinfoMessage(scip, NULL, "quadterm %d non-propagable square, using %p\n", i, (void*)sqrexpr);
#endif
            }
            else
            {
               /* we have expr1 * other_expr or other_expr * expr1; know that expr1 is non propagable, but to decide if
                * we want the bounds of expr1 or of the product expr1 * other_expr (or other_expr * expr1), we have to
                * decide whether other_expr is also non propagable; due to the way we sort bilinear terms (by
                * frequency), we can deduce that other_expr doesn't appear anywhere else (i.e. is non propagable) if the
                * product is of the form expr1 * other_expr; however, if we see other_expr * expr1 we need to find
                * other_expr and check whether it is propagable
                */
               SCIP_CONSEXPR_EXPR* expr1;
               SCIP_CONSEXPR_EXPR* prodexpr;

               assert(nbilin == 1);
               SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[0], &expr1, NULL, NULL, NULL, &prodexpr);

               if( expr1 == argexpr )
               {
                  SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, prodexpr, FALSE, TRUE, FALSE, FALSE) );

#ifdef DEBUG_DETECT
                  SCIPinfoMessage(scip, NULL, "quadterm %d non-propagable product, using %p\n", i, (void*)prodexpr);
#endif
               }
               else
               {
                  int j;
                  /* check if other_expr is propagable in which case we need the bounds of expr1; otherwise we just need
                   * the bounds of the product and this will be (or was) registered when the loop takes us to the
                   * quadexpr other_expr.
                   * TODO this should be done faster, maybe store pos1 in bilinexprterm or store quadexprterm's in bilinexprterm
                   */
                  for( j = 0; j < nquadexprs; ++j )
                  {
                     SCIP_CONSEXPR_EXPR* exprj;
                     SCIPgetConsExprQuadraticQuadTermData(quaddata, j, &exprj, NULL, NULL, NULL, NULL, NULL);
                     if( expr1 == exprj )
                     {
                        if( isPropagableTerm(quaddata, j) )
                        {
                           SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, argexpr, FALSE, TRUE, FALSE, FALSE) );
#ifdef DEBUG_DETECT
                           SCIPinfoMessage(scip, NULL, "quadterm %d non-propagable alien product, using %p\n", i, (void*)argexpr);
#endif
                        }
                        break;
                     }
                  }
               }
            }
         }
      }
   }

   /* check if we are going to separate or not */
   nlexprdata->curvature = SCIP_EXPRCURV_UNKNOWN;

   /* for now, we do not care about separation if it is not required */
   if( (*enforcing & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == SCIP_CONSEXPR_EXPRENFO_SEPABOTH )
   {
      SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate\n", (void*)expr);
      return SCIP_OKAY;
   }

   assert(SCIPgetStage(scip) >= SCIP_STAGE_INITSOLVE);  /* separation should only be required in (init)solving stage */

   /* if nlhdlr_convex handles convex quadratic, then we don't (and vice versa)
    * TODO nlhdlr_convex is supposed to take this over permanently, but for now I keep both possibilities
    */
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/expr/nlhdlr/convex/cvxquadratic", &nlhdlrconvexdoesquadratic) );
   if( nlhdlrconvexdoesquadratic )
   {
      SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate, nlhldr_convex will separate\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* check if we can do something more: check curvature of quadratic function stored in nlexprdata
    * this is currently only used to decide whether we want to separate, so it can be skipped if in presolve
    */
   SCIPdebugMsg(scip, "checking convexity of expr %p\n", (void*)expr);
   SCIP_CALL( SCIPgetConsExprQuadraticCurvature(scip, quaddata, &nlexprdata->curvature, NULL) );

   /* if convex and non propagable -> quadratic is of the form sum x_i^2 + sum y_i -> we do not want to handle it
    * TODO: we can check if it is of this form before we allocate all the data
    */
   if( (nlexprdata->curvature == SCIP_EXPRCURV_CONVEX || nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE) && ! propagable ) /*lint !e845*/
   {
      assert(*participating == SCIP_CONSEXPR_EXPRENFO_NONE);
      SCIP_CALL( nlhdlrfreeExprDataQuadratic(scip, nlhdlr, expr, nlhdlrexprdata) );
      return SCIP_OKAY;
   }

   if( nlexprdata->curvature == SCIP_EXPRCURV_CONVEX )
   {
      SCIPdebugMsg(scip, "expr %p is convex when replacing factors of bilinear terms, bases of squares and every other term by their aux vars\n",
            (void*)expr);

      /* we will estimate the expression from below, that is, handle expr <= auxvar */
      *participating |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
      *enforcing |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   }
   else if( nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE )
   {
      SCIPdebugMsg(scip, "expr %p is concave when replacing factors of bilinear terms, bases of squares and every other term by their aux vars\n",
            (void*)expr);

      /* we will estimate the expression from above, that is, handle expr >= auxvar */
      *participating |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
      *enforcing |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   }

   /* we only need auxiliary variables if we are going to separate */
   if( *participating & SCIP_CONSEXPR_EXPRENFO_SEPABOTH )
   {
      SCIP_CONSEXPR_EXPR** linexprs;
      int nquadexprs;
      int nlinexprs;
      int i;

      SCIPgetConsExprQuadraticData(quaddata, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL);

      for( i = 0; i < nlinexprs; ++i ) /* expressions appearing linearly */
      {
         SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, linexprs[i], TRUE, FALSE, FALSE, FALSE) );
      }
      for( i = 0; i < nquadexprs; ++i ) /* expressions appearing quadratically */
      {
         SCIP_CONSEXPR_EXPR* quadexpr;
         SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &quadexpr, NULL, NULL, NULL, NULL, NULL);
         SCIP_CALL( SCIPregisterConsExprExprUsage(scip, conshdlr, quadexpr, TRUE, FALSE, FALSE, FALSE) );
      }

      SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate and separate\n", (void*)expr);
   }
   else
   {
      SCIPdebugMsg(scip, "expr %p is quadratic and propagable -> propagate only\n", (void*)expr);
   }

   if( SCIPareConsExprQuadraticExprsVariables(quaddata) )
   {
      SCIPsetConsExprExprCurvature(expr, nlexprdata->curvature);
      SCIPdebugMsg(scip, "expr is %s in the original variables\n", nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE ? "concave" : "convex");
   }

   return SCIP_OKAY;
}

/** nonlinear handler auxiliary evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxQuadratic)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(auxvalue != NULL);

   /* this handler can also handle quadratic expressions whose curvature is unknown or indefinite, since it can
    * propagate them, but it does not separate these
    * we then cannot evaluate w.r.t. auxvars, so we return the value of the expression instead
    */
   if( nlhdlrexprdata->curvature == SCIP_EXPRCURV_UNKNOWN )
   {
      *auxvalue = SCIPgetConsExprExprValue(expr);
      return SCIP_OKAY;
   }

   *auxvalue = SCIPevalConsExprQuadraticAux(scip, nlhdlrexprdata->quaddata, sol);

   return SCIP_OKAY;
}

/** nonlinear handler estimation callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateQuadratic)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_QUADEXPR* quaddata;
   SCIP_CONSEXPR_EXPR** linexprs;
   SCIP_Real* lincoefs;
   int nbilinexprterms;
   int nquadexprs;
   int nlinexprs;
   int j;

   SCIP_Real constant;
   SCIP_Real coef;
   SCIP_Real coef2;
   SCIP_ROWPREP* rowprep;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowpreps != NULL);
   assert(success != NULL);

   /* this handler can also handle quadratic expressions whose curvature is unknown or indefinite, since it can
    * propagate them, but it does not separate these; however, we should not be called to estimate on indefinite quadratics
    */
   assert(nlhdlrexprdata->curvature != SCIP_EXPRCURV_UNKNOWN);
   assert(!overestimate || nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE);
   assert( overestimate || nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX);

#ifndef NDEBUG
   /* check that quaddata hasn't changed (or at least the pointer to it) */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, conshdlr, expr, &quaddata) );
   assert(quaddata == nlhdlrexprdata->quaddata);
#endif
   quaddata = nlhdlrexprdata->quaddata;

   *success = FALSE;
   *addedbranchscores = FALSE;

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );

   SCIPgetConsExprQuadraticData(quaddata, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, &nbilinexprterms);

   /*
    * compute estimator: quadfun(sol) + \nabla quadfun(sol) (x - sol)
    */

   /* constant */
   SCIPaddRowprepConstant(rowprep, constant);

   /* handle purely linear variables */
   for( j = 0; j < nlinexprs; ++j )
   {
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprAuxVar(linexprs[j]), lincoefs[j]) );
   }

   /* quadratic variables */
   *success = TRUE;
   for( j = 0; j < nquadexprs; ++j )
   {
      SCIP_VAR* var;
      SCIP_CONSEXPR_EXPR* qexpr;
      SCIP_Real sqrcoef;
      SCIP_Real lincoef;
      int nadjbilin;

      SCIPgetConsExprQuadraticQuadTermData(quaddata, j, &qexpr, &lincoef, &sqrcoef, &nadjbilin, NULL, NULL);
      assert(qexpr != NULL);

      var = SCIPgetConsExprExprAuxVar(qexpr);
      assert(var != NULL);

      /* initialize coefficients to linear coefficients of quadratic variables */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, lincoef) );

      /* add linearization of square term */
      coef = 0.0;
      constant = 0.0;
      SCIPaddSquareLinearization(scip, sqrcoef, SCIPgetSolVal(scip, sol, var),
         nadjbilin == 0 && SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS, &coef, &constant, success);
      if( !*success )
      {
         SCIPfreeRowprep(scip, &rowprep);
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef) );
      SCIPaddRowprepConstant(rowprep, constant);
   }

   /* add linearization of bilinear terms */
   for( j = 0; j < nbilinexprterms; ++j )
   {
      SCIP_CONSEXPR_EXPR* qexpr1;
      SCIP_CONSEXPR_EXPR* qexpr2;
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_Real qcoef;

      SCIPgetConsExprQuadraticBilinTermData(quaddata, j, &qexpr1, &qexpr2, &qcoef, NULL, NULL);

      var1 = SCIPgetConsExprExprAuxVar(qexpr1);
      var2 = SCIPgetConsExprExprAuxVar(qexpr2);
      assert(var1 != NULL);
      assert(var2 != NULL);

      coef = 0.0;
      coef2 = 0.0;
      constant = 0.0;
      SCIPaddBilinLinearization(scip, qcoef, SCIPgetSolVal(scip, sol, var1), SCIPgetSolVal(scip, sol,
         var2), &coef, &coef2, &constant, success);
      if( !*success )
      {
         SCIPfreeRowprep(scip, &rowprep);
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var1, coef) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var2, coef2) );
      SCIPaddRowprepConstant(rowprep, constant);
   }

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   rowprep->local = FALSE;

   SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%sestimate_quadratic%p_%s%d",
      overestimate ? "over" : "under",
      (void*)expr,
      sol != NULL ? "sol" : "lp",
      sol != NULL ? SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** nonlinear handler forward propagation callback
 *
 * This method should solve the problem
 * max/min quad expression over box constraints
 * However, this problem is difficult so we are satisfied with a proxy.
 * Interval arithmetic suffices when no variable appears twice, however this is seldom the case, so we try
 * to take care of the dependency problem to some extent:
 * Let P_l = \{i : expr_l expr_i is a bilinear expr\}.
 * 1. partition the quadratic expression as sum of quadratic functions
 * \sum_l q_l
 * where q_l = a_l expr_l^2 + c_l expr_l + \sum_{i \in P_l} b_il expr_i expr_l
 * 2. build interval quadratic functions, i.e, a x^2 + b x where b is an interval, i.e.,
 * a_l expr_l^2 + [\sum_{i \in P_l} b_il expr_i + c_l] expr_l
 * 3. compute \min and \max { a x^2 + b x : x \in [x] } for each interval quadratic, i.e.,
 * \min and \max a_l expr_l^2 + expr_l [\sum_{i \in P_l} b_il expr_i + c_l] : expr_l \in [expr_l]
 *
 * Notes:
 * 1. The l-th quadratic expr (expressions that appear quadratically) is associated with q_l
 * 2. nlhdlrdata->quadactivities[l] is the activity of q_l as computed in the description above.
 * 3. The q_l of a quadratic term might be empty, in which case nlhdlrdata->quadactivities[l] is [0,0].
 * For example, consider x^2 + x*y. There are two quadratic expressions, 'x' and 'y'.
 * The q associated to 'x' is 'x^2 + x*y', while the q associated to 'y' is empty.
 * Thus, nlhdlrdata->quadactivities[1] is [0,0] in this case.
 * The logic is to avoid considering the term 'x*y' twice.
 *
 * @note The order matters! If expr_i * expr_l is a term in of the quadratic, then i is *not* in P_l
 */
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalQuadratic)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_QUADEXPR* quaddata;
   SCIP_CONSEXPR_EXPR** linexprs;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   int nquadexprs;
   int nlinexprs;

   assert(scip != NULL);
   assert(expr != NULL);

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->quadactivities != NULL);

   SCIPdebugMsg(scip, "Interval evaluation of quadratic expr\n");

#ifndef NDEBUG
   /* check that quaddata hasn't changed (or at least the pointer to it) */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, SCIPfindConshdlr(scip, "expr"), expr, &quaddata) );
   assert(quaddata == nlhdlrexprdata->quaddata);
#endif
   quaddata = nlhdlrexprdata->quaddata;

   SCIPgetConsExprQuadraticData(quaddata, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, NULL);

   /*
    * compute activity of linear part, if some linear term has changed
    */
   {
      int i;

      SCIPdebugMsg(scip, "Computing activity of linear part\n");

      SCIPintervalSet(&nlhdlrexprdata->linactivity, constant);
      for( i = 0; i < nlinexprs; ++i )
      {
         SCIP_INTERVAL linterminterval;

         linterminterval = SCIPgetConsExprExprActivity(scip, linexprs[i]);
         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, linterminterval) )
         {
            SCIPdebugMsg(scip, "Activity of linear part is empty due to child %d\n", i);
            SCIPintervalSetEmpty(interval);
            return SCIP_OKAY;
         }
         SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &linterminterval, linterminterval, lincoefs[i]);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &nlhdlrexprdata->linactivity, nlhdlrexprdata->linactivity, linterminterval);
      }

      SCIPdebugMsg(scip, "Activity of linear part is [%g, %g]\n", nlhdlrexprdata->linactivity.inf,
            nlhdlrexprdata->linactivity.sup);
   }

   /*
    * compute activity of quadratic part
    */
   {
      int i;

      SCIPdebugMsg(scip, "Computing activity of quadratic part\n");

      nlhdlrexprdata->nneginfinityquadact = 0;
      nlhdlrexprdata->nposinfinityquadact = 0;
      nlhdlrexprdata->minquadfiniteact = 0.0;
      nlhdlrexprdata->maxquadfiniteact = 0.0;
      SCIPintervalSet(&nlhdlrexprdata->quadactivity, 0.0);

      for( i = 0; i < nquadexprs; ++i )
      {
         SCIP_Real quadlb;
         SCIP_Real quadub;
         SCIP_CONSEXPR_EXPR* qexpr;
         SCIP_Real lincoef;
         SCIP_Real sqrcoef;
         int nadjbilin;
         int* adjbilin;
         SCIP_CONSEXPR_EXPR* sqrexpr;

         SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin, &sqrexpr);

         if( !isPropagableTerm(quaddata, i) )
         {
            /* term is not propagable, i.e., the exprs involved in term only appear once; thus use the activity of the
             * quadratic term directly and not the activity of the exprs involed in the term. See also documentation of
             * DETECT
             */
            SCIP_INTERVAL tmp;

            assert(lincoef == 0.0);

            if( sqrcoef != 0.0 )
            {
               assert(sqrexpr != NULL);
               assert(nadjbilin == 0);

               tmp = SCIPgetConsExprExprActivity(scip, sqrexpr);
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, tmp) )
               {
                  SCIPintervalSetEmpty(interval);
                  return SCIP_OKAY;
               }

               SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &tmp, tmp, sqrcoef);
               quadlb = tmp.inf;
               quadub = tmp.sup;

#ifdef DEBUG_PROP
               SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term %g <expr>, where <expr> is: ", sqrcoef);
               SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), sqrexpr, NULL) );
#endif
            }
            else
            {
               SCIP_CONSEXPR_EXPR* expr1;
               SCIP_CONSEXPR_EXPR* prodexpr;
               SCIP_Real prodcoef;

               assert(nadjbilin == 1);
               SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[0], &expr1, NULL, &prodcoef, NULL, &prodexpr);

               if( expr1 == qexpr )
               {
                  /* the quadratic expression expr1 appears only as expr1 * expr2, so its 'q' is expr1 * expr2 */
                  tmp = SCIPgetConsExprExprActivity(scip, prodexpr);
                  if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, tmp) )
                  {
                     SCIPintervalSetEmpty(interval);
                     return SCIP_OKAY;
                  }

                  SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &tmp, tmp, prodcoef);
                  quadlb = tmp.inf;
                  quadub = tmp.sup;

#ifdef DEBUG_PROP
                  SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term %g <expr>, where <expr> is: ", prodcoef);
                  SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), prodexpr, NULL) );
#endif
               }
               else
               {
                  /* the quadratic expression expr1 appears as expr2 * expr1, thus its 'q' is empty, see also the Notes
                   * in the documentation of the function
                   */
                  SCIPintervalSet(&nlhdlrexprdata->quadactivities[i], 0.0);
                  continue;
               }
            }
         }
         else
         {
            int j;
            SCIP_INTERVAL b;

            SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin, NULL);

            if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPgetConsExprExprActivity(scip, qexpr)) )
            {
               SCIPintervalSetEmpty(interval);
               return SCIP_OKAY;
            }

            /* b = [c_l] */
            SCIPintervalSet(&b, lincoef);
#ifdef DEBUG_PROP
            SCIPinfoMessage(scip, NULL, "b := %g\n", lincoef);
#endif
            for( j = 0; j < nadjbilin; ++j )
            {
               SCIP_INTERVAL bterm;
               SCIP_CONSEXPR_EXPR* expr1;
               SCIP_CONSEXPR_EXPR* expr2;
               SCIP_Real bilincoef;

               SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[j], &expr1, &expr2, &bilincoef, NULL, NULL);

               if( expr1 != qexpr )
                  continue;

               bterm = SCIPgetConsExprExprActivity(scip, expr2);
               if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, bterm) )
               {
                  SCIPintervalSetEmpty(interval);
                  return SCIP_OKAY;
               }

               /* b += [b_jl * expr_j] for j \in P_l */
               SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, bterm, bilincoef);
               SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);

#ifdef DEBUG_PROP
               SCIPinfoMessage(scip, NULL, "b += %g * [expr2], where <expr2> is: ", bilincoef);
               SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), expr2, NULL) );
               SCIPinfoMessage(scip, NULL, " [%g,%g]\n", SCIPgetConsExprExprActivity(scip, expr2).inf, SCIPgetConsExprExprActivity(scip, expr2).sup);
#endif
            }

            /* TODO: under which assumptions do we know that we just need to compute min or max? its probably the locks that give some information here */
            quadub = SCIPintervalQuadUpperBound(SCIP_INTERVAL_INFINITY, sqrcoef, b,
               SCIPgetConsExprExprActivity(scip, qexpr));

            /* TODO: implement SCIPintervalQuadLowerBound */
            {
               SCIP_INTERVAL minusb;
               SCIPintervalSetBounds(&minusb, -SCIPintervalGetSup(b), -SCIPintervalGetInf(b));

               quadlb = -SCIPintervalQuadUpperBound(SCIP_INTERVAL_INFINITY, -sqrcoef, minusb,
                  SCIPgetConsExprExprActivity(scip, qexpr));
            }

#ifdef DEBUG_PROP
            SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term %g <expr>^2 + [%g,%g] <expr>, where <expr> is: ", sqrcoef, b.inf, b.sup);
            SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), qexpr, NULL) );
#endif
         }
#ifdef DEBUG_PROP
         SCIPinfoMessage(scip, NULL, " -> [%g, %g]\n", quadlb, quadub);
#endif

         SCIPintervalSetBounds(&nlhdlrexprdata->quadactivities[i], quadlb, quadub);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &nlhdlrexprdata->quadactivity, nlhdlrexprdata->quadactivity, nlhdlrexprdata->quadactivities[i]);

         /* get number of +/-infinity contributions and compute finite activity */
         if( quadlb <= -SCIP_INTERVAL_INFINITY )
            nlhdlrexprdata->nneginfinityquadact++;
         else
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeDownwards();

            nlhdlrexprdata->minquadfiniteact += quadlb;

            SCIPintervalSetRoundingMode(roundmode);
         }
         if( quadub >= SCIP_INTERVAL_INFINITY )
            nlhdlrexprdata->nposinfinityquadact++;
         else
         {
            SCIP_ROUNDMODE roundmode;

            roundmode = SCIPintervalGetRoundingMode();
            SCIPintervalSetRoundingModeUpwards();

            nlhdlrexprdata->maxquadfiniteact += quadub;

            SCIPintervalSetRoundingMode(roundmode);
         }
      }

      SCIPdebugMsg(scip, "Activity of quadratic part is [%g, %g]\n", nlhdlrexprdata->quadactivity.inf, nlhdlrexprdata->quadactivity.sup);
   }

   /* interval evaluation is linear activity + quadactivity */
   SCIPintervalAdd(SCIP_INTERVAL_INFINITY, interval, nlhdlrexprdata->linactivity,  nlhdlrexprdata->quadactivity);

   nlhdlrexprdata->activitiestag = SCIPgetConsExprCurBoundsTag(SCIPfindConshdlr(scip, "expr"));

   return SCIP_OKAY;
}

/** returns max of a/x - c*x for x in {x1, x2} with x1, x2 > 0 */
static
SCIP_Real computeMaxBoundaryForBilinearProp(
   SCIP_Real a,
   SCIP_Real c,
   SCIP_Real x1,
   SCIP_Real x2
   )
{
   SCIP_Real cneg;
   SCIP_Real cand1;
   SCIP_Real cand2;
   SCIP_ROUNDMODE roundmode;

   assert(x1 > 0.0);
   assert(x2 > 0.0);

   cneg = SCIPintervalNegateReal(c);

   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();
   cand1 = a/x1 + cneg*x1;
   cand2 = a/x2 + cneg*x2;
   SCIPintervalSetRoundingMode(roundmode);

   return MAX(cand1, cand2);
}

/** returns max of a/x - c*x for x in dom; it assumes that dom is contained in (0, +inf) */
static
SCIP_Real computeMaxForBilinearProp(
   SCIP_Real a,
   SCIP_Real c,
   SCIP_INTERVAL dom
   )
{
   SCIP_ROUNDMODE roundmode;
   SCIP_INTERVAL argmax;
   SCIP_Real negunresmax;
   SCIP_Real boundarymax;
   assert(dom.inf > 0);

   /* if a >= 0, then the function is convex which means the maximum is at one of the boundaries
    *
    * if c = 0, then the function is monotone which means the maximum is also at one of the boundaries
    *
    * if a < 0, then the function is concave. The function then has a maximum if and only if there is a point with derivative 0,
    * that is, iff -a/x^2 - c = 0 has a solution; i.e. if -a/c >= 0, i.e. (using a<0 and c != 0), c > 0.
    * Otherwise (that is, c<0), the maximum is at one of the boundaries.
    */
   if( a >= 0.0 || c <= 0.0 )
      return computeMaxBoundaryForBilinearProp(a, c, dom.inf, dom.sup);

   /* now, the (unrestricted) maximum is at sqrt(-a/c).
    * if the argmax is not in the interior of dom then the solution is at a boundary, too
    * we check this by computing an interval that contains sqrt(-a/c) first
    */
   SCIPintervalSet(&argmax, -a);
   SCIPintervalDivScalar(SCIP_INTERVAL_INFINITY, &argmax, argmax, c);
   SCIPintervalSquareRoot(SCIP_INTERVAL_INFINITY, &argmax, argmax);

   /* if the interval containing sqrt(-a/c) does not intersect with the interior of dom, then
    * the (restricted) maximum is at a boundary (we could even say at which boundary, but that doesn't save much)
    */
   if( argmax.sup <= dom.inf || argmax.inf >= dom.sup )
      return computeMaxBoundaryForBilinearProp(a, c, dom.inf, dom.sup);

   /* the maximum at sqrt(-a/c) is -2*sqrt(-a*c), so we compute an upper bound for that by computing a lower bound for 2*sqrt(-a*c) */
   roundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();
   negunresmax = 2.0*SCIPnextafter(sqrt(SCIPintervalNegateReal(a)*c), 0.0);
   SCIPintervalSetRoundingMode(roundmode);

   /* if the interval containing sqrt(-a/c) is contained in dom, then we can return -negunresmax */
   if( argmax.inf >= dom.inf && argmax.sup <= dom.sup )
      return -negunresmax;

   /* now what is left is the case where we cannot say for sure whether sqrt(-a/c) is contained in dom or not
    * so we are conservative and return the max of both cases, i.e.,
    * the max of the upper bounds on -2*sqrt(-a*c), a/dom.inf-c*dom.inf, a/dom.sup-c*dom.sup.
    */
   boundarymax = computeMaxBoundaryForBilinearProp(a, c, dom.inf, dom.sup);
   return MAX(boundarymax, -negunresmax);
}

/** computes the range of rhs/x - coef * x for x in exprdom; this is used for the propagation of bilinear terms
 *
 * If 0 is in the exprdom, we set range to R (even though this is not quite correct, it is correct for the
 * intended use of the function).
 * TODO: maybe check before calling it whether 0 is in the domain and then just avoid calling it
 *
 * If rhs is [A,B] and x > 0, then we want the min of A/x - coef*x and max of B/x - coef * x for x in [exprdom].
 * If rhs is [A,B] and x < 0, then we want the min of B/x - coef*x and max of A/x - coef * x for x in [exprdom].
 * However, this is the same as min of -B/x + coef*x and max of -A/x + coef * x for x in -[exprdom].
 * Thus, we can reduce to x > 0 always by multiplying [exprdom], rhs, and coef by -1.
 */
static
void computeRangeForBilinearProp(
   SCIP_INTERVAL         exprdom,            /**< expression for which to solve */
   SCIP_Real             coef,               /**< expression for which to solve */
   SCIP_INTERVAL         rhs,
   SCIP_INTERVAL*        range
   )
{
   SCIP_Real max;
   SCIP_Real min;

   if( exprdom.inf <= 0.0 && 0.0 <= exprdom.sup )
   {
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, range);
      return;
   }

   /* reduce to positive case */
   if( exprdom.sup < 0 )
   {
      SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &exprdom, exprdom, -1.0);
      SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &rhs, rhs, -1.0);
      coef *= -1.0;
   }
   assert(exprdom.inf > 0.0);

   /* compute maximum and minimum */
   max = computeMaxForBilinearProp(rhs.sup, coef, exprdom);
   min = -computeMaxForBilinearProp(-rhs.inf, -coef, exprdom);

   /* set interval */
   SCIPintervalSetBounds(range, min, max);
}

/** nonlinear handler reverse propagation callback
 * @note: the implemented technique is a proxy for solving the OBBT problem min/max{ x_i : quad expr in [quad expr] }
 * and as such can be improved.
 *
 * input:
 *  - scip : SCIP main data structure
 *  - nlhdlr : nonlinear handler
 *  - expr : expression
 *  - nlhdlrexprdata : expression specific data of the nonlinear handler
 *  - reversepropqueue : expression queue in reverse propagation, to be passed on to SCIPtightenConsExprExprInterval
 *  - infeasible: buffer to store whether an expression's bounds were propagated to an empty interval
 *  - nreductions : buffer to store the number of interval reductions of all children
 *  - force : force tightening even if it is below the bound strengthening tolerance
 */
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropQuadratic)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_QUADEXPR* quaddata;
   SCIP_CONSEXPR_EXPR** linexprs;
   SCIP_CONSEXPR_EXPR** bilinexprs; /* TODO: should this be stored in the nlhdlr expr data? */
   SCIP_Real* bilincoefs;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   int nquadexprs;
   int nlinexprs;

   SCIP_INTERVAL rhs;
   SCIP_INTERVAL quadactivity;
   int i;

   SCIPdebugMsg(scip, "Reverse propagation of quadratic expr\n");

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);
   assert(nreductions != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->quadactivities != NULL);
#ifndef NDEBUG
   /* check that quaddata hasn't changed (or at least the pointer to it) */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, SCIPfindConshdlr(scip, "expr"), expr, &quaddata) );
   assert(quaddata == nlhdlrexprdata->quaddata);
#endif
   quaddata = nlhdlrexprdata->quaddata;

   *nreductions = 0;

   /* not possible to conclude finite bounds if the interval of the expression is [-inf,inf] */
   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, SCIPgetConsExprExprActivity(scip, expr)) )
   {
      SCIPdebugMsg(scip, "expr's range is R -> cannot reverse propagate\n");
      return SCIP_OKAY;
   }

   /* ensure that partial activities as stored in nlhdlrexprdata are uptodate
    * if the activity stored in expr is more recent than the partial activities stored in this nlhdlrexprdata,
    * then we should reevaluate the partial activities
    */
   if( SCIPgetConsExprExprActivityTag(expr) > nlhdlrexprdata->activitiestag )
   {
      SCIP_CALL( nlhdlrIntevalQuadratic(scip, nlhdlr, expr, nlhdlrexprdata, &quadactivity, NULL, NULL) );
   }

   SCIPgetConsExprQuadraticData(quaddata, &constant, &nlinexprs, &linexprs, &lincoefs, &nquadexprs, NULL);

   /* propagate linear part in rhs = expr's interval - quadratic activity; first, reconstruct the quadratic activity */
   SCIPintervalSetBounds(&quadactivity,
         nlhdlrexprdata->nneginfinityquadact > 0 ? -SCIP_INTERVAL_INFINITY : nlhdlrexprdata->minquadfiniteact,
         nlhdlrexprdata->nposinfinityquadact > 0 ?  SCIP_INTERVAL_INFINITY : nlhdlrexprdata->maxquadfiniteact);

   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, SCIPgetConsExprExprActivity(scip, expr), quadactivity);
   SCIP_CALL( SCIPreverseConsExprExprPropagateWeightedSum(scip, conshdlr, nlinexprs,
            linexprs, lincoefs, constant, rhs, reversepropqueue, infeasible, nreductions, force) );

   /* stop if we find infeasibility */
   if( *infeasible )
      return SCIP_OKAY;

   /* propagate quadratic part in expr's interval - linear activity, where linear activity was computed in INTEVAL.
    * The idea is basically to write interval quadratics for each expr and then solve for expr.
    *
    * One way of achieving this is:
    * - for each expression expr_i, write the quadratic expression as a_i expr^2_i + expr_i ( \sum_{j \in J_i} b_ij
    *   expr_j + c_i ) + quadratic expression in expr_k for k \neq i
    * - compute the interval b = [\sum_{j \in J_i} b_ij expr_j + c_i], where J_i are all the indices j such that the
    *   bilinear expression expr_i expr_j appears
    * - use some technique (like the one in nlhdlrIntevalQuadratic), to evaluate the activity of rest_i = [quadratic
    *   expression in expr_k for k \neq i].
    * - solve a_i expr_i^2 + b expr_i \in rhs_i := [expr activity] - rest_i
    *
    * However, this might be expensive, especially computing rest_i. Hence, we implement a simpler version.
    * - we use the same partition as in nlhdlrIntevalQuadratic for the bilinear terms. This way, b = [\sum_{j \in P_i}
    *   b_ij expr_j + c_i], where P_i is the set of indices j such that expr_i * expr_j appears in that order
    * - we evaluate the activity of rest_i as sum_{k \neq i} [\min q_k, \max q_k] where q_k = a_k expr_k^2 + [\sum_{j
    *   \in P_k} b_jk expr_j + c_k] expr_k. The intervals [\min q_k, \max q_k] were already computed in
    *   nlhdlrIntevalQuadratic, so we just reuse them.
    *
    * A downside of the above is that we might not deduce any bounds for variables that appear less often. For example,
    * consider x^2 + x * y + x * z + y * z + z. This quadratic gets partitioned as (x^2 + x*y + x*z) + (z*y + z). The
    * first parenthesis is interpreted as a function of x, while the second one as a function of z.
    * To also get bounds on y, after reverse propagating x in x^2 + x*y + x*z \in rhs, we rewrite this as y + z \in rhs/x -
    * x and propagate the y + z).
    * In general, after reverse propagating expr_i, we consider
    *   \sum_{j \in J_i} b_ij expr_j in ([expr activity] - quadratic expression in expr_k for k \neq i - c_i) / expr_i - a_i expr_i,
    * compute an interval for the right hand side (see computeRangeForBilinearProp) and use that to propagate the
    * linear sum on the left hand side.
    *
    * Note: this last step generalizes a technique that appeared in the classic cons_quadratic.
    * The idea of that technique was to to borrow a bilinear term expr_k expr_l when propagating expr_l and the quadratic
    * function for expr_k was simple enough.
    * Since in P_l we only consider the indices of expressions that appear multiplying expr_l as _second_ factor, we
    * would lose the bilinear terms expr_k * expr_l, which contributes to the dependency problem.
    * The problem is that the contribution of b_kl * expr_k * expr_l to rest_i is not just [b_kl * expr_k * expr_l], but
    * rather quadactivities[k] (= max/min of a_k expr_k^2 + expr_k * [c_k + sum_i \in P_k b_ki expr_i]).
    * Thus, we _cannot_ just substract [b_kl * expr_k * expr_l] from rest_i.
    * But, if expr_k only appears as expr_k * expr_l, then  quadactivities[k] = [b_kl * expr_k * expr_l]. So this
    * case was handled in old cons_quadratic.
    *
    *
    * TODO: handle simple cases
    * TODO: identify early when there is nothing to be gain
    */
   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, SCIPgetConsExprExprActivity(scip, expr), nlhdlrexprdata->linactivity);
   SCIP_CALL( SCIPallocBufferArray(scip, &bilinexprs, nquadexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bilincoefs, nquadexprs) );

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_INTERVAL rhs_i;
      SCIP_INTERVAL rest_i;
      SCIP_CONSEXPR_EXPR* qexpr;
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      int nadjbilin;
      int* adjbilin;
      SCIP_CONSEXPR_EXPR* sqrexpr;

      SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin, &sqrexpr);

      /* rhs_i = rhs - rest_i.
       * to compute rest_i = [\sum_{k \neq i} q_k] we just have to substract
       * the activity of q_i from quadactivity; however, care must be taken about infinities;
       * if [q_i].sup = +infinity and there is = 1 contributing +infinity -> rest_i.sup = maxquadfiniteact
       * if [q_i].sup = +infinity and there is > 1 contributing +infinity -> rest_i.sup = +infinity
       * if [q_i].sup = finite and there is > 0 contributing +infinity -> rest_i.sup = +infinity
       * if [q_i].sup = finite and there is = 0 contributing +infinity -> rest_i.sup = maxquadfiniteact - [q_i].sup
       *
       * the same holds when replacing sup with inf, + with - and max(quadfiniteact) with min(...)
       */
      /* compute rest_i.sup */
      if( SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]) < SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nposinfinityquadact == 0 )
      {
         SCIP_ROUNDMODE roundmode;

         roundmode = SCIPintervalGetRoundingMode();
         SCIPintervalSetRoundingModeUpwards();
         rest_i.sup = nlhdlrexprdata->maxquadfiniteact - SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]);

         SCIPintervalSetRoundingMode(roundmode);
      }
      else if( SCIPintervalGetSup(nlhdlrexprdata->quadactivities[i]) >= SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nposinfinityquadact == 1 )
         rest_i.sup = nlhdlrexprdata->maxquadfiniteact;
      else
         rest_i.sup = SCIP_INTERVAL_INFINITY;

      /* compute rest_i.inf */
      if( SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]) > -SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nneginfinityquadact == 0 )
      {
         SCIP_ROUNDMODE roundmode;

         roundmode = SCIPintervalGetRoundingMode();
         SCIPintervalSetRoundingModeDownwards();
         rest_i.inf = nlhdlrexprdata->minquadfiniteact - SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]);

         SCIPintervalSetRoundingMode(roundmode);
      }
      else if( SCIPintervalGetInf(nlhdlrexprdata->quadactivities[i]) <= -SCIP_INTERVAL_INFINITY &&
         nlhdlrexprdata->nneginfinityquadact == 1 )
         rest_i.inf = nlhdlrexprdata->minquadfiniteact;
      else
         rest_i.inf = -SCIP_INTERVAL_INFINITY;

#if 0 /* I (SV) added the following in cons_quadratic to fix/workaround some bug. Maybe we'll need this here, too? */
      /* FIXME in theory, rest_i should not be empty here
       * what we tried to do here is to remove the contribution of the i'th bilinear term (=bilinterm) to [minquadactivity,maxquadactivity] from rhs
       * however, quadactivity is computed differently (as x*(a1*y1+...+an*yn)) than q_i (a*ak*yk) and since interval arithmetics do overestimation,
       * it can happen that q_i is actually slightly larger than quadactivity, which results in rest_i being (slightly) empty
       * a proper fix could be to compute the quadactivity also as x*a1*y1+...+x*an*yn if sqrcoef=0, but due to taking
       * also infinite bounds into account, this complicates the code even further
       * instead, I'll just work around this by turning an empty rest_i into a small non-empty one
       */
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rest_i) )
      {
         assert(SCIPisSumRelEQ(scip, rest_i.inf, rest_i.sup));
         SCIPswapReals(&rest_i.inf, &rest_i.sup);
      }
#endif
      assert(!SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, rest_i));

      /* compute rhs_i */
      SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs_i, rhs, rest_i);

      if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, rhs_i) )
         continue;


      /* try to propagate */
      if( !isPropagableTerm(quaddata, i) )
      {
         assert(lincoef == 0.0);

         if( sqrcoef != 0.0 )
         {
            assert(sqrexpr != NULL);
            assert(nadjbilin == 0);

            /* solve sqrcoef sqrexpr in rhs_i */
            SCIP_CALL( propagateBoundsLinExpr(scip, conshdlr, sqrexpr, sqrcoef, rhs_i, reversepropqueue, infeasible, nreductions, force) );
         }
         else
         {
            /* qexpr only appears in a term of the form qexpr * other_expr (or other_expr * qexpr); we only care about
             * getting bounds for the product, thus we will compute these bounds when qexpr appears as qexpr *
             * other_expr; note that if it appears as other_expr * qexpr, then when we process other_expr bounds for the
             * product will be computed
             * TODO: we can actually avoid computing rhs_i in the case that qexpr is not propagable and it appears as
             * other_expr * qexpr
             */
            SCIP_CONSEXPR_EXPR* expr1;
            SCIP_CONSEXPR_EXPR* prodexpr;
            SCIP_Real prodcoef;

            assert(nadjbilin == 1);
            SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[0], &expr1, NULL, &prodcoef, NULL, &prodexpr);

            if( expr1 == qexpr )
            {
               /* solve prodcoef prodexpr in rhs_i */
               SCIP_CALL( propagateBoundsLinExpr(scip, conshdlr, prodexpr, prodcoef, rhs_i, reversepropqueue, infeasible, nreductions, force) );
            }
         }
      }
      else
      {
         SCIP_INTERVAL b;
         SCIP_CONSEXPR_EXPR* expr1 = NULL;
         SCIP_CONSEXPR_EXPR* expr2 = NULL;
         SCIP_Real bilincoef = 0.0;
         int nbilin = 0;
         int pos2 = 0;
         int j;

         /* set b to [c_l] */
         SCIPintervalSet(&b, lincoef);

         /* add [\sum_{j \in P_l} b_lj expr_j + c_l] into b */
         for( j = 0; j < nadjbilin; ++j )
         {
            SCIP_INTERVAL bterm;

            SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[j], &expr1, &expr2, &bilincoef, &pos2, NULL);

            if( expr1 != qexpr )
               continue;

            /* b += [b_lj * expr_j] for j \in P_l */
            SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, SCIPgetConsExprExprActivity(scip, expr2), bilincoef);
            SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);

            /* remember b_lj and expr_j to propagate them too */
            bilinexprs[nbilin] = expr2;
            bilincoefs[nbilin] = bilincoef;
            nbilin++;
         }

         /* solve a_i expr_i^2 + b expr_i in rhs_i */
         SCIP_CALL( propagateBoundsQuadExpr(scip, conshdlr, qexpr, sqrcoef, b, rhs_i, reversepropqueue, infeasible, nreductions, force) );

         if( nbilin > 0 && !*infeasible )
         {
            /* if 0 is not in [expr_i], then propagate bilincoefs^T bilinexpr in rhs_i/expr_i - a_i expr_i - c_i */
            SCIP_INTERVAL bilinrhs;

            /* compute bilinrhs := [rhs_i/expr_i - a_i expr_i] */
            computeRangeForBilinearProp(SCIPgetConsExprExprActivity(scip, qexpr), sqrcoef, rhs_i, &bilinrhs);

            if( !SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, bilinrhs) )
            {
               int nreds;

               /* propagate \sum_{j \in P_i} b_ij expr_j + c_i in bilinrhs */
               SCIP_CALL( SCIPreverseConsExprExprPropagateWeightedSum(scip, conshdlr, nbilin,
                        bilinexprs, bilincoefs, lincoef, bilinrhs, reversepropqueue, infeasible, &nreds, force) );

               /* TODO FIXME: we are overestimating of the number of reductions: an expr might be tightened many times! */
               *nreductions += nreds;
            }
         }
      }

      /* stop if we find infeasibility */
      if( *infeasible )
         break;
   }

   SCIPfreeBufferArray(scip, &bilincoefs);
   SCIPfreeBufferArray(scip, &bilinexprs);

   return SCIP_OKAY;
}

/** nonlinear handler copy callback
 *
 * the method includes the nonlinear handler into a expression constraint handler
 *
 * This method is usually called when doing a copy of an expression constraint handler.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrcopyHdlrQuadratic)
{  /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrQuadratic(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes quadratic nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectQuadratic, nlhdlrEvalAuxQuadratic, NULL) );

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrcopyHdlrQuadratic);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataQuadratic);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateQuadratic, NULL);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalQuadratic, nlhdlrReversepropQuadratic);

   return SCIP_OKAY;
}
