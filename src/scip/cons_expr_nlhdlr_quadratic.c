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

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "quadratic"
#define NLHDLR_DESC         "handler for quadratic expressions"
#define NLHDLR_PRIORITY     100

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_CONSEXPR_QUADEXPR* quaddata;         /**< data of quadratic form as stored in expr (stored here again for convenient access) */

   SCIP_EXPRCURV         curvature;          /**< curvature of the quadratic representation of the expression */

   SCIP_INTERVAL         linactivity;        /**< activity of linear part */

   /* activities of quadratic parts as defined in nlhdlrIntervalQuadratic */
   SCIP_Real             minquadfiniteact;   /**< minimum activity of quadratic part where only terms with finite min
                                               activity contribute */
   SCIP_Real             maxquadfiniteact;   /**< maximum activity of quadratic part where only terms with finite max
                                               activity contribute */
   int                   nneginfinityquadact;/**< number of quadratic terms contributing -infinity to activity */
   int                   nposinfinityquadact;/**< number of quadratic terms contributing +infinity to activity */
   SCIP_INTERVAL*        quadactivities;     /**< activity of each quadratic term as defined in nlhdlrIntervalQuadratic */
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

      SCIPgetConsExprQuadraticQuadTermData(quaddata, i, NULL, &lincoef, &sqrcoef, &nadjbilin, NULL);

      if( (lincoef != 0.0) + (sqrcoef != 0.0) + nadjbilin >= 2 )  /*lint !e514*/ /* actually MIN(2, nadjbilin), but we check >= 2 */
         return TRUE;
   }

   return FALSE;
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
 * For propagation, we store the quadratic in our data structure in the following way: We count how often a variable
 * appears. Then, a bilinear product expr_i * expr_j is stored as expr_i * expr_j if # expr_i appears > # expr_j
 * appears. When # expr_i appears == # expr_j appears, it then it will be stored as expr_i * expr_j if and only if
 * expr_i < expr_j, where '<' is the expression order (see Ordering Rules in cons_expr.c documentation).
 * Heuristically, this should be useful for propagation. The intuition is that by factoring out the variable that
 * appears most often we should be able to take care of the dependency problem better.
 *
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
   SCIP_Bool separate = FALSE;
   SCIP_Bool nlhdlrconvexdoesquadratic;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* don't check if enforcement is already ensured */
   if( *enforcedbelow && *enforcedabove )
      return SCIP_OKAY;

   /* if it is not a sum of at least two terms, it is not interesting */
   /* TODO: constraints of the form l<= x*y <= r ? */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(expr) < 2 )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "Nlhdlr quadratic detecting expr %p aka", (void*)expr);
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "Have to enforce: Below? %s. Above? %s\n", *enforcedbelow ? "no" : "yes", *enforcedabove ? "no" : "yes");
#endif
   SCIPdebugMsg(scip, "checking if expr %p is a propagable quadratic\n", (void*)expr);

   /* check whether expression is quadratic (a sum with at least one square or bilinear term) */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, conshdlr, expr, &quaddata) );

   /* not quadratic -> nothing for us */
   if( quaddata == NULL )
      return SCIP_OKAY;

   propagable = isPropagable(quaddata);

   /* if we are not propagable and are in presolving, return */
   if( !propagable && SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIPdebugMsg(scip, "expr %p is not propagable and in presolving -> abort detect\n");
      return SCIP_OKAY;
   }

   /* XXX: right now we only generate cut for convex quadratics; however, if a convex quadratic is not propagable, we do
    * not want to handle it; in other words, if it is not propagable and convex -> we don't want to do anything, and if
    * it is not propaable and nonconvex -> we can't do anything; thus, there is no point of storing the quadratic and we
    * just return here. However, if the handler is extended later to do something with nonconvex quadratics one should
    * delete this
    */
   if( !propagable )
      return SCIP_OKAY;

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
      int nquadexprs;

      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_INTEVAL | SCIP_CONSEXPR_EXPRENFO_REVERSEPROP;

      SCIPgetConsExprQuadraticData(quaddata, NULL, NULL, NULL, NULL, &nquadexprs, NULL);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlexprdata->quadactivities, nquadexprs) );
   }

   /* check if we are going to separate or not */
   nlexprdata->curvature = SCIP_EXPRCURV_UNKNOWN;

   /* for now, we do not care about separation if we are not solving */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);

   /* if nlhdlr_convex handles convex quadratic, then we don't (and vice versa)
    * TODO nlhdlr_convex is supposed to take this over permanently, but for now I keep both possibilities
    */
   SCIP_CALL( SCIPgetBoolParam(scip, "constraints/expr/nlhdlr/convex/cvxquadratic", &nlhdlrconvexdoesquadratic) );
   if( nlhdlrconvexdoesquadratic )
      return SCIP_OKAY;

   /* check if we can do something more: check curvature of quadratic function stored in nlexprdata
    * this is currently only used to decide whether we want to separate, so it can be skipped if in presolve
    */
   SCIPdebugMsg(scip, "checking convexity of expr %p\n", (void*)expr);
   SCIP_CALL( SCIPgetConsExprQuadraticCurvature(scip, quaddata, &nlexprdata->curvature) );

   /* if convex and non propagable -> quadratic is of the form sum x_i^2 + sum y_i -> we do not want to handle it
    * TODO: we can check if it is of this form before we allocate all the data
    */
   if( (nlexprdata->curvature == SCIP_EXPRCURV_CONVEX || nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE) && ! propagable ) /*lint !e845*/
   {
      *success = FALSE;
      SCIP_CALL( nlhdlrfreeExprDataQuadratic(scip, nlhdlr, expr, nlhdlrexprdata) );
      return SCIP_OKAY;
   }

   if( nlexprdata->curvature == SCIP_EXPRCURV_CONVEX )
   {
      SCIPdebugMsg(scip, "expr %p is convex when replacing factors of bilinear terms, bases of squares and every other term by their aux vars\n",
            (void*)expr);

      /* we will estimate the expression from below, that is handle expr <= auxvar */
      *enforcedbelow = TRUE;
      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
      separate = TRUE;
   }
   else if( nlexprdata->curvature == SCIP_EXPRCURV_CONCAVE )
   {
      SCIPdebugMsg(scip, "expr %p is concave when replacing factors of bilinear terms, bases of squares and every other term by their aux vars\n",
            (void*)expr);

      /* we will estimate the expression from above, that is handle expr >= auxvar */
      *enforcedabove = TRUE;
      *success = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
      separate = TRUE;
   }

   /* we only need auxiliary variables if we are going to separate */
   if( separate )
   {
      SCIP_CONSEXPR_EXPR** linexprs;
      int nquadexprs;
      int nlinexprs;
      int i;

      SCIPgetConsExprQuadraticData(quaddata, NULL, &nlinexprs, &linexprs, NULL, &nquadexprs, NULL);

      for( i = 0; i < nlinexprs; ++i ) /* expressions appearing linearly */
      {
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, linexprs[i], NULL) );
      }
      for( i = 0; i < nquadexprs; ++i ) /* expressions appearing quadratically */
      {
         SCIP_CONSEXPR_EXPR* quadexpr;
         SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &quadexpr, NULL, NULL, NULL, NULL);
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, quadexpr, NULL) );
      }
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

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

#ifndef NDEBUG
   /* check that quaddata hasn't changed (or at least the pointer to it) */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, conshdlr, expr, &quaddata) );
   assert(quaddata == nlhdlrexprdata->quaddata);
#endif
   quaddata = nlhdlrexprdata->quaddata;

   *success = FALSE;
   *addedbranchscores = FALSE;

   /* this handler can also handle quadratic expressions whose curvature is unknown or indefinite, since it can
    * propagate them, but it does not separate these
    */
   if( nlhdlrexprdata->curvature == SCIP_EXPRCURV_UNKNOWN )
   {
      SCIPdebugMsg(scip, "not estimating due to unknown curvature\n");
      return SCIP_OKAY;
   }

   /* if estimating on non-convex side, then do nothing */
   if( ( overestimate && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && nlhdlrexprdata->curvature == SCIP_EXPRCURV_CONCAVE) )
   {
      SCIPdebugMsg(scip, "not estimating on nonconvex side (overestimate=%d, curv=%s)\n", overestimate, SCIPexprcurvGetName(nlhdlrexprdata->curvature));
      return SCIP_OKAY;
   }

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

      SCIPgetConsExprQuadraticQuadTermData(quaddata, j, &qexpr, &lincoef, &sqrcoef, &nadjbilin, NULL);
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
         return SCIP_OKAY;

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

      SCIPgetConsExprQuadraticBilinTermData(quaddata, j, &qexpr1, &qexpr2, &qcoef, NULL);

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
         return SCIP_OKAY;

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var1, coef) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var2, coef2) );
      SCIPaddRowprepConstant(rowprep, constant);
   }

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   rowprep->local = FALSE;

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
         int j;
         SCIP_INTERVAL b;
         SCIP_Real quadub;
         SCIP_Real quadlb;
         SCIP_CONSEXPR_EXPR* qexpr;
         SCIP_Real lincoef;
         SCIP_Real sqrcoef;
         int nadjbilin;
         int* adjbilin;

         SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin);

         if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, SCIPgetConsExprExprActivity(scip, qexpr)) )
         {
            SCIPintervalSetEmpty(interval);
            return SCIP_OKAY;
         }

         /* b = [c_l] */
         SCIPintervalSet(&b, lincoef);
         for( j = 0; j < nadjbilin; ++j )
         {
            SCIP_INTERVAL bterm;
            SCIP_CONSEXPR_EXPR* expr1;
            SCIP_CONSEXPR_EXPR* expr2;
            SCIP_Real bilincoef;

            SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[j], &expr1, &expr2, &bilincoef, NULL);

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
            SCIPinfoMessage(scip, NULL, "b += %g * [expr2], where <expr2> is:", bilincoef);
            SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), expr2, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
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
         SCIPinfoMessage(scip, NULL, "Computing activity for quadratic term a <expr>^2 + b <expr>, where <expr> is:");
         SCIP_CALL( SCIPprintConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), qexpr, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
         SCIPinfoMessage(scip, NULL, "a = %g, b = [%g, %g] and activity [%g, %g]\n", sqrcoef, b.inf, b.sup, quadlb, quadub);
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
      return SCIP_OKAY;

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
   SCIP_CALL( SCIPreverseConsExprExprPropagateWeightedSum(scip, conshdlr, nlinexprs, linexprs, lincoefs, constant, rhs,
            reversepropqueue, infeasible, nreductions, force) );

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
    * - use some technique (like the one in nlhdlrIntervalQuadratic), to evaluate the activity of rest_i = [quadratic
    *   expression in expr_k for k \neq i].
    * - solve a_i expr_i^2 + b expr_i \in rhs_i := [expr activity] - rest_i
    *
    * However, this might be expensive, especially computing rest_i. Hence, we implement a simpler version.
    * - we use the same partition as in nlhdlrIntervalQuadratic for the bilinear terms. This way, b = [\sum_{j \in P_i}
    *   b_ij expr_j + c_i], where P_i is the set of indices j such that expr_i * expr_j appears in that order
    * - we evaluate the activity of rest_i as sum_{k \neq i} [\min q_k, \max q_k] where q_k = a_k expr_k^2 + [\sum_{j
    *   \in P_k} b_jk expr_j + c_k] expr_k. The intervals [\min q_k, \max q_k] were already computed in
    *   nlhdlrIntervalQuadratic, so we just reuse them.
    *
    * We treat one special case:
    * - if there is a term of the form expr_l * expr_k and this is the only term that contains expr_l, then we not only
    *   propagate expr_l but also expr_k.
    *   Note that due to the above partition, the quadratic that we obtain when processing expr_k is 0*expr_k^2 +
    *   0*expr_k, thus it yields nothing.
    * - Note that rest_l = rest_k, so while processing expr_l, we only need to compute the b_k to propagate expr_k.
    *   b_k is just [b_lk * expr_l].
    *
    * Note: this implements a technique that appeared in the classic cons_quadratic.
    * The idea of the technique was to to borrow a bilinear term expr_k expr_l when propagating expr_l and the quadratic
    * function for expr_k was simple enough.
    * Since in P_l we only consider the indices of expressions that appear multiplying expr_l as _second_ factor, we
    * would lose the bilinear terms expr_k * expr_l, which contributes to the dependency problem.
    * The problem is that the contribution of b_kl * expr_k * expr_l to rest_i is not just [b_kl * expr_k * expr_l], but
    * rather quadactivities[k] (= max/min of a_k expr_k^2 + expr_k * [c_k + sum_i \in P_k b_ki expr_i]).
    * Thus, we _cannot_ just substract [b_kl * expr_k * expr_l] from rest_i.
    * But, if expr_k only appears as expr_k * expr_l, then  quadactivities[k] = [b_kl * expr_k * expr_l]. So this
    * case was handled in old cons_quadratic.
    * However, with our ordering of terms, this reduced to the special case. Indeed, assume that expr_k * expr_l is the
    * only appearance of expr_k. This means its frequency is 1 and so the frequency of expr_l has to be 1, as otherwise
    * it would be expr_l * expr_k.
    *
    *
    * TODO: handle simple cases
    * TODO: identify early when there is nothing to be gain
    */
   SCIPintervalSub(SCIP_INTERVAL_INFINITY, &rhs, SCIPgetConsExprExprActivity(scip, expr), nlhdlrexprdata->linactivity);

   for( i = 0; i < nquadexprs; ++i )
   {
      int j;
      SCIP_INTERVAL b;
      SCIP_INTERVAL rhs_i;
      SCIP_INTERVAL rest_i;
      SCIP_CONSEXPR_EXPR* qexpr;
      SCIP_Real lincoef;
      SCIP_Real sqrcoef;
      int nadjbilin;
      int* adjbilin;

      SCIP_CONSEXPR_EXPR* expr1 = NULL;
      SCIP_CONSEXPR_EXPR* expr2 = NULL;
      SCIP_Real bilincoef = 0.0;
      int pos2 = 0;

      SCIPgetConsExprQuadraticQuadTermData(quaddata, i, &qexpr, &lincoef, &sqrcoef, &nadjbilin, &adjbilin);

      /* set b to [c_l] */
      SCIPintervalSet(&b, lincoef);

      /* compute [\sum_{j \in P_l} b_lj expr_j + c_l] in b*/
      for( j = 0; j < nadjbilin; ++j )
      {
         SCIP_INTERVAL bterm;

         SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[j], &expr1, &expr2, &bilincoef, &pos2);

         if( expr1 != qexpr )
            continue;

         /* b += [b_lj * expr_j] for j \in P_l */
         SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bterm, SCIPgetConsExprExprActivity(scip, expr2), bilincoef);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, &b, b, bterm);
      }

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

      /* solve a_i expr_i^2 + b expr_i = rhs_i */
      SCIP_CALL( propagateBoundsQuadExpr(scip, conshdlr, qexpr, sqrcoef, b, rhs_i, reversepropqueue, infeasible, nreductions, force) );

      /* stop if we find infeasibility */
      if( *infeasible )
         return SCIP_OKAY;

      /* handle special case: check if the quadratic expr is of the form expr_i * expr_k and expr_i appears only once
       * (if nadjbilin == 1, then expr1, expr2, bilincoef, pos2 are still set to SCIPgetConsExprQuadraticBilinTermData(quaddata, adjbilin[0], ...))
       */
      if( lincoef == 0.0 && sqrcoef == 0.0 && nadjbilin == 1 && expr1 == qexpr )
      {
         /* expr_k should also only appear in expr_i * expr_k */
         SCIP_Real sqrcoef2;
#ifndef NDEBUG
         SCIP_Real lincoef2;
         int nadjbilin2;
         int* adjbilin2;

         SCIPgetConsExprQuadraticQuadTermData(quaddata, pos2, NULL, &lincoef2, &sqrcoef2, &nadjbilin2, &adjbilin2);
         assert(lincoef2 == 0.0 && sqrcoef2 == 0.0 && nadjbilin2 == 1 && adjbilin[0] == adjbilin2[0]);
#else
         SCIPgetConsExprQuadraticQuadTermData(quaddata, pos2, NULL, NULL, &sqrcoef2, NULL, NULL);
#endif

         /* propagate expr_k; rhs_k is equal to rhs_i; but b must change now */
         SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &b, SCIPgetConsExprExprActivity(scip, expr1), bilincoef);

         /* solve a_k expr_k^2 + b expr_k = rhs_k */
         SCIP_CALL( propagateBoundsQuadExpr(scip, conshdlr, expr2, sqrcoef2, b, rhs_i, reversepropqueue, infeasible, nreductions, force) );
      }
   }

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

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY,
            nlhdlrDetectQuadratic, nlhdlrEvalAuxQuadratic, NULL) );

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrcopyHdlrQuadratic);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataQuadratic);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateQuadratic, NULL);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalQuadratic, nlhdlrReversepropQuadratic);

   return SCIP_OKAY;
}
