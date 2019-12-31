/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_convex.c
 * @brief  nonlinear handlers for convex and concave expressions
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 *
 * TODO curvature information that have been computed during the detection
 *      of other nonlinear handler can not be used right now
 * TODO convex: perturb reference point if separation fails due to too large numbers
 * TODO convex: if univariate integer, then do secant on 2 nearest integers instead of tangent
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_convex.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_iterator.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_sum.h"

/* fundamental nonlinear handler properties */
#define CONVEX_NLHDLR_NAME     "convex"
#define CONVEX_NLHDLR_DESC     "handler that identifies and estimates convex expressions"
#define CONVEX_NLHDLR_PRIORITY 50

#define CONCAVE_NLHDLR_NAME     "concave"
#define CONCAVE_NLHDLR_DESC     "handler that identifies and estimates concave expressions"
#define CONCAVE_NLHDLR_PRIORITY 40

#define DEFAULT_DETECTSUM      FALSE
#define DEFAULT_PREFEREXTENDED TRUE
#define DEFAULT_CVXSIGNOMIAL   TRUE
#define DEFAULT_CVXPRODCOMP    TRUE
#define DEFAULT_HANDLETRIVIAL  FALSE

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_CONSEXPR_EXPR*   nlexpr;             /**< expression (copy) for which this nlhdlr estimates */
   SCIP_HASHMAP*         nlexpr2origexpr;    /**< mapping of our copied expression to original expression */

   int                   nleafs;             /**< number of distinct leafs of nlexpr, i.e., number of distinct (auxiliary) variables handled */
   SCIP_CONSEXPR_EXPR**  leafexprs;          /**< distinct leaf expressions (excluding value-expressions), thus variables */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool             isnlhdlrconvex;     /**< whether this data is used for the convex nlhdlr (TRUE) or the concave one (FALSE) */
   SCIP_SOL*             vpevalsol;          /**< solution used when evaluating vertex-polyhedral function in facet computation */

   /* parameters */
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a sum */
   SCIP_Bool             preferextended;     /**< whether to prefer extended formulations */

   /* advanced parameters (maybe remove some day) */
   SCIP_Bool             cvxsignomial;       /**< whether to use convexity check on signomials */
   SCIP_Bool             cvxprodcomp;        /**< whether to use convexity check on product composition f(h)*h */
   SCIP_Bool             handletrivial;      /**< whether to handle trivial expressions, i.e., those where all children are variables */
};

/** data struct to be be passed on to vertexpoly-evalfunction (see SCIPcomputeFacetVertexPolyhedral) */
typedef struct
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata;
   SCIP_SOL*                     vpevalsol;
   SCIP*                         scip;
   SCIP_CONSHDLR*                conshdlr;
} VERTEXPOLYFUN_EVALDATA;

/** stack used in constructExpr to store expressions that need to be investigated ("to do list") */
typedef struct
{
   SCIP_CONSEXPR_EXPR**  stack;              /**< stack elements */
   int                   stacksize;          /**< allocated space (in number of pointers) */
   int                   stackpos;           /**< position of top element of stack */
} EXPRSTACK;

#define DECL_CURVCHECK(x) SCIP_RETCODE x( \
   SCIP*                 scip,               /**< SCIP data structure */ \
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */ \
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< nlhdlr-expr to check */ \
   EXPRSTACK*            stack,              /**< stack where to add generated leafs */ \
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original expression */ \
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata,     /**< data of nlhdlr */ \
   SCIP_Bool*            success             /**< whether we found something */ \
   )


/*
 * static methods
 */

/** create nlhdlr-expression
 *
 * does not create children, i.e., assumes that this will be a leaf
 */
static
SCIP_RETCODE nlhdlrExprCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from copied to original expression */
   SCIP_CONSEXPR_EXPR**  nlhdlrexpr,         /**< buffer to store created expr */
   SCIP_CONSEXPR_EXPR*   origexpr,           /**< original expression to be copied */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
)
{
   assert(scip != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(nlhdlrexpr != NULL);
   assert(origexpr != NULL);

   if( SCIPgetConsExprExprNChildren(origexpr) == 0 )
   {
      /* for leaves, do not copy */
      *nlhdlrexpr = origexpr;
      SCIPcaptureConsExprExpr(*nlhdlrexpr);
      if( !SCIPhashmapExists(nlexpr2origexpr, (void*)*nlhdlrexpr) )
      {
         SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)*nlhdlrexpr, (void*)origexpr) );
      }
      return SCIP_OKAY;
   }

   /* create copy of expression, but without children */
   SCIP_CALL( SCIPduplicateConsExprExpr(scip, conshdlr, origexpr, nlhdlrexpr, FALSE) );
   assert(*nlhdlrexpr != NULL);  /* copies within the same SCIP must always work */

   /* store the curvature we want to get in the curvature flag of the copied expression
    * it's a bit of a misuse, but once we are done with everything, this is actually correct
    */
   SCIPsetConsExprExprCurvature(*nlhdlrexpr, curv);

   /* remember which the original expression was */
   SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)*nlhdlrexpr, (void*)origexpr) );

   return SCIP_OKAY;
}

/** expand nlhdlr-expression by adding children according to original expression */
static
SCIP_RETCODE nlhdlrExprGrowChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from copied to original expression */
   SCIP_CONSEXPR_EXPR*   nlhdlrexpr,         /**< expression for which to create children */
   SCIP_EXPRCURV*        childrencurv        /**< curvature required for children, or NULL if to set to UNKNOWN */
   )
{
   SCIP_CONSEXPR_EXPR* origexpr;
   SCIP_CONSEXPR_EXPR* child;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(nlhdlrexpr) == 0);

   origexpr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlhdlrexpr);

   nchildren = SCIPgetConsExprExprNChildren(origexpr);
   if( nchildren == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprCreate(scip, conshdlr, nlexpr2origexpr, &child, SCIPgetConsExprExprChildren(origexpr)[i],
         childrencurv != NULL ? childrencurv[i] : SCIP_EXPRCURV_UNKNOWN) );
      SCIP_CALL( SCIPappendConsExprExpr(scip, nlhdlrexpr, child) );
      /* append captures child, so we can release the capture from nlhdlrExprCreate */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );
   }

   assert(SCIPgetConsExprExprNChildren(nlhdlrexpr) == SCIPgetConsExprExprNChildren(origexpr));

   return SCIP_OKAY;
}

static
SCIP_DECL_VERTEXPOLYFUN(nlhdlrExprEvalConcave)
{
   VERTEXPOLYFUN_EVALDATA* evaldata = (VERTEXPOLYFUN_EVALDATA*)funcdata;
   int i;

   assert(args != NULL);
   assert(nargs == evaldata->nlhdlrexprdata->nleafs);
   assert(evaldata != NULL);

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(evaldata->scip, "eval vertexpolyfun at\n");
#endif
   for( i = 0; i < nargs; ++i )
   {
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMsg(evaldata->scip, "  <%s> = %g\n", SCIPvarGetName(SCIPgetConsExprExprVarVar(evaldata->nlhdlrexprdata->leafexprs[i])), args[i]);
#endif
      SCIP_CALL_ABORT( SCIPsetSolVal(evaldata->scip, evaldata->vpevalsol, SCIPgetConsExprExprVarVar(evaldata->nlhdlrexprdata->leafexprs[i]), args[i]) );
   }

   SCIP_CALL_ABORT( SCIPevalConsExprExpr(evaldata->scip, evaldata->conshdlr, evaldata->nlhdlrexprdata->nlexpr, evaldata->vpevalsol, 0) );

   return SCIPgetConsExprExprValue(evaldata->nlhdlrexprdata->nlexpr);
}

static
SCIP_RETCODE exprstackInit(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRSTACK*            exprstack,          /**< stack to initialize */
   int                   initsize            /**< initial size */
   )
{
   assert(scip != NULL);
   assert(exprstack != NULL);
   assert(initsize > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &exprstack->stack, initsize) );
   exprstack->stacksize = initsize;
   exprstack->stackpos = -1;

   return SCIP_OKAY;
}

static
void exprstackFree(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRSTACK*            exprstack           /**< free expression stack */
   )
{
   assert(scip != NULL);
   assert(exprstack != NULL);

   SCIPfreeBufferArray(scip, &exprstack->stack);
}

static
SCIP_RETCODE exprstackPush(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRSTACK*            exprstack,          /**< expression stack */
   int                   nexprs,             /**< number of expressions to push */
   SCIP_CONSEXPR_EXPR**  exprs               /**< expressions to push */
   )
{
   assert(scip != NULL);
   assert(exprstack != NULL);

   if( nexprs == 0 )
      return SCIP_OKAY;

   assert(exprs != NULL);

   if( exprstack->stackpos+1 + nexprs > exprstack->stacksize )
   {
      exprstack->stacksize = SCIPcalcMemGrowSize(scip, exprstack->stackpos+1 + nexprs);
      SCIP_CALL( SCIPreallocBufferArray(scip, &exprstack->stack, exprstack->stacksize) );
   }

   memcpy(exprstack->stack + (exprstack->stackpos+1), exprs, nexprs * sizeof(SCIP_CONSEXPR_EXPR*));
   exprstack->stackpos += nexprs;

   return SCIP_OKAY;
}

static
SCIP_CONSEXPR_EXPR* exprstackPop(
   EXPRSTACK*            exprstack           /**< expression stack */
   )
{
   assert(exprstack != NULL);
   assert(exprstack->stackpos >= 0);

   return exprstack->stack[exprstack->stackpos--];
}

static
SCIP_Bool exprstackIsEmpty(
   EXPRSTACK*            exprstack           /**< expression stack */
   )
{
   assert(exprstack != NULL);

   return exprstack->stackpos < 0;
}

/** looks whether top of given expression looks like a signomial that can have a given curvature
 * e.g., sqrt(x)*sqrt(y) is convex if x,y >= 0 and x and y are convex
 * unfortunately, doesn't work for tls, because i) it's originally sqrt(x*y), and ii) it is expanded into some sqrt(z*y+y)
 * but works for cvxnonsep_nsig
 */
static
DECL_CURVCHECK(curvCheckSignomial)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_Real* exponents;
   SCIP_INTERVAL* bounds;
   SCIP_EXPRCURV* curv;
   int nfactors;
   int i;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( !nlhdlrdata->cvxsignomial )
      return SCIP_OKAY;

   if( SCIPgetConsExprExprHdlr(nlexpr) != SCIPgetConsExprExprHdlrProduct(conshdlr) )
      return SCIP_OKAY;

   expr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
   assert(expr != NULL);

   nfactors = SCIPgetConsExprExprNChildren(expr);
   if( nfactors <= 1 )  /* boooring */
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &curv, nfactors) );

   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);

      if( SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         exponents[i] = 1.0;
         bounds[i] = SCIPgetConsExprExprActivity(scip, child);
      }
      else
      {
         exponents[i] = SCIPgetConsExprExprPowExponent(child);
         bounds[i] = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(child)[0]);
      }
   }

   if( !SCIPexprcurvMonomialInv(SCIPexprcurvMultiply(SCIPgetConsExprExprProductCoef(expr), SCIPgetConsExprExprCurvature(nlexpr)), nfactors, exponents, bounds, curv) )
      goto TERMINATE;

   /* add immediate children to nlexpr
    * some entries in curv actually apply to arguments of pow's, will correct this next
    */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, curv) );
   assert(SCIPgetConsExprExprNChildren(nlexpr) == nfactors);

   /* put children that are not power on stack
    * grow child for children that are power and put this child on stack
    * if preferextended, then require children with more than one child to be linear
    * unless they are linear, an auxvar will be introduced for them and thus they will be handled as var here
    */
   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPgetConsExprExprChildren(nlexpr)[i];
      assert(child != NULL);

      if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, child, &curv[i]) );
         assert(SCIPgetConsExprExprNChildren(child) == 1);
         child = SCIPgetConsExprExprChildren(child)[0];
      }
      assert(SCIPgetConsExprExprNChildren(child) == 0);

      if( nlhdlrdata->preferextended && SCIPgetConsExprExprNChildren(child) > 1 )
      {
         SCIPsetConsExprExprCurvature(child, SCIP_EXPRCURV_LINEAR);
#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "Extendedform: Require linearity for ");
         SCIPprintConsExprExpr(scip, conshdlr, child, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
#endif
      }

      SCIP_CALL( exprstackPush(scip, stack, 1, &child) );
   }

   *success = TRUE;

TERMINATE:
   SCIPfreeBufferArray(scip, &curv);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &exponents);

   return SCIP_OKAY;
}

/** looks for f(c*h(x)+d)*h(x) * constant-factor
 *
 * Assume h is univariate:
 * - First derivative is f'(c h + d) c h' h + f(c h + d) h'.
 * - Second derivative is f''(c h + d) c h' c h' h + f'(c h + d) (c h'' h + c h' h') + f'(c h + d) c h' h' + f(c h + d) h''
 *   = f''(c h + d) c^2 h'^2 h + f'(c h + d) c h'' h + 2 f'(c h + d) c h'^2 + f(c h + d) h''
 *   Remove always positive factors: f''(c h + d) h, f'(c h + d) c h'' h, f'(c h + d) c, f(c h + d) h''
 *   For convexity we want all these terms to be nonnegative. For concavity we want all of them to be nonpositive.
 *   Note, that in each term either f'(c h + d) and c occur, or none of them.
 * - Thus, f(c h(x) + d)h(x) is convex if c*f is monotonically increasing (c f' >= 0) and either
 *   - f is convex (f'' >= 0) and h is nonnegative (h >= 0) and h is convex (h'' >= 0) and [f is nonnegative (f >= 0) or h is linear (h''=0)], or
 *   - f is concave (f'' <= 0) and h is nonpositive (h <= 0) and h is concave (h'' <= 0) and [f is nonpositive (f <= 0) or h is linear (h''=0)]
 *   and f(c h(x) + d)h(x) is concave if c*f is monotonically decreasing (c f' <= 0) and either
 *   - f is convex (f'' >= 0) and h is nonpositive (h <= 0) and h is concave (h'' <= 0) and [f is nonnegative (f >= 0) or h is linear (h''=0)], or
 *   - f is concave (f'' <= 0) and h is nonnegative (h >= 0) and h is convex (h'' >= 0) and [f is nonpositive (f <= 0) or h is linear (h''=0)]
 *
 * This should hold also for multivariate and linear h, as things are invariant under linear transformations.
 * Similar to signomial, I'll assume that this will also hold for other multivariate h (someone has a formal proof?).
 */
static
DECL_CURVCHECK(curvCheckProductComposite)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* f;
   SCIP_CONSEXPR_EXPR* h = NULL;
   SCIP_Real c = 0.0;
   SCIP_CONSEXPR_EXPR* ch = NULL; /* c * h */
   SCIP_INTERVAL fbounds;
   SCIP_INTERVAL hbounds;
   SCIP_MONOTONE fmonotonicity;
   SCIP_EXPRCURV desiredcurv;
   SCIP_EXPRCURV hcurv;
   SCIP_EXPRCURV dummy;
   int fidx;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( !nlhdlrdata->cvxprodcomp )
      return SCIP_OKAY;

   if( SCIPgetConsExprExprHdlr(nlexpr) != SCIPgetConsExprExprHdlrProduct(conshdlr) )
      return SCIP_OKAY;

   expr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
   assert(expr != NULL);

   if( SCIPgetConsExprExprNChildren(expr) != 2 )
      return SCIP_OKAY;

   /* check whether we have f(c * h(x)) * h(x) or h(x) * f(c * h(x)) */
   for( fidx = 0; fidx <= 1; ++fidx )
   {
      f = SCIPgetConsExprExprChildren(expr)[fidx];

      if( SCIPgetConsExprExprNChildren(f) != 1 )
         continue;

      ch = SCIPgetConsExprExprChildren(f)[0];
      c = 1.0;
      h = ch;

      /* check whether ch is of the form c*h(x), then switch h to child ch */
      if( SCIPgetConsExprExprHdlr(ch) == SCIPgetConsExprExprHdlrSum(conshdlr) && SCIPgetConsExprExprNChildren(ch) == 1 )
      {
         c = SCIPgetConsExprExprSumCoefs(ch)[0];
         h = SCIPgetConsExprExprChildren(ch)[0];
         assert(c != 1.0 || SCIPgetConsExprExprSumConstant(ch) != 0.0);  /* we could handle this, but it should have been simplified away */
      }

#ifndef NLHDLR_CONVEX_UNITTEST
      /* can assume that duplicate subexpressions have been identified and comparing pointer is sufficient */
      if( SCIPgetConsExprExprChildren(expr)[1-fidx] == h )
#else
      /* called from unittest -> duplicate subexpressions were not identified -> compare more expensively */
      if( SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr)[1-fidx], h) == 0 )
#endif
         break;
   }
   if( fidx == 2 )
      return SCIP_OKAY;

#ifdef SCIP_MORE_DEBUG
   SCIPinfoMessage(scip, NULL, "f(c*h+d)*h with f = %s, c = %g, d = %g, h = ", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(f)), c, h != ch ? SCIPgetConsExprExprSumConstant(ch) : 0.0);
   SCIPprintConsExprExpr(scip, conshdlr, h, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   assert(c != 0.0);

   fbounds = SCIPgetConsExprExprActivity(scip, f);
   hbounds = SCIPgetConsExprExprActivity(scip, h);

   /* if h has mixed sign, then cannot conclude anything */
   if( hbounds.inf < 0.0 && hbounds.sup > 0.0 )
      return SCIP_OKAY;

   fmonotonicity = SCIPgetConsExprExprMonotonicity(scip, f, 0);

   /* if f is not monotone, then cannot conclude anything */
   if( fmonotonicity == SCIP_MONOTONE_UNKNOWN )
      return SCIP_OKAY;

   /* curvature we want to achieve (negate if product has negative coef) */
   desiredcurv = SCIPexprcurvMultiply(SCIPgetConsExprExprProductCoef(nlexpr), SCIPgetConsExprExprCurvature(nlexpr));

   /* now check the conditions as stated above */
   if( desiredcurv == SCIP_EXPRCURV_CONVEX )
   {
      /* f(c h(x)+d)h(x) is convex if c*f is monotonically increasing (c f' >= 0) and either
      *   - f is convex (f'' >= 0) and h is nonnegative (h >= 0) and h is convex (h'' >= 0) and [f is nonnegative (f >= 0) or h is linear (h''=0)], or
      *   - f is concave (f'' <= 0) and h is nonpositive (h <= 0) and h is concave (h'' <= 0) and [f is nonpositive (f <= 0) or h is linear (h''=0)]
      *  as the curvature requirements on f are on f only and not the composition f(h), we can ignore the requirements returned by SCIPcurvatureConsExprExprHdlr (last arg)
      */
      if( (c > 0.0 && fmonotonicity != SCIP_MONOTONE_INC) || (c < 0.0 && fmonotonicity != SCIP_MONOTONE_DEC) )
         return SCIP_OKAY;

      /* check whether f can be convex (h>=0) or concave (h<=0), resp., and derive requirements for h */
      if( hbounds.inf >= 0 )
      {
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, f, SCIP_EXPRCURV_CONVEX, success, &dummy) );

         /* now h also needs to be convex; and if f < 0, then h actually needs to be linear */
         if( fbounds.inf < 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONVEX;
      }
      else
      {
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, f, SCIP_EXPRCURV_CONCAVE, success, &dummy) );

         /* now h also needs to be concave; and if f > 0, then h actually needs to be linear */
         if( fbounds.sup > 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONCAVE;
      }

   }
   else
   {
      /* f(c h(x)+d)*h(x) is concave if c*f is monotonically decreasing (c f' <= 0) and either
      *   - f is convex (f'' >= 0) and h is nonpositive (h <= 0) and h is concave (h'' <= 0) and [f is nonnegative (f >= 0) or h is linear (h''=0)], or
      *   - f is concave (f'' <= 0) and h is nonnegative (h >= 0) and h is convex (h'' >= 0) and [f is nonpositive (f <= 0) or h is linear (h''=0)]
      *  as the curvature requirements on f are on f only and not the composition f(h), we can ignore the requirements returned by SCIPcurvatureConsExprExprHdlr (last arg)
      */
      if( (c > 0.0 && fmonotonicity != SCIP_MONOTONE_DEC) || (c < 0.0 && fmonotonicity != SCIP_MONOTONE_INC) )
         return SCIP_OKAY;

      /* check whether f can be convex (h<=0) or concave (h>=0), resp., and derive requirements for h */
      if( hbounds.sup <= 0 )
      {
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, f, SCIP_EXPRCURV_CONVEX, success, &dummy) );

         /* now h also needs to be concave; and if f < 0, then h actually needs to be linear */
         if( fbounds.inf < 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONCAVE;
      }
      else
      {
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, f, SCIP_EXPRCURV_CONCAVE, success, &dummy) );

         /* now h also needs to be convex; and if f > 0, then h actually needs to be linear */
         if( fbounds.sup > 0.0 )
            hcurv = SCIP_EXPRCURV_LINEAR;
         else
            hcurv = SCIP_EXPRCURV_CONVEX;
      }
   }

   if( !*success )
      return SCIP_OKAY;

   /* add immediate children (f and ch) to nlexpr; we set required curvature for h further below */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, NULL) );
   assert(SCIPgetConsExprExprNChildren(nlexpr) == 2);

   /* copy of f (and h) should have same child position in nlexpr as f (and h) has on expr (resp) */
   assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)SCIPgetConsExprExprChildren(nlexpr)[fidx]) == (void*)f);
#ifndef NLHDLR_CONVEX_UNITTEST
   assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)SCIPgetConsExprExprChildren(nlexpr)[1-fidx]) == (void*)h);
#endif
   /* push this h onto stack for further checking */
   SCIP_CALL( exprstackPush(scip, stack, 1, &(SCIPgetConsExprExprChildren(nlexpr)[1-fidx])) );

   /* h-child of product should have curvature hcurv */
   SCIPsetConsExprExprCurvature(SCIPgetConsExprExprChildren(nlexpr)[1-fidx], hcurv);

   if( h != ch )
   {
      /* add copy of ch as child to copy of f */
      SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, SCIPgetConsExprExprChildren(nlexpr)[fidx], NULL) );
      assert(SCIPgetConsExprExprNChildren(SCIPgetConsExprExprChildren(nlexpr)[fidx]) == 1);
      assert(SCIPhashmapGetImage(nlexpr2origexpr, (void*)SCIPgetConsExprExprChildren(SCIPgetConsExprExprChildren(nlexpr)[fidx])[0]) == (void*)ch);

      /* add copy of h (created above as child of product) as child in copy of ch */
      SCIP_CALL( SCIPappendConsExprExpr(scip,
         SCIPgetConsExprExprChildren(SCIPgetConsExprExprChildren(nlexpr)[fidx])[0] /* copy of ch */,
         SCIPgetConsExprExprChildren(nlexpr)[1-fidx] /* copy of h */) );
   }
   else
   {
      /* add copy of h (created above as child of product) as child in copy of f */
      SCIP_CALL( SCIPappendConsExprExpr(scip,
         SCIPgetConsExprExprChildren(nlexpr)[fidx] /* copy of f */,
         SCIPgetConsExprExprChildren(nlexpr)[1-fidx] /* copy of h */) );
   }

   return SCIP_OKAY;
}

/** use expression handlers curvature callback to check whether given curvature can be achieved */
static
DECL_CURVCHECK(curvCheckExprhdlr)
{
   SCIP_CONSEXPR_EXPR* origexpr;
   int nchildren;
   SCIP_EXPRCURV* childcurv;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   origexpr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, nlexpr);
   assert(origexpr != NULL);
   nchildren = SCIPgetConsExprExprNChildren(origexpr);

   if( nchildren == 0 )
   {
      /* if originally no children, then should be var or value, which should have every curvature, so should always be success */
      SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, origexpr, SCIPgetConsExprExprCurvature(nlexpr), success, NULL) );
      assert(*success);

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, nchildren) );

   /* check whether and under which conditions origexpr can have desired curvature */
   SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, origexpr, SCIPgetConsExprExprCurvature(nlexpr), success, childcurv) );
#ifdef SCIP_MORE_DEBUG
   SCIPprintConsExprExpr(scip, conshdlr, origexpr, NULL);
   SCIPinfoMessage(scip, NULL, " is %s? %d\n", SCIPexprcurvGetName(SCIPgetConsExprExprCurvature(nlexpr)), *success);
#endif
   if( !*success )
      goto TERMINATE;

   /* if origexpr can have curvature curv, then don't treat it as leaf, but include its children */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, childcurv) );
   assert(SCIPgetConsExprExprChildren(nlexpr) != NULL);
   assert(SCIPgetConsExprExprNChildren(nlexpr) == nchildren);

   /* If more than one child and we prefer extended formulations, then require all children to be linear.
    * Unless they are, auxvars will be introduced and they will be handles as variables, which can be an advantage in the context of extended formulations.
    */
   if( nchildren > 1 && nlhdlrdata->preferextended )
   {
      int i;
      for( i = 0; i < nchildren; ++i )
         SCIPsetConsExprExprCurvature(SCIPgetConsExprExprChildren(nlexpr)[i], SCIP_EXPRCURV_LINEAR);
#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "require linearity for children of ");
      SCIPprintConsExprExpr(scip, conshdlr, origexpr, NULL);
      SCIPinfoMessage(scip, NULL, "\n");
#endif
   }

   /* add children expressions to to-do list (stack) */
   SCIP_CALL( exprstackPush(scip, stack, nchildren, SCIPgetConsExprExprChildren(nlexpr)) );

TERMINATE:
   SCIPfreeBufferArray(scip, &childcurv);

   return SCIP_OKAY;
}

/** curvature check and expression-growing methods
 * some day this could be plugins added by users at runtime, but for now we have a fixed list here
 * NOTE: curvCheckExprhdlr should be last
 */
static DECL_CURVCHECK((*CURVCHECKS[])) = { curvCheckProductComposite, curvCheckSignomial, curvCheckExprhdlr };
/** number of curvcheck methods */
static const int NCURVCHECKS = sizeof(CURVCHECKS) / sizeof(void*);

/** checks whether expression is a sum with more than one child and each child being a variable */
static
SCIP_Bool exprIsMultivarLinear(
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression to check */
   )
{
   int nchildren;
   int c;

   assert(conshdlr != NULL);
   assert(expr != NULL);

   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      return FALSE;

   nchildren = SCIPgetConsExprExprNChildren(expr);
   if( nchildren <= 1 )
      return FALSE;

   for( c = 0; c < nchildren; ++c )
      if( SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[c]) != SCIPgetConsExprExprHdlrVar(conshdlr) )
         return FALSE;

   return TRUE;
}

/** construct a subexpression (as nlhdlr-expression) of maximal size that has a given curvature
 *
 * If the curvature cannot be achieved for an expression in the original expression graph,
 * then this expression becomes a leaf in the nlhdlr-expression.
 *
 * Sets *rootnlexpr to NULL if failed.
 */
static
SCIP_RETCODE constructExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata,     /**< nonlinear handler data */
   SCIP_CONSEXPR_EXPR**  rootnlexpr,         /**< buffer to store created expression */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original expression */
   int*                  nleafs,             /**< number of leafs in constructed expression */
   SCIP_CONSEXPR_EXPR*   rootexpr,           /**< expression */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
   )
{
   SCIP_CONSEXPR_EXPR* nlexpr;
   EXPRSTACK stack; /* to do list: expressions where to check whether they can have the desired curvature when taking their children into account */
   int oldstackpos;

   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(rootnlexpr != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(nleafs != NULL);
   assert(rootexpr != NULL);
   assert(curv == SCIP_EXPRCURV_CONVEX || curv == SCIP_EXPRCURV_CONCAVE);

   /* create root expression */
   SCIP_CALL( nlhdlrExprCreate(scip, conshdlr, nlexpr2origexpr, rootnlexpr, rootexpr, curv) );

   *nleafs = 0;

   SCIP_CALL( exprstackInit(scip, &stack, 20) );
   SCIP_CALL( exprstackPush(scip, &stack, 1, rootnlexpr) );
   while( !exprstackIsEmpty(&stack) )
   {
      /* take expression from stack */
      nlexpr = exprstackPop(&stack);
      assert(nlexpr != NULL);
      assert(SCIPgetConsExprExprNChildren(nlexpr) == 0);

      oldstackpos = stack.stackpos;
      if( nlhdlrdata->isnlhdlrconvex && !SCIPhasConsExprExprHdlrBwdiff(SCIPgetConsExprExprHdlr(nlexpr)) )
      {
         /* if bwdiff is not implemented, then we could not generate cuts in the convex nlhdlr, so "stop" (treat nlexpr as variable) */
      }
      else if( !nlhdlrdata->isnlhdlrconvex && exprIsMultivarLinear(conshdlr, (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr)) )
      {
         /* if we are in the concave handler, we would like to treat linear multivariate subexpressions by a new auxvar always,
          * e.g., handle log(x+y) as log(z), z=x+y, because the estimation problem will be smaller then without making the estimator worse
          * (cons_nonlinear does this, too)
          * TODO this check isn't sufficient to also handle sums that become linear after we add auxvars for some children
          */
#ifdef SCIP_MORE_DEBUG
         SCIPprintConsExprExpr(scip, conshdlr, SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr), NULL);
         SCIPinfoMessage(scip, NULL, "... is a multivariate linear sum that we'll treat as auxvar\n");
#endif
      }
      else if( SCIPgetConsExprExprCurvature(nlexpr) != SCIP_EXPRCURV_UNKNOWN )
      {
         SCIP_Bool success;
         int method;

         /* try through curvature check methods until one succeeds */
         for( method = 0; method < NCURVCHECKS; ++method )
         {
            SCIP_CALL( CURVCHECKS[method](scip, conshdlr, nlexpr, &stack, nlexpr2origexpr, nlhdlrdata, &success) );
            if( success )
               break;
         }
      }
      else
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, NULL) );

         /* add children expressions, if any, to to-do list (stack) */
         SCIP_CALL( exprstackPush(scip, &stack, SCIPgetConsExprExprNChildren(nlexpr), SCIPgetConsExprExprChildren(nlexpr)) );
      }
      assert(stack.stackpos >= oldstackpos);  /* none of the methods above should have removed something from the stack */

      /* if nothing was added, then none of the successors of nlexpr were added to the stack
       * this is either because nlexpr was already a variable or value expressions, thus a leaf,
       * or because the desired curvature could not be achieved, so it will be handled as variables, thus a leaf
       */
      if( stack.stackpos == oldstackpos )
         ++*nleafs;
   }

   exprstackFree(scip, &stack);

   if( *rootnlexpr != NULL )
   {
      SCIP_Bool istrivial = TRUE;

      /* if handletrivial is enabled, then only require that rootnlexpr itself has required curvature (so has children; see below) and
       * that we are not a trivial sum  (because the previous implementation of this nlhdlr didn't allow this, either)
       */
      if( !nlhdlrdata->handletrivial || SCIPgetConsExprExprHdlr(*rootnlexpr) == SCIPgetConsExprExprHdlrSum(conshdlr) )
      {
         /* if all children do not have children, i.e., are variables, or will be replaced by auxvars, then free
          * also if rootnlexpr has no children, then free
          */
         int i;
         for( i = 0; i < SCIPgetConsExprExprNChildren(*rootnlexpr); ++i )
         {
            if( SCIPgetConsExprExprNChildren(SCIPgetConsExprExprChildren(*rootnlexpr)[i]) > 0 )
            {
               istrivial = FALSE;
               break;
            }
         }
      }
      else if( SCIPgetConsExprExprNChildren(*rootnlexpr) > 0 )  /* if handletrivial, then just require children */
            istrivial = FALSE;

      if( istrivial )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, rootnlexpr) );
      }
   }

   return SCIP_OKAY;
}

/** collect (non-value) leaf expressions and ensure that they correspond to a variable (original or auxiliary)
 *
 * For children where we could not achieve the desired curvature, introduce an auxvar and replace the child by a var-expression that points to this auxvar.
 * Collect all leaf expressions (if not a value-expression) and index them.
 */
static
SCIP_RETCODE collectLeafs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< nlhdlr-expr */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original */
   SCIP_HASHMAP*         leaf2index,         /**< mapping from leaf to index */
   int*                  nindices            /**< number of indices */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;

   assert(nlexpr != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(leaf2index != NULL);
   assert(nindices != NULL);

   assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);
   assert(SCIPgetConsExprExprChildren(nlexpr) != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );
   SCIP_CALL( SCIPexpriteratorInit(it, nlexpr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_CONSEXPRITERATOR_VISITINGCHILD);

   for( nlexpr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); nlexpr = SCIPexpriteratorGetNext(it) ) /*lint !e441*/
   {
      SCIP_CONSEXPR_EXPR* child;

      assert(nlexpr != NULL);

      /* check whether to-be-visited child needs to be replaced by a new expression (representing the auxvar) */
      child = SCIPexpriteratorGetChildExprDFS(it);
      if( SCIPgetConsExprExprNChildren(child) == 0 )
      {
         SCIP_CONSEXPR_EXPR* origexpr;

         origexpr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)child);
         assert(origexpr != NULL);

         if( SCIPgetConsExprExprNChildren(origexpr) > 0 )
         {
            SCIP_CONSEXPR_EXPR* newchild;
            int childidx;
            SCIP_VAR* var;

            /* having a child that had children in original but not in copy means that we could not achieve the desired curvature
             * thus, replace by a new child that points to the auxvar of the original expression
             */
            SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, origexpr, &var) );
            assert(var != NULL);
            SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &newchild, var) );  /* this captures newchild once */

            childidx = SCIPexpriteratorGetChildIdxDFS(it);
            SCIP_CALL( SCIPreplaceConsExprExprChild(scip, nlexpr, childidx, newchild) );  /* this captures newchild again */

            /* do not remove child->origexpr from hashmap, as child may appear again due to common subexprs (created by curvCheckProductComposite, for example)
             * if it doesn't reappear, though, but the memory address is reused, we need to make sure it points to the right origexpr
             */
            /* SCIP_CALL( SCIPhashmapRemove(nlexpr2origexpr, (void*)child) ); */
            SCIPhashmapSetImage(nlexpr2origexpr, (void*)newchild, (void*)origexpr);

            if( !SCIPhashmapExists(leaf2index, (void*)newchild) )
            {
               /* new leaf -> new index and remember in hashmap */
               SCIP_CALL( SCIPhashmapInsertInt(leaf2index, (void*)newchild, (*nindices)++) );
            }

            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &newchild) );  /* because it was captured by both create and replace */
         }
         else if( SCIPisConsExprExprVar(child) )
         {
            /* if variable, then add to hashmap, if not already there */
            if( !SCIPhashmapExists(leaf2index, (void*)child) )
            {
               SCIP_CALL( SCIPhashmapInsertInt(leaf2index, (void*)child, (*nindices)++) );
            }
         }
         /* else: it's probably a value-expression, nothing to do */
      }
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** creates nonlinear handler expression data structure */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nlhdlr expression data */
   SCIP_CONSEXPR_EXPR*   expr,               /**< original expression */
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< our copy of expression */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping of expression copy to original */
   int                   nleafs              /**< number of leafs as counted by constructExpr */
   )
{
   SCIP_HASHMAP* leaf2index;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);
   assert(nlexpr != NULL);
   assert(nlexpr2origexpr != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->nlexpr = nlexpr;
   (*nlhdlrexprdata)->nlexpr2origexpr = nlexpr2origexpr;

   /* make sure there are auxvars and collect all variables */
   SCIP_CALL( SCIPhashmapCreate(&leaf2index, SCIPblkmem(scip), nleafs) );
   (*nlhdlrexprdata)->nleafs = 0;  /* we start a new count, this time skipping value-expressions */
   SCIP_CALL( collectLeafs(scip, conshdlr, nlexpr, nlexpr2origexpr, leaf2index, &(*nlhdlrexprdata)->nleafs) );
   assert((*nlhdlrexprdata)->nleafs <= nleafs);  /* we should not have seen more leafs now than in constructExpr */

   /* assemble auxvars array */
   assert((*nlhdlrexprdata)->nleafs > 0);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->leafexprs, (*nlhdlrexprdata)->nleafs) );
   for( i = 0; i < SCIPhashmapGetNEntries(leaf2index); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_CONSEXPR_EXPR* leaf;
      int idx;

      entry = SCIPhashmapGetEntry(leaf2index, i);
      if( entry == NULL )
         continue;

      leaf = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetOrigin(entry);
      assert(leaf != NULL);
      assert(SCIPgetConsExprExprAuxVar(leaf) != NULL);

      idx = SCIPhashmapEntryGetImageInt(entry);
      assert(idx >= 0);
      assert(idx < (*nlhdlrexprdata)->nleafs);

      (*nlhdlrexprdata)->leafexprs[idx] = leaf;

      SCIPdebugMsg(scip, "leaf %d: <%s>\n", idx, SCIPvarGetName(SCIPgetConsExprExprAuxVar(leaf)));
   }

   SCIPhashmapFree(&leaf2index);

#ifdef SCIP_DEBUG
   SCIPprintConsExprExpr(scip, conshdlr, nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " (%p) is handled as %s\n", SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr), SCIPexprcurvGetName(SCIPgetConsExprExprCurvature(nlexpr)));
#endif

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrfreeHdlrDataConvexConcave)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(*nlhdlrdata != NULL);
   assert((*nlhdlrdata)->vpevalsol == NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataConvexConcave)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->leafexprs, (*nlhdlrexprdata)->nleafs);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*nlhdlrexprdata)->nlexpr) );
   SCIPhashmapFree(&(*nlhdlrexprdata)->nlexpr2origexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_EXPR* nlexpr = NULL;
   SCIP_HASHMAP* nlexpr2origexpr;
   int nleafs = 0;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* we currently cannot contribute in presolve */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);
   assert(nlhdlrdata->isnlhdlrconvex);

   /* ignore sums if > 1 children
    * NOTE: this means that for something like 1+f(x), even if f is a trivial convex expression, we would handle 1+f(x)
    * with this nlhdlr, instead of formulating this as 1+z and handling z=f(x) with the default nlhdlr, i.e., the exprhdlr
    * today, I prefer handling this here, as it avoid introducing an extra auxiliary variable
    */
   if( !nlhdlrdata->detectsum && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) && SCIPgetConsExprExprNChildren(expr) > 1 )
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPgetConsExprExprNChildren(expr) == 0 )
      return SCIP_OKAY;

   /* initialize mapping from copied expression to original one
    * 20 is not a bad estimate for the size of convex subexpressions that we can usually discover
    * when expressions will be allowed to store "user"data, we could get rid of this hashmap (TODO)
    */
   SCIP_CALL( SCIPhashmapCreate(&nlexpr2origexpr, SCIPblkmem(scip), 20) );

   if( !*enforcedbelow )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr, SCIP_EXPRCURV_CONVEX) );
      if( nlexpr != NULL )
      {
         assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);  /* should not be trivial */

         *enforcedbelow = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr <= auxvar\n", (void*)expr);
      }
      else
      {
         SCIP_CALL( SCIPhashmapRemoveAll(nlexpr2origexpr) );
      }
   }

   if( !*enforcedabove && nlexpr == NULL )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr, SCIP_EXPRCURV_CONCAVE) );
      if( nlexpr != NULL )
      {
         assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);  /* should not be trivial */

         *enforcedabove = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   assert(*success || nlexpr == NULL);
   if( !*success )
   {
      SCIPhashmapFree(&nlexpr2origexpr);
      return SCIP_OKAY;
   }

   /* store variable expressions into the expression data of the nonlinear handler */
   SCIP_CALL( createNlhdlrExprData(scip, conshdlr, nlhdlrexprdata, expr, nlexpr, nlexpr2origexpr, nleafs) );

   return SCIP_OKAY;
}

/** auxiliary evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxConvexConcave)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(auxvalue != NULL);

   SCIP_CALL( SCIPevalConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), nlhdlrexprdata->nlexpr, sol, 0) );
   *auxvalue = SCIPgetConsExprExprValue(nlhdlrexprdata->nlexpr);

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* nlexpr;
   SCIP_EXPRCURV curvature;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlexpr) == expr);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;
   *addedbranchscores = FALSE;

   /* if estimating on non-convex side, then do nothing */
   curvature = SCIPgetConsExprExprCurvature(nlexpr);
   assert(curvature == SCIP_EXPRCURV_CONVEX || curvature == SCIP_EXPRCURV_CONCAVE);
   if( ( overestimate && curvature == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && curvature == SCIP_EXPRCURV_CONCAVE) )
      return SCIP_OKAY;

   /* we can skip eval as nlhdlrEvalAux should have been called for same solution before */
   /* SCIP_CALL( nlhdlrExprEval(scip, nlexpr, sol) ); */
   assert(auxvalue == SCIPgetConsExprExprValue(nlexpr)); /* given value (originally from nlhdlrEvalAuxConvexConcave) should coincide with the one stored in nlexpr */  /*lint !e777*/
   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(auxvalue)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", auxvalue, (void*)expr);
      return SCIP_OKAY;
   }

   /* compute gradient (TODO: this also reevaluates (soltag=0), which shouldn't be necessary) */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, nlexpr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(nlexpr) == SCIP_INVALID ) /*lint !e777*/
   {
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* add gradient underestimator to rowprep: first contribution of each variable, (x - sol) \nabla f(sol) */
   *success = TRUE;
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real deriv;
      SCIP_Real varval;

      var = SCIPgetConsExprExprAuxVar(nlhdlrexprdata->leafexprs[i]);
      assert(var != NULL);

      deriv = SCIPgetConsExprExprPartialDiff(scip, conshdlr, nlexpr, var);
      if( deriv == SCIP_INVALID ) /*lint !e777*/
      {
         SCIPdebugMsg(scip, "gradient evaluation error for component %d of %p\n", i, (void*)expr);
         *success = FALSE;
         break;
      }

      varval = SCIPgetSolVal(scip, sol, var);

      SCIPdebugMsg(scip, "add %g * (<%s> - %g) to rowprep\n", deriv, SCIPvarGetName(var), varval);

      /* add deriv * (var - varval) to rowprep */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, deriv) );
      SCIPaddRowprepConstant(rowprep, -deriv * varval);
   }

   /* next add f(sol) */
   SCIPaddRowprepConstant(rowprep, auxvalue);
   rowprep->local = FALSE;

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%sestimate_convex%p_%s%d",
      overestimate ? "over" : "under",
      (void*)expr,
      sol != NULL ? "sol" : "lp",
      sol != NULL ? SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrConvex)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), CONVEX_NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrConvex(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes convex nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   nlhdlrdata->isnlhdlrconvex = TRUE;
   nlhdlrdata->vpevalsol = NULL;

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, CONVEX_NLHDLR_NAME, CONVEX_NLHDLR_DESC, CONVEX_NLHDLR_PRIORITY, nlhdlrDetectConvex, nlhdlrEvalAuxConvexConcave, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONVEX_NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONVEX_NLHDLR_NAME "/preferextended",
      "whether to prefer extended formulations",
      &nlhdlrdata->preferextended, FALSE, DEFAULT_PREFEREXTENDED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONVEX_NLHDLR_NAME "/cvxsignomial",
      "whether to use convexity check on signomials",
      &nlhdlrdata->cvxsignomial, TRUE, DEFAULT_CVXSIGNOMIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONVEX_NLHDLR_NAME "/cvxprodcomp",
      "whether to use convexity check on product composition f(h)*h",
      &nlhdlrdata->cvxprodcomp, TRUE, DEFAULT_CVXPRODCOMP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONVEX_NLHDLR_NAME "/handletrivial",
      "whether to also handle trivial convex expressions",
      &nlhdlrdata->handletrivial, TRUE, DEFAULT_HANDLETRIVIAL, NULL, NULL) );

   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrfreeHdlrDataConvexConcave);
   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvexConcave);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateConvex, NULL);

   return SCIP_OKAY;
}






static
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitConcave)
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->vpevalsol != NULL )
   {
      SCIPfreeSol(scip, &nlhdlrdata->vpevalsol);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectConcave)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_EXPR* nlexpr = NULL;
   SCIP_HASHMAP* nlexpr2origexpr;
   int nleafs = 0;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* we currently cannot contribute in presolve */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);
   assert(!nlhdlrdata->isnlhdlrconvex);

   /* ignore sums if > 1 children for now
    * - for f(x) + g(y), i.e., distinct variables, it's better to handle f(x) and g(y) separately, as this keeps the estimation problem smaller and doesn't make the estimators worse
    * - for f(x) + g(x), i.e., same variables, it could actually be better to handle them jointly, because we might get tighter estimators (?)
    * - but we have no simple check which situation we are in (could well be something in between), so I'm going for the first way by default for now
    */
   if( !nlhdlrdata->detectsum && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) && SCIPgetConsExprExprNChildren(expr) > 1 )
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPgetConsExprExprNChildren(expr) == 0 )
      return SCIP_OKAY;

   /* initialize mapping from copied expression to original one
    * 20 is not a bad estimate for the size of concave subexpressions that we can usually discover
    * when expressions will be allowed to store "user"data, we could get rid of this hashmap (TODO)
    */
   SCIP_CALL( SCIPhashmapCreate(&nlexpr2origexpr, SCIPblkmem(scip), 20) );

   if( !*enforcedbelow )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr, SCIP_EXPRCURV_CONCAVE) );

      if( nlexpr != NULL && nleafs > SCIP_MAXVERTEXPOLYDIM )
      {
         SCIPdebugMsg(scip, "Too many variables (%d) in constructed expression. Will not be able to estimate. Rejecting.\n", nleafs);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &nlexpr) );
      }

      if( nlexpr != NULL )
      {
         assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);  /* should not be trivial */

         *enforcedbelow = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr <= auxvar\n", (void*)expr);
      }
      else
      {
         SCIP_CALL( SCIPhashmapRemoveAll(nlexpr2origexpr) );
      }
   }

   if( !*enforcedabove && nlexpr == NULL )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, nlexpr2origexpr, &nleafs, expr, SCIP_EXPRCURV_CONVEX) );

      if( nlexpr != NULL && nleafs > 14 )
      {
         SCIPdebugMsg(scip, "Too many variables (%d) in constructed expression. Will not be able to estimate. Rejecting.\n", nleafs);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &nlexpr) );
      }

      if( nlexpr != NULL )
      {
         assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);  /* should not be trivial */

         *enforcedabove = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   assert(*success || nlexpr == NULL);
   if( !*success )
   {
      SCIPhashmapFree(&nlexpr2origexpr);
      return SCIP_OKAY;
   }

   /* store variable expressions into the expression data of the nonlinear handler */
   SCIP_CALL( createNlhdlrExprData(scip, conshdlr, nlhdlrexprdata, expr, nlexpr, nlexpr2origexpr, nleafs) );

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateConcave)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_EXPR* nlexpr;
   SCIP_EXPRCURV curvature;
   VERTEXPOLYFUN_EVALDATA evaldata;
   SCIP_Real* xstar;
   SCIP_Real* box;
   SCIP_Real facetconstant;
   SCIP_VAR* var;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlexpr) == expr);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;
   *addedbranchscores = FALSE;

   /* if estimating on non-concave side, then do nothing */
   curvature = SCIPgetConsExprExprCurvature(nlexpr);
   assert(curvature == SCIP_EXPRCURV_CONVEX || curvature == SCIP_EXPRCURV_CONCAVE);
   if( ( overestimate && curvature == SCIP_EXPRCURV_CONCAVE) ||
       (!overestimate && curvature == SCIP_EXPRCURV_CONVEX) )
      return SCIP_OKAY;

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "%sestimate concave expression ", overestimate ? "over" : "under");
   SCIPprintConsExprExpr(scip, conshdlr, nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " at point\n");
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      var = SCIPgetConsExprExprVarVar(nlhdlrexprdata->leafexprs[i]);
      assert(var != NULL);

      SCIPinfoMessage(scip, NULL, "  <%s> = %g [%g,%g]\n", SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
   }
#endif

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( nlhdlrdata->vpevalsol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, &nlhdlrdata->vpevalsol, NULL) );
   }

   evaldata.nlhdlrexprdata = nlhdlrexprdata;
   evaldata.vpevalsol = nlhdlrdata->vpevalsol;
   evaldata.scip = scip;
   evaldata.conshdlr = conshdlr;

   SCIP_CALL( SCIPallocBufferArray(scip, &xstar, nlhdlrexprdata->nleafs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &box, 2*nlhdlrexprdata->nleafs) );

   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
   {
      var = SCIPgetConsExprExprVarVar(nlhdlrexprdata->leafexprs[i]);
      assert(var != NULL);

      box[2*i] = SCIPvarGetLbLocal(var);
      if( SCIPisInfinity(scip, -box[2*i]) )
      {
         SCIPdebugMsg(scip, "lower bound at -infinity, no estimate possible\n");
         goto TERMINATE;
      }

      box[2*i+1] = SCIPvarGetUbLocal(var);
      if( SCIPisInfinity(scip, box[2*i+1]) )
      {
         SCIPdebugMsg(scip, "upper bound at +infinity, no estimate possible\n");
         goto TERMINATE;
      }

      xstar[i] = SCIPgetSolVal(scip, sol, var);
      assert(xstar[i] != SCIP_INVALID);
   }

   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nlhdlrexprdata->nleafs) );

   SCIP_CALL( SCIPcomputeFacetVertexPolyhedral(scip, conshdlr, overestimate, nlhdlrExprEvalConcave, (void*)&evaldata,
      xstar, box, nlhdlrexprdata->nleafs, targetvalue, success, rowprep->coefs, &facetconstant) );

   if( !*success )
   {
      SCIPdebugMsg(scip, "failed to compute facet of convex hull\n");
      goto TERMINATE;
   }

   rowprep->local = TRUE;
   rowprep->side = -facetconstant;
   rowprep->nvars = nlhdlrexprdata->nleafs;
   for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
      rowprep->vars[i] = SCIPgetConsExprExprVarVar(nlhdlrexprdata->leafexprs[i]);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "computed estimator: ");
   SCIPprintRowprep(scip, rowprep, NULL);
#endif

 TERMINATE:
   if( addbranchscores )
   {
      SCIP_Real violation;

      /* check how much is the violation on the side that we estimate */
      if( auxvalue == SCIP_INVALID ) /*lint !e777*/
      {
         /* if cannot evaluate, then always branch */
         violation = SCIPinfinity(scip);
      }
      else
      {
         SCIP_Real auxval;

         /* get value of auxiliary variable of this expression */
         assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
         auxval = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));

         /* compute the violation
          * if we underestimate, then we enforce expr <= auxval, so violation is (positive part of) auxvalue - auxval
          * if we overestimate,  then we enforce expr >= auxval, so violation is (positive part of) auxval - auxvalue
          */
         if( !overestimate )
            violation = MAX(0.0, auxvalue - auxval);
         else
            violation = MAX(0.0, auxval - auxvalue);
      }
      assert(violation >= 0.0);

      /* TODO should/could do something more elaborate as in cons_expr_product */
      for( i = 0; i < nlhdlrexprdata->nleafs; ++i )
         SCIPaddConsExprExprBranchScore(scip, conshdlr, SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, nlhdlrexprdata->leafexprs[i]), REALABS(violation));

      *addedbranchscores = TRUE;
   }

   SCIPfreeBufferArrayNull(scip, &box);
   SCIPfreeBufferArrayNull(scip, &xstar);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrConcave)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), CONCAVE_NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrConcave(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes concave nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrConcave(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   nlhdlrdata->isnlhdlrconvex = FALSE;
   nlhdlrdata->vpevalsol = NULL;

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, CONCAVE_NLHDLR_NAME, CONCAVE_NLHDLR_DESC, CONCAVE_NLHDLR_PRIORITY, nlhdlrDetectConcave, nlhdlrEvalAuxConvexConcave, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONCAVE_NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   /*SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONCAVE_NLHDLR_NAME "/preferextended",
      "whether to prefer extended formulations",
      &nlhdlrdata->preferextended, FALSE, DEFAULT_PREFEREXTENDED, NULL, NULL) );*/
   /* "extended" formulations of a concave expressions can give worse estimators */
   nlhdlrdata->preferextended = FALSE;

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONCAVE_NLHDLR_NAME "/cvxsignomial",
      "whether to use convexity check on signomials",
      &nlhdlrdata->cvxsignomial, TRUE, DEFAULT_CVXSIGNOMIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONCAVE_NLHDLR_NAME "/cvxprodcomp",
      "whether to use convexity check on product composition f(h)*h",
      &nlhdlrdata->cvxprodcomp, TRUE, DEFAULT_CVXPRODCOMP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" CONCAVE_NLHDLR_NAME "/handletrivial",
      "whether to also handle trivial convex expressions",
      &nlhdlrdata->handletrivial, TRUE, DEFAULT_HANDLETRIVIAL, NULL, NULL) );

   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrfreeHdlrDataConvexConcave);
   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConcave);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvexConcave);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateConcave, NULL);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, NULL, nlhdlrExitConcave);

   return SCIP_OKAY;
}
