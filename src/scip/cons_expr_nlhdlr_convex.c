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
 * @brief  nonlinear handler for convex expressions
 * @author Benjamin Mueller
 *
 * TODO curvature information that have been computed during the detection
 *      of other nonlinear handler can not be used right now
 *
 * TODO perturb reference point if separation fails due to too large numbers
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_convex.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_iterator.h"
#include "scip/cons_expr_var.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "convex"
#define NLHDLR_DESC         "convex handler for expressions"
#define NLHDLR_PRIORITY     50

#define DETECTSUM    FALSE

/*
 * Data structures
 */

/** expression tree */
typedef struct NlHdlr_Expr NLHDLR_EXPR;
struct NlHdlr_Expr
{
   SCIP_CONSEXPR_EXPR* origexpr;             /**< the original expression that is represented */
   NLHDLR_EXPR**       children;             /**< nlhdlr-expressions for children, or NULL if no children or to use auxvar */

   NLHDLR_EXPR*        parent;               /**< parent in expression tree TODO: remove as unused? */

   SCIP_EXPRCURV       curv;                 /**< required curvature */
   int                 depth;                /**< distance from root, root as distance 1 */

   SCIP_Real           val;                  /**< value from last eval */
   SCIP_Real           deriv;                /**< partial derivative w.r.t. root from last gradient cut computation */
};

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   NLHDLR_EXPR*          nlexpr;             /**< subexpression for which this nlhdlr estimates */
};

/*
 * static methods
 */

/** create nlhdlr-expression
 *
 * does not create children, i.e., assumes that this will be a leaf
 */
static
SCIP_RETCODE createNlHdlrExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR**         nlhdlrexpr,         /**< buffer to store created nlhdlr-expr */
   SCIP_CONSEXPR_EXPR*   origexpr,           /**< original expression to be represented */
   NLHDLR_EXPR*          parent,             /**< parent nlhdlr-expr, or NULL if root */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
)
{
   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);
   assert(origexpr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexpr) );
   (*nlhdlrexpr)->origexpr = origexpr;
   (*nlhdlrexpr)->children = NULL;
   (*nlhdlrexpr)->parent = parent;
   (*nlhdlrexpr)->curv = curv;
   (*nlhdlrexpr)->depth = parent != NULL ? parent->depth + 1 : 0;

   SCIPcaptureConsExprExpr(origexpr);

   return SCIP_OKAY;
}

/** expand nlhdlr-expression by adding children according to original expression */
static
SCIP_RETCODE growChildrenNlHdlrExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlhdlrexpr,         /**< nlhdlr-expr for which to create children */
   SCIP_EXPRCURV*        childrencurv        /**< curvature required for children, or NULL if to set to UNKNOWN */
   )
{
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);

   nchildren = SCIPgetConsExprExprNChildren(nlhdlrexpr->origexpr);
   if( nchildren == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexpr->children, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( createNlHdlrExpr(scip, &nlhdlrexpr->children[i], SCIPgetConsExprExprChildren(nlhdlrexpr->origexpr)[i], nlhdlrexpr,
         childrencurv != NULL ? childrencurv[i] : SCIP_EXPRCURV_UNKNOWN) );
   }

   return SCIP_OKAY;
}

/** free nlhdlr-expression, incl children */
static
void freeNlHdlrExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR**         nlhdlrexpr          /**< pointer to nlhdlr-expr to be freed */
)
{
   assert(nlhdlrexpr != NULL);

   if( *nlhdlrexpr == NULL )
      return;

   assert((*nlhdlrexpr)->origexpr != NULL);

   if( (*nlhdlrexpr)->children != NULL )
   {
      int nchildren = SCIPgetConsExprExprNChildren((*nlhdlrexpr)->origexpr);
      int i;

      for( i = 0; i < nchildren; ++i )
         freeNlHdlrExpr(scip, &(*nlhdlrexpr)->children[i]);

      SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexpr)->children, nchildren);
   }

   SCIPreleaseConsExprExpr(scip, &(*nlhdlrexpr)->origexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexpr);

   assert(*nlhdlrexpr == NULL);
}

/** ensure that for every leaf of a nlhdlr-expression there is an auxiliary variable */
static
SCIP_RETCODE ensureVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr              /**< pointer to nlhdlr-expr to be freed */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);

   if( nlexpr->children == NULL )
   {
      assert(SCIPgetConsExprExprAuxVar(nlexpr->origexpr) != NULL);  /* TODO remove again when we allow to introduce auxvars */
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexpr->origexpr, NULL) );
      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( ensureVariables(scip, conshdlr, nlexpr->children[i]) );
   }

   return SCIP_OKAY;
}

/** evaluate nlhdlr-expression for given solution */
static
SCIP_RETCODE evalNlHdlrExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL */
   )
{
   SCIP_Real* childrenvals;
   int nchildren;
   int i;

   assert(nlexpr != NULL);

   if( nlexpr->children == NULL )
   {
      nlexpr->val = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(nlexpr->origexpr));
      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   SCIP_CALL( SCIPallocBufferArray(scip, &childrenvals, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( evalNlHdlrExpr(scip, nlexpr->children[i], sol) );
      if( nlexpr->children[i]->val == SCIP_INVALID )
      {
         nlexpr->val = SCIP_INVALID;
         goto TERMINATE;
      }
      childrenvals[i] = nlexpr->children[i]->val;
   }

   SCIP_CALL( SCIPevalConsExprExprHdlr(scip, nlexpr->origexpr, &nlexpr->val, childrenvals, sol) );

TERMINATE:
   SCIPfreeBufferArray(scip, &childrenvals);

   return SCIP_OKAY;
}

/* differentiate nlhdlr-expression and store contribution to gradient cut in rowprep */
static
SCIP_RETCODE gradientCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_ROWPREP*         rowprep,            /**< row where to add coefficients and constants */
   SCIP_Bool*            success             /**< set to FALSE if gradient evaluation error */
   )
{
   SCIP_Real* childrenvals;
   int nchildren;
   int i;

   assert(nlexpr != NULL);

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   if( nchildren == 0 )
   {
      SCIP_VAR* var;

      if( SCIPgetConsExprExprHdlr(nlexpr->origexpr) == SCIPgetConsExprExprHdlrValue(conshdlr) )
         return SCIP_OKAY;

      var = SCIPgetConsExprExprAuxVar(nlexpr->origexpr);

      SCIPdebugMsg(scip, "add %g * (<%s> - %g) to rowprep\n", nlexpr->deriv, SCIPvarGetName(var), nlexpr->val);

      /* add deriv * (var - varval) to rowprep */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, nlexpr->deriv) );
      SCIPaddRowprepConstant(rowprep, -nlexpr->deriv * nlexpr->val);

      return SCIP_OKAY;
   }

   /* assemble children values */
   SCIP_CALL( SCIPallocBufferArray(scip, &childrenvals, nchildren) );
   for( i = 0; i < nchildren; ++i )
      childrenvals[i] = nlexpr->children[i]->val;

   /* differentiate w.r.t. each child */
   for( i = 0; i < nchildren; ++i )
   {
      /* compute partial derivative w.r.t. child i */
      SCIP_CALL( SCIPbwdiffConsExprExprHdlr(scip, nlexpr->origexpr, i, &nlexpr->children[i]->deriv, childrenvals, nlexpr->val) );
      if( nlexpr->children[i]->deriv == SCIP_INVALID )
      {
         *success = FALSE;
         goto TERMINATE;
      }
      nlexpr->children[i]->deriv *= nlexpr->deriv;

      /* push derivatives further down */
      SCIP_CALL( gradientCut(scip, conshdlr, nlexpr->children[i], rowprep, success) );
      if( !*success )
         goto TERMINATE;
   }

TERMINATE:
   SCIPfreeBufferArray(scip, &childrenvals);

   return SCIP_OKAY;
}

/** register violation as branchscore in leafs of nlhdlr-expression */
static
void branchscoreNlhdlrExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   unsigned int          brscoretag,         /**< branchscore tag */
   SCIP_Real             violation           /**< violation */
   )
{
   int nchildren;
   int c;

   if( nlexpr->children == NULL )
   {
      SCIPaddConsExprExprBranchScore(scip, nlexpr->origexpr, brscoretag, violation);
      return;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   for( c = 0; c < nchildren; ++c )
      branchscoreNlhdlrExpr(scip, nlexpr->children[c], brscoretag, violation);
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
   NLHDLR_EXPR**         rootnlexpr,         /**< buffer to store created expression */
   int*                  maxdepth,           /**< buffer to store created depth of expression */
   SCIP_CONSEXPR_EXPR*   rootexpr,           /**< expression */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
   )
{
   NLHDLR_EXPR* nlexpr;
   SCIP_EXPRCURV* childcurv;
   int childcurvsize;
   int nchildren;
   SCIP_Bool success;

   /* to do list: expressions where to check whether they can have the desired curvature when taking their children into account */
   NLHDLR_EXPR** stack;
   int stacksize;
   int stackpos;

   assert(scip != NULL);
   assert(rootnlexpr != NULL);
   assert(maxdepth != NULL);
   assert(rootexpr != NULL);
   assert(curv == SCIP_EXPRCURV_CONVEX || curv == SCIP_EXPRCURV_CONCAVE);

   childcurvsize = SCIPgetConsExprExprNChildren(rootexpr);
   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, childcurvsize) );

   /* create root nlhdlr-expression */
   SCIP_CALL( createNlHdlrExpr(scip, rootnlexpr, rootexpr, NULL, curv) );
   *maxdepth = 0;

   stacksize = 20;
   SCIP_CALL( SCIPallocBufferArray(scip, &stack, stacksize) );

   stack[0] = *rootnlexpr;
   stackpos = 0;
   while( stackpos >= 0 )
   {
      /* take expression from stack */
      nlexpr = stack[stackpos--];
      assert(nlexpr != NULL);
      assert(nlexpr->origexpr != NULL);

      *maxdepth = MAX(*maxdepth, nlexpr->depth);
      nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);

      /* TODO if sum with > N children, then always turn into variable?
       * especially if we have a big linear term, this would help to consider less vars for vertex-polyhedral funcs
       */

      /* TODO if bwdiff not implemented, then proceed further? */

      if( childcurvsize < nchildren )
      {
         childcurvsize = SCIPcalcMemGrowSize(scip, nchildren);
         SCIP_CALL( SCIPreallocBufferArray(scip, &childcurv, childcurvsize) );
      }

      success = TRUE;
      if( nlexpr->curv != SCIP_EXPRCURV_UNKNOWN )
      {
         /* check whether and under which conditions nlexpr->origexpr can have desired curvature */
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, nlexpr->origexpr, nlexpr->curv, &success, childcurv) );
         /* SCIPprintConsExprExpr(scip, conshdlr, nlexpr->origexpr, NULL);
         SCIPinfoMessage(scip, NULL, " is %s? %d\n", SCIPexprcurvGetName(nlexpr->curv), success); */
         if( success )
         {
            /* if origexpr can have curvature curv, then don't treat it as leaf, but include its children */
            SCIP_CALL( growChildrenNlHdlrExpr(scip, nlexpr, childcurv) );
         }
         else if( nlexpr == *rootnlexpr )
         {
            /* if there is no way to ensure curv for root expression, then we failed */
            freeNlHdlrExpr(scip, rootnlexpr);
            break;
         }
         else
         {
            /* for now, also require that desired curvature can be achieved w.r.t. original variables, i.e., without introducing auxvars (TODO: remove again) */
            freeNlHdlrExpr(scip, rootnlexpr);
            break;
         }
      }
      else if( nchildren > 0 )
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         SCIP_CALL( growChildrenNlHdlrExpr(scip, nlexpr, NULL) );
      }

      if( success && nchildren > 0 )
      {
         /* add children expressions to to do list (stack) */
         assert(nlexpr->children != NULL);

         if( stackpos+1 + nchildren > stacksize )
         {
            stacksize = SCIPcalcMemGrowSize(scip, stackpos + nchildren);
            SCIP_CALL( SCIPreallocBufferArray(scip, &stack, stacksize) );
         }
         memcpy(stack + (stackpos+1), nlexpr->children, nchildren * sizeof(NLHDLR_EXPR*));
         stackpos += nchildren;
      }
   }

   SCIPfreeBufferArray(scip, &stack);
   SCIPfreeBufferArray(scip, &childcurv);

   return SCIP_OKAY;
}

#if 0
static
SCIP_RETCODE constructExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   rootexpr,           /**< expression */
   SCIP_EXPRCURV         curv,               /**< curvature to achieve */
   SCIP_Bool*            success             /**< whether we were successful in enforcing this curvature */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_CONSEXPRITERATOR_USERDATA itdata;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_EXPRCURV* childcurv;
   int childcurvsize;
   int maxdepth;
   SCIP_Bool localsuccess;
   SCIP_CONSEXPR_EXPR** leafexprs = NULL;
   int nleafexprs = 0;
   int leafexprssize = 0;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(curv == SCIP_EXPRCURV_CONVEX || curv == SCIP_EXPRCURV_CONCAVE);
   assert(success != NULL);

   childcurvsize = SCIPgetConsExprExprNChildren(expr);
   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, childcurvsize) );

   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );

   /* allow revisit exprs since we may require them both convex and concave FIXME collecting leafexprs doesn't work this way! */
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, TRUE) );
   itdata.intvals[0] = (int)curv;  /* desired curvature */
   itdata.intvals[1] = 0;          /* depth in exprtree (w.r.t. root expression) */
   SCIPexpriteratorSetCurrentUserData(it, itdata);

   maxdepth = 0;
   expr = rootexpr;
   while( !SCIPexpriteratorIsEnd(it) )
   {
      itdata = SCIPexpriteratorGetCurrentUserData(it);

      maxdepth = MAX(maxdepth, itdata.intvals[1]);

      /* TODO if sum with > N children, then always turn into variable?
       * especially if we have a big linear term, this would help to consider less vars for vertex-polyhedral funcs
       */

      if( childcurvsize < SCIPgetConsExprExprNChildren(expr) )
      {
         childcurvsize = SCIPcalcMemGrowSize(scip, SCIPgetConsExprExprNChildren(expr));
         SCIP_CALL( SCIPreallocBufferArray(scip, &childcurv, childcurvsize) );
      }

      /* check whether and under which conditions expr can have desired curvature (itdata[0]) */
      if( (SCIP_EXPRCURV)itdata.intvals[0] != SCIP_EXPRCURV_UNKNOWN )
      {
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, expr, (SCIP_EXPRCURV)itdata.intvals[0], &localsuccess, childcurv) );
      }
      else
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
            childcurv[c] = SCIP_EXPRCURV_UNKNOWN;
      }

      if( localsuccess && SCIPgetConsExprExprNChildren(expr) > 0 )
      {
         /* store curvature conditions on children and depth in children */
         SCIP_CONSEXPRITERATOR_USERDATA itdatachild;

         itdatachild.intvals[1] = itdata.intvals[1] + 1;
         for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
         {
            SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[c];

            itdatachild.intvals[0] = (int)childcurv[c];
            SCIPexpriteratorSetExprUserData(it, child, itdatachild);
         }

         expr = SCIPexpriteratorGetNext(it);
      }
      else
      {
         /* either expr cannot have the desired curvature, or we are at a variable (or constant?)
          * so remember that expr should be handled as if a variable (i.e., add auxvar if not a var)
          */
         if( nleafexprs + 1 > leafexprssize )
         {
            leafexprssize = SCIPcalcMemGrowSize(scip, nleafexprs+1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &leafexprs, leafexprssize) );
         }
         leafexprs[nleafexprs++] = expr;

         /* skip remaining subtree */
         expr = SCIPexpriteratorSkipDFS(it);
      }
   }

   if( maxdepth <= 1 )
   {
      /* maxdepth = 0: root expression couldn't have the desired curvature
       * maxdepth = 1: only has desired curvature if all immediate children are (made) variables: better let exprhdlr handle then
       */
      *success = FALSE;
   }
   else
   {
      *success = TRUE;
      /* TODO ensure auxvars for leafexprs, store in nlhdlr exprdata */
   }


   SCIPfreeBufferArray(scip, &childcurv);
   SCIPfreeBufferArrayNull(scip, &leafexprs);

   SCIPexpriteratorFree(&it);


   return SCIP_OKAY;
}
#endif

/** creates nonlinear handler expression data structure */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nlhdlr expression data */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   NLHDLR_EXPR*          nlexpr              /**< nlhdlr expression */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);
   assert(nlexpr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->nlexpr = nlexpr;

#ifdef SCIP_DEBUG
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, " is handled as %s\n", SCIPexprcurvGetName(nlexpr->curv));
#endif

   /* make sure there are auxvars */
   SCIP_CALL( ensureVariables(scip, conshdlr, nlexpr) );

   return SCIP_OKAY;
}

/** frees nonlinear handler expression data structure */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to free nlhdlr expression data */
   )
{
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   freeNlHdlrExpr(scip, &(*nlhdlrexprdata)->nlexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataConvex)
{  /*lint --e{715}*/
   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );

   return SCIP_OKAY;
}

/** the detection assumes that the curvature information of the expression has been computed already */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectConvex)
{ /*lint --e{715}*/
   NLHDLR_EXPR* nlexpr = NULL;
   int depth;

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

   /* ignore sums */
   if( !DETECTSUM && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) ) /*lint !e506 !e774*/
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPgetConsExprExprNChildren(expr) == 0 )
      return SCIP_OKAY;

   /* TODO we are also interested in handling the concave side */

   if( !*enforcedbelow )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, &nlexpr, &depth, expr, SCIP_EXPRCURV_CONVEX) );
      assert(nlexpr == NULL || depth >= 1);
      if( depth < 1 )  /* TODO change back to <= 1, i.e. free if only immediate children, but no grand-daugthers */
      {
         freeNlHdlrExpr(scip, &nlexpr);
      }
      else if( nlexpr != NULL )
      {
         *enforcedbelow = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr <= auxvar\n", (void*)expr);
      }
   }

   if( !*enforcedabove && nlexpr == NULL )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, &nlexpr, &depth, expr, SCIP_EXPRCURV_CONCAVE) );
      assert(nlexpr == NULL || depth >= 1);
      if( depth < 1 )  /* TODO change back to <= 1 */
      {
         freeNlHdlrExpr(scip, &nlexpr);
      }
      else if( nlexpr != NULL )
      {
         *enforcedabove = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   /* store variable expressions into the expression data of the nonlinear handler */
   if( *success )
   {
      SCIP_CALL( createNlhdlrExprData(scip, conshdlr, nlhdlrexprdata, expr, nlexpr) );
   }
   else
   {
      assert(nlexpr == NULL);
   }

   return SCIP_OKAY;
}

/** auxiliary evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxConvex)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(nlhdlrexprdata->nlexpr->origexpr == expr);
   assert(auxvalue != NULL);

   SCIP_CALL( evalNlHdlrExpr(scip, nlhdlrexprdata->nlexpr, sol) );
   *auxvalue = nlhdlrexprdata->nlexpr->val;

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateConvex)
{ /*lint --e{715}*/
   NLHDLR_EXPR* nlexpr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(nlexpr->origexpr == expr);
   assert(nlexpr->curv == SCIP_EXPRCURV_CONVEX || nlexpr->curv == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* if estimating on non-convex side, then do nothing
    * TODO we are vertex-polyhedral and so should compute something
    */
   if( ( overestimate && nlexpr->curv == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && nlexpr->curv == SCIP_EXPRCURV_CONCAVE) )
      return SCIP_OKAY;

   /* TODO we can probably skip this as nlhdlrEvalAux was called before */
   SCIP_CALL( evalNlHdlrExpr(scip, nlexpr, sol) );
   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(nlexpr->val)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", nlexpr->val, (void*)expr);
      return SCIP_OKAY;
   }
   assert(auxvalue == nlexpr->val); /* given value (originally from nlhdlrEvalAuxConvex) should coincide with expression value (so we could skip eval) */  /*lint !e777*/

   /* compute gradient cut, add to rowprep */
   nlexpr->deriv = 1.0;
   *success = TRUE;
   SCIP_CALL( gradientCut(scip, conshdlr, nlexpr, rowprep, success) );

   if( !*success )
      return SCIP_OKAY;

   SCIPaddRowprepConstant(rowprep, nlexpr->val);
   rowprep->local = FALSE;

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%sestimate_convex%p_%s%d",
      overestimate ? "over" : "under",
      (void*)expr,
      sol != NULL ? "sol" : "lp",
      sol != NULL ? SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscoreConvex)
{ /*lint --e{715}*/
   SCIP_Real violation;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(auxvalue == SCIPgetConsExprExprValue(expr)); /* given auxvalue should have been computed by nlhdlrEvalAuxConvex */  /*lint !e777*/
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* we separate only convex functions here, so there should be little use for branching
    * if violations are small or there are numerical issues, then we will not have generated a cut, though
    * in that case, we will still branch, that is, register branchscores for all depending var exprs
    */

   /* compute violation */
   if( auxvalue == SCIP_INVALID ) /*lint !e777*/
      violation = SCIPinfinity(scip); /* evaluation error -> we should branch */
   else if( SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX  )
      violation = auxvalue - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - auxvalue;

   /* if violation is not on the side that we need to enforce, then no need for branching */
   if( violation <= 0.0 )
      return SCIP_OKAY;

   /* TODO try to figure out which variables appear linear and skip them here */
   branchscoreNlhdlrExpr(scip, nlhdlrexprdata->nlexpr, brscoretag, violation);

   *success = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrConvex)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

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

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectConvex, nlhdlrEvalAuxConvex, NULL) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvex);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateConvex, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscoreConvex);

   return SCIP_OKAY;
}
