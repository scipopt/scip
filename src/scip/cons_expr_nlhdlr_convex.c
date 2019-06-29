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

#define DEFAULT_DETECTSUM   FALSE

/*
 * Data structures
 */

/** expression tree */
typedef struct NlHdlr_Expr NLHDLR_EXPR;
struct NlHdlr_Expr
{
   SCIP_CONSEXPR_EXPR*   origexpr;           /**< the original expression that is represented */
   NLHDLR_EXPR**         children;           /**< nlhdlr-expressions for children, or NULL if no children or to use auxvar */
   int                   varidx;             /**< if expression without children, then index of auxvar in nlhdlr exprdata (auxvarexprs) */

   NLHDLR_EXPR*          parent;             /**< parent in expression tree TODO: remove as unused? */

   SCIP_EXPRCURV         curv;               /**< required curvature */
   int                   depth;              /**< distance from root, root as distance 1 */

   SCIP_Real*            childrenval;        /**< buffer to store value of children */
   SCIP_Real             val;                /**< value from last eval */
   SCIP_Real             deriv;              /**< partial derivative w.r.t. root from last gradient cut computation */
};

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   NLHDLR_EXPR*          nlexpr;             /**< subexpression for which this nlhdlr estimates */

   int                   nauxvars;           /**< number of distinct (auxiliary) variables handled */
   SCIP_CONSEXPR_EXPR**  auxvarexprs;        /**< original expressions which auxiliary variable we use */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a sum */
};

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
   (*nlhdlrexpr)->childrenval = NULL;
   (*nlhdlrexpr)->parent = parent;
   (*nlhdlrexpr)->curv = curv;
   (*nlhdlrexpr)->depth = parent != NULL ? parent->depth + 1 : 0;

   SCIPcaptureConsExprExpr(origexpr);

   return SCIP_OKAY;
}

/** expand nlhdlr-expression by adding children according to original expression */
static
SCIP_RETCODE nlhdlrExprGrowChildren(
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
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexpr->childrenval, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprCreate(scip, &nlhdlrexpr->children[i], SCIPgetConsExprExprChildren(nlhdlrexpr->origexpr)[i], nlhdlrexpr,
         childrencurv != NULL ? childrencurv[i] : SCIP_EXPRCURV_UNKNOWN) );
   }

   return SCIP_OKAY;
}

/** free nlhdlr-expression, incl children */
static
void nlhdlrExprFree(
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
         nlhdlrExprFree(scip, &(*nlhdlrexpr)->children[i]);

      SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexpr)->childrenval, nchildren);
      SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexpr)->children, nchildren);
   }

   SCIPreleaseConsExprExpr(scip, &(*nlhdlrexpr)->origexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexpr);

   assert(*nlhdlrexpr == NULL);
}

/** evaluate nlhdlr-expression for given solution */
static
SCIP_RETCODE nlhdlrExprEval(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);

   if( nlexpr->children == NULL )
   {
      nlexpr->val = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(nlexpr->origexpr));
      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprEval(scip, nlexpr->children[i], sol) );
      if( nlexpr->children[i]->val == SCIP_INVALID )
      {
         nlexpr->val = SCIP_INVALID;
         return SCIP_OKAY;
      }
      nlexpr->childrenval[i] = nlexpr->children[i]->val;
   }

   SCIP_CALL( SCIPevalConsExprExprHdlr(scip, nlexpr->origexpr, &nlexpr->val, nlexpr->childrenval, sol) );

   return SCIP_OKAY;
}

/* differentiate nlhdlr-expression and store contribution to gradient cut in rowprep */
static
SCIP_RETCODE nlhdlrExprGradientCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_ROWPREP*         rowprep,            /**< row where to add coefficients and constants */
   SCIP_Bool*            success             /**< set to FALSE if gradient evaluation error */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);
   assert(success != NULL);

   *success = TRUE;

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

   /* assemble children values
    * should still be uptodate from last eval
   for( i = 0; i < nchildren; ++i )
      nlexpr->childrenval[i] = nlexpr->children[i]->val;
    */

   /* differentiate w.r.t. each child */
   for( i = 0; i < nchildren; ++i )
   {
      /* compute partial derivative w.r.t. child i */
      SCIP_CALL( SCIPbwdiffConsExprExprHdlr(scip, nlexpr->origexpr, i, &nlexpr->children[i]->deriv, nlexpr->childrenval, nlexpr->val) );
      if( nlexpr->children[i]->deriv == SCIP_INVALID )
      {
         *success = FALSE;
         break;
      }
      nlexpr->children[i]->deriv *= nlexpr->deriv;

      /* push derivatives further down */
      SCIP_CALL( nlhdlrExprGradientCut(scip, conshdlr, nlexpr->children[i], rowprep, success) );
      if( !*success )
         break;
   }

   return SCIP_OKAY;
}

/** register violation as branchscore in leafs of nlhdlr-expression */
static
void nlhdlrExprBranchscore(
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
      nlhdlrExprBranchscore(scip, nlexpr->children[c], brscoretag, violation);
}

/** collect the original expressions which are leafs in the nlexpr, thus which auxvar we use, and assign an index to them
 *
 * Also ensure that for every leaf there is an auxiliary variable.
 */
static
SCIP_RETCODE collectAuxVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expr */
   SCIP_HASHMAP*         origexpr2index,     /**< mapping from original expression to index */
   int*                  nindices            /**< number of indices */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);
   assert(origexpr2index != NULL);
   assert(nindices != NULL);

   if( nlexpr->children == NULL )
   {
      /* TODO skip for value expressions */

      assert(SCIPgetConsExprExprAuxVar(nlexpr->origexpr) != NULL);  /* TODO remove again when we allow to introduce auxvars */
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexpr->origexpr, NULL) );

      if( SCIPhashmapExists(origexpr2index, (void*)nlexpr->origexpr) )
      {
         nlexpr->varidx = SCIPhashmapGetImageInt(origexpr2index, (void*)nlexpr->origexpr);
      }
      else
      {
         nlexpr->varidx = (*nindices)++;
         SCIP_CALL( SCIPhashmapInsertInt(origexpr2index, (void*)nlexpr->origexpr, nlexpr->varidx) );
      }

      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( collectAuxVarExprs(scip, conshdlr, nlexpr->children[i], origexpr2index, nindices) );
   }

   return SCIP_OKAY;
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
   SCIP_CALL( nlhdlrExprCreate(scip, rootnlexpr, rootexpr, NULL, curv) );
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
            SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr, childcurv) );
         }
         else if( nlexpr == *rootnlexpr )
         {
            /* if there is no way to ensure curv for root expression, then we failed */
            nlhdlrExprFree(scip, rootnlexpr);
            break;
         }
         else
         {
            /* for now, also require that desired curvature can be achieved w.r.t. original variables, i.e., without introducing auxvars (TODO: remove again) */
            nlhdlrExprFree(scip, rootnlexpr);
            break;
         }
      }
      else if( nchildren > 0 )
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr, NULL) );
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
   SCIP_HASHMAP* origexpr2index;
   int i;

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

   /* make sure there are auxvars and collect them */
   SCIP_CALL( SCIPhashmapCreate(&origexpr2index, SCIPblkmem(scip), 10) );   /* TODO replace 10 by number of leafs, which we could count in constructExpr */
   (*nlhdlrexprdata)->nauxvars = 0;
   SCIP_CALL( collectAuxVarExprs(scip, conshdlr, nlexpr, origexpr2index, &(*nlhdlrexprdata)->nauxvars) );

   /* assemble auxvarexprs array */
   assert((*nlhdlrexprdata)->nauxvars > 0);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->auxvarexprs, (*nlhdlrexprdata)->nauxvars) );
   for( i = 0; i < SCIPhashmapGetNEntries(origexpr2index); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_CONSEXPR_EXPR* origexpr;
      int idx;

      entry = SCIPhashmapGetEntry(origexpr2index, i);
      if( entry == NULL )
         continue;

      origexpr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetOrigin(entry);
      assert(expr != NULL);

      idx = SCIPhashmapEntryGetImageInt(entry);
      assert(idx >= 0);
      assert(idx < (*nlhdlrexprdata)->nauxvars);

      (*nlhdlrexprdata)->auxvarexprs[idx] = origexpr;

      SCIPdebugMsg(scip, "auxvar %d: <%s>\n", idx, SCIPvarGetName(SCIPgetConsExprExprAuxVar(origexpr)));
   }

   SCIPhashmapFree(&origexpr2index);

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrfreeHdlrDataConvex)
{
   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(*nlhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataConvex)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->auxvarexprs, (*nlhdlrexprdata)->nauxvars);
   nlhdlrExprFree(scip, &(*nlhdlrexprdata)->nlexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** the detection assumes that the curvature information of the expression has been computed already */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
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

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* ignore sums */
   if( !nlhdlrdata->detectsum && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) ) /*lint !e506 !e774*/
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPgetConsExprExprNChildren(expr) == 0 )
      return SCIP_OKAY;

   /* TODO we are also interested in handling the concave side? (-> nlhdlr_vertexpolyhedral?) */

   if( !*enforcedbelow )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, &nlexpr, &depth, expr, SCIP_EXPRCURV_CONVEX) );
      assert(nlexpr == NULL || depth >= 1);
      if( depth < 1 )  /* TODO change back to <= 1, i.e. free if only immediate children, but no grand-daugthers */
      {
         nlhdlrExprFree(scip, &nlexpr);
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
         nlhdlrExprFree(scip, &nlexpr);
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

   SCIP_CALL( nlhdlrExprEval(scip, nlhdlrexprdata->nlexpr, sol) );
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
   SCIP_CALL( nlhdlrExprEval(scip, nlexpr, sol) );
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
   SCIP_CALL( nlhdlrExprGradientCut(scip, conshdlr, nlexpr, rowprep, success) );

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
   nlhdlrExprBranchscore(scip, nlhdlrexprdata->nlexpr, brscoretag, violation);

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
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectConvex, nlhdlrEvalAuxConvex, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrfreeHdlrDataConvex);
   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvex);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateConvex, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscoreConvex);

   return SCIP_OKAY;
}
