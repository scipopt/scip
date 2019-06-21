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

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_CONSEXPR_EXPR**  varexprs;           /**< all dependent variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
   int                   varexprssize;       /**< size of varexprs array */
};

/*
 * static methods
 */

/** creates nonlinear handler expression data structure */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to store nlhdlr expression data */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);

   /* compute the number of unique variable expressions */
   SCIP_CALL( SCIPgetConsExprExprNVars(scip, conshdlr, expr, &nvars) );
   assert(nvars > 0);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->varexprs, nvars) );
   (*nlhdlrexprdata)->varexprssize = nvars;

   /* collect all variable expressions that are contained in expr (the function also captures all variable expressions) */
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, (*nlhdlrexprdata)->varexprs, &(*nlhdlrexprdata)->nvarexprs) );
   assert((*nlhdlrexprdata)->nvarexprs > 0);
   assert((*nlhdlrexprdata)->nvarexprs == nvars);

   return SCIP_OKAY;
}

/** frees nonlinear handler expression data structure */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to free nlhdlr expression data */
   )
{
   int i;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   /* release variable expressions */
   for( i = 0; i < (*nlhdlrexprdata)->nvarexprs; ++i )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*nlhdlrexprdata)->varexprs[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->varexprs, (*nlhdlrexprdata)->varexprssize);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

static
SCIP_RETCODE enforceCurvature(
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
   SCIP_Bool success2;

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

   /* TODO we are also interested in handling the concave side */

   if( !*enforcedbelow )
   {
      SCIP_CALL( enforceCurvature(scip, conshdlr, expr, SCIP_EXPRCURV_CONVEX, &success2) );

      if( success2 )
      {
         *enforcedbelow = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr <= auxvar\n", (void*)expr);
      }
   }

   if( !*enforcedabove )
   {
      SCIP_CALL( enforceCurvature(scip, conshdlr, expr, SCIP_EXPRCURV_CONCAVE, &success2) );

      if( success2 )
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
      SCIP_CALL( createNlhdlrExprData(scip, conshdlr, expr, nlhdlrexprdata) );
   }

   return SCIP_OKAY;
}

/** auxiliary evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxConvex)
{ /*lint --e{715}*/
   assert(auxvalue != NULL);
   assert(expr != NULL);

   /* currently this nlhdlr does not introduce auxiliary variables,
    * so we can return the value of the expression in the original variables
    */
   *auxvalue = SCIPgetConsExprExprValue(expr);

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateConvex)
{ /*lint --e{715}*/
   SCIP_Real constant;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varexprs != NULL);
   assert(nlhdlrexprdata->nvarexprs > 0);
   assert(SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX || SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* if estimating on non-convex side, then do nothing
    * TODO we are vertex-polyhedral and so should compute something
    */
   if( ( overestimate && SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONCAVE) )
      return SCIP_OKAY;

   /* compute gradient */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(expr) == SCIP_INVALID ) /*lint !e777*/
   {
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* get g(x*) */
   constant = SCIPgetConsExprExprValue(expr);
   assert(auxvalue == constant); /* given value (originally from nlhdlrEvalAuxConvex) should coincide with expression value */  /*lint !e777*/

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(constant)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", constant, (void*)expr);
      return SCIP_OKAY;
   }

   /* compute gradient cut */
   rowprep->local = FALSE;
   for( i = 0; i < nlhdlrexprdata->nvarexprs; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real derivative;
      SCIP_Real val;

      assert(nlhdlrexprdata->varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(nlhdlrexprdata->varexprs[i]));

      /* get the variable of the variable expression */
      var = SCIPgetConsExprExprVarVar(nlhdlrexprdata->varexprs[i]);
      assert(var != NULL);

      /* get solution value */
      val = SCIPgetSolVal(scip, sol, var);

      /* avoid overhead of SCIPgetConsExprExprPartialDiff by accessing the derivative directly */
      derivative = SCIPgetConsExprExprDerivative(nlhdlrexprdata->varexprs[i]);
      assert(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var) == derivative); /*lint !e777*/

      /* evaluation error or too large values -> skip */
      if( SCIPisInfinity(scip, REALABS(derivative * val)) )
      {
         SCIPdebugMsg(scip, "evaluation error / too large values (%g %g) for %s in %p\n", derivative, val,
            SCIPvarGetName(var), (void*)expr);
         return SCIP_OKAY;
      }

      /* - grad(g(x*))_i x*_i */
      constant -= derivative * val;

      /* grad(g(x*))_i x_i */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, derivative) );
   }

   /* add constant */
   SCIPaddRowprepConstant(rowprep, constant);

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%sestimate_convex%p_%s%d",
      overestimate ? "over" : "under",
      (void*)expr,
      sol != NULL ? "sol" : "lp",
      sol != NULL ? SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));

   *success = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscoreConvex)
{ /*lint --e{715}*/
   SCIP_Real violation;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX || SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(auxvalue == SCIPgetConsExprExprValue(expr)); /* given auxvalue should have been computed by nlhdlrEvalAuxConvex */  /*lint !e777*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varexprs != NULL);
   assert(nlhdlrexprdata->nvarexprs > 0);
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
   for( i = 0; i < nlhdlrexprdata->nvarexprs; ++i )
   {
      assert(nlhdlrexprdata->varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(nlhdlrexprdata->varexprs[i]));

      SCIPaddConsExprExprBranchScore(scip, nlhdlrexprdata->varexprs[i], brscoretag, violation);
   }

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
