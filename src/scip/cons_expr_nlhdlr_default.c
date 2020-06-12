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

/**@file   cons_expr_nlhdlr_default.c
 * @brief  default nonlinear handler that calls expression handler methods
 * @author Stefan Vigerske
 *
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_default.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_iterator.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "default"
#define NLHDLR_DESC               "default handler for expressions"
#define NLHDLR_DETECTPRIORITY     0
#define NLHDLR_ENFOPRIORITY       0

/** evaluates an expression w.r.t. the values in the auxiliary variables */
static
SCIP_RETCODE evalExprInAux(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be evaluated */
   SCIP_Real*            val,                /**< buffer to store value of expression */
   SCIP_SOL*             sol                 /**< solution to be evaluated */
)
{
   SCIP_Real* childvals;
   SCIP_VAR* childvar;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(val != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) > 0);

   SCIP_CALL( SCIPallocBufferArray(scip, &childvals, SCIPgetConsExprExprNChildren(expr)) );

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      childvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[c]);
      /* there should be an auxiliary variable, because we created them in detect for every child if we said that we will separate
       * at the moment, EVALAUX should only be called for nlhdlr for which we said that we will separate
       * if that changes, then we should handle this here, e.g., via *val = SCIPgetConsExprExprValue(expr); break;
       */
      assert(childvar != NULL);

      childvals[c] = SCIPgetSolVal(scip, sol, childvar);
   }

   SCIP_CALL( SCIPevalConsExprExprHdlr(scip, expr, val, childvals, sol) );

   SCIPfreeBufferArray(scip, &childvals);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_Bool estimateusesactivity = FALSE;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcing != NULL);
   assert(participating != NULL);
   assert(nlhdlrexprdata != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);

   if( (*enforcing & SCIP_CONSEXPR_EXPRENFO_ACTIVITY) == 0 )
   {
      /* having reverseprop but no inteval is something that we don't support at the moment for simplicity */
      assert(!SCIPhasConsExprExprHdlrReverseProp(exprhdlr) || SCIPhasConsExprExprHdlrIntEval(exprhdlr));

      /* participate in inteval and/or reverseprop if that is not yet provided in enforcing and we have inteval */
      if( SCIPhasConsExprExprHdlrIntEval(exprhdlr) )
         *participating = SCIP_CONSEXPR_EXPRENFO_ACTIVITY;
   }

   /* participate in sepa if exprhdlr for expr has an estimate callback and sepa below or above is still missing */
   if( ((*enforcing & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) != SCIP_CONSEXPR_EXPRENFO_SEPABOTH) && SCIPhasConsExprExprHdlrEstimate(exprhdlr) )
   {
      /* communicate back that the nlhdlr will provide the separation on the currently missing sides */
      if( (*enforcing & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) == 0 )
         *participating |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;

      if( (*enforcing & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) == 0 )
         *participating |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   }

   if( !*participating )
      return SCIP_OKAY;

   /* since this is the default handler, we enforce where we participate */
   *enforcing |= *participating;

   /* increment activity usage counter and create auxiliary variables if necessary
    * if separating, first guess whether we will use activities in estimate
    * we assume that the exprhdlr will use activity on all children iff we are estimating on a nonconvex side
    * TODO it would be better to request this information directly from the exprhdlr than inferring it from curvature, but with the currently available exprhdlr that wouldn't make a difference
    */
   if( (*participating & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) != 0 )
   {
      SCIP_EXPRCURV* childcurv;
      SCIP_EXPRCURV requiredcurv;

      /* allocate memory to store the required curvature of the children (though we don't use it) */
      SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, SCIPgetConsExprExprNChildren(expr)) );

      /* check whether the expression is convex, if sepabelow, and whether the expression is concave, if sepaabove */
      requiredcurv  = ((*participating & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) ? SCIP_EXPRCURV_CONVEX : SCIP_EXPRCURV_UNKNOWN);
      requiredcurv |= ((*participating & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) ? SCIP_EXPRCURV_CONCAVE : SCIP_EXPRCURV_UNKNOWN);
      assert(requiredcurv != SCIP_EXPRCURV_UNKNOWN);  /* because we have sepabelow or sepaabove */

      SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, expr, requiredcurv, &estimateusesactivity, childcurv) );

      /* free memory */
      SCIPfreeBufferArray(scip, &childcurv);
   }

   /* indicate enforcement methods required in children:
    * - if separating, make sure that (auxiliary) variables exist
    * - if separation requires curvature, then increment activityusage count with usedforsepa == TRUE
    * - if activity computation, then increment activityusage count with usedforprop == TRUE
    */
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      /* todo skip this for value-expressions? would then need update in evalExprInAux, too */
      if( (*participating & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) != 0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], NULL) );
      }

      SCIP_CALL( SCIPincrementConsExprExprNActivityUses(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c],
         *participating & SCIP_CONSEXPR_EXPRENFO_ACTIVITY, estimateusesactivity) );
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxDefault)
{ /*lint --e{715}*/
   assert(expr != NULL);
   assert(auxvalue != NULL);

   SCIP_CALL( evalExprInAux(scip, expr, auxvalue, sol) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaDefault)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* call the separation initialization callback of the expression handler */
   SCIP_CALL( SCIPinitsepaConsExprExprHdlr(scip, conshdlr, cons, expr, overestimate, underestimate, infeasible) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateDefault)
{ /*lint --e{715}*/
   SCIP_Real constant;
   SCIP_Bool* branchcand = NULL;
   int nchildren;
   int c;
   SCIP_ROWPREP* rowprep;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(rowpreps != NULL);
   assert(success != NULL);

   *addedbranchscores = FALSE;

   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, TRUE) );

   nchildren = SCIPgetConsExprExprNChildren(expr);

   /* make sure enough space is available in rowprep arrays */
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, nchildren) );
   assert(rowprep->varssize >= nchildren);

   /* we need to pass a branchcand array to exprhdlr's estimate also if not asked to add branching scores */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchcand, nchildren) );
   for( c = 0; c < nchildren; ++c )
      branchcand[c] = TRUE;

   /* call the estimation callback of the expression handler */
   SCIP_CALL( SCIPestimateConsExprExprHdlr(scip, conshdlr, expr, sol, overestimate, targetvalue, rowprep->coefs, &constant, &rowprep->local, success, branchcand) );

   if( *success )
   {
      int i;

      /* add variables to rowprep */
      rowprep->nvars = SCIPgetConsExprExprNChildren(expr);
      for( i = 0; i < rowprep->nvars; ++i )
      {
         rowprep->vars[i] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[i]);
         assert(rowprep->vars[i] != NULL);
      }

      rowprep->side = -constant;

      SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, 0, rowprep) );

      (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%sestimate_%s%p_%s%d",
         overestimate ? "over" : "under",
         SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)),
         (void*)expr,
         sol != NULL ? "sol" : "lp",
         sol != NULL ? SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));
   }
   else
   {
      SCIPfreeRowprep(scip, &rowprep);
   }

   if( addbranchscores )
   {
      SCIP_Real violation;

#ifndef BRSCORE_ABSVIOL
      SCIP_CALL( SCIPgetConsExprExprRelAuxViolation(scip, conshdlr, expr, auxvalue, sol, &violation, NULL, NULL) );
#else
      SCIP_CALL( SCIPgetConsExprExprAbsAuxViolation(scip, conshdlr, expr, auxvalue, sol, &violation, NULL, NULL) );
#endif
      assert(violation > 0.0);  /* there should be a violation if we were called to enforce */

      if( nchildren == 1 )
      {
         if( branchcand[0] )
         {
            SCIP_CALL( SCIPaddConsExprExprsViolScore(scip, conshdlr, SCIPgetConsExprExprChildren(expr), 1, violation, sol, addedbranchscores) );
         }
      }
      else
      {
         SCIP_CONSEXPR_EXPR** exprs;
         int nexprs = 0;

         /* get list of those children that have the branchcand-flag set */
         SCIP_CALL( SCIPallocBufferArray(scip, &exprs, nchildren) );

         for( c = 0; c < nchildren; ++c )
            if( branchcand[c] )
               exprs[nexprs++] = SCIPgetConsExprExprChildren(expr)[c];

         SCIP_CALL( SCIPaddConsExprExprsViolScore(scip, conshdlr, exprs, nexprs, violation, sol, addedbranchscores) );

         SCIPfreeBufferArray(scip, &exprs);
      }

      if( *addedbranchscores )
      {
         /* count this branchscore as belonging to the exprhdlr, too
          * thus, it will be counted for the default nlhdlr, but also for this exprhdlr
          */
         SCIPincrementConsExprExprHdlrNBranchScore(SCIPgetConsExprExprHdlr(expr));
      }
   }

   SCIPfreeBufferArray(scip, &branchcand);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaDefault)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* call the separation deinitialization callback of the expression handler */
   SCIP_CALL( SCIPexitsepaConsExprExprHdlr(scip, expr) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalDefault)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* call the interval evaluation callback of the expression handler */
   SCIP_CALL( SCIPintevalConsExprExprHdlr(scip, expr, interval, intevalvar, global, intevalvardata) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropDefault)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* call the reverse propagation callback of the expression handler */
   SCIP_CALL( SCIPreversepropConsExprExprHdlr(scip, conshdlr, expr, reversepropqueue, infeasible, nreductions, force) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrDefault)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrDefault(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes default nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectDefault, nlhdlrEvalAuxDefault, NULL) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrDefault);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, nlhdlrInitSepaDefault, NULL, nlhdlrEstimateDefault, nlhdlrExitSepaDefault);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalDefault, nlhdlrReversepropDefault);

   return SCIP_OKAY;
}
