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
      assert(childvar != NULL); /* because we created auxvars in detect for every child */

      childvals[c] = SCIPgetSolVal(scip, sol, childvar);
   }

   SCIP_CALL( SCIPevalConsExprExprHdlr(scip, expr, val, childvals, sol) );

   SCIPfreeBufferArray(scip, &childvals);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectDefault)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRENFO_METHOD mymethods;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* we currently do not get active in presolve, the core will call the exprhdlr directly */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   mymethods = SCIP_CONSEXPR_EXPRENFO_NONE;

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   assert(exprhdlr != NULL);

   /* return interval evaluation possibility if exprhdlr for expr has a inteval callback and no one already provides (a good) inteval */
   if( SCIPhasConsExprExprHdlrIntEval(exprhdlr) && (*enforcemethods & SCIP_CONSEXPR_EXPRENFO_INTEVAL) == 0 )
   {
      mymethods |= SCIP_CONSEXPR_EXPRENFO_INTEVAL;
      *success = TRUE;
   }

   /* return reverse propagation possibility if exprhdlr for expr has a reverseprop callback and no one already provides (a good) reverseprop */
   if( SCIPhasConsExprExprHdlrReverseProp(exprhdlr) && (*enforcemethods & SCIP_CONSEXPR_EXPRENFO_REVERSEPROP) == 0 )
   {
      /* one could claim that reverse propagation is sufficient for enforcement, but separation is probably stronger
       * so, not setting enforcedbelow/above to TRUE here for now
       */
      mymethods |= SCIP_CONSEXPR_EXPRENFO_REVERSEPROP;
      *success = TRUE;
   }

   /* notify children if that we will need their activity for domain propagation */
   if( (mymethods & (SCIP_CONSEXPR_EXPRENFO_INTEVAL | SCIP_CONSEXPR_EXPRENFO_REVERSEPROP)) != 0 )
      for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
      {
         SCIP_CALL( SCIPincrementConsExprExprNActivityUses(scip, conshdlr,
            SCIPgetConsExprExprChildren(expr)[c], TRUE, FALSE) );
      }


   /* return sepa possibility if exprhdlr for expr has an estimate callback and enforcement is not ensured already */
   if( SCIPhasConsExprExprHdlrEstimate(exprhdlr) && (!*enforcedbelow || !*enforcedabove) )
   {
      /* make sure that an (auxiliary) variable exists for every child */
      for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
      {
         /* todo skip this for value-expressions? would then need update in evalExprInAux, too */
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], NULL) );
      }

      /* communicate back what the nlhdlr will do
       * - it will enforce via estimation/separation on those sides that are not enforced yet
       * - it needs to be called for this expression (success = TRUE)
       */
      if( !*enforcedbelow )
      {
         mymethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *enforcedbelow = TRUE;
         *success = TRUE;
      }

      if( !*enforcedabove )
      {
         mymethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *enforcedabove = TRUE;
         *success = TRUE;
      }
   }
#if 0 /* TODO branching method needs to distinguish whether we do separation (thus added auxvar) or only propagate (no auxvar) */
   else if( (!*enforcedbelow || !*enforcedabove) &&
      (mymethods & SCIP_CONSEXPR_EXPRENFO_INTEVAL) != 0 &&
      (mymethods & SCIP_CONSEXPR_EXPRENFO_REVERSEPROP) != 0 )
   {
      /* return branching score possibility if enforcement is not ensured yet, but we provide propagation,
       * since propagation and branching should be sufficient for enforcement, too
       */
      mymethods |= SCIP_CONSEXPR_EXPRENFO_BRANCHSCORE;
      *enforcedbelow = TRUE;
      *enforcedabove = TRUE;
      *success = TRUE;
   }
#endif

   /* it does not makes much sense to advertise a brscore callback if we do not also enforce via separation or propagation */

   if( *success )
   {
      SCIP_Bool sepabelow;
      SCIP_Bool sepaabove;

      /* remember in the nlhdlr exprdata (pointer) which methods we advertised */
      *nlhdlrexprdata = (SCIP_CONSEXPR_NLHDLREXPRDATA*)(size_t)mymethods;
      /* augment mymethods in enforcemethods */
      *enforcemethods |= mymethods;

      /* since this is the default handler, it should always enforce, even if none of the stronger enforcement methods (estimate/separate) are available
       * this allows to handle value-expressions, for example
       */
      *enforcedbelow = TRUE;
      *enforcedabove = TRUE;

      sepabelow = (mymethods & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) != 0;
      sepaabove = (mymethods & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) != 0;

      /* update ndomainuses counter for the children if under- or overestimators will be computed */
      if( sepabelow || sepaabove )
      {
         SCIP_EXPRCURV* childcurv;
         SCIP_Bool isconvex;
         SCIP_Bool isconcave;

         /* allocate memory to store the required curvature of the children */
         SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, SCIPgetConsExprExprNChildren(expr)) );

         /* check whether the expression is convex and concave
          *
          * TODO add a method that computes the curvature of an expression when the children are considered to be variables
          */
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, expr, SCIP_EXPRCURV_CONVEX, &isconvex, childcurv) );
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, expr, SCIP_EXPRCURV_CONCAVE, &isconcave, childcurv) );

         /* use the curvature to decide whether bounds on the children are used to refine under- or overestimates */
         if( (sepabelow && !isconvex) || (sepaabove && !isconcave) )
         {
            for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
            {
               SCIP_CALL( SCIPincrementConsExprExprNActivityUses(scip, conshdlr,
                  SCIPgetConsExprExprChildren(expr)[c], FALSE, TRUE) );
            }
         }

         /* free memory */
         SCIPfreeBufferArray(scip, &childcurv);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxDefault)
{ /*lint --e{715}*/
   assert(expr != NULL);
   assert(auxvalue != NULL);

   if( ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == 0 )
   {
      /* if we did not say that we separated, then we did not introduce auxvars
       * in that case, return the expression value, though it is a bit odd that we are still called
       */
      *auxvalue = SCIPgetConsExprExprValue(expr);

      return SCIP_OKAY;
   }

   SCIP_CALL( evalExprInAux(scip, expr, auxvalue, sol) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaDefault)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);

   /* if we will not separate, then don't call initsepa */
   if( ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == 0 )
      return SCIP_OKAY;

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

   /* if we did not say that we will separate on this side, then stand by it */
   if( !overestimate && ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   if(  overestimate && ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) == 0 )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

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

   /* if we have not separated, then don't call exitsepa */
   if( ((SCIP_CONSEXPR_EXPRENFO_METHOD)(size_t)nlhdlrexprdata & SCIP_CONSEXPR_EXPRENFO_SEPABOTH) == 0 )
      return SCIP_OKAY;

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
