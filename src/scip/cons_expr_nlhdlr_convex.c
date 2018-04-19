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

#define SCIP_PRIVATE_ROWPREP

#include "scip/cons_quadratic.h" /* for SCIP_ROWPREP */
#include "scip/cons_expr_nlhdlr_convex.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "convex"
#define NLHDLR_DESC         "convex handler for expressions"
#define NLHDLR_PRIORITY     50

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
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to store nlhdlr expression data */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);

   /* compute a decent upper bound on the number of unique variable expressions */
   SCIP_CALL( SCIPgetConsExprExprNVars(scip, expr, &nvars) );
   assert(nvars > 0);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->varexprs, nvars) );
   (*nlhdlrexprdata)->varexprssize = nvars;

   /* collect all variable expressions that are contained in expr (the function also captures all variable expressions) */
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, expr, (*nlhdlrexprdata)->varexprs, &(*nlhdlrexprdata)->nvarexprs) );
   assert((*nlhdlrexprdata)->nvarexprs > 0);

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
   SCIP_EXPRCURV curvature;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;
   curvature = SCIPgetConsExprExprCurvature(expr);

   /* check whether expression is nonlinear, convex or concave, and is not handled by another nonlinear handler */
   if( curvature == SCIP_EXPRCURV_CONVEX && !*enforcedabove )
   {
      *enforcedabove = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
      *success = TRUE;

      SCIPdebugMsg(scip, "detected convexity of %p\n", (void*)expr);
   }
   else if( curvature == SCIP_EXPRCURV_CONCAVE && !*enforcedbelow )
   {
      *enforcedbelow = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
      *success = TRUE;

      SCIPdebugMsg(scip, "detected concavity of %p\n", (void*)expr);
   }

   /* store variable expressions into the expression data of the nonlinear handler */
   if( *success )
   {
      SCIP_CALL( createNlhdlrExprData(scip, expr, nlhdlrexprdata) );
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaConvex)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_SIDETYPE type;
   SCIP_Real constant;
   SCIP_Real violation;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varexprs != NULL);
   assert(nlhdlrexprdata->nvarexprs > 0);
   assert(result != NULL);
   assert(ncuts != NULL);
   assert(SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX || SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);

   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;
   type = SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT;

   /* compute gradient */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(expr) == SCIP_INVALID ) /*lint !e777*/
   {
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* compute g(x*) */
   constant = SCIPgetConsExprExprValue(expr);

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(constant)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", constant, (void*)expr);
      return SCIP_OKAY;
   }

   /* compute violation */
   if( type == SCIP_SIDETYPE_RIGHT )
      violation = constant - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - constant;

   /* expression is not violated -> skip */
   if( violation < minviolation )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "violation of %p = %.12f\n", (void*)expr, violation);

   /* compute gradient cut */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, type, FALSE) );

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
         goto CLEANUP;
      }

      /* - grad(g(x*))_i x*_i */
      constant -= derivative * val;

      /* grad(g(x*))_i x_i */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, derivative) );
   }

   /* add auxiliary variable */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprAuxVar(expr), -1.0) );

   /* add constant */
   SCIPaddRowprepConstant(rowprep, constant);

   /* check build cut and check violation */
   {
      SCIP_ROW* row;
      SCIP_Bool infeasible;
      SCIP_Bool success;

      SCIPmergeRowprepTerms(scip, rowprep);

      /* improve coefficients */
      SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, minviolation, NULL, &success) );

      if( !success )
         goto CLEANUP;

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, conshdlr) );
      SCIP_CALL( SCIPaddRow(scip, row, TRUE, &infeasible) );
      (*ncuts)++;

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

CLEANUP:
   SCIPfreeRowprep(scip, &rowprep);

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
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varexprs != NULL);
   assert(nlhdlrexprdata->nvarexprs > 0);
   assert(success != NULL);

   /* we separate only convex functions here, so there should be little use for branching
    * if violations are small or there are numerical issues, then we will not have generated a cut, though
    * in that case, we will still branch, that is, register branchscores for all depending var exprs
    */

   /* compute violation */
   if( SCIPgetConsExprExprValue(expr) == SCIP_INVALID )
      violation = SCIPinfinity(scip); /* evaluation error -> we should branch */
   else if( SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX  )
      violation = SCIPgetConsExprExprValue(expr) - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - SCIPgetConsExprExprValue(expr);

   *success = TRUE;

   /* if no violation, then no need for branching
    * TODO doing this only if violation > epsilon is correct, or better do this for any violation > 0?
    */
   if( !SCIPisPositive(scip, violation) )
      return SCIP_OKAY;

   for( i = 0; i < nlhdlrexprdata->nvarexprs; ++i )
   {
      assert(nlhdlrexprdata->varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(nlhdlrexprdata->varexprs[i]));

      SCIPaddConsExprExprBranchScore(scip, nlhdlrexprdata->varexprs[i], brscoretag, violation);
   }

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

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectConvex, NULL) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvex);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, nlhdlrSepaConvex, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscoreConvex);

   return SCIP_OKAY;
}
