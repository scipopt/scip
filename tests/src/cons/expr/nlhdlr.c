/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr.c
 * @brief  tests basic nonlinear handler methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_var.h"

#include "include/scip_test.h"

static SCIP* testscip;
static SCIP_VAR* x;
static SCIP_VAR* y;

struct SCIP_ConsExpr_NlHdlrData
{
   SCIP_Bool             initialized;        /**< whether handler has been initialized and not yet de-initialized */
};

/** compact storage for variables and coefficients in a bivariate quadratic term that is either convex or concave */
struct SCIP_ConsExpr_NlHdlrExprData
{
   SCIP_VAR*             varx;               /**< first variable */
   SCIP_VAR*             vary;               /**< second variable */
   SCIP_Real             xcoef;              /**< coefficient of first variable linear term */
   SCIP_Real             ycoef;              /**< coefficient of second variable linear term */
   SCIP_Real             xycoef;             /**< coefficient of bilinear term */
   SCIP_Real             xxcoef;             /**< coefficient of first variable square term */
   SCIP_Real             yycoef;             /**< coefficient of second variable square term */
   SCIP_Real             constant;           /**< constant term */
   SCIP_Bool             convex;             /**< whether convex or concave */
};

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   const char* input1 = "[expr] <test>: (<x>+<y>-0)^2 + (<y>-0)^2 <= 1.5;";
   const char* input2 = "[expr] <test>: (<x>+<y>-1)^2 + (<y>-1)^2 <= 1.0;";

   /* create scip with all plugins */
   SCIP_CALL( SCIPcreate(&testscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(testscip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(testscip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(testscip, &x, "x", 0.0, 5.0, -1.5, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(testscip, &y, "y", 0.0, 5.0, -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(testscip, x) );
   SCIP_CALL( SCIPaddVar(testscip, y) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(testscip, &consexpr, input1,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(testscip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(testscip, &consexpr) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(testscip, &consexpr, input2,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(testscip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(testscip, &consexpr) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(testscip, &x) );
   SCIP_CALL( SCIPreleaseVar(testscip, &y) );
   SCIP_CALL( SCIPfree(&testscip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(freeHdlrData)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrdata != NULL);
   assert(*nlhdlrdata != NULL);
   cr_assert(!(*nlhdlrdata)->initialized, "nlhdlr should not yet be initialized or be de-initialized when freeing");

   SCIPfreeMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(freeExprData)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);

   SCIPfreeMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRINIT(initHdlr)
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);

   nlhdlrdata = SCIPgetConsExprNlHdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   cr_assert(!nlhdlrdata->initialized, "nlhdlr cannot be initialized already");

   nlhdlrdata->initialized = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLREXIT(exitHldr)
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(nlhdlr != NULL);

   nlhdlrdata = SCIPgetConsExprNlHdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   cr_assert(nlhdlrdata->initialized, "nlhdlr must have been initialized");

   nlhdlrdata->initialized = FALSE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(detectHdlr)
{
   SCIP_CONSEXPR_EXPRHDLR* powhdlr;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSEXPR_NLHDLREXPRDATA exprdata;
   int c;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* only look at sum expressions */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      return SCIP_OKAY;

   powhdlr = SCIPfindConsExprExprHdlr(conshdlr, "pow");
   cr_assert(powhdlr != NULL, "pow hdlr should be present");

   BMSclearMemory(&exprdata);
   exprdata.constant = SCIPgetConsExprExprSumConstant(expr);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];
      assert(child != NULL);

      if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrVar(conshdlr) )
      {
         SCIP_VAR* var;

         var = SCIPgetConsExprExprVarVar(child);
         assert(var != NULL);

         if( var == exprdata.varx )
         {
            exprdata.xcoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( var == exprdata.vary )
         {
            exprdata.ycoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            exprdata.varx = var;
            assert(exprdata.xcoef == 0.0);
            exprdata.xcoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.vary == NULL )
         {
            exprdata.vary = var;
            assert(exprdata.ycoef == 0.0);
            exprdata.ycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else
         {
            /* more than two variables -> give up */
            return SCIP_OKAY;
         }
      }
      else if( SCIPgetConsExprExprHdlr(child) == powhdlr )
      {
         SCIP_VAR* var;

         if( SCIPgetConsExprExprPowExponent(child) != 2.0 )
            return SCIP_OKAY;  /* only allow for exponent 2 (square terms) */

         assert(SCIPgetConsExprExprNChildren(child) == 1);
         child = SCIPgetConsExprExprChildren(child)[0];

         if( SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrVar(conshdlr) )
            return SCIP_OKAY;  /* only allow for variable as base of power */

         var = SCIPgetConsExprExprVarVar(child);
         if( var == exprdata.varx )
         {
            exprdata.xxcoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( var == exprdata.vary )
         {
            exprdata.yycoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            assert(exprdata.xxcoef == 0.0);
            exprdata.varx = var;
            exprdata.xxcoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.vary == NULL )
         {
            assert(exprdata.yycoef == 0.0);
            exprdata.vary = var;
            exprdata.yycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else
         {
            /* more than two variables -> give up */
            return SCIP_OKAY;
         }
      }
      else if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;

         if( SCIPgetConsExprExprNChildren(child) != 2 )
            return SCIP_OKAY; /* only allow for bilinear term */

         if( SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(child)[0]) != SCIPgetConsExprExprHdlrVar(conshdlr) )
            return SCIP_OKAY; /* only allow for variables in bilinear term */
         if( SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(child)[1]) != SCIPgetConsExprExprHdlrVar(conshdlr) )
            return SCIP_OKAY; /* only allow for variables in bilinear term */

         var1 = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0]);
         var2 = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[1]);
         assert(var1 != var2);

         if( (var1 == exprdata.varx && var2 == exprdata.vary) || (var1 == exprdata.vary && var2 == exprdata.varx)  )
         {
            exprdata.xycoef += SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( ((var1 == exprdata.varx) || (var2 == exprdata.varx)) && exprdata.vary == NULL )
         {
            assert(exprdata.xycoef == 0.0);
            exprdata.vary = (var1 == exprdata.varx) ? var2 : var1;
            exprdata.xycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else if( exprdata.varx == NULL )
         {
            assert(exprdata.xycoef == 0.0);
            assert(exprdata.vary == NULL);
            exprdata.varx = var1;
            exprdata.vary = var2;
            exprdata.xycoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
         else
         {
            /* more than two variables */
            return SCIP_OKAY;
         }
      }
      else
      {
         /* unknown expression type */
         return SCIP_OKAY;
      }
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, " -> %gx%+gy%+gxy%+gx^2%+gy^2%+g (x=<%s>, y=<%s>)\n", exprdata.xcoef, exprdata.ycoef, exprdata.xycoef, exprdata.xxcoef, exprdata.yycoef, exprdata.constant, SCIPvarGetName(exprdata.varx), SCIPvarGetName(exprdata.vary));
#endif

   /* separable function is not of interest (for this unittest) */
   if( exprdata.xycoef == 0.0 )
      return SCIP_OKAY;

   /* check if convex or concave */
   if( exprdata.xxcoef >= 0.0 && exprdata.yycoef >= 0.0 && SCIPisGE(scip, 4.0 * exprdata.xxcoef * exprdata.yycoef, exprdata.xycoef*exprdata.xycoef) )
      exprdata.convex = TRUE;
   else if( exprdata.xxcoef <= 0.0 && exprdata.yycoef <= 0.0 && SCIPisGE(scip, 4.0 * exprdata.xxcoef * exprdata.yycoef, exprdata.xycoef*exprdata.xycoef) )
      exprdata.convex = FALSE;
   else
      return SCIP_OKAY; /* indefinite */

   *success = TRUE;
   SCIP_CALL( SCIPduplicateMemory(scip, nlhdlrexprdata, &exprdata) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(sepaHdlr)
{
   SCIP_VAR* auxvar;
   SCIP_ROW* cut = NULL;
   SCIP_Bool overestimate;
   SCIP_Bool infeasible;
   SCIP_Real violation;
   SCIP_Real xval;
   SCIP_Real yval;
   SCIP_Real xcoef;
   SCIP_Real ycoef;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrexprdata != NULL);

   *result = SCIP_DIDNOTFIND;
   *ncuts = 0;

   /* get auxiliary variable */
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);

   xval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->varx);
   yval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vary);

   /* compute the violation; this determines whether we need to over- or underestimate */
   violation = nlhdlrexprdata->constant;
   violation += nlhdlrexprdata->xcoef * xval + nlhdlrexprdata->ycoef * yval;
   violation += nlhdlrexprdata->xxcoef * xval*xval + nlhdlrexprdata->yycoef * yval*yval;
   violation += nlhdlrexprdata->xycoef * xval*yval;
   violation -= SCIPgetSolVal(scip, sol, auxvar);

   /* check if there is a violation */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   /* determine if we need to under- or overestimate */
   overestimate = SCIPisLT(scip, violation, 0.0);

   /* adjust the reference point */
   xval = MIN(MAX(xval, SCIPvarGetLbLocal(nlhdlrexprdata->varx)), SCIPvarGetUbLocal(nlhdlrexprdata->varx));
   yval = MIN(MAX(yval, SCIPvarGetLbLocal(nlhdlrexprdata->vary)), SCIPvarGetUbLocal(nlhdlrexprdata->vary));

   if( nlhdlrexprdata->convex == !overestimate )
   {
      /* convex side -> linearize
       * constant + xcoef * x + ycoef * y + xxcoef*xval*xval + yycoef*yval*yval + xycoef*xval*yval
       *   + (2*xxcoef*xval + xycoef*yval) * (x - xval) + (2*yycoef*yval + xycoef*xval) * (y - yval)
       * = (xcoef + 2*xxcoef*xval + xycoef*yval) * x + (ycoef * 2*yycoef*yval + xycoef*xval) * y
       *   - xxcoef*xval*xval - yycoef*yval*yval - xycoef*xval*yval + constant
       */
      xcoef = nlhdlrexprdata->xcoef + 2.0 * nlhdlrexprdata->xxcoef * xval + nlhdlrexprdata->xycoef * yval;
      ycoef = nlhdlrexprdata->ycoef + 2.0 * nlhdlrexprdata->yycoef * yval + nlhdlrexprdata->xycoef * xval;
      if( !overestimate )
      {
         lhs = -SCIPinfinity(scip);
         rhs = nlhdlrexprdata->xxcoef * xval * xval + nlhdlrexprdata->yycoef * yval * yval + nlhdlrexprdata->xycoef * xval * yval - nlhdlrexprdata->constant;
      }
      else
      {
         lhs = nlhdlrexprdata->xxcoef * xval * xval + nlhdlrexprdata->yycoef * yval * yval + nlhdlrexprdata->xycoef * xval * yval - nlhdlrexprdata->constant;
         rhs = SCIPinfinity(scip);
      }

      SCIP_CALL( SCIPcreateRowCons(scip, &cut, conshdlr, "testhdlrcut_cvx", 0, NULL, NULL, lhs, rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, nlhdlrexprdata->varx, xcoef) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, nlhdlrexprdata->vary, ycoef) );
      SCIP_CALL( SCIPaddVarToRow(scip, cut, auxvar, -1.0) );

      /* SCIP_CALL( SCIPprintRow(scip, cut, NULL) ); */
   }
   else
   {
      /* todo? */
   }

   if( cut == NULL )
      return SCIP_OKAY;

   if( nlhdlrexprdata->convex == !overestimate )
   {
      /* if separated on convex side, then violation of cut in sol should be same as violation of cons in sol */
      assert(SCIPisRelEQ(scip, -SCIPgetRowSolFeasibility(scip, cut, sol), violation));
   }

   /* check whether its violation and numerical properties are ok (and maybe improve) */
   SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cut, sol, minviolation) );

   if( cut == NULL )
      return SCIP_OKAY;

   assert(-SCIPgetRowSolFeasibility(scip, cut, sol) >= minviolation);

   /* add cut */
   SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
   *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   *ncuts += 1;

#ifdef SCIP_DEBUG
   if( *result == SCIP_CUTOFF )
   {
      SCIPdebugMsg(scip, "add cut makes node infeasible!\n");
   }
   else
   {
      SCIPdebugMsg(scip, "add cut with violation %e\n", -SCIPgetRowSolFeasibility(scip, cut, sol));
   }
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   /* release cut */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(copyHdlr)
{
   SCIP_CONSEXPR_NLHDLR* targetnlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourceconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   cr_assert(strcmp(SCIPgetConsExprNlHdlrName(sourcenlhdlr), "testhdlr") == 0, "source nlhdlr is not testhdlr");

   SCIP_CALL( SCIPallocClearMemory(testscip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeConsExprNlHdlrBasic(targetscip, targetconsexprhdlr, &targetnlhdlr,
      SCIPgetConsExprNlHdlrName(sourcenlhdlr), SCIPgetConsExprNlHdlrDesc(sourcenlhdlr), SCIPgetConsExprNlHdlrPriority(sourcenlhdlr), detectHdlr, nlhdlrdata) );
   SCIPsetConsExprNlHdlrFreeHdlrData(targetscip, targetnlhdlr, freeHdlrData);
   SCIPsetConsExprNlHdlrFreeExprData(targetscip, targetnlhdlr, freeExprData);
   SCIPsetConsExprNlHdlrCopyHdlr(targetscip, targetnlhdlr, copyHdlr);
   SCIPsetConsExprNlHdlrInitExit(targetscip, targetnlhdlr, initHdlr, exitHldr);
   SCIPsetConsExprNlHdlrSepa(targetscip, targetnlhdlr, sepaHdlr);

   return SCIP_OKAY;
}

Test(conshdlr, nlhdlr, .init = setup, .fini = teardown,
   .description = "test basic functionality of nonlinear handler of the cons_expr constraint handler."
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(testscip, "expr");
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPallocClearMemory(testscip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeConsExprNlHdlrBasic(testscip, conshdlr, &nlhdlr, "testhdlr", "tests nonlinear handler functionality", 0, detectHdlr, nlhdlrdata) );

   SCIPsetConsExprNlHdlrFreeHdlrData(testscip, nlhdlr, freeHdlrData);
   SCIPsetConsExprNlHdlrFreeExprData(testscip, nlhdlr, freeExprData);
   SCIPsetConsExprNlHdlrCopyHdlr(testscip, nlhdlr, copyHdlr);
   SCIPsetConsExprNlHdlrInitExit(testscip, nlhdlr, initHdlr, exitHldr);
   SCIPsetConsExprNlHdlrSepa(testscip, nlhdlr, sepaHdlr);

   SCIP_CALL( SCIPsetIntParam(testscip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   /* SCIP_CALL( SCIPsetRealParam(testscip, "limits/gap", 1e-6) ); */
/*
   SCIP_CALL( SCIPpresolve(testscip) );
   SCIP_CALL( SCIPprintOrigProblem(testscip, NULL, NULL, FALSE) );
   SCIP_CALL( SCIPprintTransProblem(testscip, NULL, NULL, FALSE) );
*/
   SCIP_CALL( SCIPsolve(testscip) );
   /* SCIP_CALL( SCIPprintBestSol(testscip, NULL, TRUE) ); */

   cr_assert(SCIPgetStatus(testscip) == SCIP_STATUS_OPTIMAL, "not solved to optimality");
   cr_assert(SCIPisFeasEQ(testscip, SCIPgetPrimalbound(testscip), -1.93649230212515), "optimal value not correct, expected -1.93649230212515, but got %.20g", SCIPgetPrimalbound(testscip));
   cr_assert(SCIPgetNNodes(testscip) <= 1, "convex NLP should be solved without branching, but took %" SCIP_LONGINT_FORMAT " nodes", SCIPgetNNodes(testscip));
}
