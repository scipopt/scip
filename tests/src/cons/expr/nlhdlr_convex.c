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

/**@file   nlhdlr_convex.c
 * @brief  tests convex nonlinear handler
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr.c"

#define NLHDLR_CONVEX_UNITTEST
#include "scip/cons_expr_nlhdlr_convex.c"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x_1;
static SCIP_VAR* x_2;
static SCIP_VAR* x_3;

static SCIP_CONSHDLR* conshdlr;
static SCIP_CONSEXPR_NLHDLR* nlhdlr = NULL;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   int h;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get perspective handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* get nlhdlr */
   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
   {
      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), CONVEX_NLHDLR_NAME) == 0 )
      {
         nlhdlr = conshdlrdata->nlhdlrs[h];
         break;
      }
   }

   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x_1, "x1", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x_2, "x2", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x_3, "x3", 1.0, 4.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x_1) );
   SCIP_CALL( SCIPaddVar(scip, x_2) );
   SCIP_CALL( SCIPaddVar(scip, x_3) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   nlhdlrExitConvex(scip, nlhdlr);

   SCIP_CALL( SCIPreleaseVar(scip, &x_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x_2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x_3) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/** given a string for f(x) and its curvature, run nlhdlr_convex detect on f(x) = 0 and see whether that gives correct flags */
static
SCIP_RETCODE detect(
   const char*           exprstr,
   SCIP_EXPRCURV         exprrootcurv,
   SCIP_Bool             simplify
)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* oexpr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_INTERVAL activity;
   SCIP_Bool changed;
   SCIP_Bool infeas;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_CONS* cons;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)exprstr, NULL, &oexpr) );
   if( simplify )
   {
      SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, oexpr, &expr, &changed, &infeas) );
      cr_expect(!infeas);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &oexpr) );
   }
   else
      expr = oexpr;
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, 0.0, 0.0)  );

   SCIPprintCons(scip, cons, NULL);
   SCIPinfoMessage(scip, NULL, " and %s\n", SCIPexprcurvGetName(exprrootcurv));

   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &activity, FALSE) );

   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_NONE;
   enforcebelow = FALSE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( SCIPdetectConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   if( (exprrootcurv & SCIP_EXPRCURV_CONVEX) != 0 )
   {
      cr_expect(success);
      cr_expect(enforcebelow);
      cr_expect((provided & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) != 0);
   }

   if( (exprrootcurv & SCIP_EXPRCURV_CONCAVE) != 0 )
   {
      cr_expect(success);
      cr_expect(enforceabove);
      cr_expect((provided & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) != 0);
   }

   if( success )
   {
      cr_assert_not_null(nlhdlrexprdata);
      SCIP_CALL( nlhdlrfreeExprDataConvexConcave(scip, nlhdlr, expr, &nlhdlrexprdata) );
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/* tests detection of convex/concave subexpressions */
Test(nlhdlrconvex, detect, .init = setup, .fini = teardown)
{
   detect("exp(exp(<x1>))", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(exp(log(<x1>)))", SCIP_EXPRCURV_CONVEX, FALSE);

   /* use signomial rule */
   detect("<x1>^0.5*<x2>^0.5)", SCIP_EXPRCURV_CONCAVE, FALSE);
   detect("<x1>^0.5*(exp(<x2>))^0.5", SCIP_EXPRCURV_CONCAVE, FALSE);  /* here an auxvar will be introduced for exp() */
   detect("<x1>^0.5*(log(<x2>))^0.5", SCIP_EXPRCURV_CONCAVE, FALSE);  /* here no auxvar will be introduced for exp() */
   detect("<x1>^(-0.5)*<x2>^(-0.5))", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>^3*<x2>^(-1)*<x3>^(-1)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>^3*<x2>^(-1)*log(<x3>)^(-1)", SCIP_EXPRCURV_CONVEX, FALSE);

   /* product-composition */
   detect("<x1>*exp(<x1>)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>*exp(<x1>+1.0)", SCIP_EXPRCURV_CONVEX, TRUE);
   detect("<x1>*exp(2*<x1>)", SCIP_EXPRCURV_CONVEX, TRUE);
   detect("(<x1>+<x2>)*exp(<x1>+<x2>)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(<x1>)*<x1>", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(<x1>^2)*<x1>^2", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("exp(2*<x1>^2)*<x1>^2", SCIP_EXPRCURV_CONVEX, TRUE);
   detect("log(4-<x1>)*<x1>", SCIP_EXPRCURV_CONCAVE, TRUE);   /* similar to arki0017 */

   /* quadratic */
   detect("exp(<x1>^2+<x2>^2)", SCIP_EXPRCURV_CONVEX, FALSE);
   detect("<x1>^2+2*<x1>*<x2>+<x2>^2", SCIP_EXPRCURV_CONVEX, TRUE);
}

/** given a string for f(x) and its curvature, run nlhdlr_convex detect on f(x) = 0 and estimate on the enforced side and check whether the estimator is as expected */
static
SCIP_RETCODE estimate(
   const char*           exprstr,
   SCIP_Bool             simplify,
   SCIP_SOL*             sol,
   SCIP_Real             x1coef_expected,
   SCIP_Real             x2coef_expected,
   SCIP_Real             x3coef_expected,
   SCIP_Real             constant_expected
)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* oexpr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_INTERVAL activity;
   SCIP_Real auxvalue;
   SCIP_Real targetvalue;
   SCIP_ROWPREP* rowprep;
   SCIP_Bool changed;
   SCIP_Bool infeas;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_Bool addedbranchscores;
   SCIP_CONS* cons;
   SCIP_Real x1coef = 0.0;
   SCIP_Real x2coef = 0.0;
   SCIP_Real x3coef = 0.0;
   SCIP_Real constant;
   int i;
   SCIP_PTRARRAY* rowpreps;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)exprstr, NULL, &oexpr) );
   if( simplify )
   {
      SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, oexpr, &expr, &changed, &infeas) );
      cr_expect(!infeas);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &oexpr) );
   }
   else
      expr = oexpr;
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, 0.0, 0.0)  );

   SCIPprintCons(scip, cons, NULL);
   SCIPinfoMessage(scip, NULL, " at x1=%g x2=%g x3=%g\n",
      SCIPgetSolVal(scip, sol, x_1), SCIPgetSolVal(scip, sol, x_2), SCIPgetSolVal(scip, sol, x_3));

   SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, expr, &activity, FALSE) );

   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_NONE;
   enforcebelow = FALSE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( SCIPdetectConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   cr_expect(success);
   cr_expect(!enforcebelow || (provided & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) != 0);
   cr_expect(!enforceabove || (provided & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) != 0);
   cr_expect(enforcebelow || enforceabove);

   /* estimate on enforced side */
   targetvalue = enforceabove ? SCIPinfinity(scip) : -SCIPinfinity(scip);
   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr, expr, nlhdlrexprdata, &auxvalue, sol) );
   SCIP_CALL( SCIPestimateConsExprNlhdlr(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, enforceabove,
         targetvalue, rowpreps, &success, FALSE, &addedbranchscores) );

   cr_assert(success);
   cr_expect(SCIPgetPtrarrayMinIdx(scip, rowpreps) == 0);
   cr_expect(SCIPgetPtrarrayMaxIdx(scip, rowpreps) == 0);
   rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, 0);
   cr_assert(!rowprep->local);  /* nlhdlr should have set rowprep->local to FALSE */

   SCIPmergeRowprepTerms(scip, rowprep);
   constant = -rowprep->side;
   for( i = 0; i < rowprep->nvars; ++i )
   {
      if( rowprep->vars[i] == x_1 )
         x1coef = rowprep->coefs[i];
      else if( rowprep->vars[i] == x_2 )
         x2coef = rowprep->coefs[i];
      else if( rowprep->vars[i] == x_3 )
         x3coef = rowprep->coefs[i];
      else
      {
         SCIPerrorMessage("unexpected variable in rowprep");
         cr_assert(0);
      }
   }

   cr_expect_float_eq(x1coef, x1coef_expected, 1e-9, "x1 coef wrong. Expected %g, but got %g", x1coef_expected, x1coef);
   cr_expect_float_eq(x2coef, x2coef_expected, 1e-9, "x2 coef wrong. Expected %g, but got %g", x2coef_expected, x2coef);
   cr_expect_float_eq(x3coef, x3coef_expected, 1e-9, "x3 coef wrong. Expected %g, but got %g", x3coef_expected, x3coef);
   cr_expect_float_eq(constant, constant_expected, 1e-9, "Constant wrong. Expected %g, but got %g", constant_expected, constant);

   SCIPfreeRowprep(scip, &rowprep);

   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   cr_assert_not_null(nlhdlrexprdata);
   SCIP_CALL( nlhdlrfreeExprDataConvexConcave(scip, nlhdlr, expr, &nlhdlrexprdata) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

static
void gradest_univariate(
   SCIP_Real x,
   SCIP_Real fx,
   SCIP_Real fp,
   SCIP_Real* xcoef,
   SCIP_Real* constant
   )
{
   *xcoef = fp;
   *constant = fx - fp * x;
}

static
void gradest_bivariate(
   SCIP_Real x,
   SCIP_Real y,
   SCIP_Real fxy,
   SCIP_Real fpx,
   SCIP_Real fpy,
   SCIP_Real* xcoef,
   SCIP_Real* ycoef,
   SCIP_Real* constant
   )
{
   *xcoef = fpx;
   *ycoef = fpy;
   *constant = fxy - fpx * x - fpy * y;
}

static
void secantest(
   SCIP_Real left,
   SCIP_Real fleft,
   SCIP_Real right,
   SCIP_Real fright,
   SCIP_Real* xcoef,
   SCIP_Real* constant
   )
{
   *xcoef = (fright - fleft) / (right - left);
   *constant = fleft - *xcoef * left;
}

/* tests detection of convex/concave subexpressions */
Test(nlhdlrconvex, estimate, .init = setup, .fini = teardown)
{
   SCIP_Real x1coef = 0.0;
   SCIP_Real x2coef = 0.0;
   SCIP_Real x3coef = 0.0;
   SCIP_Real constant = 0.0;
   SCIP_SOL* sol;

   SCIPcreateSol(scip, &sol, NULL);

   /* convex exp(exp(x)) at x=2 */
   SCIPsetSolVal(scip, sol, x_1, 2.0);
   gradest_univariate(2.0, exp(exp(2.0)), exp(2.0)*exp(exp(2.0)), &x1coef, &constant);
   estimate("exp(exp(<x1>))", FALSE, sol, x1coef, 0.0, 0.0, constant);

   /* concave sqrt(x1*x2) at x1=4, x2=9 */
   SCIPsetSolVal(scip, sol, x_1, 4.0);
   SCIPsetSolVal(scip, sol, x_2, 9.0);
   gradest_bivariate(4.0, 9.0, 2.0*3.0, 0.5*3.0/2.0, 0.5*2.0/3.0, &x1coef, &x2coef, &constant);
   estimate("<x1>^0.5*<x2>^0.5)", FALSE, sol, x1coef, x2coef, 0.0, constant);

   /* convex x3*exp(x3) with x3 being integer */
   SCIPsetSolVal(scip, sol, x_3, 2.5);
   secantest(2.0, 2.0*exp(2.0), 3.0, 3.0*exp(3.0), &x3coef, &constant);
   estimate("<x3>*exp(<x3>)", FALSE, sol, 0.0, 0.0, x3coef, constant);

   SCIPfreeSol(scip, &sol);
}
