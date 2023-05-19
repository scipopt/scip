/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_convex.c
 * @brief  tests convex nonlinear handler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"
#include "scip/nlpi_ipopt.h" /* to check whether LAPACK is around */

#define NLHDLR_CONVEX_UNITTEST
#include "scip/nlhdlr_convex.c"


/*
 * TEST
 */

static SCIP* scip;
static SCIP_VAR* x_1;
static SCIP_VAR* x_2;
static SCIP_VAR* x_3;

static SCIP_CONSHDLR* conshdlr;
static SCIP_NLHDLR* nlhdlr = NULL;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);

   /* get nlhdlr */
   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, CONVEX_NLHDLR_NAME);
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE);
   SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);

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

   SCIP_CALL( SCIPclearCuts(scip) );

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
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* oexpr;
   SCIP_EXPR* expr;
   SCIP_Bool changed;
   SCIP_Bool infeas;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_CONS* cons;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &oexpr, (char*)exprstr, NULL, NULL, NULL) );
   if( simplify )
   {
      SCIP_CALL( SCIPsimplifyExpr(scip, oexpr, &expr, &changed, &infeas, NULL, NULL) );
      cr_expect(!infeas);
      SCIP_CALL( SCIPreleaseExpr(scip, &oexpr) );
   }
   else
      expr = oexpr;
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, 0.0, 0.0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIPprintCons(scip, cons, NULL);
   SCIPinfoMessage(scip, NULL, " and %s\n", SCIPexprcurvGetName(exprrootcurv));

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectConvex(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(enforcing, participating);
   if( (exprrootcurv & SCIP_EXPRCURV_CONVEX) != 0 )
   {
      cr_expect(enforcing & SCIP_NLHDLR_METHOD_SEPABELOW);
   }

   if( (exprrootcurv & SCIP_EXPRCURV_CONCAVE) != 0 )
   {
      cr_expect(enforcing & SCIP_NLHDLR_METHOD_SEPAABOVE);
   }

   if( participating != SCIP_NLHDLR_METHOD_NONE )
   {
      cr_assert_not_null(nlhdlrexprdata);
      SCIP_CALL( nlhdlrfreeExprDataConvexConcave(scip, nlhdlr, expr, &nlhdlrexprdata) );
   }

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

/* tests detection of convex/concave subexpressions */
Test(nlhdlrconvex, detect, .init = setup, .fini = teardown)
{
   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/convex/extendedform", FALSE) );

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
   if( SCIPisIpoptAvailableIpopt() )
      detect("<x1>^2+2*<x1>*<x2>+<x2>^2", SCIP_EXPRCURV_CONVEX, TRUE);

   /* assumeconvex */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/nonlinear/assumeconvex", TRUE) );
   detect("<x1>*<x2>-<x3>^2", SCIP_EXPRCURV_CONVEX, FALSE);
   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/convex/handletrivial", TRUE) );
   detect("<x1>*<x2>", SCIP_EXPRCURV_CONVEX, FALSE);
}

/** test detection for block-decomposable quadratic */
Test(nlhdlrconvex, detectquad, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* oexpr;
   SCIP_EXPR* expr;
   SCIP_Bool changed;
   SCIP_Bool infeas;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_CONS* cons;
   SCIP_Bool isquadratic;
   int nquadexprs;
   int i;

   /* need Lapack for convexity check to succeed at the moment */
   if( !SCIPisIpoptAvailableIpopt() )
      return;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &oexpr, "<x1>^2+2*<x1>*<x2>+<x2>^2+5*<x3>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, oexpr, &expr, &changed, &infeas, NULL, NULL) );
   cr_expect(!infeas);
   SCIP_CALL( SCIPreleaseExpr(scip, &oexpr) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, 0.0, 0.0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIPprintCons(scip, cons, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectConvex(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );

   cr_expect_eq(enforcing, participating);
   cr_expect(enforcing == SCIP_NLHDLR_METHOD_SEPABELOW);
   cr_assert_not_null(nlhdlrexprdata);

   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );
   cr_expect(isquadratic);
   SCIPexprGetQuadraticData(expr, NULL, NULL, NULL, NULL, &nquadexprs, NULL, NULL, NULL);
   cr_expect(nquadexprs == 3);

   for( i = 0; i < nquadexprs; ++i )
   {
      SCIP_EXPR* varexpr;
      SCIP_EXPR* sqrexpr;
      SCIP_VAR* var;

      SCIPexprGetQuadraticQuadTerm(expr, i, &varexpr, NULL, NULL, NULL, NULL, &sqrexpr);
      cr_assert(SCIPisExprVar(scip, varexpr));
      cr_assert_not_null(sqrexpr);

      var = SCIPgetVarExprVar(varexpr);
      cr_expect(var != x_3 || SCIPgetExprNAuxvarUsesNonlinear(sqrexpr) == 1);
      cr_expect(var == x_3 || SCIPgetExprNAuxvarUsesNonlinear(sqrexpr) == 0);
   }

   SCIP_CALL( nlhdlrfreeExprDataConvexConcave(scip, nlhdlr, expr, &nlhdlrexprdata) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
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
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* oexpr;
   SCIP_EXPR* expr;
   SCIP_Real auxvalue;
   SCIP_Real targetvalue;
   SCIP_ROWPREP* rowprep;
   SCIP_Bool changed;
   SCIP_Bool infeas;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
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
   SCIP_CALL( SCIPparseExpr(scip, &oexpr, (char*)exprstr, NULL, NULL, NULL) );
   if( simplify )
   {
      SCIP_CALL( SCIPsimplifyExpr(scip, oexpr, &expr, &changed, &infeas, NULL, NULL) );
      cr_expect(!infeas);
      SCIP_CALL( SCIPreleaseExpr(scip, &oexpr) );
   }
   else
      expr = oexpr;
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, 0.0, 0.0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);

   SCIPprintCons(scip, cons, NULL);
   SCIPinfoMessage(scip, NULL, " at x1=%g x2=%g x3=%g\n",
      SCIPgetSolVal(scip, sol, x_1), SCIPgetSolVal(scip, sol, x_2), SCIPgetSolVal(scip, sol, x_3));

   /* detect */
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );

   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectConvex(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );
   cr_expect_eq(enforcing, participating);
   cr_expect((participating & SCIP_NLHDLR_METHOD_SEPABOTH) != 0);
   enforceabove = (participating & SCIP_NLHDLR_METHOD_SEPAABOVE) != 0;

   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeas) );
   SCIP_CALL( nlhdlrInitSepaConvex(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, TRUE, TRUE, &infeas) );

   /* estimate on enforced side */
   targetvalue = enforceabove ? SCIPinfinity(scip) : -SCIPinfinity(scip);
   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );
   SCIP_CALL( nlhdlrEvalAuxConvexConcave(scip, nlhdlr, expr, nlhdlrexprdata, &auxvalue, sol) );
   SCIP_CALL( nlhdlrEstimateConvex(scip, conshdlr, nlhdlr, expr, nlhdlrexprdata, sol, auxvalue, enforceabove,
         targetvalue, FALSE, rowpreps, &success, &addedbranchscores) );

   cr_assert(success);
   cr_expect(SCIPgetPtrarrayMinIdx(scip, rowpreps) == 0);
   cr_expect(SCIPgetPtrarrayMaxIdx(scip, rowpreps) == 0);
   rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps, 0);
   cr_assert(!SCIProwprepIsLocal(rowprep));  /* nlhdlr should have set rowprep->local to FALSE */

   SCIPmergeRowprepTerms(scip, rowprep);
   constant = -SCIProwprepGetSide(rowprep);
   for( i = 0; i < SCIProwprepGetNVars(rowprep); ++i )
   {
      if( SCIProwprepGetVars(rowprep)[i] == x_1 )
         x1coef = SCIProwprepGetCoefs(rowprep)[i];
      else if( SCIProwprepGetVars(rowprep)[i] == x_2 )
         x2coef = SCIProwprepGetCoefs(rowprep)[i];
      else if( SCIProwprepGetVars(rowprep)[i] == x_3 )
         x3coef = SCIProwprepGetCoefs(rowprep)[i];
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

   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/convex/extendedform", FALSE) );

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
