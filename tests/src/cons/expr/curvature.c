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

/**@file   curvature.c
 * @brief  tests curvature expression handler callbacks
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"

#include "include/scip_test.h"


static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 1.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -4.0, -3.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 5.0, 6.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, w) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(curvature, .init = setup, .fini = teardown);

/**< auxiliary function for creating an expression and checking its curvature */
static
SCIP_RETCODE checkCurvature(
   const char*          input,              /**< input creating an expression */
   const char*          exprhdlrname,       /**< target expression handler name */
   SCIP_EXPRCURV        expectedcur         /**< expected curvature */
   )
{
   SCIP_CONSEXPR_EXPR* origexpr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* create and print expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &origexpr), SCIP_OKAY);

   /* simplify expression */
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, origexpr, &expr, &changed, &infeasible) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &origexpr) );

   /* print simplified expression */
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* check name of the corresponding expression handler */
   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   cr_assert(exprhdlr != NULL);
   cr_expect(strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), exprhdlrname) == 0, "expect nonlinear handler %s, got %s\n",
      exprhdlrname, SCIPgetConsExprExprHdlrName(exprhdlr));

   /* compute curvature */
   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );

   /* check curvature */
   cr_expect(SCIPgetConsExprExprCurvature(expr) == expectedcur, "expect %d, got %d", expectedcur, SCIPgetConsExprExprCurvature(expr));

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/* check for absolute expression */
Test(curvature, absolute)
{
   SCIP_CALL( checkCurvature("abs(<x>[C])", "abs", SCIP_EXPRCURV_CONVEX) );
}

/* check for cosine expression */
Test(curvature, cosine)
{
   SCIP_CALL(checkCurvature("cos(<x>[C])", "cos", SCIP_EXPRCURV_UNKNOWN));
   SCIP_CALL(checkCurvature("cos(<y>[C])", "cos", SCIP_EXPRCURV_UNKNOWN));
   SCIP_CALL(checkCurvature("cos(<z>[C])", "cos", SCIP_EXPRCURV_CONVEX));
   SCIP_CALL(checkCurvature("cos(<w>[C])", "cos", SCIP_EXPRCURV_CONCAVE));
}

/* check for exponential expression */
Test(curvature, exponential)
{
   SCIP_CALL( checkCurvature("exp(<x>[C])", "exp", SCIP_EXPRCURV_CONVEX) );
}

/* check for logarithm expression */
Test(curvature, logarithm)
{
   SCIP_CALL( checkCurvature("log(<x>[C])", "log", SCIP_EXPRCURV_CONCAVE) );
}

/* check for power expression */
Test(curvature, power)
{
   SCIP_CALL( checkCurvature("(<x>[C])^2", "pow", SCIP_EXPRCURV_CONVEX) );

   /* 0 is contained in the interior of x -> neither convex nor concave */
   SCIP_CALL( checkCurvature("(<x>[C])^(-1)", "pow", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("(<x>[C])^(-2)", "pow", SCIP_EXPRCURV_UNKNOWN) );

   /* check cases for y > 0 */
   SCIP_CALL( checkCurvature("(<y>[C])^(-1)", "pow", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("(<y>[C])^(-2.5)", "pow", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("(<y>[C])^(0.5)", "pow", SCIP_EXPRCURV_CONCAVE) );

   /* check cases for z < 0 */
   SCIP_CALL( checkCurvature("(<z>[C])^(-1)", "pow", SCIP_EXPRCURV_CONCAVE) );
   SCIP_CALL( checkCurvature("(<z>[C])^(-2)", "pow", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("(<z>[C])^(-3)", "pow", SCIP_EXPRCURV_CONCAVE) );
}

/* check for signpower expression */
Test(curvature, signpower)
{
   /* 0 is contained in the interior of x -> neither convex nor concave */
   SCIP_CALL( checkCurvature("signpower(<x>[C],2)", "signpower", SCIP_EXPRCURV_UNKNOWN) );

   /* check cases for y > 0 */
   SCIP_CALL( checkCurvature("signpower(<y>[C],2)", "signpower", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("signpower(<y>[C],2.5)", "signpower", SCIP_EXPRCURV_CONVEX) );

   /* check cases for z < 0 */
   SCIP_CALL( checkCurvature("signpower(<z>[C],2)", "signpower", SCIP_EXPRCURV_CONCAVE) );
   SCIP_CALL( checkCurvature("signpower(<z>[C],2.5)", "signpower", SCIP_EXPRCURV_CONCAVE) );
}

/* check for product expression */
Test(curvature, product)
{
   SCIP_CALL( checkCurvature("(<x>[C] * <y>[C])", "prod", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("(<x>[C] * <y>[C] * <z>[C])", "prod", SCIP_EXPRCURV_UNKNOWN) );
}

/* check for sine expression */
Test(curvature, sine)
{
   SCIP_CALL( checkCurvature("sin(<x>[C])", "sin", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("sin(<y>[C])", "sin", SCIP_EXPRCURV_CONCAVE) );
   SCIP_CALL( checkCurvature("sin(<z>[C])", "sin", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("sin(<w>[C])", "sin", SCIP_EXPRCURV_CONVEX) );

}

/* check for sum expression */
Test(curvature, sum)
{
   /* sum of linear expressions -> linear */
   SCIP_CALL( checkCurvature("<x>[C] + <y>[C] + <z>[C]", "sum", SCIP_EXPRCURV_LINEAR) );

   /* sum of convex expressions -> convex */
   SCIP_CALL( checkCurvature("(<x>[C])^2 + <y>[C] + <z>[C]", "sum", SCIP_EXPRCURV_CONVEX) );

   /* sum of concave expressions -> concave */
   SCIP_CALL( checkCurvature("<x>[C] + <y>[C]^(0.5) + <z>[C]^(-1)", "sum", SCIP_EXPRCURV_CONCAVE) );

   /* sum of concave and convex expressions -> unknown */
   SCIP_CALL( checkCurvature("<x>[C]^2 + <y>[C]^(0.5) + <z>[C]^(-1)", "sum", SCIP_EXPRCURV_UNKNOWN) );

   /* -1 * convex = concave */
   SCIP_CALL( checkCurvature("-1 * <x>[C]^2", "sum", SCIP_EXPRCURV_CONCAVE) );

   /* -1 * concave = convex */
   SCIP_CALL( checkCurvature("-1 * <z>[C]^(-1)", "sum", SCIP_EXPRCURV_CONVEX) );
}

/* check for value expression */
Test(curvature, value)
{
   SCIP_CALL( checkCurvature("5.2", "val", SCIP_EXPRCURV_LINEAR) );
}

/* check for variable expression */
Test(curvature, variable)
{
   SCIP_CALL( checkCurvature("<x>[C]", "var", SCIP_EXPRCURV_LINEAR) );
}

/* check curvature in the constraint data and in the nonlinear rows */
Test(curvature, cons_and_nlrows)
{
   const char* inputs[3] = {
      "[expr] <c1>: (<y>[C] + <z>[C])^2 <= 12;",
      "[expr] <c2>: -(<x>[C] + <y>[C])^2 + <x>[C] + <y>[C] + <z>[C] == -2;",
      "[expr] <c3>: <x>[C] * <y>[C] * <z>[C] - <x>[C] <= 1;"
      };
   SCIP_EXPRCURV targetcurvs[3] = {SCIP_EXPRCURV_CONVEX, SCIP_EXPRCURV_CONCAVE, SCIP_EXPRCURV_UNKNOWN};
   int ninputs = 3;
   int i;

   /* disable presolving */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* create, add, and release expression constraints */
   for( i = 0; i < ninputs; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Bool success;

      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[i],
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      cr_expect(success);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to the solving stage; this should have triggered CONSINITSOL */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );
   cr_assert(SCIPgetNConss(scip) == ninputs);
   cr_assert(SCIPgetNNLPNlRows(scip) == ninputs);

   for( i = 0; i < ninputs; ++i )
   {
      SCIP_CONS* cons;
      SCIP_NLROW* nlrow;

      cons = SCIPgetConss(scip)[i];
      assert(cons != NULL);

      /* check curvature that is stored in the constraint data */
      cr_expect(SCIPgetCurvatureConsExpr(scip, cons) == targetcurvs[i], "for cons %d (%s): expected %d got %d", i,
         SCIPconsGetName(cons), targetcurvs[i], SCIPgetCurvatureConsExpr(scip, cons));

      /* check curvature that is stored in the nonlinear row */
      nlrow = SCIPgetNLPNlRows(scip)[i];
      assert(nlrow != NULL);
      cr_expect(SCIPnlrowGetCurvature(nlrow) == targetcurvs[i]);
   }
}