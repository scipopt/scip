/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   curvature.c
 * @brief  tests curvature expression handler callbacks
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_VAR* w;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

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

   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory leak!!");
}


/** auxiliary function for creating an expression and checking its curvature */
static
SCIP_RETCODE checkCurvature(
   const char*           input,              /**< input creating an expression */
   const char*           exprhdlrname,       /**< target expression handler name */
   SCIP_EXPRCURV         expectedcur         /**< expected curvature */
   )
{
   SCIP_EXPR* origexpr;
   SCIP_EXPR* expr;
   SCIP_EXPRHDLR* exprhdlr;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* create and print expression */
   SOFT_ASSERT_EQUAL(SCIPparseExpr(scip, &origexpr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* simplify expression */
   SCIP_CALL( SCIPsimplifyExpr(scip, origexpr, &expr, &changed, &infeasible, NULL, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &origexpr) );

   /* print simplified expression */
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* check name of the corresponding expression handler */
   exprhdlr = SCIPexprGetHdlr(expr);
   TEST_ASSERT(exprhdlr != NULL);
   SOFT_ASSERT(strcmp(SCIPexprhdlrGetName(exprhdlr), exprhdlrname) == 0, "expect expression handler %s, got %s\n",
      exprhdlrname, SCIPexprhdlrGetName(exprhdlr));

   /* compute curvature */
   SCIP_CALL( SCIPcomputeExprCurvature(scip, expr) );

   /* check curvature */
   SOFT_ASSERT(SCIPexprGetCurvature(expr) == expectedcur, "expect %s, got %s", SCIPexprcurvGetName(expectedcur), SCIPexprcurvGetName(SCIPexprGetCurvature(expr)));

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   return SCIP_OKAY;
}

/* check for absolute expression */
void test_curvature_absolute(void)
{
   SCIP_CALL( checkCurvature("abs(<x>[C])", "abs", SCIP_EXPRCURV_CONVEX) );
}

/* check for cosine expression */
void test_curvature_cosine(void)
{
   SCIP_CALL(checkCurvature("cos(<x>[C])", "cos", SCIP_EXPRCURV_UNKNOWN));
   SCIP_CALL(checkCurvature("cos(<y>[C])", "cos", SCIP_EXPRCURV_UNKNOWN));
   SCIP_CALL(checkCurvature("cos(<z>[C])", "cos", SCIP_EXPRCURV_CONVEX));
   SCIP_CALL(checkCurvature("cos(<w>[C])", "cos", SCIP_EXPRCURV_CONCAVE));
}

/* check for exponential expression */
void test_curvature_exponential(void)
{
   SCIP_CALL( checkCurvature("exp(<x>[C])", "exp", SCIP_EXPRCURV_CONVEX) );
}

/* check for logarithm expression */
void test_curvature_logarithm(void)
{
   SCIP_CALL( checkCurvature("log(<x>[C])", "log", SCIP_EXPRCURV_CONCAVE) );
}

/* check for power expression */
void test_curvature_power(void)
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
void test_curvature_signpower(void)
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
void test_curvature_product(void)
{
   SCIP_CALL( checkCurvature("(<x>[C] * <y>[C])", "prod", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("(<x>[C] * <y>[C] * <z>[C])", "prod", SCIP_EXPRCURV_UNKNOWN) );
}

/* check for sine expression */
void test_curvature_sine(void)
{
   SCIP_CALL( checkCurvature("sin(<x>[C])", "sin", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("sin(<y>[C])", "sin", SCIP_EXPRCURV_CONCAVE) );
   SCIP_CALL( checkCurvature("sin(<z>[C])", "sin", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("sin(<w>[C])", "sin", SCIP_EXPRCURV_CONVEX) );
}

/* check for sum expression */
void test_curvature_sum(void)
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
void test_curvature_value(void)
{
   SCIP_CALL( checkCurvature("5.2", "val", SCIP_EXPRCURV_LINEAR) );
}

/* check for variable expression */
void test_curvature_variable(void)
{
   SCIP_CALL( checkCurvature("<x>[C]", "var", SCIP_EXPRCURV_LINEAR) );
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_curvature_absolute
);
   RUN_TEST(test_curvature_cosine);
   RUN_TEST(test_curvature_exponential);
   RUN_TEST(test_curvature_logarithm);
   RUN_TEST(test_curvature_power);
   RUN_TEST(test_curvature_signpower);
   RUN_TEST(test_curvature_product);
   RUN_TEST(test_curvature_sine);
   RUN_TEST(test_curvature_sum);
   RUN_TEST(test_curvature_value);
   RUN_TEST(test_curvature_variable);
   return UNITY_END();
}
