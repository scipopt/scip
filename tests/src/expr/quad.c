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

/**@file   quad.c
 * @brief  tests quadratic representation of expression
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/struct_expr.h"
#include "scip/nlpi_ipopt.h" /* to check whether LAPACK is around */
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;

/* creates scip, problem, includes expression handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 1.49, 1.51, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -0.9, 0.7, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   BMScheckEmptyMemory();
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* detects x^2 + x as quadratic expression */
void test_quad_detectandfree1(void)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_QUADEXPR_QUADTERM quad;
   SCIP_EXPRCURV curv;
   SCIP_Bool isquadratic;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;
   SCIP_VAR* var;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 + <x>", NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   SOFT_ASSERT(changed);
   SOFT_ASSERT_NOT(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &expr, 1, &changed) );

   /* detect */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );
   SOFT_ASSERT(isquadratic);

   SOFT_ASSERT_NOT_NULL(expr->quaddata);
   SOFT_ASSERT_EQUAL(expr->quaddata->nlinexprs, 0, "Expecting 0 linear expr, got %d\n", expr->quaddata->nlinexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nquadexprs, 1, "Expecting 1 quadratic terms, got %d\n", expr->quaddata->nquadexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nbilinexprterms, 0, "Expecting 0 bilinear terms, got %d\n", expr->quaddata->nbilinexprterms);
   SOFT_ASSERT(expr->quaddata->allexprsarevars);

   quad = expr->quaddata->quadexprterms[0];
   TEST_ASSERT_NOT_NULL(quad.expr);
   var = SCIPgetVarExprVar(quad.expr);
   SOFT_ASSERT_EQUAL(var, x, "Expecting var %s in quad term, got %s\n", SCIPvarGetName(x), SCIPvarGetName(var));
   SOFT_ASSERT_EQUAL(1.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 1.0, quad.lincoef);
   SOFT_ASSERT_EQUAL(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 1.0, quad.sqrcoef);

   SCIP_CALL( SCIPcomputeExprQuadraticCurvature(scip, expr, &curv, NULL, FALSE) );
   if( expr->quaddata->curvaturechecked )  /* currently check only with Ipopt */
      TEST_ASSERT_EQUAL(curv, SCIP_EXPRCURV_CONVEX);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* detect x^2 + 2*x cos(y x^2) + cos(y x^2)^2 as convex quadratic expression
 * simplify yields x^2 + 2 * x cos(x^2 y) + cos(x^2 y)^2 <= 1 --> should detect x^2 + 2 x * w + w^2
 */
void test_quad_detectandfree2(void)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_EXPR* cosexpr;
   SCIP_QUADEXPR_QUADTERM quad;
   SCIP_QUADEXPR_BILINTERM bilin;
   SCIP_EXPRCURV curv;
   SCIP_Bool infeasible;
   SCIP_Bool changed;
   SCIP_Bool isquadratic;

   /* create expression, simplify it and find common subexpressions*/
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 + 2 * <x> * cos(<y> * <x>^2) + cos(<y> * <x>^2)^2", NULL, NULL, NULL) );

   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   TEST_ASSERT(!infeasible);
   SCIPreleaseExpr(scip, &expr);
   expr = simplified;
   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &expr, 1, &changed) );

   /* get cosine expression */
   TEST_ASSERT_EQUAL(SCIPexprGetNChildren(expr), 3);
   cosexpr = SCIPexprGetChildren(expr)[1]; /*  x * cos(x^2 y) */
   cosexpr = SCIPexprGetChildren(cosexpr)[1]; /* cos(x^2 y) */
   TEST_ASSERT_EQUAL_STRING(SCIPexprhdlrGetName(SCIPexprGetHdlr(cosexpr)), "cos", "expecting cos got %s\n",
         SCIPexprhdlrGetName(SCIPexprGetHdlr(cosexpr)));
   /* detect */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );
   SOFT_ASSERT(isquadratic);

   SOFT_ASSERT_EQUAL(expr->quaddata->nlinexprs, 0, "Expecting 0 linear vars, got %d\n", expr->quaddata->nlinexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nquadexprs, 2, "Expecting 2 quadratic terms, got %d\n", expr->quaddata->nquadexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nbilinexprterms, 1, "Expecting 1 bilinear terms, got %d\n", expr->quaddata->nbilinexprterms);
   SOFT_ASSERT(!expr->quaddata->allexprsarevars);

   /* x var */
   quad = expr->quaddata->quadexprterms[0];
   TEST_ASSERT_NOT_NULL(quad.expr);
   SOFT_ASSERT_EQUAL(x, SCIPgetVarExprVar(quad.expr), "Expecting var %s in quad term, got %s\n",
         SCIPvarGetName(x), SCIPvarGetName(SCIPgetVarExprVar(quad.expr)));
   SOFT_ASSERT_EQUAL(0.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 0.0, quad.lincoef);
   SOFT_ASSERT_EQUAL(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 1.0, quad.sqrcoef);

   /* expr cos(x^2 y) is quadratic */
   quad = expr->quaddata->quadexprterms[1];
   TEST_ASSERT_NOT_NULL(quad.expr);
   SOFT_ASSERT_EQUAL(cosexpr, quad.expr);
   SOFT_ASSERT_EQUAL(0.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 0.0, quad.lincoef);
   SOFT_ASSERT_EQUAL(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 0.0, quad.sqrcoef);

   bilin = expr->quaddata->bilinexprterms[0];
   TEST_ASSERT_NOT_NULL(bilin.expr1);
   TEST_ASSERT_NOT_NULL(bilin.expr2);
   SOFT_ASSERT_EQUAL(SCIPgetVarExprVar(bilin.expr1), x, "Expecting expr's auxvar %s in bilin term, got %s\n",
         SCIPvarGetName(x), SCIPvarGetName(SCIPgetVarExprVar(bilin.expr1)));
   SOFT_ASSERT_EQUAL(bilin.expr2, cosexpr);
   SOFT_ASSERT_EQUAL(2.0, bilin.coef, "Expecting bilinear coef of %g, got %g\n", 2.0, bilin.coef);
   SOFT_ASSERT_EQUAL(bilin.pos2, 1);  /* because quaddata->quadexprterms[1].expr == cosexpr */

   SCIP_CALL( SCIPcomputeExprQuadraticCurvature(scip, expr, &curv, NULL, FALSE) );
   if( expr->quaddata->curvaturechecked )  /* currently check only with Ipopt */
      TEST_ASSERT_EQUAL(curv, SCIP_EXPRCURV_CONVEX);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* detects x^2 + 2*x*y + y^2 + y*z - z^2 as quadratic expression and convexity in (x,y) */
void test_quad_detectandfree3(void)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_EXPRCURV curv;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;
   SCIP_Bool isquadratic;
   SCIP_HASHMAP* assumevarfixed;

   /* create expression and simplify it */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 + 2*<x>*<y> + <y>^2 + <y>*<z> - <z>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   SOFT_ASSERT(changed);
   SOFT_ASSERT_NOT(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &expr, 1, &changed) );

   /* detect */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );
   SOFT_ASSERT(isquadratic);

   SOFT_ASSERT_EQUAL(expr->quaddata->nlinexprs, 0, "Expecting 0 linear expr, got %d\n", expr->quaddata->nlinexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nquadexprs, 3, "Expecting 3 quadratic terms, got %d\n", expr->quaddata->nquadexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nbilinexprterms, 2, "Expecting 2 bilinear terms, got %d\n", expr->quaddata->nbilinexprterms);
   SOFT_ASSERT(expr->quaddata->allexprsarevars);

   /* we are indefinite (aka "unknown") */
   SCIP_CALL( SCIPcomputeExprQuadraticCurvature(scip, expr, &curv, NULL, FALSE) );
   if( expr->quaddata->curvaturechecked )  /* currently check only with Ipopt */
      TEST_ASSERT_EQUAL(curv, SCIP_EXPRCURV_UNKNOWN);

   /* check again, but now assume that z will be fixed; then we should be convex */
   SCIP_CALL( SCIPhashmapCreate(&assumevarfixed, SCIPblkmem(scip), 1) );
   SCIP_CALL( SCIPhashmapInsert(assumevarfixed, (void*)z, NULL) );

   SCIP_CALL( SCIPcomputeExprQuadraticCurvature(scip, expr, &curv, assumevarfixed, FALSE) );
   if( SCIPisIpoptAvailableIpopt() )  /* currently check only with Ipopt; curvaturechecked isn't set if assumevarfixed is set */
      TEST_ASSERT_EQUAL(curv, SCIP_EXPRCURV_CONVEX);

   SCIPhashmapFree(&assumevarfixed);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* detects x^2 + 0*x*y + y^2 as quadratic expression */
void test_quad_detectandfree4(void)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* simplified;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;
   SCIP_Bool isquadratic;

   /* create expression and simplify it */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x>^2 + <x>*<y> + <y>^2 + <z>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr, &simplified, &changed, &infeasible, NULL, NULL) );
   SOFT_ASSERT(!changed);
   SOFT_ASSERT_NOT(infeasible);
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = simplified;
   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &expr, 1, &changed) );

   /* set coef of x*y and of z^2 to 0 */
   TEST_ASSERT_EQUAL(SCIPexprGetNChildren(expr), 4);
   SCIPgetCoefsExprSum(expr)[1] = 0.0;
   SCIPgetCoefsExprSum(expr)[3] = 0.0;

   /* detect */
   SCIP_CALL( SCIPcheckExprQuadratic(scip, expr, &isquadratic) );
   SOFT_ASSERT(isquadratic);

   SOFT_ASSERT_EQUAL(expr->quaddata->nlinexprs, 0, "Expecting 0 linear expr, got %d\n", expr->quaddata->nlinexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nquadexprs, 2, "Expecting 2 quadratic terms, got %d\n", expr->quaddata->nquadexprs);
   SOFT_ASSERT_EQUAL(expr->quaddata->nbilinexprterms, 0, "Expecting 0 bilinear terms, got %d\n", expr->quaddata->nbilinexprterms);
   SOFT_ASSERT(expr->quaddata->allexprsarevars);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_quad_detectandfree1
);
   RUN_TEST(test_quad_detectandfree2);
   RUN_TEST(test_quad_detectandfree3);
   RUN_TEST(test_quad_detectandfree4);
   return UNITY_END();
}
