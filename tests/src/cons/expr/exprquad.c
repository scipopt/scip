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

/**@file   exprquad.c
 * @brief  tests quadratic representation of expression
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

/*
 * TEST
 */

#include "scip/cons_expr.h"
#include "scip/struct_cons_expr.h"
#include "scip/cons_expr.c"  /* so can call canonicalizeConstraints */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;

static SCIP_CONSHDLR* conshdlr;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get quadratic handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

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
   //BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* detects x^2 + x as quadratic expression */
Test(exprquad, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONSEXPR_QUADEXPR* quaddata = NULL;
   SCIP_CONSEXPR_QUADEXPRTERM quad;
   SCIP_EXPRCURV curv;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;
   SCIP_VAR* var;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>^2 + <x>", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed, &infeasible) );
   cr_expect(changed);
   cr_expect_not(infeasible);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;

   /* detect */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, conshdlr, expr, &quaddata) );
   cr_assert_not_null(quaddata);

   cr_expect_eq(quaddata->nlinexprs, 0, "Expecting 0 linear expr, got %d\n", quaddata->nlinexprs);
   cr_expect_eq(quaddata->nquadexprs, 1, "Expecting 1 quadratic terms, got %d\n", quaddata->nquadexprs);
   cr_expect_eq(quaddata->nbilinexprterms, 0, "Expecting 0 bilinear terms, got %d\n", quaddata->nbilinexprterms);

   quad = quaddata->quadexprterms[0];
   cr_assert_not_null(quad.expr);
   var = SCIPgetConsExprExprAuxVar(quad.expr);
   fprintf(stderr, "x = %s, quad.expr's auxvar %s\n", SCIPvarGetName(x), SCIPvarGetName(var));
   cr_expect_eq(var, x, "Expecting var %s in quad term, got %s\n", SCIPvarGetName(x), SCIPvarGetName(var));
   cr_expect_eq(1.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 1.0, quad.lincoef);
   cr_expect_eq(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 1.0, quad.sqrcoef);

   SCIP_CALL( SCIPgetConsExprQuadraticCurvature(scip, quaddata, &curv) );
   cr_assert_eq(curv, SCIP_EXPRCURV_CONVEX);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* detects x^2 + 2*x cos(y x^2) + cos(y x^2)^2 <= 1 as convex quadratic expression:
 * simplify yields x^2 + 2 * x cos(x^2 y) + cos(x^2 y)^2 <= 1 --> should detect x^2 + 2 x * w + w^2
 */
Test(exprquad, detectandfree2, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* cosexpr;
   SCIP_CONSEXPR_QUADEXPR* quaddata = NULL;
   SCIP_CONSEXPR_QUADEXPRTERM quad;
   SCIP_CONSEXPR_BILINEXPRTERM bilin;
   SCIP_CONS* cons;
   SCIP_EXPRCURV curv;
   SCIP_Bool infeasible;
   SCIP_Bool success;

   /* create expression, simplify it and find common subexpressions*/
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <x>^2 + 2 * <x> * cos(<y> * <x>^2) + cos(<y> * <x>^2)^2 <= 1", TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   success = FALSE;
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_assert(!infeasible);

   /* get expr and work with it */
   expr = SCIPgetExprConsExpr(scip, cons);

   /* get cosine expression */
   cr_assert_eq(SCIPgetConsExprExprNChildren(expr), 3);
   cosexpr = SCIPgetConsExprExprChildren(expr)[1]; /*  x * cos(x^2 y) */
   cosexpr = SCIPgetConsExprExprChildren(cosexpr)[1]; /* cos(x^2 y) */
   cr_assert_str_eq(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(cosexpr)), "cos", "expecting cos got %s\n",
         SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(cosexpr)));
   /* detect */
   SCIP_CALL( SCIPgetConsExprQuadratic(scip, conshdlr, expr, &quaddata) );
   cr_assert_not_null(quaddata);

   cr_expect_eq(quaddata->nlinexprs, 0, "Expecting 0 linear vars, got %d\n", quaddata->nlinexprs);
   cr_expect_eq(quaddata->nquadexprs, 2, "Expecting 2 quadratic terms, got %d\n", quaddata->nquadexprs);
   cr_expect_eq(quaddata->nbilinexprterms, 1, "Expecting 1 bilinear terms, got %d\n", quaddata->nbilinexprterms);

   /* x var */
   quad = quaddata->quadexprterms[0];
   cr_assert_not_null(quad.expr);
   cr_expect_eq(x, SCIPgetConsExprExprAuxVar(quad.expr), "Expecting var %s in quad term, got %s\n",
         SCIPvarGetName(x), SCIPvarGetName(SCIPgetConsExprExprAuxVar(quad.expr)));
   cr_expect_eq(0.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 0.0, quad.lincoef);
   cr_expect_eq(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 1.0, quad.sqrcoef);

   /* expr cos(x^2 y) is quadratic */
   quad = quaddata->quadexprterms[1];
   cr_assert_not_null(quad.expr);
   cr_expect_eq(cosexpr, quad.expr);
   cr_expect_eq(0.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 0.0, quad.lincoef);
   cr_expect_eq(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 0.0, quad.sqrcoef);

   bilin = quaddata->bilinexprterms[0];
   cr_assert_not_null(bilin.expr1);
   cr_assert_not_null(bilin.expr2);
   cr_expect_eq(SCIPgetConsExprExprAuxVar(bilin.expr1), x, "Expecting expr's auxvar %s in bilin term, got %s\n",
         SCIPvarGetName(x), SCIPvarGetName(SCIPgetConsExprExprAuxVar(bilin.expr1)));
   cr_expect_eq(bilin.expr2, cosexpr);
   cr_expect_eq(2.0, bilin.coef, "Expecting bilinear coef of %g, got %g\n", 2.0, bilin.coef);

   SCIP_CALL( SCIPgetConsExprQuadraticCurvature(scip, quaddata, &curv) );
   cr_assert_eq(curv, SCIP_EXPRCURV_CONVEX);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
