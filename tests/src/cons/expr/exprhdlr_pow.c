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

/**@file   exprhdlr_pow.c
 * @brief  tests expression handler functions of xzy an expression
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;

/* creates scip, problem, includes expression constraint handler, creates  and adds variables */
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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/* test suite */
TestSuite(pow, .init = setup, .fini = teardown);

/*
 * TESTS
 */

Test(pow, creation, .description = "Tests the expression creation.")
{
   /* TODO */
}

Test(pow, print, .description = "Tests the expression printing function.")
{
   /* TODO */
}

Test(pow, parse, .description = "Tests the expression parsing.")
{
   /* TODO */
}

Test(pow, eval, .description = "Tests the expression evaluation.")
{
   /* TODO */
}

Test(pow, inteval, .description = "Tests the expression interval evaluation.")
{
   /* TODO */
}

Test(pow, derivative, .description = "Tests the expression derivation.")
{
   /* TODO */
}

Test(pow, hash, .description = "Tests the expression hash.")
{
   /* TODO */
}

Test(pow, simplify, .description = "Tests the expression simplification.")
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* binary variable to positive exponent */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x>^17.43", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed, &infeasible) );

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_assert_eq(SCIPgetConsExprExprHdlr(expr), SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_assert_eq(SCIPgetConsExprExprHdlr(simplified), SCIPgetConsExprExprHdlrVar(conshdlr));
   cr_assert_eq(x, SCIPgetConsExprExprVarVar(simplified));
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplified) );

   /* binary variable to negative exponent -> nothing happens */
   changed = FALSE;
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x>^(-1.43)", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed, &infeasible) );

   cr_expect_not(changed);
   cr_assert_not(infeasible);
   cr_assert_eq(SCIPgetConsExprExprHdlr(expr), SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   cr_assert_eq(SCIPgetConsExprExprHdlr(simplified), SCIPfindConsExprExprHdlr(conshdlr, "pow"));
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplified) );
}
