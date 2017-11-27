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

/**@file   nlhdlr_quadratic.c
 * @brief  tests quadratic nonlinear handler methods
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

/* XXX: need the consdata struct because we don't have getNlhdlrs or findNlhdlrs; I don't add those function because I'm unsure
 * we actually need them
 */
#include "scip/cons_expr.c"

#include "scip/cons_expr_nlhdlr_quadratic.c"

/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static SCIP_CONSHDLR* conshdlr;
static SCIP_CONSEXPR_NLHDLR* nlhdlr = NULL;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   int h;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get quadratic handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* get nlhdlr */
   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), "quadratic") == 0 )
      {
         nlhdlr = conshdlrdata->nlhdlrs[h];
         break;
      }
   cr_assert_not_null(nlhdlr);


   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z ", -1.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* go to TRANSFORMED stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, TRUE) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   //BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* detects x^2 + x as quadratic expression */
Test(nlhdlrquadratic, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_Bool success;


   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>^2 + <x>", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, expr, &simplified) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;

   /* detect */
   SCIP_CALL( detectHdlrQuadratic(scip, conshdlr, nlhdlr, expr, &success, &nlhdlrexprdata) );
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nlinvars, 0, "Expecting 0 linear vars, got %d\n", nlhdlrexprdata->nlinvars);
   cr_expect_eq(nlhdlrexprdata->nquadvars, 1, "Expecting 1 quadratic terms, got %d\n", nlhdlrexprdata->nquadvars);
   cr_expect_eq(nlhdlrexprdata->nbilinterms, 0, "Expecting 0 bilinear terms, got %d\n", nlhdlrexprdata->nbilinterms);

   SCIP_QUADVARTERM quad;
   quad = nlhdlrexprdata->quadvarterms[0];
   cr_assert_not_null(quad.var);
   fprintf(stderr, "x = %p, quad.var %p\n", x, quad.var);
   cr_expect_eq(quad.var, x, "Expecting var %s in quad term, got %s\n", SCIPvarGetName(x), SCIPvarGetName(quad.var));
   cr_expect_eq(1.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 1.0, quad.lincoef);
   cr_expect_eq(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 1.0, quad.sqrcoef);

   /* register nlhdlr info in expr and free */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &expr->nlhdlrs, &nlhdlr, 1) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &expr->nlhdlrsexprdata, &nlhdlrexprdata, 1) );
   expr->nnlhdlrs = 1;

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* detects x^2 + 2*x exp(y x^2) <= 1 as quadratic expression:
 * simplify yields x^2 + 2 * x exp(x^2 y) <= 1 --> should detect x^2 + 2 x * w
 */
Test(nlhdlrquadratic, detectandfree2, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* expexpr;
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeasible;


   /* create expression and simplify it; introduce auxiliary variables */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <x>^2 + 2 * <x> * exp(<y> * <x>^2) <= 1", TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   success = FALSE;
   SCIP_CALL( simplifyConstraints(scip, &cons, 1, &success) );
   cr_assert(success);

   infeasible = FALSE;
   SCIP_CALL( createAuxVars(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   /* get expr and work with it */
   expr = SCIPgetExprConsExpr(scip, cons);

   /* get exp expression */
   cr_assert_eq(SCIPgetConsExprExprNChildren(expr), 2);
   expexpr = SCIPgetConsExprExprChildren(expr)[1]; /* pow < product since x < exp(x^2 y) (var < anything) */
   expexpr = SCIPgetConsExprExprChildren(expexpr)[1]; /* var < anything */
   cr_assert_str_eq(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expexpr)), "exp", "expecting exp got %s\n",
         SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expexpr)));

   /* detect */
   SCIP_CALL( detectHdlrQuadratic(scip, conshdlr, nlhdlr, expr, &success, &nlhdlrexprdata) );
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nlinvars, 0, "Expecting 0 linear vars, got %d\n", nlhdlrexprdata->nlinvars);
   cr_expect_eq(nlhdlrexprdata->nquadvars, 2, "Expecting 2 quadratic terms, got %d\n", nlhdlrexprdata->nquadvars);
   cr_expect_eq(nlhdlrexprdata->nbilinterms, 1, "Expecting 1 bilinear terms, got %d\n", nlhdlrexprdata->nbilinterms);

   /* x var */
   SCIP_QUADVARTERM quad;
   quad = nlhdlrexprdata->quadvarterms[0];
   cr_assert_not_null(quad.var);
   cr_expect_eq(x, quad.var, "Expecting var %s in quad term, got %s\n", SCIPvarGetName(x), SCIPvarGetName(quad.var));
   cr_expect_eq(0.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 0.0, quad.lincoef);
   cr_expect_eq(1.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 1.0, quad.sqrcoef);

   /* auxiliary var for exp(x^2 y) */
   quad = nlhdlrexprdata->quadvarterms[1];
   cr_assert_not_null(quad.var);
   cr_expect_eq(SCIPgetConsExprExprLinearizationVar(expexpr), quad.var, "Expecting var %s in quad term, got %s\n",
         SCIPvarGetName(SCIPgetConsExprExprLinearizationVar(expexpr)), SCIPvarGetName(quad.var));
   cr_expect_eq(0.0, quad.lincoef, "Expecting lincoef %g in quad term, got %g\n", 0.0, quad.lincoef);
   cr_expect_eq(0.0, quad.sqrcoef, "Expecting sqrcoef %g in quad term, got %g\n", 0.0, quad.sqrcoef);

   SCIP_BILINTERM bilin;
   bilin = nlhdlrexprdata->bilinterms[0];
   cr_assert_not_null(bilin.var1);
   cr_assert_not_null(bilin.var2);
   cr_expect_eq(bilin.var1, x, "Expecting var %s in bilin term, got %s\n", SCIPvarGetName(x), SCIPvarGetName(bilin.var1));
   cr_expect_eq(bilin.var2, SCIPgetConsExprExprLinearizationVar(expexpr), "Expecting var %s in bilin term, got %s\n",
         SCIPvarGetName(SCIPgetConsExprExprLinearizationVar(expexpr)), SCIPvarGetName(bilin.var2));
   cr_expect_eq(2.0, bilin.coef, "Expecting bilinear coef of %g, got %g\n", 2.0, bilin.coef);

   /* register nlhdlr info in expr and free */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &expr->nlhdlrs, &nlhdlr, 1) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &expr->nlhdlrsexprdata, &nlhdlrexprdata, 1) );
   expr->nnlhdlrs = 1;

   /* XXX: why does release cons doesn't free the auxvars (of the transformed constraint?)? */
   SCIP_CALL( freeAuxVars(scip, conshdlr, &cons, 1) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
