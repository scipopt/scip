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

/**@file   exprhdlr_sin.c
 * @brief  tests expression handler functions of xzy an expression
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_var.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONSEXPR_EXPR* sinexpr;
static SCIP_CONSEXPR_EXPR* xexpr;
static SCIP_CONSEXPR_EXPR* yexpr;

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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create variable and sine expressions */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &sinexpr, xexpr) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sinexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/* test suite */
TestSuite(sin, .init = setup, .fini = teardown);

/*
 * TESTS
 */

Test(sin, creation, .description = "Tests the expression creation.")
{
   SCIP_CONSEXPR_EXPR* expr;

   /* create sine expression */
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr, xexpr) );

   cr_assert(expr != NULL);
   cr_expect(SCIPgetConsExprExprNChildren(expr) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(sin, print, .description = "Tests the expression printing function.")
{
   /* TODO */
}

Test(sin, parse, .description = "Tests the expression parsing.")
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "sin(<x>[C])";

   /* create expression for product of -5, sin(x), and sin(y) */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );

   cr_expect(expr != NULL);
   cr_expect(SCIPgetConsExprExprNChildren(expr) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(sin, eval, .description = "Tests the expression evaluation.")
{
   /* get random real number between -10 and 10 */
   SCIP_RANDNUMGEN* rndgen;
   SCIP_CALL( SCIPcreateRandom(scip, &rndgen, 1) );
   SCIP_Real randnum = SCIPrandomGetReal(rndgen, -10.0, 10.0);

   /* create 5 different testvalues and corresponding expected results */
   SCIP_Real testvalues[5] = {0.0, M_PI, 2.5*M_PI, -0.5*M_PI, randnum};
   SCIP_Real results[5] = {0.0, 0.0, 1.0, -1.0, SIN(randnum)};

   for(int i = 0; i < 5; i++)
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(sinexpr), results[i]));
   }

   SCIPfreeRandom(scip, &rndgen);
}

Test(sin, inteval, .description = "Tests the expression interval evaluation.")
{
   /* create 5 sets of bounds to cover different cases */
   SCIP_Real testlowerbounds[5] = {0.0, -2.7, M_PI, 0.5*M_PI, 1.5};
   SCIP_Real testupperbounds[5] = {2*M_PI, 1.9, 6.0, 1.5*M_PI, 4.0};
   SCIP_Real reslowerbounds[5] = {-1.0, -1.0, -1.0, -1.0, SIN(4.0)};
   SCIP_Real resupperbounds[5] = {1.0, 1.0, 0.0, 1.0, 1.0};

   for(int i = 0; i < 5; i++)
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, testlowerbounds[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, testupperbounds[i]) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, sinexpr, sol, 0) );
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, sinexpr, FALSE, 0, SCIPepsilon(scip)) );

      SCIP_INTERVAL interval = SCIPgetConsExprExprInterval(sinexpr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), reslowerbounds[i]));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), resupperbounds[i]));
   }

}

Test(sin, derivative, .description = "Tests the expression derivation.")
{
   /* TODO */
}

Test(sin, hash, .description = "Tests the expression hash.")
{
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;
   SCIP_CONSEXPR_EXPR* expr3;
   unsigned int hashkey1;
   unsigned int hashkey2;
   unsigned int hashkey3;

   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr1, xexpr) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr2, xexpr) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr3, yexpr) );

   SCIP_CALL( SCIPgetConsExprExprHashkey(scip, expr1, &hashkey1) );
   SCIP_CALL( SCIPgetConsExprExprHashkey(scip, expr2, &hashkey2) );
   SCIP_CALL( SCIPgetConsExprExprHashkey(scip, expr3, &hashkey3) );

   cr_expect(hashkey1 != 0);
   cr_expect(hashkey2 != 0);
   cr_expect(hashkey3 != 0);
   cr_expect(hashkey1 == hashkey2);
   cr_expect(hashkey1 != hashkey3);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr3) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );
}

Test(sin, simplify, .description = "Tests the expression simplification.")
{
   /* TODO */
}
