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

/**@file   exprhdlr_sin.c
 * @brief  tests expression handler functions of sine expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_sin.h"
#include "scip/cons_expr_var.h"
#include <scip/cons_expr_value.h>

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONSEXPR_EXPR* sinexpr;
static SCIP_CONSEXPR_EXPR* xexpr;
static SCIP_CONSEXPR_EXPR* yexpr;
static SCIP_RANDNUMGEN* rndgen;


/* creates scip, problem, includes expression constraint handler, creates  and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* disable relaxing variable bounds in activity evaluation */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/expr/varboundrelax", 'n') );

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

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &rndgen, 1, TRUE) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free random number generator */
   SCIPfreeRandom(scip, &rndgen);

   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sinexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
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

Test(sin, parse, .description = "Tests the expression parsing.")
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "sin(<x>[C])";

   /* create sine expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );

   cr_assert(expr != NULL);
   cr_expect(SCIPgetConsExprExprNChildren(expr) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(sin, eval, .description = "Tests the expression evaluation.")
{
   int i;
   SCIP_Real randnum;
   SCIP_Real testvalues[5] = {0.0, M_PI, 2.5*M_PI, -0.5*M_PI, -2.0*M_PI};
   SCIP_Real results[5] = {0.0, 0.0, 1.0, -1.0, 0.0};

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(sinexpr), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i )
   {
      randnum = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(sinexpr), SIN(randnum)));
   }
}

Test(sin, inteval, .description = "Tests the expression interval evaluation.")
{
   int i;
   SCIP_INTERVAL interval;
   SCIP_Real rndlb[5];
   SCIP_Real rndub[5];
   SCIP_Real rndreslb[5];
   SCIP_Real rndresub[5];

   /* pick 5 special cases with well known results */
   SCIP_Real detlb[5] = {0.0, -0.5*M_PI, -M_PI, 2.0*M_PI, 0.5*M_PI};
   SCIP_Real detub[5] = {2.0*M_PI, 0.5*M_PI, -0.5*M_PI, 3.0*M_PI, 1.5*M_PI};
   SCIP_Real detreslb[5] = {-1.0, -1.0, -1.0, 0.0, -1.0};
   SCIP_Real detresub[5] = {1.0, 1.0, 0.0, 1.0, 1.0};

   /* create 5 random cases within specific bounds that have non-trivial results */
   rndlb[0]    = SCIPrandomGetReal(rndgen, 0.0, 0.5*M_PI);
   rndub[0]    = SCIPrandomGetReal(rndgen, 0.5*M_PI, M_PI);
   rndreslb[0] = MIN(SIN(rndlb[0]), SIN(rndub[0]));
   rndresub[0] = 1.0;

   rndlb[1]    = SCIPrandomGetReal(rndgen, 0.0, 0.5*M_PI);
   rndub[1]    = SCIPrandomGetReal(rndgen, M_PI, 1.5*M_PI);
   rndreslb[1] = SIN(rndub[1]);
   rndresub[1] = 1.0;

   rndlb[2]    = SCIPrandomGetReal(rndgen, 0.5*M_PI, M_PI);
   rndub[2]    = SCIPrandomGetReal(rndgen, M_PI, 1.5*M_PI);
   rndreslb[2] = SIN(rndub[2]);
   rndresub[2] = SIN(rndlb[2]);

   rndlb[3]    = SCIPrandomGetReal(rndgen, M_PI, 1.5*M_PI);
   rndub[3]    = SCIPrandomGetReal(rndgen, 1.5*M_PI, 2.0*M_PI);
   rndreslb[3] = -1.0;
   rndresub[3] = MAX(SIN(rndlb[3]), SIN(rndub[3]));

   rndlb[4]    = SCIPrandomGetReal(rndgen, 0.0, 0.5*M_PI);
   rndub[4]    = SCIPrandomGetReal(rndgen, 1.5*M_PI, 2.0*M_PI);
   rndreslb[4] = -1.0;
   rndresub[4] = 1.0;

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, detlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, detub[i]) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, sinexpr, sol, 0) );
      SCIPincrementConsExprCurBoundsTag(conshdlr, TRUE);
      SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, sinexpr, &interval, FALSE) );

      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), detreslb[i]));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), detresub[i]));
   }

   /* random part */
   for ( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, rndlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, rndub[i]) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, sinexpr, sol, 0) );
      SCIPincrementConsExprCurBoundsTag(conshdlr, TRUE);
      SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, sinexpr, &interval, FALSE) );

      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), rndreslb[i]));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), rndresub[i]));
   }
}

Test(sin, derivative, .description = "Tests the expression derivation.")
{
   int i;
   SCIP_Real randnum;
   SCIP_Real testvalues[5] = {0.0, M_PI, 2.0*M_PI, 2.5*M_PI, -0.5*M_PI};
   SCIP_Real results[5] = {1.0, -1.0, 1.0, 0.0, 0.0};

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, sinexpr, sol, 0) );

      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, sinexpr, x), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i)
   {
      randnum = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );
      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, sinexpr, x), COS(randnum)));
   }
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

   SCIP_CALL( SCIPgetConsExprExprHash(scip, expr1, &hashkey1) );
   SCIP_CALL( SCIPgetConsExprExprHash(scip, expr2, &hashkey2) );
   SCIP_CALL( SCIPgetConsExprExprHash(scip, expr3, &hashkey3) );

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
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;
   SCIP_CONSEXPR_EXPR* expr3;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* expr1 = <5.0>, expr2 = sin(<5.0>), expr3 is buffer for simplification */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr1, 5.0) );
   SCIP_CALL( SCIPcreateConsExprExprSin(scip, conshdlr, &expr2, expr1) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr2, &expr3, &changed, &infeasible) );
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr2, sol, 0) );

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPgetConsExprExprHdlr(expr3) == SCIPgetConsExprExprHdlrValue(conshdlr));
   cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr2), SIN(5.0)));

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr3) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );
}
