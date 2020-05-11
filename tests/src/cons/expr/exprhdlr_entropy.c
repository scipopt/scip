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

/**@file   exprhdlr_entropy.c
 * @brief  tests expression handler functions of entropy an expression
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_entropy.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_CONSEXPR_EXPR* entropyexpr;
static SCIP_CONSEXPR_EXPR* prodexpr;          /* xlogx as product expression */
static SCIP_CONSEXPR_EXPR* negprodexpr;       /* -xlogx as product expression */
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

   /* create variable and entropy expressions */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &entropyexpr, xexpr) );
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x>[C]*log(<x>[C])", NULL, &prodexpr) );
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "-<x>[C]*log(<x>[C])", NULL, &negprodexpr) );

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
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &negprodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &entropyexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/* test suite */
TestSuite(entropy, .init = setup, .fini = teardown);

/*
 * TESTS
 */

Test(entropy, creation, .description = "Tests the expression creation.")
{
   SCIP_CONSEXPR_EXPR* expr;

   /* create entropy expression */
   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr, xexpr) );

   cr_assert(expr != NULL);
   cr_expect(SCIPgetConsExprExprNChildren(expr) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(entropy, parse, .description = "Tests the expression parsing.")
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "entropy(<x>[C])";

   /* create entropy expression */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr) );

   cr_assert(expr != NULL);
   cr_expect(SCIPgetConsExprExprNChildren(expr) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(entropy, eval, .description = "Tests the expression evaluation.")
{
   SCIP_Real testvalues[4] = {-1.0, 0.0, exp(-1.0), 1.0};
   SCIP_Real results[4] = {SCIP_INVALID, 0.0, exp(-1.0), 0.0};
   SCIP_Real randnum;
   int i;

   /* deterministic part */
   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, entropyexpr, sol, 0) );

      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprValue(entropyexpr), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i )
   {
      randnum = SCIPrandomGetReal(rndgen, 0.0, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );

      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, entropyexpr, sol, 0) );
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprValue(entropyexpr), -randnum * log(randnum)));

      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, negprodexpr, sol, 0) );
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprValue(negprodexpr), -randnum * log(randnum)));
   }
}

Test(entropy, inteval, .description = "Tests the expression interval evaluation.")
{
   SCIP_INTERVAL intervalEntropy;
   /*FIXME SCIP_INTERVAL intervalProd;*/
   SCIP_Real rndlb[4];
   SCIP_Real rndub[4];
   SCIP_Real rndreslb[4];
   SCIP_Real rndresub[4];
   int i;

   /* pick 5 special cases with well known results */
   SCIP_Real detlb[4] = {0.0, 0.0, exp(-1.0), 1.0};
   SCIP_Real detub[4] = {exp(-1.0), 1.0, 1.0, 2.0};
   SCIP_Real detreslb[4] = {0.0, 0.0, 0.0, -2.0 * log(2.0)};
   SCIP_Real detresub[4] = {exp(-1.0), exp(-1.0), exp(-1.0), 0.0};

   /* create 5 random cases within specific bounds that have non-trivial resulsts */
   rndlb[0]    = SCIPrandomGetReal(rndgen, 0.0, exp(-1.0));
   rndub[0]    = SCIPrandomGetReal(rndgen, exp(-1.0), 1.0);
   rndreslb[0] = MIN(-rndlb[0] * log(rndlb[0]), -rndub[0] * log(rndub[0]));
   rndresub[0] = exp(-1.0);

   rndlb[1]    = SCIPrandomGetReal(rndgen, 0.0, exp(-1.0));
   rndub[1]    = SCIPrandomGetReal(rndgen, 1.0, 10.0);
   rndreslb[1] = -rndub[1] * log(rndub[1]);
   rndresub[1] = exp(-1.0);

   rndlb[2]    = SCIPrandomGetReal(rndgen, exp(-1.0), 1.0);
   rndub[2]    = SCIPrandomGetReal(rndgen, 1.0, 10.0);
   rndreslb[2] = -rndub[2] * log(rndub[2]);
   rndresub[2] = -rndlb[2] * log(rndlb[2]);

   rndlb[3]    = SCIPrandomGetReal(rndgen, 1.0, 10.0);
   rndub[3]    = SCIPrandomGetReal(rndgen, rndlb[3], 10.0);
   rndreslb[3] = -rndub[3] * log(rndub[3]);
   rndresub[3] = -rndlb[3] * log(rndlb[3]);

   /* deterministic part */
   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, detlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, detub[i]) );
      SCIPincrementConsExprCurBoundsTag(conshdlr, TRUE);
      SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, entropyexpr, &intervalEntropy, FALSE, FALSE) );

      cr_expect(SCIPisEQ(scip, intervalEntropy.inf, detreslb[i]));
      cr_expect(SCIPisEQ(scip, intervalEntropy.sup, detresub[i]));

      /* FIXME what was this code meant to do?
      intervalProd = SCIPgetConsExprExprInterval(entropyexpr);
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, negprodexpr, 0, NULL, NULL) );
      cr_expect(SCIPisLE(scip, intervalEntropy.inf, intervalProd.inf));
      cr_expect(SCIPisGE(scip, intervalEntropy.sup, intervalProd.sup));
      */
   }

   /* random part */
   for ( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, rndlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, rndub[i]) );
      SCIPincrementConsExprCurBoundsTag(conshdlr, TRUE);
      SCIP_CALL( SCIPevalConsExprExprActivity(scip, conshdlr, entropyexpr, &intervalEntropy, FALSE, FALSE) );

      cr_expect(SCIPisEQ(scip, intervalEntropy.inf, rndreslb[i]));
      cr_expect(SCIPisEQ(scip, intervalEntropy.sup, rndresub[i]));

      /* FIXME what was this code meant to do?
      intervalProd = SCIPgetConsExprExprInterval(entropyexpr);
      SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, negprodexpr, 0, NULL, NULL) );
      cr_expect(SCIPisLE(scip, intervalEntropy.inf, intervalProd.inf));
      cr_expect(SCIPisGE(scip, intervalEntropy.sup, intervalProd.sup));
      */
   }
}

Test(entropy, derivative, .description = "Tests the expression derivation.")
{
   SCIP_Real testvalues[5] = {-1.0, 0.0, exp(-1.0), 1.0, 2.0};
   SCIP_Real results[5] = {SCIP_INVALID, SCIP_INVALID, 0.0, -1.0, -1.0 - log(2.0)};
   SCIP_Real randnum;
   int i;

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, entropyexpr, sol, 0) );

      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, entropyexpr, x), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i)
   {
      randnum = SCIPrandomGetReal(rndgen, 1e-12, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );

      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, entropyexpr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, entropyexpr, x), -1.0 - log(randnum)));

      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, negprodexpr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, negprodexpr, x), -1.0 - log(randnum)));
   }
}

Test(entropy, hash, .description = "Tests the expression hash.")
{
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;
   SCIP_CONSEXPR_EXPR* expr3;
   unsigned int hashkey1;
   unsigned int hashkey2;
   unsigned int hashkey3;

   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr1, xexpr) );
   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr2, xexpr) );
   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr3, yexpr) );

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

Test(entropy, simplify, .description = "Tests the expression simplification.")
{
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;
   SCIP_CONSEXPR_EXPR* expr3;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* expr1 = <5.0>, expr2 = entropy(<5.0>), expr3 is buffer for simplification */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr1, 5.0));
   SCIP_CALL( SCIPcreateConsExprExprEntropy(scip, conshdlr, &expr2, expr1));
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr2, &expr3, &changed, &infeasible));
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr2, sol, 0) );

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPgetConsExprExprHdlr(expr3) == SCIPgetConsExprExprHdlrValue(conshdlr));
   cr_expect(SCIPisFeasEQ(scip, SCIPgetConsExprExprValue(expr2), -5.0 * log(5.0)));

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr3) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );

   /* test if product of x and log(x) is transformed to sum of entropy expression
    * expr1 is buffer for simplification and expr2 is used to store children
    */
   changed = FALSE;
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, prodexpr, &expr1, &changed, &infeasible));

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPgetConsExprExprHdlr(expr1) == SCIPgetConsExprExprHdlrSum(conshdlr));
   cr_expect(SCIPgetConsExprExprNChildren(expr1) == 1);
   cr_expect(SCIPgetConsExprExprSumCoefs(expr1)[0] == -1.0);

   expr2 = SCIPgetConsExprExprChildren(expr1)[0];
   cr_expect(SCIPgetConsExprExprHdlr(expr2) == SCIPfindConsExprExprHdlr(conshdlr, "entropy"));
   cr_expect(SCIPgetConsExprExprNChildren(expr2) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr2)[0] == xexpr);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );


   /* test if product of -x and log(x) is transformed to entropy expression
    * expr1 is buffer for simplification
    */
   changed = FALSE;
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, negprodexpr, &expr1, &changed, &infeasible));

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPgetConsExprExprHdlr(expr1) == SCIPfindConsExprExprHdlr(conshdlr, "entropy"));
   cr_expect(SCIPgetConsExprExprNChildren(expr1) == 1);
   cr_expect(SCIPgetConsExprExprChildren(expr1)[0] == xexpr);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );
}
