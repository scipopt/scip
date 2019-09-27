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

/**@file   derivative.c
 * @brief  tests expression derivative callbacks
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_pow.h"

#include "include/scip_test.h"


static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_RANDNUMGEN* rndgen;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &rndgen, 1, TRUE) );
}

static
void teardown(void)
{
   SCIPfreeRandom(scip, &rndgen);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(derivative, .init = setup, .fini = teardown);

Test(derivative, value)
{
   SCIP_CONSEXPR_EXPR* expr;

   /* create expression */
   SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, &expr, 1.0) );

   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 0.0);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, var)
{
   SCIP_CONSEXPR_EXPR* expr;

   /* create expression */
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &expr, x) );

   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 1.0);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, sum)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "5*<x>[C] -2.2 * <y>[I]";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 5.0);
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, y), -2.2);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, product)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "<x>[C] * <y>[I] * <z>[I]";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.3) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, -3.5) );

   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 4.0 * 3.5);
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, y), -2.3 * 3.5);
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, z), -2.3 * 4.0);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, pow)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "<x>[C]^(-3)";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.3) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), -3.0 * pow(2.3, -4));

   /* try an undefined point */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), SCIP_INVALID);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, signpow)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "signpower(<x>[C],2.5)";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   /* try positive */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.3) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 2.5 * pow(2.3, 1.5));

   /* try 0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 0.0);

   /* try negative */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.42) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 2.5 * pow(1.42, 1.5));

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, exp)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "exp(<x>[C])";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), exp(-2.0));

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, log)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "log(<x>[C])";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 0.5);

   /* try an undefined point */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), SCIP_INVALID);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, abs)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "abs(<x>[C])";

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -5.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), -1.0);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 5.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 1.0);

   /* for x = 0 the gradient should can be everything in [-1,1] */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x) >= -1.0);
   cr_expect(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x) <= 1.0);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, quadratic)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "3.0*<x>[C]^2 -4.0*<y>[C] * <x>[C] + <z>[C] * <x>[C] -7.0*<z>[C]";
   int i;

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   for( i = 0; i < 100; ++i )
   {
      SCIP_Real xval, yval, zval;

      xval = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      yval = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      zval = SCIPrandomGetReal(rndgen, -10.0, 10.0);

      SCIP_CALL( SCIPsetSolVal(scip, sol, x, xval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, yval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, zval) );

      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), 6.0 * xval - 4.0 * yval + zval));
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, y), -4.0 * xval));
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, z), xval - 7.0));
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, complex1)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "<x>[C] * <y>[C] * exp(<x>[C]^2 + <z>[C]^2 - 5)";
   int i;

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   for( i = 0; i < 100; ++i )
   {
      SCIP_Real xval, yval, zval;

      xval = SCIPrandomGetReal(rndgen, -2.0, 2.0);
      yval = SCIPrandomGetReal(rndgen, -2.0, 2.0);
      zval = SCIPrandomGetReal(rndgen, -2.0, 2.0);

      SCIP_CALL( SCIPsetSolVal(scip, sol, x, xval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, yval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, zval) );

      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), (yval + 2*pow(xval,2) * yval) * exp(pow(xval,2) + pow(zval,2) -5)));
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, y), xval * exp(pow(xval,2) + pow(zval,2) - 5)));
      cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, z), 2*xval*yval*zval*exp(pow(xval,2) + pow(zval,2) - 5)));
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

Test(derivative, complex2)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "-2*<x>[C] * (<x>[C]^2 - 9)^(-2)";
   int i;

   /* create expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   for( i = 0; i < 100; ++i )
   {
      SCIP_Real xval;

      xval = SCIPrandomGetReal(rndgen, -2.0, 2.0);

      SCIP_CALL( SCIPsetSolVal(scip, sol, x, xval) );

      SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

      if( xval == 3.0 )
         cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), SCIP_INVALID));
      else
         cr_expect(SCIPisEQ(scip, SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), (6*pow(xval,2) + 18) / pow(pow(xval,2)-9, 3)));
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

static const int K = 40;
Test(performance, rosenbrock)
{
   SCIP_VAR* vars[K];
   SCIP_CONSEXPR_EXPR* exprx[K];
   SCIP_CONSEXPR_EXPR* rosenbrock;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   for( int i = 0; i < K; ++i )
   {
      char name[100];

      snprintf(name, 100, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   SCIP_CONSEXPR_EXPR* sum;
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sum, 0, NULL, NULL, 0.0) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprx[0], vars[0]) );
   for( int i = 1; i < K; ++i )
   {
      SCIP_CONSEXPR_EXPR* sqr;
      SCIP_CONSEXPR_EXPR* pow1;
      SCIP_CONSEXPR_EXPR* pow2;
      SCIP_CONSEXPR_EXPR* diff1;
      SCIP_CONSEXPR_EXPR* diff2;

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &exprx[i], vars[i]) );
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &sqr, exprx[i-1], 2.0) );

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &diff1, 0, NULL, NULL, 0.0) );
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &diff2, 0, NULL, NULL, 0.0) );

      /* diff1 = x_i - x_{i-1}^2 */
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, diff1, exprx[i], 1.0) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, diff1, sqr, -1.0) );

      /* pow1 = (x_i - x_{i-1}^2)^2 */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &pow1, diff1, 2.0) );

      /* diff2 = 1 - x_{i-1} */
      SCIPsetConsExprExprSumConstant(diff2, 1.0);
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, diff2, exprx[i-1], -1.0) );

      /* pow2 = (1 - x_{i-1})^2 */
      SCIP_CALL( SCIPcreateConsExprExprPow(scip, conshdlr, &pow2, diff2, 2.0) );

      /* sum += 100(x_i - x_{i-1}^2)^2 + (1 - x_{i-1})^2 */
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sum, pow1, 100.0) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sum, pow2, 1.0) );

      /* release */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sqr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &pow1) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &pow2) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &diff1) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &diff2) );
   }

   /* simplify is incredibly slow for large k! */
   //SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, sum, &rosenbrock) );
   //SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sum) );
   rosenbrock = sum;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( int i = 0; i < K; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], (SCIP_Real)i) );
   }

   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, rosenbrock, sol, 0) );

   for( int i = 1; i < K - 1; ++i )
   {
      SCIP_Real expected;
      SCIP_Real partial;

      expected = 200.0 * ( i - (i - 1.0) * (i - 1.0) ) - 400.0 * i * (i + 1 - i * i) - 2.0 * (1.0 - i);
      partial = SCIPgetConsExprExprPartialDiff(scip, conshdlr, rosenbrock, vars[i]);
      if( fabs(partial - expected) > 1e-5 )
         printf("Something terribly wrong: got %f, expected %f (error %f)\n", partial, expected, fabs(partial - expected));
   }

   for( int i = 0; i < K; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], 0.0) );
   }

   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, rosenbrock, sol, 0) );

   for( int i = 0; i < K - 1; ++i )
   {
      if( fabs(SCIPgetConsExprExprPartialDiff(scip, conshdlr, rosenbrock, vars[i]) - (-2.0)) > 1e-5 )
      {
         printf("Something terribly wrong: got %.10f", SCIPgetConsExprExprPartialDiff(scip, conshdlr, rosenbrock, vars[i]));
         printf(" error is %.10f\n", fabs(SCIPgetConsExprExprPartialDiff(scip, conshdlr, rosenbrock, vars[i]) - (-2.0)));
      }
   }

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &rosenbrock) );

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   for( int i = 0; i < K; ++i )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exprx[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}
