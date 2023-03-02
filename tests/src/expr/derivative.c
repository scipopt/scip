/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   derivative.c
 * @brief  tests expression derivative callbacks
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_value.h"
#include "scip/expr_var.h"
#include "scip/expr_sum.h"
#include "scip/expr_pow.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_RANDNUMGEN* rndgen;


/* give derivative in expr that belongs to var */
static
SCIP_Real getPartialDiff(
   SCIP_EXPR* expr,
   SCIP_VAR*           var
   )
{
   SCIP_EXPRITER* it;
   SCIP_Real deriv = 0.0;

   SCIP_CALL_ABORT( SCIPcreateExpriter(scip, &it) );

   /* we need to sum-up the derivative value from all expr that represent variable var; but only visit each expr once */
   for( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      if( SCIPisExprVar(scip, expr) && SCIPgetVarExprVar(expr) == var )
         deriv += SCIPexprGetDerivative(expr);

   SCIPfreeExpriter(&it);

   return deriv;
}


static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

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
   SCIP_EXPR* expr;

   /* create expression */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 0.0);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, var)
{
   SCIP_EXPR* expr;

   /* create expression */
   SCIP_CALL( SCIPcreateExprVar(scip, &expr, x, NULL, NULL) );

   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 1.0);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, sum)
{
   SCIP_EXPR* expr;
   const char* input = "5*<x>[C] -2.2 * <y>[I]";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 5.0);
   cr_expect_eq(getPartialDiff(expr, y), -2.2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, product)
{
   SCIP_EXPR* expr;
   const char* input = "<x>[C] * <y>[I] * <z>[I]";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.3) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, -3.5) );

   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 4.0 * 3.5);
   cr_expect_eq(getPartialDiff(expr, y), -2.3 * 3.5);
   cr_expect_eq(getPartialDiff(expr, z), -2.3 * 4.0);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, pow)
{
   SCIP_EXPR* expr;
   const char* input = "<x>[C]^(-3)";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.3) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), -3.0 * pow(2.3, -4));

   /* try an undefined point */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(SCIPexprGetDerivative(expr), SCIP_INVALID);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, signpow)
{
   SCIP_EXPR* expr;
   const char* input = "signpower(<x>[C],2.5)";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* try positive */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.3) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 2.5 * pow(2.3, 1.5));

   /* try 0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 0.0);

   /* try negative */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.42) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 2.5 * pow(1.42, 1.5));

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, exp)
{
   SCIP_EXPR* expr;
   const char* input = "exp(<x>[C])";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), exp(-2.0));

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, log)
{
   SCIP_EXPR* expr;
   const char* input = "log(<x>[C])";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 0.5);

   /* try an undefined point */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(SCIPexprGetDerivative(expr), SCIP_INVALID);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, abs)
{
   SCIP_EXPR* expr;
   const char* input = "abs(<x>[C])";

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -5.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), -1.0);

   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 5.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect_eq(getPartialDiff(expr, x), 1.0);

   /* for x = 0 the gradient should can be everything in [-1,1] */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
   cr_expect(getPartialDiff(expr, x) >= -1.0);
   cr_expect(getPartialDiff(expr, x) <= 1.0);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, quadratic)
{
   SCIP_EXPR* expr;
   const char* input = "3.0*<x>[C]^2 -4.0*<y>[C] * <x>[C] + <z>[C] * <x>[C] -7.0*<z>[C]";
   int i;

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   for( i = 0; i < 100; ++i )
   {
      SCIP_Real xval, yval, zval;

      xval = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      yval = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      zval = SCIPrandomGetReal(rndgen, -10.0, 10.0);

      SCIP_CALL( SCIPsetSolVal(scip, sol, x, xval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, yval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, zval) );

      SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
      cr_expect(SCIPisEQ(scip, getPartialDiff(expr, x), 6.0 * xval - 4.0 * yval + zval));
      cr_expect(SCIPisEQ(scip, getPartialDiff(expr, y), -4.0 * xval));
      cr_expect(SCIPisEQ(scip, getPartialDiff(expr, z), xval - 7.0));
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, complex1)
{
   SCIP_EXPR* expr;
   const char* input = "<x>[C] * <y>[C] * exp(<x>[C]^2 + <z>[C]^2 - 5)";
   int i;

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   for( i = 0; i < 100; ++i )
   {
      SCIP_Real xval, yval, zval;

      xval = SCIPrandomGetReal(rndgen, -2.0, 2.0);
      yval = SCIPrandomGetReal(rndgen, -2.0, 2.0);
      zval = SCIPrandomGetReal(rndgen, -2.0, 2.0);

      SCIP_CALL( SCIPsetSolVal(scip, sol, x, xval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, yval) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, z, zval) );

      SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );
      cr_expect(SCIPisEQ(scip, getPartialDiff(expr, x), (yval + 2*pow(xval,2) * yval) * exp(pow(xval,2) + pow(zval,2) -5)));
      cr_expect(SCIPisEQ(scip, getPartialDiff(expr, y), xval * exp(pow(xval,2) + pow(zval,2) - 5)));
      cr_expect(SCIPisEQ(scip, getPartialDiff(expr, z), 2*xval*yval*zval*exp(pow(xval,2) + pow(zval,2) - 5)));
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(derivative, complex2)
{
   SCIP_EXPR* expr;
   const char* input = "-2*<x>[C] * (<x>[C]^2 - 9)^(-2)";
   int i;

   /* create expression */
   cr_expect_eq(SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   for( i = 0; i < 100; ++i )
   {
      SCIP_Real xval;

      xval = SCIPrandomGetReal(rndgen, -2.0, 2.0);

      SCIP_CALL( SCIPsetSolVal(scip, sol, x, xval) );

      SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, 0) );

      if( xval == 3.0 )
         cr_expect(SCIPisEQ(scip, getPartialDiff(expr, x), SCIP_INVALID));
      else
         cr_expect(SCIPisEQ(scip, getPartialDiff(expr, x), (6*pow(xval,2) + 18) / pow(pow(xval,2)-9, 3)), "expected %g, got %g", (6*pow(xval,2) + 18) / pow(pow(xval,2)-9, 3), getPartialDiff(expr, x));
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

#define K 40
Test(performance, rosenbrock)
{
   SCIP_VAR* vars[K];
   SCIP_EXPR* exprx[K];
   SCIP_EXPR* rosenbrock;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   for( int i = 0; i < K; ++i )
   {
      char name[100];

      snprintf(name, 100, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   SCIP_EXPR* sum;
   SCIP_CALL( SCIPcreateExprSum(scip, &sum, 0, NULL, NULL, 0.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &exprx[0], vars[0], NULL, NULL) );
   for( int i = 1; i < K; ++i )
   {
      SCIP_EXPR* sqr;
      SCIP_EXPR* pow1;
      SCIP_EXPR* pow2;
      SCIP_EXPR* diff1;
      SCIP_EXPR* diff2;

      SCIP_CALL( SCIPcreateExprVar(scip, &exprx[i], vars[i], NULL, NULL) );
      SCIP_CALL( SCIPcreateExprPow(scip, &sqr, exprx[i-1], 2.0, NULL, NULL) );

      SCIP_CALL( SCIPcreateExprSum(scip, &diff1, 0, NULL, NULL, 0.0, NULL, NULL) );
      SCIP_CALL( SCIPcreateExprSum(scip, &diff2, 0, NULL, NULL, 0.0, NULL, NULL) );

      /* diff1 = x_i - x_{i-1}^2 */
      SCIP_CALL( SCIPappendExprSumExpr(scip, diff1, exprx[i], 1.0) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, diff1, sqr, -1.0) );

      /* pow1 = (x_i - x_{i-1}^2)^2 */
      SCIP_CALL( SCIPcreateExprPow(scip, &pow1, diff1, 2.0, NULL, NULL) );

      /* diff2 = 1 - x_{i-1} */
      SCIPsetConstantExprSum(diff2, 1.0);
      SCIP_CALL( SCIPappendExprSumExpr(scip, diff2, exprx[i-1], -1.0) );

      /* pow2 = (1 - x_{i-1})^2 */
      SCIP_CALL( SCIPcreateExprPow(scip, &pow2, diff2, 2.0, NULL, NULL) );

      /* sum += 100(x_i - x_{i-1}^2)^2 + (1 - x_{i-1})^2 */
      SCIP_CALL( SCIPappendExprSumExpr(scip, sum, pow1, 100.0) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, sum, pow2, 1.0) );

      /* release */
      SCIP_CALL( SCIPreleaseExpr(scip, &sqr) );
      SCIP_CALL( SCIPreleaseExpr(scip, &pow1) );
      SCIP_CALL( SCIPreleaseExpr(scip, &pow2) );
      SCIP_CALL( SCIPreleaseExpr(scip, &diff1) );
      SCIP_CALL( SCIPreleaseExpr(scip, &diff2) );
   }

   /* simplify is incredibly slow for large k! */
   //SCIP_CALL( SCIPsimplifyExpr(scip, conshdlr, sum, &rosenbrock) );
   //SCIP_CALL( SCIPreleaseExpr(scip, &sum) );
   rosenbrock = sum;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( int i = 0; i < K; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], (SCIP_Real)i) );
   }

   SCIP_CALL( SCIPevalExprGradient(scip, rosenbrock, sol, 0) );

   for( int i = 1; i < K - 1; ++i )
   {
      SCIP_Real expected;
      SCIP_Real partial;

      expected = 200.0 * ( i - (i - 1.0) * (i - 1.0) ) - 400.0 * i * (i + 1 - i * i) - 2.0 * (1.0 - i);
      partial = getPartialDiff(rosenbrock, vars[i]);
      if( fabs(partial - expected) > 1e-5 )
         printf("Something terribly wrong: got %f, expected %f (error %f)\n", partial, expected, fabs(partial - expected));
   }

   for( int i = 0; i < K; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], 0.0) );
   }

   SCIP_CALL( SCIPevalExprGradient(scip, rosenbrock, sol, 0) );

   for( int i = 0; i < K - 1; ++i )
   {
      if( fabs(getPartialDiff(rosenbrock, vars[i]) - (-2.0)) > 1e-5 )
      {
         printf("Something terribly wrong: got %.10f", getPartialDiff(rosenbrock, vars[i]));
         printf(" error is %.10f\n", fabs(getPartialDiff(rosenbrock, vars[i]) - (-2.0)));
      }
   }

   SCIP_CALL( SCIPreleaseExpr(scip, &rosenbrock) );

   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   for( int i = 0; i < K; ++i )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &exprx[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}
