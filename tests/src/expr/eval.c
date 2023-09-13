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

/**@file   eval.c
 * @brief  tests expression's pointwise and interval evaluation
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/expr_sum.h"
#include "scip/expr_product.h"
#include "scip/expr_pow.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* turn off log() and pow() assuming arguments to a away from zero
    * good for solving, but complicates unittesting
    */
   SCIP_CALL( SCIPsetRealParam(scip, "expr/log/minzerodistance", 0.0) );
   SCIP_CALL( SCIPsetRealParam(scip, "expr/pow/minzerodistance", 0.0) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1.0e-7) );
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

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

TestSuite(evalexpr, .init = setup, .fini = teardown);

/** TESTS **/
Test(evalexpr, absolute, .description = "Tests expression evaluation for absolute expressions.")
{
   int i;
   SCIP_EXPR* expr;
   SCIP_INTERVAL interval;
   const char* inputs[2] = {"abs(<x>[C]) + abs(<x>[C])",
      "abs(abs(abs(<x>[C]) + abs(<y>[C])) * abs(<x>[C])^3 * abs(<y>[C]))"};

   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)inputs[0], NULL, NULL, NULL) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = -10; i <= 10; ++i )
   {
      /* evaluate expression at i */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr),  2.0 * ABS(i)));

      /* evaluate expression at interval [-|i|, i^2]*/
      SCIP_CALL( SCIPchgVarLb(scip, x, -ABS(i)) );
      SCIP_CALL( SCIPchgVarUb(scip, x, i*i) );
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
      interval = SCIPexprGetActivity(expr);

      cr_expect_float_eq(SCIPintervalGetInf(interval), 0.0, SCIPepsilon(scip), "interval.inf = %g", interval.inf);
      cr_expect_float_eq(SCIPintervalGetSup(interval), 2.0 * i*i, SCIPepsilon(scip), "interval.sup = %g, expected %g", interval.sup, 2.0*i*i);
   }
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* test a more complicated expression |x + y| |x|^3 |y|*/
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)inputs[1], NULL, NULL, NULL) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = -10; i <= 10; ++i )
   {
      /* evaluate expression at i */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr), 2.0 * pow(ABS(i),5)));

      /* evaluate expression at interval [-|i|, |i|]^2 */
      SCIP_CALL( SCIPchgVarLb(scip, x, -ABS(i)) );
      SCIP_CALL( SCIPchgVarUb(scip, x, ABS(i)) );
      SCIP_CALL( SCIPchgVarLb(scip, y, -ABS(i)) );
      SCIP_CALL( SCIPchgVarUb(scip, y, ABS(i)) );
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
      interval = SCIPexprGetActivity(expr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 0.0));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2.0 * pow(ABS(i),5)));
   }
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(evalexpr, exponential, .description = "Tests expression evaluation for exponential expressions.")
{
   SCIP_EXPR* expr;
   SCIP_INTERVAL interval;
   const char* inputs[2] = {"exp(<x>[C]) + exp(<x>[C])", "exp(exp(<x>[C])) * exp(<y>[C])^2"};
   int i;

   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)inputs[0], NULL, NULL, NULL) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = -10; i <= 10; ++i )
   {
      /* evaluate expression at i */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr), exp(i) + exp(i)));

      /* evaluate expression at interval [i, i + 1/(|i| + 1)] */
      SCIP_CALL( SCIPchgVarLb(scip, x, i) );
      SCIP_CALL( SCIPchgVarUb(scip, x, i + 1.0 / (ABS(i) + 1)) );
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
      interval = SCIPexprGetActivity(expr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 2*exp(i)));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2*exp(i + 1.0 / (ABS(i) + 1))));
   }
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* complicated exponential expression e^(2y + e^x) */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)inputs[1], NULL, NULL, NULL) );
   SCIPinfoMessage(scip, NULL, "testing expression: ");
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* evaluate expression for different points */
   for( i = 1; i <= 10; ++i )
   {
      /* evaluate expression at (1/i, i) */
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) 1.0 / i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i) );
      SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr), exp(exp(1.0 / i)) * exp(2*i)));

      /* evaluate expression at interval [-1/i, 1/i] x [i, i + 1/i] */
      SCIP_CALL( SCIPchgVarLb(scip, x, -1.0 / i) );
      SCIP_CALL( SCIPchgVarUb(scip, x,  1.0 / i) );
      SCIP_CALL( SCIPchgVarLb(scip, y, i) );
      SCIP_CALL( SCIPchgVarUb(scip, y, i + 1.0 / i) );
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
      interval = SCIPexprGetActivity(expr);
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), exp(exp(-1.0 / i)) * exp(2*i)));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), exp(exp(1.0 / i)) * exp(2*i + 2.0 / i)));
   }
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(evalexpr, logarithm, .description = "Tests expression evaluation for logarithmic expressions.")
{
      SCIP_EXPR* expr;
      SCIP_INTERVAL interval;
      const char* inputs[2] = {"log(<x>[C]) + log(<x>[C])", "log(log(exp(<x>[C]) * exp(<y>[C])))"};
      int i;

      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)inputs[0], NULL, NULL, NULL) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = -10; i <= 10; ++i )
      {
         SCIP_Real xlb, xub;

         /* evaluate expression at i */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );

         if( i <= 0 )
            cr_expect_eq(SCIPexprGetEvalValue(expr), SCIP_INVALID);
         else
            cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr), log(i) + log(i)));

         /* evaluate expression at interval [i, i + 1/(|i|+1)] */
         xlb = i;
         xub = i + 1.0 / (ABS(i) + 1);
         SCIP_CALL( SCIPchgVarLb(scip, x, xlb) );
         SCIP_CALL( SCIPchgVarUb(scip, x, xub) );
         SCIP_CALL( SCIPevalExprActivity(scip, expr) );
         interval = SCIPexprGetActivity(expr);

         /* interval is empty if both bounds are non-positive */
         if( xub <= 0 )
            cr_expect(SCIPintervalIsEmpty(SCIPinfinity(scip), interval));
         else
         {
            cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), 2*log(xub)));

            if( xlb <= 0 )
               cr_expect(SCIPisInfinity(scip, -SCIPintervalGetInf(interval)));
            else
               cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), 2*log(xlb)));
         }
      }
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

      /* complicated logarithmic expression: log( log( exp(x) * exp(y) ) ) = log( x + y ) */
      SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)inputs[1], NULL, NULL, NULL) );
      SCIPinfoMessage(scip, NULL, "testing expression: ");
      SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      /* evaluate expression for different points */
      for( i = 1; i <= 10; ++i )
      {
         /* evaluate expression at (i, i + 1) */
         SCIP_CALL( SCIPsetSolVal(scip, sol, x, (SCIP_Real) i) );
         SCIP_CALL( SCIPsetSolVal(scip, sol, y, (SCIP_Real) i + 1) );
         SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
         cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr), log(2*i + 1) ));

         /* evaluate expression at intervar [1/i, 2/i] x [3/i, 4/i] */
         SCIP_CALL( SCIPchgVarLb(scip, x,  1.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, x,  2.0 / i) );
         SCIP_CALL( SCIPchgVarLb(scip, y,  3.0 / i) );
         SCIP_CALL( SCIPchgVarUb(scip, y,  4.0 / i) );
         SCIP_CALL( SCIPevalExprActivity(scip, expr) );
         interval = SCIPexprGetActivity(expr);
         cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), log(1.0 / i + 3.0 / i)));
         cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), log(2.0 / i + 4.0 / i)));
      }

      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(evalexpr, power, .description = "Tests expression evaluation for power expressions.")
{
   SCIP_EXPR* expr;
   SCIP_EXPR* xexpr;
   SCIP_INTERVAL interval;

   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );

   /* create expressions x^2.5 */
   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &expr, xexpr, 2.5, NULL, NULL) );

   /* evaluate expression at 2.0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
   cr_assert(SCIPisEQ(scip, SCIPexprGetEvalValue(expr), pow(2.0, 2.5)));

   /* evaluate expression for an undefined point: -1.0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
   cr_assert(SCIPexprGetEvalValue(expr) == SCIP_INVALID);

   /* evaluate expression at interval [1.0, 3.0] */
   SCIP_CALL( SCIPchgVarLb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   interval = SCIPexprGetActivity(expr);
   cr_assert(SCIPisEQ(scip, interval.inf, 1.0));
   cr_assert(SCIPisEQ(scip, interval.sup, pow(3.0, 2.5)));

   /* evaluate expression at an undefined interval [-2.0, -1.0] => resulting interval should be empty */
   SCIP_CALL( SCIPchgVarLb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   interval = SCIPexprGetActivity(expr);
   cr_assert(SCIPintervalIsEmpty(SCIPinfinity(scip), interval));

   /* free expressions */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
}

Test(evalexpr, signpower, .description = "Tests expression evaluation for signpower expressions.")
{
   SCIP_EXPR* expr;
   SCIP_EXPR* xexpr;
   SCIP_INTERVAL interval;

   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );

   /* create expressions signpower(x,2.5) */
   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSignpower(scip, &expr, xexpr, 2.5, NULL, NULL) );

   /* evaluate expression at 2.0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
   cr_assert(SCIPisEQ(scip, SCIPexprGetEvalValue(expr), pow(2.0, 2.5)));

   /* evaluate expression at -1.0 */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPevalExpr(scip, expr, sol, 0) );
   cr_assert(SCIPisEQ(scip, SCIPexprGetEvalValue(expr), -1.0));

   /* evaluate expression at interval [-1.0, 3.0] */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 3.0) );
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   interval = SCIPexprGetActivity(expr);
   cr_assert(SCIPisEQ(scip, interval.inf, -1.0));
   cr_assert(SCIPisEQ(scip, interval.sup, pow(3.0, 2.5)));

   /* evaluate expression at an interval [-2.0, -1.0] */
   SCIP_CALL( SCIPchgVarLb(scip, x, -2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   SCIP_CALL( SCIPevalExprActivity(scip, expr) );
   interval = SCIPexprGetActivity(expr);
   cr_assert(SCIPisEQ(scip, interval.inf, -pow(2.0, 2.5)));
   cr_assert(SCIPisEQ(scip, interval.sup, -1.0));

   /* free expressions */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
}

/* creates expression for f(x,y) = 0.5 * ( (x^2*y^(-1)*5^(-4))^2 * (2*x + 1)^(-1) ) */
static
SCIP_RETCODE createExpr(
   SCIP_EXPR**  xexpr,              /**< pointer to store variable expression */
   SCIP_EXPR**  yexpr,              /**< pointer to store variable expression */
   SCIP_EXPR**  const_expr,         /**< pointer to store constant expression */
   SCIP_EXPR**  prodexpr,           /**< pointer to store product expression */
   SCIP_EXPR**  sumexpr,            /**< pointer to store sum expression */
   SCIP_EXPR**  mainexpr            /**< pointer to store full expression */
   )
{
   SCIP_EXPR* exprs[] = {NULL, NULL, NULL};
   SCIP_Real coef = 2.0;

   /* create variable and constant expressions */
   SCIP_CALL( SCIPcreateExprVar(scip, xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprValue(scip, const_expr, 5.0, NULL, NULL) );

   /* create sum expression */
   SCIP_CALL( SCIPcreateExprSum(scip, sumexpr, 1, xexpr, &coef, 1.0, NULL, NULL) );  /* 2*x+1 */

   /* create product expression */
   SCIP_CALL( SCIPcreateExprPow(scip, &exprs[0], *xexpr,  2.0, NULL, NULL) );  /* x^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &exprs[1], *yexpr, -1.0, NULL, NULL) );  /* y^(-1) */
   SCIP_CALL( SCIPcreateExprPow(scip, &exprs[2], *const_expr, -4.0, NULL, NULL) );  /* 5^(-4) */
   SCIP_CALL( SCIPcreateExprProduct(scip, prodexpr, 3, exprs, 1.0, NULL, NULL) );  /* x^2*y^(-1)*5^(-4) */

   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[2]) );

   /* create main expression */
   SCIP_CALL( SCIPcreateExprPow(scip, &exprs[0], *prodexpr, 2.0, NULL, NULL) );  /* (x^2*y^(-1)*5^(-4))^2 */
   SCIP_CALL( SCIPcreateExprPow(scip, &exprs[1], *sumexpr, -1.0, NULL, NULL) );  /* (2*x + 1)^(-1) */
   SCIP_CALL( SCIPcreateExprProduct(scip, mainexpr, 2, exprs, 0.5, NULL, NULL) );

   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[0]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[1]) );

   return SCIP_OKAY;
}

/* helper function to check evaluation of expression created with createExpr() */
static
void checkExprEval(
   SCIP_EXPR*  xexpr,               /**< variable expression */
   SCIP_EXPR*  yexpr,               /**< variable expression */
   SCIP_EXPR*  const_expr,          /**< constant expression */
   SCIP_EXPR*  prodexpr,            /**< product expression */
   SCIP_EXPR*  sumexpr,             /**< sum expression */
   SCIP_EXPR*  mainexpr,            /**< full expression */
   SCIP_Real            xval,                /**< x value used for evaluation */
   SCIP_Real            yval,                /**< y value used for evaluation */
   unsigned int         tag                  /**< tag used for evaluation */
   )
{
   SCIP_Real prodval;
   SCIP_Real sumval;

   prodval = pow(xval,2)*pow(yval,-1)*pow(5,-4);
   sumval = 2*xval + 1;

   /* check values */
   cr_expect_float_eq(SCIPexprGetEvalValue(mainexpr), 0.5 * pow(prodval,2) * pow(sumval,-1), SCIPepsilon(scip));
   cr_expect_float_eq(SCIPexprGetEvalValue(sumexpr), sumval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPexprGetEvalValue(prodexpr), prodval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPexprGetEvalValue(xexpr), xval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPexprGetEvalValue(yexpr), yval, SCIPepsilon(scip));
   cr_expect_float_eq(SCIPexprGetEvalValue(const_expr), 5.0, SCIPepsilon(scip));

   /* check tags */
   cr_expect_eq(SCIPexprGetEvalTag(mainexpr), tag);
   cr_expect_eq(SCIPexprGetEvalTag(sumexpr), tag);
   cr_expect_eq(SCIPexprGetEvalTag(prodexpr), tag);
   cr_expect_eq(SCIPexprGetEvalTag(xexpr), tag);
   cr_expect_eq(SCIPexprGetEvalTag(yexpr), tag);
   cr_expect_eq(SCIPexprGetEvalTag(const_expr), tag);
}

Test(evalexpr, complicated, .description = "Tests expression evaluation for a large complicated expression.")
{
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* const_expr;
   SCIP_EXPR* prodexpr;
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* mainexpr;
   SCIP_EXPRPRINTDATA* dotdata;
   int i;

   /* create all expressions */
   SCIP_CALL( createExpr(&xexpr, &yexpr, &const_expr, &prodexpr, &sumexpr, &mainexpr) );

   /* initialize solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 4.0) );

   /* evaluate main expression, print it, and check values for sub-expressions */
   printf("evaluate and check expressions\n");
   SCIP_CALL( SCIPevalExpr(scip, mainexpr, sol, 1) );
   SCIP_CALL( SCIPprintExprDotInit(scip, &dotdata, NULL, SCIP_EXPRPRINT_EXPRSTRING | SCIP_EXPRPRINT_EVALTAG) );
   SCIP_CALL( SCIPprintExprDot(scip, dotdata, mainexpr) );
   SCIP_CALL( SCIPprintExprDotFinal(scip, &dotdata) );
   checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, 2.0, 4.0, 1);

   /* modify solution and evaluate expression with the same tag again; values should not change */
   printf("evaluate and check expressions with a modified solution but the same tag\n");
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -5.0) );
   SCIP_CALL( SCIPevalExpr(scip, mainexpr, sol, 1) );
   checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, 2.0, 4.0, 1);

   /* evaluate expression with a different tag; values should have changed */
   printf("evaluate expression with a new tag\n");
   SCIP_CALL( SCIPevalExpr(scip, mainexpr, sol, 2) );
   checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, -2.0, -5.0, 2);

   /* evaluate solution with zero tag */
   printf("evaluate expression with a zero tag\n");
   for( i = 1; i < 50; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, i*i) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, y, -5.0/i) );
      SCIP_CALL( SCIPevalExpr(scip, mainexpr, sol, 0) );
      checkExprEval(xexpr, yexpr, const_expr, prodexpr, sumexpr, mainexpr, i*i, -5.0 / i, 0);
   }

   /* mainexpr is not defined for x = -1 or y = 0; the result should be SCIP_INVALID */
   printf("evaluate expression for an undefined point\n");
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.0) );
   SCIP_CALL( SCIPevalExpr(scip, mainexpr, sol, 0) );
   cr_expect_eq(SCIPexprGetEvalValue(mainexpr), SCIP_INVALID);

   /* release all expressions */
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &const_expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &mainexpr) );
}


/* helper function to check interval evaluation of an expression */
static
void checkExprIntEval(
   SCIP_EXPR*            expr,               /**< expression to evaluate */
   SCIP_Real             targetinf,          /**< target infimum */
   SCIP_Real             targetsup,          /**< target supremum */
   SCIP_Bool             empty               /**< should the result be empty? */
   )
{
   SCIP_INTERVAL interval;

   assert(expr != NULL);

   SCIP_CALL_ABORT( SCIPevalExprActivity(scip, expr) );
   interval = SCIPexprGetActivity(expr);

   /* check if interval is and should be empty */
   cr_expect_eq(empty, SCIPintervalIsEmpty(SCIPinfinity(scip), interval));

   /* check interval */
   if( !empty )
   {
      /* check infimum */
      if( SCIPisInfinity(scip, -targetinf) )
         cr_expect(SCIPisInfinity(scip, -targetinf));
      else
         cr_expect_float_eq(SCIPintervalGetInf(interval), targetinf, SCIPepsilon(scip), "expected inf = %.15g, got %.15g", targetinf, interval.inf);

      /* check supremum */
      if( SCIPisInfinity(scip, targetsup) )
         cr_expect(SCIPisInfinity(scip, targetsup));
      else
         cr_expect_float_eq(SCIPintervalGetSup(interval), targetsup, SCIPepsilon(scip), "expected sup = %.15g, got %.15g", targetsup, interval.sup);
   }
}

/* Test expression interval evaluation method */
Test(evalexprInterval, complicated_interval, .description = "Tests expression interval evaluation for a large complicated expression.")
{
   SCIP_EXPR* xexpr;
   SCIP_EXPR* yexpr;
   SCIP_EXPR* const_expr;
   SCIP_EXPR* prodexpr;
   SCIP_EXPR* sumexpr;
   SCIP_EXPR* mainexpr;
   SCIP_EXPR* sqrtexpr;
   SCIP_Real xlb, xub, ylb, yub;
   SCIP_Real exponent;
#if 0
   int i;
   SCIP_INTERVAL interval;
#endif

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetRealParam(scip, "expr/log/minzerodistance", 0.0) );
   SCIP_CALL( SCIPsetRealParam(scip, "expr/pow/minzerodistance", 0.0) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -5.0, 5.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create different kind of expressions */
   exponent = 0.5;
   SCIP_CALL( createExpr(&xexpr, &yexpr, &const_expr, &prodexpr, &sumexpr, &mainexpr) );
   SCIP_CALL( SCIPcreateExprPow(scip, &sqrtexpr, xexpr, exponent, NULL, NULL) );

   /*
    * check interval evaluation method for constant expressions
    */
   printf("check interval-evaluation of constant expressions\n");
   checkExprIntEval(const_expr, 5.0, 5.0, FALSE);

   /*
    * check interval evaluation method for variable expressions
    */
   printf("check interval-evaluation of variable expressions\n");
   checkExprIntEval(xexpr, 0.0, 10.0, FALSE);
   checkExprIntEval(yexpr, -5.0, 5.0, FALSE);

   /*
    * check interval evaluation method for sum expression
    */
   printf("check interval-evaluation of sum expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 2.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 7.5) );
   checkExprIntEval(sumexpr, 5.0, 16.0, FALSE);
   checkExprIntEval(xexpr, 2.0, 7.5, FALSE);

   /*
    * check interval evaluation method for product expression: (x^2 / (y*5^(-4)))
    */
   printf("check interval-evaluation of product expression\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 2.0) );
   checkExprIntEval(prodexpr, pow(5,-4.0) / 8.0 , pow(5,-4.0), FALSE);
   checkExprIntEval(xexpr, 0.5, 1.0, FALSE);
   checkExprIntEval(yexpr, 1.0, 2.0, FALSE);

   /*
    * check interval-evaluation for a complicated expression: 0.5 * ( (x^2*y^(-1)*5^(-4))^2 * (2x + 1)^(-1) )
    */
   printf("check interval evaluation of a complicated expression\n");
   for( xub = 0.0; xub <= 10.0; xub += 1.0 )
      for( xlb = 0.0; xlb <= xub; xlb += 1.0 )
         for( yub = 1.0; yub <= 10.0; yub += 1.0 )
            for( ylb = 1.0; ylb <= yub; ylb += 1.0 )
            {
               SCIP_Real inf, sup, a, b;

               SCIP_CALL( SCIPchgVarLb(scip, x, xlb) );
               SCIP_CALL( SCIPchgVarUb(scip, x, xub) );
               SCIP_CALL( SCIPchgVarLb(scip, y, ylb) );
               SCIP_CALL( SCIPchgVarUb(scip, y, yub) );

               /* compute range of mainexpr */
               inf = (pow(5.0,-8) / 2.0) * MIN(pow(xlb,4), pow(xub,4)) * pow(yub, -2);
               sup = (pow(5.0,-8) / 2.0) * MAX(pow(xlb,4), pow(xub,4)) * pow(ylb, -2);
               a = MIN(1.0 / (1.0 + 2 * xlb), 1.0 / (1.0 + 2 * xub));
               b = MAX(1.0 / (1.0 + 2 * xlb), 1.0 / (1.0 + 2 * xub));
               inf *= b < 0 ? b : a;
               sup *= b < 0 ? a : b;

               /* check all expressions */
               checkExprIntEval(mainexpr, MIN(inf,sup), MAX(inf,sup), FALSE);
               checkExprIntEval(sumexpr, 2*xlb + 1.0, 2*xub + 1.0, FALSE);
               inf = MIN(xlb*xlb, xub*xub) * (1.0/yub) * pow(5,-4);
               sup = MAX(xlb*xlb, xub*xub) * (1.0/ylb) * pow(5,-4);
               checkExprIntEval(prodexpr, inf, sup, FALSE);
               checkExprIntEval(xexpr, xlb, xub, FALSE);
               checkExprIntEval(yexpr, ylb, yub, FALSE);
               checkExprIntEval(const_expr, 5.0, 5.0, FALSE);
            }

   /*
    * check if interval evaluation for 1/(1+2*x)^3 or 1/y leads to infinite bounds
    */
   printf("check interval-evaluation for expressions containing functions like 1/f(x)\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -0.5) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), TRUE);  /* (2x + 1)^(-1) for x=0 is now empty */

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE);

   SCIP_CALL( SCIPchgVarLb(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 0.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), TRUE);  /* 1/y for y=0 is now empty */

   /* (1/y)^2 should lead to [0,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   checkExprIntEval(mainexpr, 0.0, SCIPinfinity(scip), FALSE);

   /* (1/y)^2 should lead to [0,inf] but because of 1/(1+2*x)^3 we should get [-inf,inf] */
   SCIP_CALL( SCIPchgVarLb(scip, y, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarLb(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 10.0) );
   checkExprIntEval(mainexpr, -SCIPinfinity(scip), SCIPinfinity(scip), FALSE);

   /*
    * check if interval evaluation aborts for some cases like (-1)^2
    */
   printf("check interval-evaluation for undefined expressions like (-1)^2\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -1.0) );
   checkExprIntEval(sqrtexpr, 0, 0, TRUE);

   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, -0.5) );
   checkExprIntEval(sqrtexpr, 0, 0, TRUE);

   /* the result of sqrt([-1,2]) should be [0,sqrt(2)] and not an empty interval */
   SCIP_CALL( SCIPchgVarLb(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 2.0) );
   checkExprIntEval(sqrtexpr, 0, sqrt(2), FALSE);

   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 1.0) );
   checkExprIntEval(sqrtexpr, 0.0, 1.0, FALSE);

#if 0  /* TODO: what should be kept? */
   /*
    * check interval evaluation tags
    */
   printf("check interval tag behavior\n");
   SCIP_CALL( SCIPchgVarLb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUb(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarLb(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarUb(scip, y, 1.0) );

   /* do the expression store the box tag correctly? */
   for( i = 0; i < 10; ++i )
   {
      int tag = i % 3;
      checkExprIntEval(mainexpr, 0.0, 0.0, FALSE, tag);
      checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, tag);
      checkExprIntEval(sumexpr, 1.0, 1.0, FALSE, tag);
      checkExprIntEval(xexpr, 0.0, 0.0, FALSE, tag);
      checkExprIntEval(yexpr, 1.0, 1.0, FALSE, tag);
      checkExprIntEval(const_expr, 5.0, 5.0, FALSE, tag);
   }

   /* set another tag for some subexpression; result should not change */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, mainexpr, 1, NULL, NULL) );
   checkExprIntEval(mainexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(sumexpr, 1.0, 1.0, FALSE, 1);

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, prodexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, sumexpr, 2, NULL, NULL) );
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, mainexpr, 1, NULL, NULL) ); /* this should not trigger a reevaluation */
   checkExprIntEval( mainexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval( prodexpr, 0.0, 0.0, FALSE, 2);
   checkExprIntEval( sumexpr, 1.0, 1.0, FALSE, 2);

   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, mainexpr, 3, NULL, NULL) ); /* this should trigger a reevaluation */
   checkExprIntEval(mainexpr, 0.0, 0.0, FALSE, 3);
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 3);
   checkExprIntEval(sumexpr, 1.0, 1.0, FALSE, 3);

   /* manipulate evaluation interval */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, prodexpr, 1, NULL, NULL) );
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(xexpr, 0.0, 0.0, FALSE, 1);
   checkExprIntEval(yexpr, 1.0, 1.0, FALSE, 1);

   /* set bounds of x to [1,1] in xexpr */
   SCIPintervalSetBounds(&interval, 1.0, 1.0);
   SCIPsetConsExprExprEvalInterval(xexpr, &interval, 2);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, prodexpr, 2, NULL, NULL) ); /* should use [1,1] for xexpr */
   checkExprIntEval(prodexpr, pow(5.0,-4), pow(5.0,-4), FALSE, 2);
   checkExprIntEval(xexpr, 1.0, 1.0, FALSE, 2);
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, conshdlr, prodexpr, 3, NULL, NULL) ); /* should use [0,0] for xexpr */
   checkExprIntEval(prodexpr, 0.0, 0.0, FALSE, 3);
   checkExprIntEval(xexpr, 0.0, 0.0, FALSE, 3);
#endif
   /* release all expressions */
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &const_expr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &mainexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &sqrtexpr) );
}
