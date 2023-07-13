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

/**@file   exprinterpret.c
 * @brief  unit test for expression interpreter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "scip/exprinterpret.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_EXPRINT* exprint;
static SCIP_SOL* sol;

#define nvars 10
static SCIP_VAR* vars[nvars];
static SCIP_EXPR* varexprs[nvars];
static SCIP_EXPR* varidxexprs[nvars];
static SCIP_Real varvals[10][nvars];
static unsigned int soltag = 0L;

#define TOL 1e-12

static
void checkad_setup(void)
{
   char name[SCIP_MAXSTRLEN];
   int i;

   SCIP_CALL_ABORT( SCIPcreate(&scip) );
   SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );

   /* need a problem to have stat created, which is used by expr iterators */
   SCIP_CALL( SCIPcreateProbBasic(scip, "dummy") );

   for( i = 0; i < nvars; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL_ABORT( SCIPcreateVarBasic(scip, &vars[i], name, -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL_ABORT( SCIPaddVar(scip, vars[i]) );

      SCIP_CALL_ABORT( SCIPcreateExprVar(scip, &varexprs[i], vars[i], NULL, NULL) );
      SCIP_CALL_ABORT( SCIPcreateExprVaridx(scip, &varidxexprs[i], i, NULL, NULL) );
   }
   for( i = 0; i < 10; ++i )
      BMSclearMemoryArray(varvals[i], nvars);

   SCIP_CALL_ABORT( SCIPexprintCreate(scip, &exprint) );
   SCIP_CALL_ABORT( SCIPcreateSol(scip, &sol, NULL) );
}

static
void checkad_teardown(void)
{
   int i;

   SCIP_CALL_ABORT( SCIPfreeSol(scip, &sol) );
   SCIPexprintFree(scip, &exprint);

   for( i = nvars-1; i >= 0; --i )
   {
      SCIP_CALL_ABORT( SCIPreleaseExpr(scip, &varidxexprs[i]) );
      SCIP_CALL_ABORT( SCIPreleaseExpr(scip, &varexprs[i]) );
      SCIP_CALL_ABORT( SCIPreleaseVar(scip, &vars[i]) );
   }

   SCIP_CALL_ABORT( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

static
SCIP_DECL_EXPR_MAPEXPR(mapvar)
{
   int varidx;

   if( !SCIPisExprVar(scip, sourceexpr) )
      return SCIP_OKAY;

   varidx = SCIPvarGetProbindex(SCIPgetVarExprVar(sourceexpr));
   assert(varidx >= 0);
   assert(varidx < nvars);

   *targetexpr = varidxexprs[varidx];
   SCIPcaptureExpr(*targetexpr);

   return SCIP_OKAY;
}

/* checks the zeroth, first, and second derivative for an expression
 *
 * Evaluates expression, derivative, and Hessian-vector products using expression
 * callbacks and expression interpreter and check that values are the same.
 * Uses varvals[i], i<=npoints, as points where to evaluate.
 */
static
void checkAD(
   SCIP_EXPR*            expr,               /**< expression (in SCIP vars) to be checked */
   int                   dim,                /**< number of variables used in expression */
   int                   npoints             /**< number of points to check */
   )
{
   SCIP_EXPRINTDATA* exprintdata = NULL;
   SCIP_EXPR* expr2;
   SCIP_Real val;
   SCIP_Real gradient[nvars];
   SCIP_SOL* hessiandir;
   SCIP_Real hessian[nvars*nvars];
   int* hescolidx;
   int* hesrowidx;
   SCIP_Real* hesvalues;
   SCIP_Bool havehess;
   int hesnnz;
   int i, j;
   int p;

   /* if exprint cannot eval, then nothing we can check */
   if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_FUNCVALUE) )
      return;

   /* get copy of expr using varidx for var exprs */
   SCIP_CALL( SCIPduplicateExpr(scip, expr, &expr2, mapvar, NULL, NULL, NULL) );
   /* let exprint get familiar with expr */
   SCIP_CALL( SCIPexprintCompile(scip, exprint, expr2, &exprintdata) );
   /* if exprint itself has eval capabilities, then should also have this capabilities for a specific expr since every exprhdlr has to implement eval */
   cr_expect_eq(SCIPexprintGetExprCapability(scip, exprint, expr2, exprintdata) & SCIP_EXPRINTCAPABILITY_FUNCVALUE, SCIP_EXPRINTCAPABILITY_FUNCVALUE);

   for( p = 0; p < npoints; ++p )
   {
      SCIPinfoMessage(scip, NULL, "checking expression ");
      SCIPprintExpr(scip, expr, NULL);
      SCIPinfoMessage(scip, NULL, " at x =");
      for( i = 0; i < dim; ++i )
         SCIPinfoMessage(scip, NULL, " %g", varvals[p][i]);
      SCIPinfoMessage(scip, NULL, "\n");

      SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, varvals[p]) );

      SCIP_CALL( SCIPevalExpr(scip, expr, sol, ++soltag) );
      SCIP_CALL( SCIPexprintEval(scip, exprint, expr2, exprintdata, varvals[p], &val) );

      SCIPinfoMessage(scip, NULL, "  expr value: %g  exprint value: %g\n", SCIPexprGetEvalValue(expr), val);
      cr_expect(SCIPisFinite(val) == (SCIPexprGetEvalValue(expr) != SCIP_INVALID));
      if( !SCIPisFinite(val) )
         continue;
      cr_expect_float_eq(val, SCIPexprGetEvalValue(expr), TOL);


      if( !(SCIPexprintGetExprCapability(scip, exprint, expr2, exprintdata) & SCIP_EXPRINTCAPABILITY_GRADIENT) )
         continue;

      /* evaluate gradient with expr and expr2 and compare */
      SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, soltag) );
      SCIP_CALL( SCIPexprintGrad(scip, exprint, expr2, exprintdata, varvals[p], FALSE, &val, gradient) );

      cr_expect_float_eq(val, SCIPexprGetEvalValue(expr), TOL);
      for( i = 0; i < dim; ++i )
      {
         SCIPinfoMessage(scip, NULL, "  gradient[%d]: expr: %g  exprint: %g\n", i, SCIPexprGetDerivative(varexprs[i]), gradient[i]);
         /* this assumes that exprint will set all gradient values to some non-finite number if not differentiable
          * for this test, this is sufficient at the moment
          */
         cr_expect(SCIPisFinite(gradient[i]) == (SCIPexprGetDerivative(expr) != SCIP_INVALID));
         if( SCIPisFinite(gradient[i]) )
            cr_expect_float_eq(gradient[i], SCIPexprGetDerivative(varexprs[i]), TOL);
      }


      if( !(SCIPexprintGetExprCapability(scip, exprint, expr2, exprintdata) & SCIP_EXPRINTCAPABILITY_HESSIAN) )
         continue;

      /* evaluate Hessian with expr by dim times Hessian*vector-evaluation */
      SCIP_CALL( SCIPcreateSol(scip, &hessiandir, NULL) );
      for( i = 0; i < dim; ++i )
      {
         SCIP_CALL( SCIPsetSolVal(scip, hessiandir, vars[i], 1.0) );

         SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, soltag, hessiandir) );
         for( j = 0; j < dim; ++j )
            hessian[i*dim+j] = SCIPexprGetBardot(expr);

         SCIP_CALL( SCIPsetSolVal(scip, hessiandir, vars[i], 0.0) );
      }
      SCIP_CALL( SCIPfreeSol(scip, &hessiandir) );

      havehess = FALSE;
      SCIPinfoMessage(scip, NULL, "  hessian expr:\n");
      for( i = 0; i < dim; ++i )
      {
         SCIPinfoMessage(scip, NULL, "    ");
         for( j = 0; j < dim; ++j )
         {
            SCIPinfoMessage(scip, NULL, " %10g", hessian[i*dim+j]);
            if( hessian[i*dim+j] == SCIP_INVALID )
               havehess = FALSE;
         }
         SCIPinfoMessage(scip, NULL, "\n");
      }
      /* skip if no Hessian available via expr, could be an unimplemented callback
       * TODO should have this callback implemented everywhere
       */
      if( !havehess )
         continue;

      /* evaluate Hessian (sparse) with expr2 */
      SCIP_CALL( SCIPexprintHessian(scip, exprint, expr2, exprintdata, varvals[p], FALSE, &val, &hesrowidx, &hescolidx, &hesvalues, &hesnnz) );
      SCIPinfoMessage(scip, NULL, "  hessian exprint (%d nz):\n", hesnnz);
      for( i = 0; i < hesnnz; ++i )
      {
         SCIPinfoMessage(scip, NULL, "    (%d,%d) = %g", hesrowidx[i], hescolidx[i], hesvalues[i]);
         cr_assert(hesrowidx[i] >= 0);
         cr_assert(hescolidx[i] >= 0);
         cr_assert(hesrowidx[i] < dim);
         cr_assert(hescolidx[i] < dim);
         cr_assert(hescolidx[i] <= hesrowidx[i]);

         /* subtract exprint hessian from expr hessian, so we can check remaining hessian to be 0 */
         hessian[hesrowidx[i]*dim + hescolidx[i]] -= hesvalues[i];
         if( hesrowidx[i] != hescolidx[i] )
            hessian[hescolidx[i]*dim + hesrowidx[i]] -= hesvalues[i];
      }
      for( i = 0; i < dim*dim; ++i )
         cr_expect_float_eq(hessian[i], 0.0, TOL);
   }

   SCIP_CALL( SCIPexprintFreeData(scip, exprint, expr2, &exprintdata) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
}

TestSuite(checkad, .init = checkad_setup, .fini = checkad_teardown);

Test(checkad, abs)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprAbs(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0][0] = 2.0;
   varvals[1][0] = -3.0;
   checkAD(expr, 1, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, cos)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprCos(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0][0] = 0.0;
   varvals[1][0] = 1.0;
   checkAD(expr, 1, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, entropy)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0][0] = 1.0;
   varvals[1][0] = 0.0;
   varvals[2][0] = 2.0;
   checkAD(expr, 1, 3);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, exp)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprExp(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0][0] = -1.0;
   varvals[1][0] = 1.0;
   checkAD(expr, 1, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, log)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprLog(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0][0] = 0.0;
   varvals[1][0] = 0.1;
   varvals[2][0] = 1.0;
   checkAD(expr, 1, 3);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, pow)
{
   SCIP_EXPR* expr;
   SCIP_Real exponents[] = { 2.0, 3.0, 1.0, -1.0, 0.5, 1.875 };
   size_t i;

   for( i = 0; i < sizeof(exponents) / sizeof(SCIP_Real); ++i )
   {
      SCIP_CALL( SCIPcreateExprPow(scip, &expr, varexprs[0], exponents[i], NULL, NULL) );

      varvals[0][0] = 1.0;
      varvals[1][0] = 2.0;
      varvals[2][0] = -3.0;
      varvals[3][0] = 0.0;
      checkAD(expr, 1, 4);

      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }
}

Test(checkad, signpow)
{
   SCIP_EXPR* expr;
   SCIP_Real exponents[] = { 1.0, 1.875, 2.0 };
   size_t i;

   for( i = 0; i < sizeof(exponents) / sizeof(SCIP_Real); ++i )
   {
      SCIP_CALL( SCIPcreateExprSignpower(scip, &expr, varexprs[0], exponents[i], NULL, NULL) );

      varvals[0][0] = -3.0;
      varvals[1][0] = 0.0;
      varvals[2][0] = 5.0;
      checkAD(expr, 1, 3);

      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }
}

Test(checkad, product)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 2, varexprs, 5.0, NULL, NULL) );

   varvals[0][0] = 0.0;
   varvals[0][1] = 0.0;

   varvals[1][0] = 0.0;
   varvals[1][1] = 1.0;

   varvals[2][0] = 3.0;
   varvals[2][1] = 0.0;

   varvals[3][0] = 2.0;
   varvals[3][1] = -4.0;

   checkAD(expr, 2, 4);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, sin)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprSin(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0][0] = 0.0;
   varvals[1][0] = -1.0;
   checkAD(expr, 1, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, sum)
{
   SCIP_EXPR* expr;
   SCIP_Real coefs[] = { 0.0, 1.0, 2.0, 4.0 };

   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 5, varexprs, coefs, 42.0, NULL, NULL) );

   varvals[0][0] = 1.0;
   varvals[0][1] = -1.0;
   varvals[0][2] = 2.0;
   varvals[0][3] = -2.0;
   checkAD(expr, 1, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(checkad, value)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprValue(scip, &expr, 42.0, NULL, NULL) );

   checkAD(expr, 0, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

#if SCIP_DISABLED_CODE
// CppAD fails to report x^0.3y^0.7 to be non-differentiable at x=y=0, see #2590
// this test reproduces this fail
Test(checkad, issue2590)
{
   SCIP_EXPR* pows[2];
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprPow(scip, &pows[0], varexprs[0], 0.3, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprPow(scip, &pows[1], varexprs[1], 0.7, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 2, pows, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &pows[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &pows[0]) );

   varvals[0][0] = 2.0;
   varvals[0][1] = 2.0;

   varvals[1][0] = 0.0;
   varvals[1][1] = 0.0;
   checkAD(expr, 2, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
#endif

Test(checkad, quad)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* term;
   SCIP_Real coef;
   int i, j;

   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 0, NULL, NULL, 1.0, NULL, NULL) );
   /* make up some sparse matrix and two nonzero points */
   for( i = 0; i < nvars; ++i )
   {
      for( j = 0; j < nvars; ++j )
      {
         coef = (SCIP_Real)(i%4 - j%2) * (i%5 - j%3);
         if( coef == 0.0 )
            continue;
         if( i == j )
         {
            SCIP_CALL( SCIPcreateExprPow(scip, &term, varexprs[i], 2.0, NULL, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPcreateExprProduct(scip, &term, 2, (SCIP_EXPR*[2]){varexprs[i], varexprs[j]}, 1.0, NULL, NULL) );
         }
         SCIP_CALL( SCIPappendExprSumExpr(scip, expr, term, coef) );
         SCIP_CALL( SCIPreleaseExpr(scip, &term) );
      }
      varvals[0][i] = i+1.0;
      varvals[1][i] = nvars-i;
   }

   /* 3rd point is 0 */
   checkAD(expr, nvars, 3);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* https://en.wikipedia.org/wiki/Griewank_function should be there if we test AD */
Test(checkad, griewank)
{
   SCIP_EXPR* exprsum;
   SCIP_EXPR* exprprod;
   int i;

   SCIP_CALL( SCIPcreateExprSum(scip, &exprsum, 0, NULL, NULL, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &exprprod, 0, NULL, 1.0, NULL, NULL) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_EXPR* expr;
      SCIP_EXPR* expr2;
      SCIP_Real coef;

      SCIP_CALL( SCIPcreateExprPow(scip, &expr, varexprs[i], 2.0, NULL, NULL) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, exprsum, expr, 1.0/4000.0) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

      coef = 1.0/sqrt(i+1.0);
      SCIP_CALL( SCIPcreateExprSum(scip, &expr, 1, &varexprs[i], &coef, 0.0, NULL, NULL) );
      SCIP_CALL( SCIPcreateExprCos(scip, &expr2, expr, NULL, NULL) );
      SCIP_CALL( SCIPappendExprChild(scip, exprprod, expr2) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );

      varvals[0][i] = i;
   }

   SCIP_CALL( SCIPappendExprSumExpr(scip, exprsum, exprprod, -1.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprprod) );

   checkAD(exprsum, nvars, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &exprsum) );
}

/* for testing, keep these numbers down
 * but to check timing, set these numbers higher
 */
#define DIM 10      /**< dimension of quadratic coef matrix */
#define NROUNDS 1   /**< number of times to run for the same sparsity (with different sparsity pattern) */
static SCIP_EXPR* vexprs[DIM];
static SCIP_Real vvals[DIM];
static SCIP_Real grad[DIM];
static SCIP_CLOCK* cpuclock;

static
void performance_setup(void)
{
   int i;

   SCIP_CALL_ABORT( SCIPcreate(&scip) );
   SCIP_CALL_ABORT( SCIPincludeDefaultPlugins(scip) );

   /* need a problem to have stat created, which is used by expr iterators */
   SCIP_CALL_ABORT( SCIPcreateProbBasic(scip, "dummy") );

   SCIP_CALL_ABORT( SCIPexprintCreate(scip, &exprint) );

   for( i = 0; i < DIM; ++i )
   {
      SCIP_CALL_ABORT( SCIPcreateExprVaridx(scip, &vexprs[i], i, NULL, NULL) );
   }

   SCIP_CALL_ABORT( SCIPcreateCPUClock(scip, &cpuclock) );
}

static
void performance_teardown(void)
{
   int i;

   SCIP_CALL_ABORT( SCIPfreeClock(scip, &cpuclock) );

   for( i = DIM-1; i >= 0; --i )
   {
      SCIP_CALL_ABORT( SCIPreleaseExpr(scip, &vexprs[i]) );
   }

   SCIPexprintFree(scip, &exprint);

   SCIP_CALL_ABORT( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}


TestSuite(performance, .init = performance_setup, .fini = performance_teardown);

Test(performance, quad)
{
   SCIP_EXPR* expr;
   SCIP_EXPRINTDATA* exprintdata = NULL;
   SCIP_EXPR* term;
   SCIP_Real coef;
   int i, j;
   SCIP_Real val;
   int* rowidxs;
   int* colidxs;
   SCIP_Real* hessianvals;
   int nnz;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_Real fullhessian[DIM*DIM];
   int rnd;
   int sparse;

   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 42, FALSE) );

   printf("%5s %.6s %10.10s %6.6s %6.6s %6.6s %6.6s %6.6s\n", "DIM", "sparsity", "nterms", "compile", "eval", "grad", "hesssparsity", "hessian");

   for( sparse = 0; sparse <= 10; ++sparse )
   {
      SCIP_Real compiletime = 0.0;
      SCIP_Real evaltime = 0.0;
      SCIP_Real gradtime = 0.0;
      SCIP_Real hessiansparsitytime = 0.0;
      SCIP_Real hessiantime = 0.0;
      SCIP_Real sparsity = 1.0/DIM + sparse/10.0 * (1/sqrt(DIM) - 1/DIM);
      /* SCIP_Real sparsity = sparse/10.0; */
      int nterms = 0;

      for( rnd = 1; rnd <= NROUNDS; ++rnd )
      {
         BMSclearMemoryArray(fullhessian, DIM*DIM);

         SCIP_CALL( SCIPcreateExprSum(scip, &expr, 0, NULL, NULL, 1.0, NULL, NULL) );
         /* make up some sparse matrix and a nonzero point */
         for( i = 0; i < DIM; ++i )
         {
            for( j = 0; j < DIM; ++j )
            {
               if( SCIPrandomGetReal(randnumgen, 0.0, 1.0) > sparsity )
                  continue;
               coef = ABS(i - j) + 1.0;
               if( i == j )
               {
                  SCIP_CALL( SCIPcreateExprPow(scip, &term, vexprs[i], 2.0, NULL, NULL) );
               }
               else
               {
                  SCIP_CALL( SCIPcreateExprProduct(scip, &term, 2, (SCIP_EXPR*[2]){vexprs[i], vexprs[j]}, 1.0, NULL, NULL) );
               }
               SCIP_CALL( SCIPappendExprSumExpr(scip, expr, term, coef) );
               SCIP_CALL( SCIPreleaseExpr(scip, &term) );

               fullhessian[MIN(i,j)*DIM+MAX(i,j)] += (i == j) ? 2*coef : coef;
            }
            vvals[i] = i+1.0;
         }
         nterms += SCIPexprGetNChildren(expr);
         /* for( i = 0; i < DIM; ++i )
         {
            for( j = 0; j < DIM; ++j )
               printf(" %5g", fullhessian[i*DIM+j]);
            printf("\n");
          } */

         SCIPresetClock(scip, cpuclock);
         SCIPstartClock(scip, cpuclock);
         SCIP_CALL( SCIPexprintCompile(scip, exprint, expr, &exprintdata) );
         SCIPstopClock(scip, cpuclock);
         compiletime += SCIPgetClockTime(scip, cpuclock);

         if( SCIPexprintGetExprCapability(scip, exprint, expr, exprintdata) & SCIP_EXPRINTCAPABILITY_FUNCVALUE )
         {
            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintEval(scip, exprint, expr, exprintdata, vvals, &val) );
            SCIPstopClock(scip, cpuclock);
            evaltime += SCIPgetClockTime(scip, cpuclock);
         }

         if( SCIPexprintGetExprCapability(scip, exprint, expr, exprintdata) & SCIP_EXPRINTCAPABILITY_GRADIENT )
         {
            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintGrad(scip, exprint, expr, exprintdata, vvals, FALSE, &val, grad) );
            SCIPstopClock(scip, cpuclock);
            gradtime += SCIPgetClockTime(scip, cpuclock);
         }

         if( SCIPexprintGetExprCapability(scip, exprint, expr, exprintdata) & SCIP_EXPRINTCAPABILITY_HESSIAN )
         {
            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintHessianSparsity(scip, exprint, expr, exprintdata, vvals, &rowidxs, &colidxs, &nnz) );
            SCIPstopClock(scip, cpuclock);
            hessiansparsitytime += SCIPgetClockTime(scip, cpuclock);

            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintHessian(scip, exprint, expr, exprintdata, vvals, FALSE, &val, &rowidxs, &colidxs, &hessianvals, &nnz) );
            SCIPstopClock(scip, cpuclock);
            hessiantime += SCIPgetClockTime(scip, cpuclock);

            /* check that Hessian is correct */
            for( i = 0; i < nnz; ++i )
            {
               /* printf("  %d %d = %g\n", colidxs[i], rowidxs[i], hessianvals[i]); */
               fullhessian[colidxs[i]*DIM+rowidxs[i]] -= hessianvals[i];
            }
            for( i = 0; i < DIM*DIM; ++i )
               cr_assert_float_eq(fullhessian[i], 0.0, TOL);
         }

         SCIP_CALL( SCIPexprintFreeData(scip, exprint, expr, &exprintdata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
      }

      printf("%5d %6.2g %10d %6.4f %6.4f %6.4f %6.4f %6.4f\n", DIM, sparsity, nterms/rnd, compiletime/rnd, evaltime/rnd, gradtime/rnd, hessiansparsitytime/rnd, hessiantime/rnd);
   }

   SCIPfreeRandom(scip, &randnumgen);
}

Test(performance, griewank)
{
   int numvars;

   printf("%5s %6.6s %6.6s %6.6s %6.6s %6.6s\n", "nvars", "compile", "eval", "grad", "hesssparsity", "hessian");

   for( numvars = DIM/10; numvars <= DIM; numvars += DIM/10 )
   {
      SCIP_EXPR* exprsum;
      SCIP_EXPR* exprprod;
      SCIP_EXPRINTDATA* exprintdata = NULL;
      SCIP_Real compiletime = 0.0;
      SCIP_Real evaltime = 0.0;
      SCIP_Real gradtime = 0.0;
      SCIP_Real hessiansparsitytime = 0.0;
      SCIP_Real hessiantime = 0.0;
      int rnd;
      int i;

      SCIP_CALL( SCIPcreateExprSum(scip, &exprsum, 0, NULL, NULL, 1.0, NULL, NULL) );
      SCIP_CALL( SCIPcreateExprProduct(scip, &exprprod, 0, NULL, 1.0, NULL, NULL) );

      for( i = 0; i < numvars; ++i )
      {
         SCIP_EXPR* expr;
         SCIP_EXPR* expr2;
         SCIP_Real coef;

         SCIP_CALL( SCIPcreateExprPow(scip, &expr, vexprs[i], 2.0, NULL, NULL) );
         SCIP_CALL( SCIPappendExprSumExpr(scip, exprsum, expr, 1.0/4000.0) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

         coef = 1.0/sqrt(i+1.0);
         SCIP_CALL( SCIPcreateExprSum(scip, &expr, 1, &vexprs[i], &coef, 0.0, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprCos(scip, &expr2, expr, NULL, NULL) );
         SCIP_CALL( SCIPappendExprChild(scip, exprprod, expr2) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );

         vvals[i] = i;
      }

      SCIP_CALL( SCIPappendExprSumExpr(scip, exprsum, exprprod, -1.0) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprprod) );

      for( rnd = 1; rnd <= NROUNDS; ++rnd )
      {
         SCIP_Real val;
         int* rowidxs;
         int* colidxs;
         SCIP_Real* hessianvals;
         int nnz;

         SCIPresetClock(scip, cpuclock);
         SCIPstartClock(scip, cpuclock);
         SCIP_CALL( SCIPexprintCompile(scip, exprint, exprsum, &exprintdata) );
         SCIPstopClock(scip, cpuclock);
         compiletime += SCIPgetClockTime(scip, cpuclock);

         if( SCIPexprintGetExprCapability(scip, exprint, exprsum, exprintdata) & SCIP_EXPRINTCAPABILITY_FUNCVALUE )
         {
            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintEval(scip, exprint, exprsum, exprintdata, vvals, &val) );
            SCIPstopClock(scip, cpuclock);
            evaltime += SCIPgetClockTime(scip, cpuclock);
         }

         if( SCIPexprintGetExprCapability(scip, exprint, exprsum, exprintdata) & SCIP_EXPRINTCAPABILITY_GRADIENT )
         {
            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintGrad(scip, exprint, exprsum, exprintdata, vvals, FALSE, &val, grad) );
            SCIPstopClock(scip, cpuclock);
            gradtime += SCIPgetClockTime(scip, cpuclock);
         }

         if( SCIPexprintGetExprCapability(scip, exprint, exprsum, exprintdata) & SCIP_EXPRINTCAPABILITY_HESSIAN )
         {
            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintHessianSparsity(scip, exprint, exprsum, exprintdata, vvals, &rowidxs, &colidxs, &nnz) );
            SCIPstopClock(scip, cpuclock);
            hessiansparsitytime += SCIPgetClockTime(scip, cpuclock);

            SCIPresetClock(scip, cpuclock);
            SCIPstartClock(scip, cpuclock);
            SCIP_CALL( SCIPexprintHessian(scip, exprint, exprsum, exprintdata, vvals, FALSE, &val, &rowidxs, &colidxs, &hessianvals, &nnz) );
            SCIPstopClock(scip, cpuclock);
            hessiantime += SCIPgetClockTime(scip, cpuclock);
         }

         SCIP_CALL( SCIPexprintFreeData(scip, exprint, exprsum, &exprintdata) );
      }

      printf("%5d %6.4f %6.4f %6.4f %6.4f %6.4f\n", numvars, compiletime/rnd, evaltime/rnd, gradtime/rnd, hessiansparsitytime/rnd, hessiantime/rnd);
      SCIP_CALL( SCIPreleaseExpr(scip, &exprsum) );
   }
}

Test(performance, signpower)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* term;
   SCIP_EXPRINTDATA* exprintdata = NULL;
   SCIP_Real val;
   int* rowidxs;
   int* colidxs;
   SCIP_Real* hessianvals;
   int nnz;
   int i;

   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 0, NULL, NULL, 0.0, NULL, NULL) );

   for( i = 0; i < DIM; ++i )
   {
      SCIP_CALL( SCIPcreateExprSignpower(scip, &term, vexprs[i], i+1.0, NULL, NULL) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, expr, term, 1.0) );
      SCIP_CALL( SCIPreleaseExpr(scip, &term) );
   }

   SCIPresetClock(scip, cpuclock);
   SCIPstartClock(scip, cpuclock);

   SCIP_CALL( SCIPexprintCompile(scip, exprint, expr, &exprintdata) );

   if( SCIPexprintGetExprCapability(scip, exprint, expr, exprintdata) != SCIP_EXPRINTCAPABILITY_ALL )
      goto TERMINATE;

   SCIP_CALL( SCIPexprintHessianSparsity(scip, exprint, expr, exprintdata, vvals, &rowidxs, &colidxs, &nnz) );
   for( i = 0; i < NROUNDS; ++i )
   {
      for( int j = 0; j < DIM; ++j )
         vvals[j] = (j%10) * (i*j % 3 ? 1 : -1);
      SCIP_CALL( SCIPexprintEval(scip, exprint, expr, exprintdata, vvals, &val) );
      SCIP_CALL( SCIPexprintGrad(scip, exprint, expr, exprintdata, vvals, FALSE, &val, grad) );
      SCIP_CALL( SCIPexprintHessian(scip, exprint, expr, exprintdata, vvals, FALSE, &val, &rowidxs, &colidxs, &hessianvals, &nnz) );
   }

   printf("time: %g\n", SCIPgetClockTime(scip, cpuclock));

TERMINATE:
   SCIP_CALL( SCIPexprintFreeData(scip, exprint, expr, &exprintdata) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(performance, intpower)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* term;
   SCIP_EXPRINTDATA* exprintdata = NULL;
   SCIP_Real val;
   int* rowidxs;
   int* colidxs;
   SCIP_Real* hessianvals;
   int nnz;
   int i;

   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 0, NULL, NULL, 0.0, NULL, NULL) );

   for( i = 0; i < DIM; ++i )
   {
      SCIP_CALL( SCIPcreateExprPow(scip, &term, vexprs[i], i+1.0, NULL, NULL) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, expr, term, 1.0) );
      SCIP_CALL( SCIPreleaseExpr(scip, &term) );
   }

   SCIPresetClock(scip, cpuclock);
   SCIPstartClock(scip, cpuclock);

   SCIP_CALL( SCIPexprintCompile(scip, exprint, expr, &exprintdata) );

   if( SCIPexprintGetExprCapability(scip, exprint, expr, exprintdata) != SCIP_EXPRINTCAPABILITY_ALL )
      goto TERMINATE;

   SCIP_CALL( SCIPexprintHessianSparsity(scip, exprint, expr, exprintdata, vvals, &rowidxs, &colidxs, &nnz) );
   for( i = 0; i < NROUNDS; ++i )
   {
      for( int j = 0; j < DIM; ++j )
         vvals[j] = (j%10) * (i*j % 3 ? 1 : -1);
      SCIP_CALL( SCIPexprintEval(scip, exprint, expr, exprintdata, vvals, &val) );
      SCIP_CALL( SCIPexprintGrad(scip, exprint, expr, exprintdata, vvals, FALSE, &val, grad) );
      SCIP_CALL( SCIPexprintHessian(scip, exprint, expr, exprintdata, vvals, FALSE, &val, &rowidxs, &colidxs, &hessianvals, &nnz) );
   }

   printf("time: %g\n", SCIPgetClockTime(scip, cpuclock));

TERMINATE:
   SCIP_CALL( SCIPexprintFreeData(scip, exprint, expr, &exprintdata) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
