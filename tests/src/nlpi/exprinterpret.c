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
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   exprinterpret.c
 * @brief  unit test for expression interpreter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "nlpi/exprinterpret.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_EXPRINT* exprint;
static SCIP_SOL* sol;

#define nvars 10
static SCIP_VAR* vars[nvars];
static SCIP_EXPR* varexprs[nvars];
static SCIP_EXPR* varidxexprs[nvars];
static SCIP_Real varvals[nvars];
static unsigned int soltag = 0L;

#define TOL 1e-12

static
void setup(void)
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
   BMSclearMemoryArray(varvals, nvars);

   SCIP_CALL( SCIPexprintCreate(scip, &exprint) );
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPexprintFree(scip, &exprint);

   int i;
   for( i = nvars-1; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &varidxexprs[i]) );
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }

   SCIP_CALL( SCIPfree(&scip) );
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
 * Uses varvals as point where to evaluate.
 */
static
void checkAD(
   SCIP_EXPR*            expr,               /**< expression (in SCIP vars) to be checked */
   int                   dim                 /**< number of variables used in expression */
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

   SCIPinfoMessage(scip, NULL, "checking expression ");
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, " at x =");
   for( i = 0; i < dim; ++i )
      SCIPinfoMessage(scip, NULL, " %g", varvals[i]);
   SCIPinfoMessage(scip, NULL, "\n");

   /* if exprint cannot eval, then nothing we can check */
   if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_FUNCVALUE) )
      return;

   /* get copy of expr using varidx for var exprs */
   SCIP_CALL( SCIPduplicateExpr(scip, expr, &expr2, mapvar, NULL, NULL, NULL) );
   /* let exprint get familiar with expr */
   SCIP_CALL( SCIPexprintCompile(scip, exprint, expr2, &exprintdata) );

   SCIP_CALL( SCIPsetSolVals(scip, sol, nvars, vars, varvals) );

   SCIP_CALL( SCIPevalExpr(scip, expr, sol, ++soltag) );
   SCIP_CALL( SCIPexprintEval(scip, exprint, expr2, exprintdata, varvals, &val) );

   SCIPinfoMessage(scip, NULL, "  expr value: %g  exprint value: %g\n", SCIPexprGetEvalValue(expr), val);
   cr_expect(SCIPisFinite(val) == (SCIPexprGetEvalValue(expr) != SCIP_INVALID));
   if( SCIPisFinite(val) )
      cr_expect_float_eq(val, SCIPexprGetEvalValue(expr), TOL);
   else
      goto TERMINATE;


   if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_GRADIENT) )
      goto TERMINATE;

   /* evaluate gradient with expr and expr2 and compare */
   SCIP_CALL( SCIPevalExprGradient(scip, expr, sol, soltag) );
   SCIP_CALL( SCIPexprintGrad(scip, exprint, expr2, exprintdata, varvals, FALSE, &val, gradient) );

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


   if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_HESSIAN) )
      goto TERMINATE;

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
      goto TERMINATE;

   /* evaluate Hessian (sparse) with expr2 */
   SCIP_CALL( SCIPexprintHessian(scip, exprint, expr2, exprintdata, varvals, FALSE, &val, &hesrowidx, &hescolidx, &hesvalues, &hesnnz) );
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

TERMINATE:
   SCIP_CALL( SCIPexprintFreeData(scip, exprint, expr2, &exprintdata) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
}

TestSuite(exprint, .init = setup, .fini = teardown);

Test(exprint, abs)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprAbs(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0] = 2.0;
   checkAD(expr, 1);

   varvals[0] = -3.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, cos)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprCos(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0] = 0.0;
   checkAD(expr, 1);

   varvals[0] = 1.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, entropy)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr, varexprs[0], NULL, NULL) );

   // FIXME CppAD gives 0 as derivative at the moment, see FIXME there
//   varvals[0] = 0.0;
//   checkAD(expr, 1);

   varvals[0] = 1.0;
   checkAD(expr, 1);

   varvals[0] = 2.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, exp)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprExp(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0] = -1.0;
   checkAD(expr, 1);

   varvals[0] = 1.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, log)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprLog(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0] = 0.0;
   checkAD(expr, 1);

   varvals[0] = 0.1;
   checkAD(expr, 1);

   varvals[0] = 1.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, pow)
{
   SCIP_EXPR* expr;
   SCIP_Real exponents[] = { 2.0, 3.0, 1.0, -1.0, 0.5, 1.875 };
   size_t i;

   for( i = 0; i < sizeof(exponents) / sizeof(SCIP_Real); ++i )
   {
      SCIP_CALL( SCIPcreateExprPow(scip, &expr, varexprs[0], exponents[i], NULL, NULL) );

      varvals[0] = 1.0;
      checkAD(expr, 1);

      varvals[0] = 2.0;
      checkAD(expr, 1);

      // FIXME CppAD claims that x^1.875 is non-diff at x=0
      if( EPSISINT(exponents[i], 0.0) )
      {
         varvals[0] = 0.0;
         checkAD(expr, 1);
      }

      varvals[0] = -3.0;
      checkAD(expr, 1);

      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   }
}

Test(exprint, product)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 2, varexprs, 5.0, NULL, NULL) );

   varvals[0] = 0.0;
   varvals[1] = 0.0;
   checkAD(expr, 2);

   varvals[0] = 0.0;
   varvals[1] = 1.0;
   checkAD(expr, 2);

   varvals[0] = 3.0;
   varvals[1] = 0.0;
   checkAD(expr, 2);

   varvals[0] = 2.0;
   varvals[1] = -4.0;
   checkAD(expr, 2);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, sin)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprSin(scip, &expr, varexprs[0], NULL, NULL) );

   varvals[0] = 0.0;
   checkAD(expr, 1);

   varvals[0] = -1.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, sum)
{
   SCIP_EXPR* expr;
   SCIP_Real coefs[] = { 0.0, 1.0, 2.0, 4.0 };

   SCIP_CALL( SCIPcreateExprSum(scip, &expr, 5, varexprs, coefs, 42.0, NULL, NULL) );

   varvals[0] = 1.0;
   varvals[1] = -1.0;
   varvals[2] = 2.0;
   varvals[3] = -2.0;
   checkAD(expr, 1);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(exprint, value)
{
   SCIP_EXPR* expr;

   SCIP_CALL( SCIPcreateExprValue(scip, &expr, 42.0, NULL, NULL) );

   checkAD(expr, 0);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
