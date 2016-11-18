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

/**@file   derivative.c
 * @brief  tests expression derivative callbacks
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_var.h"

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
   SCIP_CALL( SCIPrandomCreate(&rndgen, SCIPblkmem(scip), 1) );
}

static
void teardown(void)
{
   SCIPrandomFree(&rndgen);
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
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), -3.0 * pow(2.3,-4));

   /* try an undefined point */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );
   cr_expect_eq(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, x), SCIP_INVALID);

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
