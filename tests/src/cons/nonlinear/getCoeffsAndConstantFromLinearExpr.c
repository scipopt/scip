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

/**@file getCoeffsAndConstantFromLinearExpr.c
 * @brief test cons_nonlinear.c getCoeffsAndConstantFromLinearExpr()
 * @author Ingmar Vierhaus
 */

#define SCIP_DEBUG

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/cons_nonlinear.c"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;

SCIP_EXPR* expr;
SCIP_Real scalar = 1.0;
SCIP_Real constant;
SCIP_Real varcoeffs[10];
SCIP_Real* varcoeffsdummy = NULL;


/* creates scip, variables, problem */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -5, 5, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -10, 10, 3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   cr_assert_eq(SCIPgetNVars(scip), 2);
   x = SCIPgetVars(scip)[0];
   y = SCIPgetVars(scip)[1];

}

static
void teardown(void)
{
   SCIPexprFreeDeep(SCIPblkmem(scip), &expr);
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(simple_expressions, .init=setup, .fini = teardown);

Test(simple_expressions, constant)
{
   SCIP_Real exprConstant = 2.0;
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_CONST, exprConstant) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffsdummy, &constant);
   cr_assert_eq(varcoeffsdummy, NULL, "varcoeffsdummy should be NULL (i.e. untouched)");
   cr_assert_eq(constant, exprConstant * scalar, "Wrong constant returned for SCIP_EXPR_CONST");
}

Test(simple_expressions, var)
{
   constant = 0.0;
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffs, &constant);
   cr_assert_eq(varcoeffs[0], scalar, "varcoeffs should be set to scalar");
   cr_assert_eq(constant, 0.0, "constant should be zero (i.e. untouched)");
}

Test(simple_expressions, mulA)
{
   constant = 10;
   scalar = 5.0;
   SCIP_Real constval = 2.0;
   SCIP_EXPR* varexpr;
   SCIP_EXPR* constexpr;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &constexpr, SCIP_EXPR_CONST, constval) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_MUL, constexpr, varexpr) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffs, &constant);
   cr_assert_eq(constant, 10.0, "constant should be untouched");
   cr_assert_eq(varcoeffs[0], scalar * constval, "varcoeffs have wrong value");
}

Test(simple_expressions, mulB)
{
   constant = 10;
   scalar = 5.0;
   SCIP_Real constval = 2.0;
   SCIP_EXPR* varexpr;
   SCIP_EXPR* constexpr;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &constexpr, SCIP_EXPR_CONST, constval) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexpr, SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_MUL, varexpr, constexpr) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffs, &constant);
   cr_assert_eq(constant, 10.0, "constant should be untouched");
   cr_assert_eq(varcoeffs[0], scalar * constval, "varcoeffs have wrong value");
}

Test(simple_expressions, plus)
{
   constant = 0.0;
   scalar = 5.0;
   SCIP_Real constval1 = 2.0;
   SCIP_Real constval2 = 3.0;
   SCIP_EXPR* constexpr1;
   SCIP_EXPR* constexpr2;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &constexpr1, SCIP_EXPR_CONST, constval1) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &constexpr2, SCIP_EXPR_CONST, constval2) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PLUS, constexpr1, constexpr2) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffsdummy, &constant);
   cr_assert_eq(constant, scalar*(constval1 + constval2) , "constant has wrong value");
   cr_assert_eq(varcoeffsdummy, NULL, "varcoeffsdummy should be NULL (i.e. untouched)");
}

Test(simple_expressions, minus)
{
   constant = 0.0;
   scalar = 5.0;
   SCIP_Real constval1 = 2.0;
   SCIP_Real constval2 = 3.0;
   SCIP_EXPR* constexpr1;
   SCIP_EXPR* constexpr2;

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &constexpr1, SCIP_EXPR_CONST, constval1) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &constexpr2, SCIP_EXPR_CONST, constval2) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_MINUS, constexpr1, constexpr2) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffsdummy, &constant);
   cr_assert_eq(constant, scalar*(constval1 - constval2) , "constant has wrong value");
   cr_assert_eq(varcoeffsdummy, NULL, "varcoeffsdummy should be NULL (i.e. untouched)");
}

Test(simple_expressions, linear)
{
   constant = 0.0;
   scalar = 5.0;
   SCIP_Real intercept = 5.0;
   SCIP_Real coefficients[2];
   coefficients[0] = 2.0;
   coefficients[1] = 3.0;
   SCIP_EXPR* varexprs[2];

   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[0], SCIP_EXPR_VARIDX, 0) );
   SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[1], SCIP_EXPR_VARIDX, 1) );
   SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &expr, 2, varexprs, coefficients, intercept) );
   getCoeffsAndConstantFromLinearExpr(expr, scalar, varcoeffs, &constant);
   cr_assert_eq(constant, scalar*intercept , "constant has wrong value");
   for( int i = 0; i < 2; ++i)
      cr_assert_eq(varcoeffs[i], scalar * coefficients[i], "varcoeffs[%i] has wrong value", i);
}
