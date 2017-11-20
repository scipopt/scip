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

/**@file   curvature.c
 * @brief  tests curvature expression handler callbacks
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"

#include "include/scip_test.h"


static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 1.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -4.0, -3.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(curvature, .init = setup, .fini = teardown);

/**< auxiliary function for creating an expression and checking its curvature */
static
SCIP_RETCODE checkCurvature(
   const char*          input,              /**< input creating an expression */
   const char*          exprhdlrname,       /**< target expression handler name */
   SCIP_EXPRCURV        expectedcur         /**< expected curvature */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* create an expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &expr), SCIP_OKAY);

   /* check name of the corresponding expression handler */
   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   cr_assert(exprhdlr != NULL);
   cr_expect(strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), exprhdlrname) == 0, "expect nonlinear handler %s, got %s\n",
      exprhdlrname, SCIPgetConsExprExprHdlrName(exprhdlr));

   /* compute curvature */
   SCIP_CALL( SCIPcomputeCurvatureExprExpr(scip, expr) );

   /* check curvature */
   cr_expect(SCIPgetCurvatureExprExpr(expr) == expectedcur, "expect %d, got %d", SCIPgetCurvatureExprExpr(expr), expectedcur);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/* check for absolute expression */
Test(curvature, absolute)
{
   SCIP_CALL( checkCurvature("abs(<x>[C])", "abs", SCIP_EXPRCURV_CONVEX) );
}

/* check for cosine expression */
Test(curvature, cosine)
{
   SCIP_CALL( checkCurvature("cos(<x>[C])", "cos", SCIP_EXPRCURV_UNKNOWN) );

   /* TODO add a test for convex and concave sine expression */
}

/* check for exponential expression */
Test(curvature, exponential)
{
   SCIP_CALL( checkCurvature("exp(<x>[C])", "exp", SCIP_EXPRCURV_CONVEX) );
}

/* check for logarithm expression */
Test(curvature, logarithm)
{
   SCIP_CALL( checkCurvature("log(<x>[C])", "log", SCIP_EXPRCURV_CONCAVE) );
}

/* check for power expression */
Test(curvature, power)
{
   SCIP_CALL( checkCurvature("(<x>[C])^2", "pow", SCIP_EXPRCURV_CONVEX) );

   /* 0 is contained in the interior of x -> neither convex nor concave */
   SCIP_CALL( checkCurvature("(<x>[C])^(-1)", "pow", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("(<x>[C])^(-2)", "pow", SCIP_EXPRCURV_UNKNOWN) );

   /* check cases for y > 0 */
   SCIP_CALL( checkCurvature("(<y>[C])^(-1)", "pow", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("(<y>[C])^(-2.5)", "pow", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("(<y>[C])^(0.5)", "pow", SCIP_EXPRCURV_CONCAVE) );

   /* check cases for z > 0 */
   SCIP_CALL( checkCurvature("(<z>[C])^(-1)", "pow", SCIP_EXPRCURV_CONCAVE) );
   SCIP_CALL( checkCurvature("(<z>[C])^(-2)", "pow", SCIP_EXPRCURV_CONVEX) );
   SCIP_CALL( checkCurvature("(<z>[C])^(-3)", "pow", SCIP_EXPRCURV_CONCAVE) );
}

/* check for product expression */
Test(curvature, product)
{
   SCIP_CALL( checkCurvature("(<x>[C] * <y>[C])", "prod", SCIP_EXPRCURV_UNKNOWN) );
   SCIP_CALL( checkCurvature("(<x>[C] * <y>[C] * <z>[C])", "prod", SCIP_EXPRCURV_UNKNOWN) );
}

/* check for sine expression */
Test(curvature, sine)
{
   SCIP_CALL( checkCurvature("sin(<x>[C])", "sin", SCIP_EXPRCURV_UNKNOWN) );

   /* TODO add a test for convex and concave sine expression */
}

/* check for sum expression */
Test(curvature, sum)
{
   /* sum of linear expressions -> linear */
   SCIP_CALL( checkCurvature("<x>[C] + <y>[C] + <z>[C]", "sum", SCIP_EXPRCURV_LINEAR) );

   /* sum of convex expressions -> convex */
   SCIP_CALL( checkCurvature("(<x>[C])^2 + <y>[C] + <z>[C]", "sum", SCIP_EXPRCURV_CONVEX) );

   /* sum of concave expressions -> concave */
   SCIP_CALL( checkCurvature("<x>[C] + <y>[C]^(0.5) + <z>[C]^(-1)", "sum", SCIP_EXPRCURV_CONCAVE) );

   /* sum of concave and convex expressions -> unknown */
   SCIP_CALL( checkCurvature("<x>[C]^2.5 + <y>[C]^(0.5) + <z>[C]^(-1)", "sum", SCIP_EXPRCURV_UNKNOWN) );
}

/* check for value expression */
Test(curvature, value)
{
   SCIP_CALL( checkCurvature("5.2", "val", SCIP_EXPRCURV_LINEAR) );
}

/* check for variable expression */
Test(curvature, variable)
{
   SCIP_CALL( checkCurvature("<x>[C]", "var", SCIP_EXPRCURV_LINEAR) );
}
