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

/**@file   integrality.c
 * @brief  tests expression handler callbacks for computing integrality information
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"

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
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -10.0, 10.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/** helper function to check integrality information */
static
SCIP_RETCODE checkIntegrality(
   const char*          input,              /**< input string to create an expression */
   SCIP_Bool            isintegral          /**< target integrality information */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* origexpr;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* create and print expression */
   cr_expect_eq(SCIPparseConsExprExpr(scip, conshdlr, (char*)input, NULL, &origexpr), SCIP_OKAY);

   /* simplify expression */
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, origexpr, &expr, &changed, &infeasible) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &origexpr) );

   /* print simplified expression */
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* compute and check integrality information */
   SCIP_CALL( SCIPcomputeConsExprExprIntegral(scip, conshdlr, expr) );
   cr_expect( SCIPisConsExprExprIntegral(expr) == isintegral);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/*
 * define the testsuite
 */

TestSuite(integrality, .init = setup, .fini = teardown);

/*
 * tests for single expressions
 */

Test(integrality, abs)
{
   SCIP_CALL( checkIntegrality("abs(<x>)", TRUE) );
   SCIP_CALL( checkIntegrality("abs(<y>)", TRUE) );
   SCIP_CALL( checkIntegrality("abs(<z>)", FALSE) );
}

Test(integrality, cos)
{
   SCIP_CALL( checkIntegrality("cos(<x>)", FALSE) );
}

Test(integrality, entropy)
{
   SCIP_CALL( checkIntegrality("entropy(<y>)", FALSE) );
}

Test(integrality, exp)
{
   SCIP_CALL( checkIntegrality("exp(<x>)", FALSE) );
}

Test(integrality, log)
{
   SCIP_CALL( checkIntegrality("log(<x>)", FALSE) );
}

Test(integrality, pow)
{
   SCIP_CALL( checkIntegrality("<x>^2", TRUE) );
   SCIP_CALL( checkIntegrality("<y>^2.2", FALSE) );
   SCIP_CALL( checkIntegrality("<y>^(-2)", FALSE) );
}

Test(integrality, signpower)
{
   SCIP_CALL( checkIntegrality("signpower(<x>,2)", TRUE) );
   SCIP_CALL( checkIntegrality("signpower(<y>,2.2)", FALSE) );
}

Test(integrality, product)
{
   SCIP_CALL( checkIntegrality("<x> * <y>", TRUE) );
   SCIP_CALL( checkIntegrality("-1.1 * <x> * <y>", FALSE) );
   SCIP_CALL( checkIntegrality("<x> * <y> * <z>", FALSE) );
}

Test(integrality, sin)
{
   SCIP_CALL( checkIntegrality("sin(<x>)", FALSE) );
}

Test(integrality, sum)
{
   SCIP_CALL( checkIntegrality("<x> + <y>", TRUE) );
   SCIP_CALL( checkIntegrality("<x> -2.2*<y>", FALSE) );
   SCIP_CALL( checkIntegrality("<x> + <y> + <z>", FALSE) );
}

Test(integrality, value)
{
   SCIP_CALL( checkIntegrality("-2.0", TRUE) );
   SCIP_CALL( checkIntegrality("3.1", FALSE) );
}

Test(integrality, var)
{
   SCIP_CALL( checkIntegrality("<x>", TRUE) );
   SCIP_CALL( checkIntegrality("<y>", TRUE) );
   SCIP_CALL( checkIntegrality("<z>", FALSE) );
}

/*
 * tests for more complex expressions
 */

Test(integrality, quadratic)
{
   SCIP_CALL( checkIntegrality("3*<x>^2 -2 * <x>*<y> + 5 * <y> + 2", TRUE) );
}

Test(integrality, polynomial)
{
   SCIP_CALL( checkIntegrality("<x>^2 * <y>^3 + <y>^2.5", FALSE) );
}
