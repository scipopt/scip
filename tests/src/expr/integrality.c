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

/**@file   integrality.c
 * @brief  tests expression handler callbacks for computing integrality information
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

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
   const char*           input,              /**< input string to create an expression */
   SCIP_Bool             isintegral          /**< target integrality information */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPR* origexpr;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* create and print expression */
   cr_expect_eq(SCIPparseExpr(scip, &origexpr, (char*)input, NULL, NULL, NULL), SCIP_OKAY);

   /* simplify expression */
   SCIP_CALL( SCIPsimplifyExpr(scip, origexpr, &expr, &changed, &infeasible, NULL, NULL) );
   SCIP_CALL( SCIPreleaseExpr(scip, &origexpr) );

   /* print simplified expression */
   SCIP_CALL( SCIPprintExpr(scip, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* compute and check integrality information */
   SCIP_CALL( SCIPcomputeExprIntegrality(scip, expr) );
   cr_expect( SCIPexprIsIntegral(expr) == isintegral);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

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
