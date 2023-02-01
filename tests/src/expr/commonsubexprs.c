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

/**@file   commonsubexprs.c
 * @brief  unit test for identifying and handling common subexpressions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;

static
void setup(void)
{
   SCIP_VAR* x;
   SCIP_VAR* y;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
}

static
void teardown(void)
{
   /* free scip and check for memory leaks */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(commonSubexpr, .init = setup, .fini = teardown);

Test(commonSubexpr, singleConstraints)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* children[2];
   SCIP_Bool replacedroot;

   SCIPinfoMessage(scip, NULL, "test single constraint\n");

   SCIP_CALL( SCIPparseExpr(scip, &children[0], "<x>^2 * <y>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPparseExpr(scip, &children[1], "<x>^2 * <y>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 2, children, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPreleaseExpr(scip, &children[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &children[0]) );

   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &expr, 1, &replacedroot) );
   cr_expect_not(replacedroot);
   cr_expect_eq(SCIPexprGetChildren(expr)[0], SCIPexprGetChildren(expr)[1]);

   /* this should not change anything */
   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &expr, 1, &replacedroot) );
   cr_expect_not(replacedroot);
   cr_expect_eq(SCIPexprGetChildren(expr)[0], SCIPexprGetChildren(expr)[1]);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(commonSubexpr, multipleExprs)
{
   SCIP_EXPR* exprs[4];
   SCIP_Bool replacedroot;

   SCIP_CALL( SCIPparseExpr(scip, &exprs[0], "<x> * <y>", NULL, NULL, NULL) );
   SCIP_CALL( SCIPparseExpr(scip, &exprs[1], "exp(<x> * <y>)", NULL, NULL, NULL) );
   SCIP_CALL( SCIPparseExpr(scip, &exprs[2], "abs(exp(<x> * <y>))", NULL, NULL, NULL) );
   SCIP_CALL( SCIPparseExpr(scip, &exprs[3], "log(abs(exp(<x> * <y>)))", NULL, NULL, NULL) );

   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, exprs, 4, &replacedroot) );
   cr_expect_not(replacedroot);

   cr_expect_eq(exprs[0], SCIPexprGetChildren(exprs[1])[0]);
   cr_expect_eq(exprs[0], SCIPexprGetChildren(SCIPexprGetChildren(exprs[2])[0])[0]);
   cr_expect_eq(exprs[0], SCIPexprGetChildren(SCIPexprGetChildren(SCIPexprGetChildren(exprs[3])[0])[0])[0]);

   cr_expect_eq(exprs[1], SCIPexprGetChildren(exprs[2])[0]);
   cr_expect_eq(exprs[1], SCIPexprGetChildren(SCIPexprGetChildren(exprs[3])[0])[0]);

   cr_expect_eq(exprs[2], SCIPexprGetChildren(exprs[3])[0]);

   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[3]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[2]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[1]) );
   SCIP_CALL( SCIPreleaseExpr(scip, &exprs[0]) );
}
