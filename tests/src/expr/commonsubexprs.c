/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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
