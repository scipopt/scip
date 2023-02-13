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

/**@file   create.c
 * @brief  tests creation of non-linear rows with nonlinear constraints
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONS* cons;
static SCIP_EXPR* expr;
static SCIP_NLROW* nlrow;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* x5;
static SCIP_VAR* tx1;
static SCIP_VAR* tx2;
static SCIP_VAR* tx3;
static SCIP_VAR* tx4;
static SCIP_VAR* tx5;
static const char* input;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "misc/usesymmetry", 0) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x5, "x5", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, x3) );
   SCIP_CALL( SCIPaddVar(scip, x4) );
   SCIP_CALL( SCIPaddVar(scip, x5) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4) );
   SCIP_CALL( SCIPreleaseVar(scip, &x5) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_expect_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* TEST SUITE */
TestSuite(test_create_nlrow, .init = setup, .fini = teardown);

Test(test_create_nlrow, general)
{
   int nvars;
   SCIP_EXPR* varexpr;
   input = "2*<x1> + 3.2*<x2> + 0.5*<x3>^3 - 4*<x4> + <x5> + 10";

   /* create expression from input string */
   SCIP_CALL( SCIPparseExpr(scip, &expr, input, NULL, NULL, NULL) );

   /* add constraint to SCIP and release it */
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "test", expr, 0.0, 100.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* goto solving stage to get NLP setup */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPgetTransformedVar(scip, x1, &tx1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x2, &tx2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x3, &tx3) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x4, &tx4) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x5, &tx5) );

   /* there should be on NLROW for the 1 constraint */
   cr_expect_eq(SCIPgetNNLPNlRows(scip), 1, "%d nlrows in NLP\n", SCIPgetNNLPNlRows(scip));
   nlrow = SCIPgetNLPNlRows(scip)[0];
   cr_assert_not_null(nlrow);

   /* check linear part */
   cr_expect_eq(SCIPnlrowGetConstant(nlrow), 10);
   cr_assert_eq(SCIPnlrowGetNLinearVars(nlrow), 4);
   cr_expect_eq(SCIPnlrowGetLinearVars(nlrow)[0], tx1);
   cr_expect_eq(SCIPnlrowGetLinearVars(nlrow)[1], tx2);
   cr_expect_eq(SCIPnlrowGetLinearVars(nlrow)[2], tx4);
   cr_expect_eq(SCIPnlrowGetLinearVars(nlrow)[3], tx5);

   /* check nonlinear part */
   cr_assert_not_null(SCIPnlrowGetExpr(nlrow));
   SCIP_CALL( SCIPgetExprVarExprs(scip, SCIPnlrowGetExpr(nlrow), &varexpr, &nvars) );
   cr_expect_eq(nvars, 1);
   cr_expect_eq(SCIPgetVarExprVar(varexpr), tx3);

   SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
}

Test(test_create_nlrow, nolin)
{
   int nvars;
   SCIP_EXPR* varexprs[5];
   SCIP_EXPR* rowexpr;
   SCIP_Bool replacedroot;
   input = "2*<x1>^2 + 3.2*<x1>*<x2> + 0.5*<x3>^3 - 4*<x4>*<x5>";

   /* create constraint from input string */
   SCIP_CALL( SCIPparseExpr(scip, &expr, input, NULL, NULL, NULL) );

   /* add constraint to SCIP and release it */
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "test", expr, 0, 1) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* goto solving stage to get NLP setup */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   /* there should be on NLROW for the 1 constraint */
   cr_expect_eq(SCIPgetNNLPNlRows(scip), 1, "%d nlrows in NLP\n", SCIPgetNNLPNlRows(scip));
   nlrow = SCIPgetNLPNlRows(scip)[0];
   cr_assert_not_null(nlrow);

   /* check linear part */
   cr_expect_eq(SCIPnlrowGetConstant(nlrow), 0);
   cr_expect_eq(SCIPnlrowGetNLinearVars(nlrow), 0);

   /* check nonlinear part */
   cr_assert_not_null(SCIPnlrowGetExpr(nlrow));

   rowexpr = SCIPnlrowGetExpr(nlrow);
   SCIP_CALL( SCIPreplaceCommonSubexpressions(scip, &rowexpr, 1, &replacedroot) );
   cr_expect_not(replacedroot);
   SCIP_CALL( SCIPgetExprVarExprs(scip, SCIPnlrowGetExpr(nlrow), varexprs, &nvars) );
   cr_expect_eq(nvars, 5, "got %d vars, but 5 expected", nvars);

   for( --nvars ; nvars >= 0; --nvars )
   {
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[nvars]) );
   }
}
