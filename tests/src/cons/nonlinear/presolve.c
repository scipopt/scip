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

/**@file   presolve.c
 * @brief  tests presolving methods of the nonlinear constraint handler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_CONSHDLR* conshdlr;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* disable heuristics to avoid that problem is solved by trivial in presolve before constraint upgrading */
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   /* disable all presolving, then reenable nonlinear presolve */
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "constraints/nonlinear/maxprerounds", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", -1) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -10.0, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
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

/*
 * define test suite
 */

TestSuite(presolve, .init = setup, .fini = teardown);

/*
 * define tests
 */

/* tests whether constraints of the form g(x) <= rhs and g(x) >= lhs are merged correctly */
Test(presolve, mergeconss)
{
   SCIP_EXPR* expr;
   SCIP_CONS* conss[3];
   SCIP_Real lhss[3] = {-SCIPinfinity(scip), 1.0, 0.0};
   SCIP_Real rhss[3] = {5.0, SCIPinfinity(scip), 4.0};
   SCIP_Bool success;
   int c;

   /* create expression for each constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x> + <y>*<z>", NULL, NULL, NULL) );

   /* create constraints, always using expr */
   for( c = 0; c < 3; ++c )
   {
      SCIP_CALL( createCons(scip, conshdlr, &conss[c], "cons", expr, lhss[c], rhss[c], FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   }

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* merge constraints */
   SCIP_CALL( presolveMergeConss(scip, conss, 3, &success) );
   cr_expect(success);
   cr_expect(SCIPgetLhsNonlinear(conss[0]) == 1.0);
   cr_expect(SCIPgetRhsNonlinear(conss[0]) == 4.0);
   cr_expect(!SCIPconsIsDeleted(conss[0]));
   cr_expect(SCIPconsIsDeleted(conss[1]));
   cr_expect(SCIPconsIsDeleted(conss[2]));

   /* release constraints */
   for( c = 0; c < 3; ++c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conss[c]) );
   }
}

/** helper method to create, add, and release a nonlinear constraint */
static
SCIP_RETCODE addCons(
   const char*           exprstr,            /**< string with the expr to parse */
   SCIP_Real             lhs,                /**< left-hand side */
   SCIP_Real             rhs                 /**< right-hand side */
   )
{
   SCIP_EXPR* expr;
   SCIP_CONS* cons;

   SCIP_CALL( SCIPparseExpr(scip, &expr, exprstr, NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "cons", expr, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** helper method to check variable types */
static
void checkTypes(
   SCIP_VARTYPE          vartypex,           /**< target variable type of x */
   SCIP_VARTYPE          vartypey,           /**< target variable type of y */
   SCIP_VARTYPE          vartypez            /**< target variable type of z */
   )
{
   cr_expect(SCIPvarGetType(SCIPvarGetTransVar(x)) == vartypex);
   cr_expect(SCIPvarGetType(SCIPvarGetTransVar(y)) == vartypey);
   cr_expect(SCIPvarGetType(SCIPvarGetTransVar(z)) == vartypez);
}

/* test for presolSingleLockedVars() */
Test(presolve, singlelockedvars1)
{
   /* -x^2 + y^2 -z^6 >= 0 implies that x and z are at their bounds */
   SCIP_CALL( addCons("+<x>^2 -<y>^2 +<z>^6 + <x>*<y> + 1/sin(20 + <y>)", 0.5, SCIPinfinity(scip)) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );
   cr_assert(SCIPgetNConss(scip) == 1);

   /* check variable types */
   checkTypes(SCIP_VARTYPE_BINARY, SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_BINARY);
}

/* test for presolSingleLockedVars() */
Test(presolve, singlelockedvars2)
{
   /* same example as singlelockedvars1 but there is another constraint that contains x and y => only z can be at its bounds */
   SCIP_CALL( addCons("+<x>^2 -<y>^2 +<z>^6 + <x>*<y> + 1/sin(20 + <y>)", 0.5, SCIPinfinity(scip)) );
   SCIP_CALL( addCons("<x> * <y> + cos(<x>)", -0.5, 0.5) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );
   cr_assert(SCIPgetNConss(scip) == 2);

   /* check variable types */
   checkTypes(SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_BINARY);
}

/* test for presolSingleLockedVars() */
Test(presolve, singlelockedvars3)
{
   /* change bounds of z */
   SCIP_CALL( SCIPchgVarLbGlobal(scip, z, -1.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, z, 2.0) );

   /* same example as singlelockedvars1; x should be binary and there should be one disjunction constraint for z */
   SCIP_CALL( addCons("+<x>^2 -<y>^2 +<z>^6 + <x>*<y> + 1/sin(20 + <y>)", 0.5, SCIPinfinity(scip)) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );
   cr_expect(SCIPgetNConss(scip) == 2);

   /* check variable types */
   checkTypes(SCIP_VARTYPE_BINARY, SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_CONTINUOUS);
}

/* test for presolveImplint()
 *
 * consider x + 2 y^2 - 3 y^3 - 4 z == 5 with x continuous, y integer, and z binary
 */
Test(presolve, implint)
{
   SCIP_Bool infeasible;

   /* change bounds of x to be [-10,10] */
   SCIP_CALL( SCIPchgVarLbGlobal(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, x, 10.0) );

   /* change variable types */
   SCIP_CALL( SCIPchgVarType(scip, x, SCIP_VARTYPE_CONTINUOUS, &infeasible) );
   SCIP_CALL( SCIPchgVarType(scip, y, SCIP_VARTYPE_INTEGER, &infeasible) );
   SCIP_CALL( SCIPchgVarType(scip, z, SCIP_VARTYPE_BINARY, &infeasible) );

   /* add nonlinear constraint */
   SCIP_CALL( addCons("<x> + 2*<y>^2 - 3*<y>^3 - 4*<z>", 5.0, 5.0) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   /* check variable types */
   checkTypes(SCIP_VARTYPE_IMPLINT, SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_BINARY);
}

/* test for presolveImplint()
 *
 * consider 0.5 x + 0.5 yz - 1.5 z == 5 with x continuous, y and z binary
 */
Test(presolve, implint2)
{
   SCIP_Bool infeasible;

   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/nonlinear/reformbinprods", FALSE) );

   /* change bounds of x to be [-10,10] */
   SCIP_CALL( SCIPchgVarLbGlobal(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, x, 10.0) );

   /* change variable types */
   SCIP_CALL( SCIPchgVarType(scip, x, SCIP_VARTYPE_CONTINUOUS, &infeasible) );
   SCIP_CALL( SCIPchgVarLbGlobal(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarType(scip, y, SCIP_VARTYPE_BINARY, &infeasible) );
   SCIP_CALL( SCIPchgVarType(scip, z, SCIP_VARTYPE_BINARY, &infeasible) );

   /* add nonlinear constraint */
   SCIP_CALL( addCons("0.5*<x> + 0.5*<y>*<z> - 1.5*<z>", 0.0, 0.0) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   /* check variable types */
   checkTypes(SCIP_VARTYPE_IMPLINT, SCIP_VARTYPE_BINARY, SCIP_VARTYPE_BINARY);
}

/* test for presolveImplint()
 *
 * consider 2 x^2 + 0.3 * y - z == 5 with x integer, y continuous, and z binary
 */
Test(presolve, implint3)
{
   SCIP_Bool infeasible;

   /* change bounds of x to be [-10,10] */
   SCIP_CALL( SCIPchgVarLbGlobal(scip, x, -10.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, x, 10.0) );

   /* change variable types */
   SCIP_CALL( SCIPchgVarType(scip, x, SCIP_VARTYPE_INTEGER, &infeasible) );
   SCIP_CALL( SCIPchgVarType(scip, y, SCIP_VARTYPE_CONTINUOUS, &infeasible) );
   SCIP_CALL( SCIPchgVarType(scip, z, SCIP_VARTYPE_BINARY, &infeasible) );

   /* add nonlinear constraint */
   SCIP_CALL( addCons("2*<x>^2 + 0.3*<y> - <z>", 5.0, 5.0) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   /* check variable types */
   checkTypes(SCIP_VARTYPE_INTEGER, SCIP_VARTYPE_CONTINUOUS, SCIP_VARTYPE_BINARY);
}

/* tests setppc upgrade */
Test(presolve, setppcupg)
{
   SCIP_Bool infeasible;

   /* change bounds of all variables */
   SCIP_CALL( SCIPchgVarLbGlobal(scip, x, 0.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, x, 1.0) );
   SCIP_CALL( SCIPchgVarLbGlobal(scip, y, 0.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, y, 1.0) );
   SCIP_CALL( SCIPchgVarObj(scip, x, -1.0) );
   SCIP_CALL( SCIPchgVarObj(scip, y, -1.0) );

   /* change variable types */
   SCIP_CALL( SCIPchgVarType(scip, x, SCIP_VARTYPE_BINARY, &infeasible) );
   SCIP_CALL( SCIPchgVarType(scip, y, SCIP_VARTYPE_BINARY, &infeasible) );

   /* add nonlinear constraint -2xy + 2x + 2y -6 = -4, which is equivalent to (x-1)(y-1) = 0 */
   SCIP_CALL( addCons("-2 * <x> * <y> + 2 * <x> + 2 * <y> - 6", -4.0, -4.0) );

   /* apply presolving */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   cr_expect_eq(SCIPconshdlrGetNUpgdConss(conshdlr), 1, "got %d upgrades, but expected one", SCIPconshdlrGetNUpgdConss(conshdlr));
}
