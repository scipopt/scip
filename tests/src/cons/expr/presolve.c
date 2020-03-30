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

/**@file   presolve.c
 * @brief  tests presolving methods of the expression constraint handler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.c"
#include "scip/cons_bounddisjunction.h"

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
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* conss[3];
   SCIP_Real lhss[3] = {-SCIPinfinity(scip), 1.0, 0.0};
   SCIP_Real rhss[3] = {5.0, SCIPinfinity(scip), 4.0};
   SCIP_Bool success;
   int c;

   /* create expression for each constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x> + <y>*<z>", NULL, &expr) );

   /* create constraints */
   for( c = 0; c < 3; ++c )
   {
      SCIP_CALL( SCIPcreateConsExprBasic(scip, &conss[c], "cons", expr, lhss[c], rhss[c]) );
   }

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   /* merge constraints */
   SCIP_CALL( presolMergeConss(scip, conss, 3, &success) );
   cr_expect(success);
   cr_expect(SCIPgetLhsConsExpr(scip, conss[0]) == 1.0);
   cr_expect(SCIPgetRhsConsExpr(scip, conss[0]) == 4.0);
   cr_expect(!SCIPconsIsDeleted(conss[0]));
   cr_expect(SCIPconsIsDeleted(conss[1]));
   cr_expect(SCIPconsIsDeleted(conss[2]));

   /* release constraints */
   for( c = 0; c < 3; ++c )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &conss[c]) );
   }
}

/** helper method to create, add, and release an expression constraint */
static
SCIP_RETCODE addCons(
   const char*           exprstr,            /**< string with the expr to parse */
   SCIP_Real             lhs,                /**< left-hand side */
   SCIP_Real             rhs                 /**< right-hand side */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons;

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, exprstr, NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "cons", expr, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

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
   /* include bound disjunction constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBounddisjunction(scip) );

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