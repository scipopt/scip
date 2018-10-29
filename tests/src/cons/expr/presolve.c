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

/**@file   presolve.c
 * @brief  tests presolving methods of the expression constraint handler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.c"

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
