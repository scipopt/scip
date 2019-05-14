/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reformbinprods.c
 * @brief  tests reformulation of products of binary variables
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;

static
void setup(void)
{
   SCIP_VAR* var;
   int i;

   /* create SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* add default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* create variables */
   for( i = 0; i < 10; ++i )
   {
      char name[SCIP_MAXSTRLEN];

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "x%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, var) );
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(reformbinprods, .init = setup, .fini = teardown);

/** tests the reformulation of binary products during presolve */
Test(reformbinprods, presolve)
{
   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "[expr] <test1>: <x0>[B] * <x1>[B] * <x2>[B] + <x0>[B] * <x1>[B] + <x0>[B] * <x1>[B] + <x0>[B] * <x1>[B] * <x2>[B] <= 1;";
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   int naddconss = 0;

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, "<x0>[B] * <x1>[B] + <x2>[B]", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "c1", expr, 1.0, 1.0) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* go to presolving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(SCIPgetNConss(scip) == 1);

   cons = SCIPgetConss(scip)[0];
   assert(cons != NULL);

   /* we need to call canonicalize here because of the expression locks */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_EXHAUSTIVE, &infeasible, NULL, &naddconss) );
   cr_expect(naddconss == 3, "expect 3 got %d", naddconss);
   cr_expect(SCIPgetNConss(scip) == 4, "expect 4 got %d", SCIPgetNConss(scip));

   /* SCIPwriteTransProblem(scip, "trans.cip", NULL, FALSE); */
}