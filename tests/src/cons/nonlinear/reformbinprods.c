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

/**@file   reformbinprods.c
 * @brief  tests reformulation of products of binary variables
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;

static
void setup(void)
{
   SCIP_VAR* var;
   int i;

   /* create SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
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

   /* change reformbinprodsand parameter */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/nonlinear/reformbinprodsand", FALSE) );
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(reformbinprods, .init = setup, .fini = teardown);

/** tests the reformulation for a single product of two binary variables */
Test(reformbinprods, presolve_single_2)
{
   SCIP_EXPR* expr;
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   int naddconss = 0;
   int nchgcoefs = 0;

   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x0>[B] * <x1>[B] + <x2>[B]", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c1", expr, 1.0, 1.0) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* go to presolving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(SCIPgetNConss(scip) == 1);

   cons = SCIPgetConss(scip)[0];
   assert(cons != NULL);

   /* call canonizalize() to replace binary products */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_EXHAUSTIVE, &infeasible, NULL, &naddconss, &nchgcoefs) );
   cr_expect(naddconss == 3, "expect 3 got %d", naddconss);
   cr_expect(SCIPgetNConss(scip) == 4, "expect 4 got %d", SCIPgetNConss(scip));

   /* SCIPwriteTransProblem(scip, "reform.cip", NULL, FALSE); */
}

/** tests the reformulation for a single product of five binary variables */
Test(reformbinprods, presolve_two)
{
   SCIP_EXPR* expr;
   SCIP_CONS* conss[2];
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   int naddconss = 0;
   int nchgcoefs = 0;

   /* create constraint x0 x1 + x2 x3 + sin(x0 x1) <= 1 */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x0>[B] * <x1>[B] + <x2>[B] * <x3>[B] * <x4>[B] + sin(<x0>[B] * <x1>[B])", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c1", expr, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* create constraint x0 x1 + x2 x3 <= 1 */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x0>[B] * <x1>[B] + <x2>[B] * <x3>[B] * <x4>[B]", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c2", expr, 0.0, 1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* go to presolving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(SCIPgetNConss(scip) == 2);

   /* note that we cannot use SCIPgetConss() directly because canonicalize() adds additional constraints */
   conss[0] = SCIPgetConss(scip)[0];
   conss[1] = SCIPgetConss(scip)[1];

   /* call canonizalize() to replace binary products; note that cannonizalize is called once in presolving to replace common subexpressions */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, 2, SCIP_PRESOLTIMING_EXHAUSTIVE, &infeasible, NULL, &naddconss, &nchgcoefs) );
   cr_expect(naddconss == 4, "expect 4 got %d", naddconss);
   cr_expect(SCIPgetNConss(scip) == 6, "expect 6 got %d", SCIPgetNConss(scip));

   /* SCIPwriteTransProblem(scip, "reform.cip", NULL, FALSE); */
}

/** tests the reformulation for a product of two variables that are contained in a clique */
Test(reformbinprods, clique)
{
   SCIP_EXPR* expr;
   SCIP_CONS* conss[2];
   SCIP_CONS* cons;
   SCIP_Bool infeasible;
   int naddconss = 0;
   int nchgcoefs = 0;
   SCIP_VAR* clique_vars[2];
   SCIP_Bool clique_vals[2];
   int nbdchgs;

   /* create constraint x0 x1 + x2 x3 + sin(x0 x1) <= 1 */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x0>[B] * <x1>[B] + <x2>[B] * <x3>[B]", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c1", expr, 0.0, 0.5) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );

   /* go to presolving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   cr_expect(SCIPgetNConss(scip) == 1, "expect 1 got %d", SCIPgetNConss(scip));
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);

   /* add a clique x0 + x1 <= 1 */
   clique_vars[0] = SCIPfindVar(scip, "x0");
   assert(clique_vars[0] != NULL);
   clique_vars[1] = SCIPfindVar(scip, "x1");
   assert(clique_vars[1] != NULL);
   clique_vals[0] = 1;
   clique_vals[1] = 1;
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /*add a clique (1-x2) + (1-x3) <= 1*/
   clique_vars[0] = SCIPfindVar(scip, "x2");
   assert(clique_vars[0] != NULL);
   clique_vars[1] = SCIPfindVar(scip, "x3");
   assert(clique_vars[1] != NULL);
   clique_vals[0] = 0;
   clique_vals[1] = 0;
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /* note that we cannot use SCIPgetConss() directly because canonicalize() adds additional constraints */
   conss[0] = SCIPgetConss(scip)[0];

   /* call canonizalize() to replace binary products; note that canonicalize is called once in presolving to replace common subexpressions */
   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, conss, 1, SCIP_PRESOLTIMING_EXHAUSTIVE, &infeasible, NULL, &naddconss, &nchgcoefs) );
   cr_expect(naddconss == 0, "expect 0 got %d", naddconss);
   cr_expect(SCIPgetNConss(scip) == 1, "expect 1 got %d", SCIPgetNConss(scip));
   cr_expect(nchgcoefs == 4, "expect 4 changed coefs, got %d", nchgcoefs);
}

/** tests the reformulation of binary quadratic expressions when factorzing variables */
Test(reformbinprods, factorize1)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* newexpr;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   int naddconss = 0;

   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x0>[B] * <x1>[B] - <x0>[B] * <x2>[B] - <x0>[B] * <x3>[B]", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c1", expr, 0.0, 0.5) );

   /* not enough terms -> nothing should happen */
   SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, cons, expr, 4, &newexpr, NULL) );
   cr_assert(newexpr == NULL);

   SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, cons, expr, 3, &newexpr, &naddconss) );
   cr_assert(newexpr != NULL);
   cr_expect(naddconss == 4);
   cr_expect(SCIPexprGetNChildren(newexpr) == 1);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(newexpr)[0]));

   /* newexpr is a sum with only one variable; the bounds of the variable correspond to the activities of the bilinear binary terms */
   var = SCIPgetVarExprVar(SCIPexprGetChildren(newexpr)[0]);
   cr_assert(var != NULL);
   cr_expect(SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
   cr_expect(SCIPvarGetLbGlobal(var) == -2.0);
   cr_expect(SCIPvarGetUbGlobal(var) == 1.0);

   SCIP_CALL( SCIPreleaseExpr(scip, &newexpr) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/** tests the reformulation of binary quadratic expressions when factorzing variables */
Test(reformbinprods, factorize2)
{
   SCIP_EXPR* expr;
   SCIP_EXPR* newexpr;
   SCIP_CONS* cons;
   SCIP_VAR* var;
   int naddconss = 0;

   /* expression is equivalent to  x0 * (x1 - x2 + x7) + x3 * (-x4 -x5) + sin(x0) */
   SCIP_CALL( SCIPparseExpr(scip, &expr, "<x0>[B] * <x1>[B] + <x0>[B] * <x7>[B] - <x4>[B] * <x3>[B] + sin(<x0>[B]) - <x0>[B] * <x2>[B] - <x3>[B] * <x5>[B]", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, "c1", expr, 0.0, 0.5) );

   SCIP_CALL( getFactorizedBinaryQuadraticExpr(scip, conshdlr, cons, expr, 2, &newexpr, &naddconss) );
   cr_assert(newexpr != NULL);
   cr_expect(naddconss == 7);
   cr_expect(SCIPexprGetNChildren(newexpr) == 3);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(newexpr)[0]));
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(newexpr)[1]));

   /* first variable represents x0 * (x1 - x2 + x7) and thus has bounds [-1,2]*/
   var = SCIPgetVarExprVar(SCIPexprGetChildren(newexpr)[0]);
   cr_expect(SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
   cr_expect(SCIPvarGetLbGlobal(var) == -1.0);
   cr_expect(SCIPvarGetUbGlobal(var) == 2.0);

   /* second variable represents  x3 * (-x4 -x5) and thus has bounds [-2,0]*/
   var = SCIPgetVarExprVar(SCIPexprGetChildren(newexpr)[1]);
   cr_expect(SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
   cr_expect(SCIPvarGetLbGlobal(var) == -2.0);
   cr_expect(SCIPvarGetUbGlobal(var) == 0.0);

   /* release memory */
   SCIP_CALL( SCIPreleaseExpr(scip, &newexpr) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
