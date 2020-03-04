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

/**@file   nlhdlr_quotient.c
 * @brief  tests quotient nonlinear handler methods
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr.c"
#include "scip/cons_expr_nlhdlr_quotient.c"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;
static SCIP_VAR* u;

static SCIP_CONSEXPR_NLHDLR* nlhdlr;
static SCIP_CONSHDLR* conshdlr;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get quadratic handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   nlhdlr = SCIPfindConsExprNlhdlr(conshdlr, "quotient");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 1.0, 5.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 2.0, 3.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -4.0, -1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &u, "u", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, u) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &u) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   //BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test suite */
TestSuite(nlhdlrsoc, .init = setup, .fini = teardown);

static
void checkData(
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata,
   SCIP_VAR*             nomvar,
   SCIP_Real             nomcoef,
   SCIP_Real             nomconst,
   SCIP_VAR*             denomvar,
   SCIP_Real             denomcoef,
   SCIP_Real             denomconst,
   SCIP_Real             constant
   )
{
   cr_expect_not_null(nlhdlrexprdata->nomvar);
   cr_expect_not_null(nlhdlrexprdata->denomvar);
   cr_expect(SCIPisEQ(scip, nomcoef, nlhdlrexprdata->nomcoef));
   cr_expect(SCIPisEQ(scip, nomconst, nlhdlrexprdata->nomconst));
   cr_expect(SCIPisEQ(scip, denomcoef, nlhdlrexprdata->denomcoef));
   cr_expect(SCIPisEQ(scip, denomconst, nlhdlrexprdata->denomconst));
   cr_expect(SCIPisEQ(scip, constant, nlhdlrexprdata->constant));
}

static
void checkCut(
   SCIP_ROW*             cut,
   SCIP_VAR**            vars,
   SCIP_Real*            vals,
   SCIP_Real             rhs,
   int                   nvars
   )
{
   int i;

   cr_assert_not_null(cut);
   cr_assert_not_null(vars);
   cr_assert_not_null(vals);

   cr_expect_eq(SCIProwGetNNonz(cut), nvars, "expected %d vars in cut, but got %d\n",
      nvars, SCIProwGetNNonz(cut));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(cut), rhs), "expected rhs = %f, but got %f\n", rhs, SCIProwGetRhs(cut));

   for( i = 0; i < nvars; ++i )
   {
      cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(cut)[i]), vars[i], "expected var%d = %s, but got %s\n",
         i + 1, SCIPvarGetName(vars[i]), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(cut)[i])));
      cr_expect_eq(SCIProwGetVals(cut)[i], vals[i], "expected val%d = %f, but got %f\n", i + 1, vals[i],
         SCIProwGetVals(cut)[i]);
   }
}

/* detects x / y */
Test(nlhdlrsoc, detectandfree1, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: <x> / <y> <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 1.0, 0.0, y, 1.0, 0.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects (4x + 1) / (-3x - 3) */
Test(nlhdlrsoc, detectandfree2, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (4*<x> + 1) / (-3*<x> - 3) <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 1.0, x, -3.0, -3.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log((4x + 3) / (x + 1)) */
Test(nlhdlrsoc, detectandfree3, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: log((4*<x> + 3) / (<x> + 1)) <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons)->children[0];

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 3.0, x, 1.0, 1.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects that (4x + 2y + 3) / (x + 1)) is invalid */
Test(nlhdlrsoc, detectandfree4, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (4*<x> + 2*<y> + 3) / (<x> + 1) <= 10",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_null(nlhdlrexprdata);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log(x) / log(y)) */
Test(nlhdlrsoc, detectandfree5, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_VAR* auxvar1;
   SCIP_VAR* auxvar2;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: log(<x>) / log(<y>) <= 10",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   auxvar1 = SCIPgetConsExprExprAuxVar(expr->children[0]);
   auxvar1 = SCIPgetConsExprExprAuxVar(expr->children[1]->children[0]);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, auxvar1, 1.0, 0.0, auxvar2, 1.0, 0.0, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects (4x + 1) / (-3x + 3) + 2 after simplification */
Test(nlhdlrsoc, detectandfree6, .description = "detects simple quotient expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: ((4*<x> + 1) / (-3*<x> - 3) + 2) <= 3",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, x, 4.0, 1.0, x, -3.0, -3.0, 2.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* tests interval evaluation for ((+/-)4x + 1) / (-3x + 3) - 2*/
Test(nlhdlrsoc, inteval, .description = "tests interval evaluation of simple quotient expression")
{
   SCIP_INTERVAL varbnds;
   SCIP_INTERVAL result;

   /* test interval including 0 in denominator*/

   varbnds.inf = 0.0;
   varbnds.sup = 2.0;

   result = intEval(varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, result));

   /* test positive denominator part for monotone increasing expression */

   varbnds.inf = 2.0;
   varbnds.sup = 9.0;

   result = intEval(varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, -5.0));
   cr_expect(SCIPisEQ(scip, result.sup, -37.0 / 24.0 - 2.0), "expected %f, but got %f\n",
      -37.0 / 24.0 - 2.0, result.sup);

   /* test negative denominator part for monotone increasing expression */

   varbnds.inf = -1.0;
   varbnds.sup = 0.0;

   result = intEval(varbnds, 4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, -2.5));
   cr_expect(SCIPisEQ(scip, result.sup, 1.0 / 3.0 - 2.0));


   /* test positive denominator part for monotone decreasing expression */

   varbnds.inf = 2.0;
   varbnds.sup = 9.0;

   result = intEval(varbnds, -4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, 35.0 / 24.0 - 2.0));
   cr_expect(SCIPisEQ(scip, result.sup, 7.0 / 3.0 - 2.0));

   /* test negative denominator part for monotone decreasing expression */

   varbnds.inf = -1.0;
   varbnds.sup = 0.0;

   result = intEval(varbnds, -4.0, 1.0, -3.0, 3.0, -2.0);

   cr_expect(SCIPisEQ(scip, result.inf, 1.0 / 3.0 - 2.0));
   cr_expect(SCIPisEQ(scip, result.sup, 5.0 / 6.0 - 2.0));

}