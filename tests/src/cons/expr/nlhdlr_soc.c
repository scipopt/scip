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

/**@file   nlhdlr_soc.c
 * @brief  tests quadratic nonlinear handler methods
 * @author Fabian Wegscheider
 *
 * @TODO: make the coefficient and variables checks independent of their order
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr.c"
#include "scip/cons_expr_nlhdlr_soc.c"


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

static SCIP_CONSEXPR_EXPR* xexpr;
static SCIP_CONSEXPR_EXPR* yexpr;
static SCIP_CONSEXPR_EXPR* wexpr;
static SCIP_CONSEXPR_EXPR* zexpr;
static SCIP_CONSEXPR_EXPR* uexpr;

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

   nlhdlr = SCIPfindConsExprNlhdlr(conshdlr, "soc");
   cr_assert_not_null(nlhdlr);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* the hacky way how conss are processed here give an assert in nlhdlrbilinearexit (because exitsepa is not called, I guess) */
   SCIP_CALL( SCIPsetBoolParam(scip, "constraints/expr/nlhdlr/bilinear/enabled", FALSE) );

   /* go to SOLVING stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -0.9, 0.7, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &u, "u", 0.91, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPaddVar(scip, u) );

   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &xexpr, x) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &yexpr, y) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &wexpr, w) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &zexpr, z) );
   SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &uexpr, u) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &uexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &zexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &wexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &xexpr) );

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
   SCIP_CONSEXPR_EXPR**  vars,
   SCIP_Real*            offsets,
   SCIP_Real*            transcoefs,
   int*                  transcoefsidx,
   int*                  nnonzeroes,
   int                   nvars,
   int                   nterms,
   int                   ntranscoefs
   )
{
   int i;

   cr_assert_not_null(nlhdlrexprdata->vars);
   cr_assert_not_null(nlhdlrexprdata->offsets);
   cr_assert_not_null(nlhdlrexprdata->transcoefs);
   cr_assert_not_null(nlhdlrexprdata->transcoefsidx);

   cr_expect_eq(nlhdlrexprdata->nvars, nvars);
   cr_expect_eq(nlhdlrexprdata->nterms, nterms);
   cr_expect_eq(nlhdlrexprdata->termbegins[nlhdlrexprdata->nterms], ntranscoefs);

   for( i = 0; i < nvars; ++i )
   {
      cr_assert_not_null(nlhdlrexprdata->vars[i]);
      cr_expect_eq(nlhdlrexprdata->vars[i], vars[i], "expected expr %d to be %p, but got %p\n",
         i + 1, vars[i], nlhdlrexprdata->vars[i]);
   }

   /* FIXME the remaining tests assume a certain order and sign and thus are not invariant to
    * valid permutations or changes in sign (multiplying a whole term by -1)
    */
   return;

   for( i = 0; i < nterms; ++i )
   {
      cr_expect(SCIPisEQ(scip, nlhdlrexprdata->offsets[i], offsets[i]), "expected offset %d to be %f, but got %f\n",
         i + 1, offsets[i], nlhdlrexprdata->offsets[i]);
   }

   for( i = 0; i < nterms; ++i )
   {
      int nnz = nlhdlrexprdata->termbegins[i + 1] - nlhdlrexprdata->termbegins[i];
      cr_expect_eq(nnz, nnonzeroes[i], "expected nnonzeroes %d to be %d, but got %d\n",
         i + 1, nnonzeroes[i], nnz);
   }

   for( i = 0; i < ntranscoefs; ++i )
   {
      cr_expect(SCIPisEQ(scip, nlhdlrexprdata->transcoefs[i], transcoefs[i]), "expected transcoef[%d] to be %f, but got %f\n",
         i + 1, transcoefs[i], nlhdlrexprdata->transcoefs[i]);
   }

   for( i = 0; i < ntranscoefs; ++i )
   {
      cr_expect_eq(nlhdlrexprdata->transcoefsidx[i], transcoefsidx[i], "expected transcoefsidx[%d] to be %d, but got %d\n",
         i + 1, transcoefsidx[i], nlhdlrexprdata->transcoefsidx[i]);
   }
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

   /* FIXME the remaining tests assume a certain order of terms and thus are not invariant to
    * valid permutations
    */
   return;

   for( i = 0; i < nvars; ++i )
   {
      cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(cut)[i]), vars[i], "expected var%d = %s, but got %s\n",
         i + 1, SCIPvarGetName(vars[i]), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(cut)[i])));
      cr_expect_eq(SCIProwGetVals(cut)[i], vals[i], "expected val%d = %f, but got %f\n", i + 1, vals[i],
         SCIProwGetVals(cut)[i]);
   }
}

/* norm <= constant shouldn't be handled by soc */
Test(nlhdlrsoc, detectandfree1, .description = "detects simple norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (<x>^2 + <y>^2 + <z>^2)^0.5 <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects ||x|| - y <= 0 as soc expression */
Test(nlhdlrsoc, detectandfree2, .description = "detects simple norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* normexpr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (<x>^2 + <y>^2 + <z>^2)^0.5 - <w> <= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);
   normexpr = SCIPgetConsExprExprChildren(expr)[0];

   /* find the nlhdlr expr data */
   for( i = 0; i < normexpr->nenfos; ++i )
   {
      if( normexpr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = normexpr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[4] = {xexpr, yexpr, zexpr, normexpr};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[4] = {1.0, 1.0, 1.0, 1.0};
   int transcoefsidx[4] = {0, 1, 2, 3};
   int nnonzeroes[4] = {1, 1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 4);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects SQRT(8 +  2*(x + 1)^2 + 4*(sin(y) - 2)^2 ) as soc expression */
Test(nlhdlrsoc, detectandfree3, .description = "detects more complex norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* normexpr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (8 + 2*(<x> + 1)^2 + 4*(sin(<y>) - 2)^2)^0.5 + 2*(<w> - 1) <= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   expr = SCIPgetExprConsExpr(scip, cons);
   normexpr = SCIPgetConsExprExprChildren(expr)[1];

   /* find the nlhdlr expr data */
   for( i = 0; i < normexpr->nenfos; ++i )
   {
      if( normexpr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = normexpr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[3] = {xexpr, normexpr->children[0]->children[2], normexpr};
   int nterms = 4;
   SCIP_Real offsets[4] = {SQRT(2.0), -4.0, SQRT(8.0) /* constant */, 0.0};
   int nnonzeroes[4] = {1, 1, 0, 1};
   SCIP_Real transcoefs[3] = {SQRT(2.0), 2.0, 1.0};
   int transcoefsidx[3] = {0, 1, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, nterms, 3);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects 2x^2 - 9y^2 + sin(z)^2 < = 0 as soc expression */
Test(nlhdlrsoc, detectandfree4, .description = "detects simple quadratic expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: 2*<x>^2 - 9*<y>^2 + sin(<z>)^2  <= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[3] = {xexpr, expr->children[2]->children[0], yexpr};
   SCIP_Real offsets[3] = {0.0, 0.0, 0.0};
   SCIP_Real transcoefs[3] = {SQRT(2.0), 1.0, 3.0};
   int transcoefsidx[3] = {0, 1, 2};
   int nnonzeroes[3] = {1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 3);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects 5 - 9cos(x)^2 + y^2 + 2sin(z)^2 <= 4 as soc expression */
Test(nlhdlrsoc, detectandfree5, .description = "detects more complication quadratic expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: 5 - 9*cos(<x>)^2 + <y>^2 + 2*sin(<z>)^2 <= 4",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[3] = {yexpr, expr->children[2]->children[0], expr->children[1]->children[0]};
   int nterms = 4;
   SCIP_Real offsets[4] = {0.0, 0.0, 1.0, 0.0};
   int nnonzeroes[4] = {1, 1, 0, 1};
   SCIP_Real transcoefs[3] = {1.0, SQRT(2.0), 3.0};
   int transcoefsidx[3] = {0, 1, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, nterms, 3);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects that -7exp(x)^2 + y^2 - 2sin(z)^2 <= 1 is not a soc expression */
Test(nlhdlrsoc, detectandfree6, .description = "detects quadratic expression that is not soc")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: - 7*cos(<x>)^2 + <y>^2 - 2*sin(<z>)^2  <= 5",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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
   cr_assert_null(nlhdlrexprdata);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects 1 + 7x^2 - 2*u*y + 2z^2 <= 0 as soc expression */
Test(nlhdlrsoc, detectandfree7, .description = "detects hyperbolic quadratic expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: 1 + 8*<x>^2 - 2*<y>*<u> + 2*<z>^2 <= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[4] = {xexpr, zexpr, yexpr, uexpr};
   int nterms = 5;
   SCIP_Real offsets[5] = {0.0, 0.0, SQRT(2.0), 0.0, 0.0}; /* hyperbolic, constant is second to last of lhs */
   SCIP_Real transcoefs[6] = {4.0, 2.0, 1.0, -1.0, 1.0, 1.0};
   int transcoefsidx[6] = {0, 1, 2, 3, 2, 3};
   int nnonzeroes[5] = {1, 1, 0, 2, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, nterms, 6);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects that 7x^2 - 2*u*y - 2z^2 <= 0 is no expression */
Test(nlhdlrsoc, detectandfree8, .description = "detects hyperbolic quadratic expression that is not soc")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: 7*<x>^2 - 2*<y>*<u> - 2*<z>^2 <= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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
   cr_assert_null(nlhdlrexprdata);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects -5 + 7cos(x)^2 - y^2 - 2sin(z)^2 >= -4 as soc expression */
Test(nlhdlrsoc, detectandfree9, .description = "detects negated quadratic expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: -5 + 9*cos(<x>)^2 - <y>^2 - 2*sin(<z>)^2 >= -4",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[3] = {yexpr, expr->children[2]->children[0], expr->children[1]->children[0]};
   int nterms = 4;
   SCIP_Real offsets[4] = {0.0, 0.0, 1.0, 0.0};
   int nnonzeroes[4] = {1, 1, 0, 1};
   SCIP_Real transcoefs[3] = {1.0, SQRT(2.0), 3.0};
   int transcoefsidx[3] = {0, 1, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, nterms, 3);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects -1 -7x^2 + 2*u*y - 2z^2 >= 0 as soc expression */
Test(nlhdlrsoc, detectandfree10, .description = "detects negated hyperbolic quadratic expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: -1 - 8*<x>^2 + 2*<y>*<u> - 2*<z>^2 >= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[4] = {xexpr, zexpr, yexpr, uexpr};
   int nterms = 5;
   SCIP_Real offsets[5] = {0.0, 0.0, SQRT(2.0), 0.0, 0.0};
   int nnonzeroes[5] = {1, 1, 0, 2, 2};
   SCIP_Real transcoefs[6] = {4.0, 2.0, 1.0, -1.0, 1.0, 1.0};
   int transcoefsidx[6] = {0, 1, 2, 3, 2, 3};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, nterms, 6);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects x^2 + 5xy + 5xz + 2y^2 + 2yz + 3z^2 + 8 <= 0 as soc expression */
Test(nlhdlrsoc, detectandfree11, .description = "detects complex quadratic constraint")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons,
         (char*) "[expr] <test>: <x>^2 + 2*<y>^2 + 3*<z>^2 + 5*<x>*<y> + 5*<x>*<z> + 2*<y>*<z> + 10*<x> + 8 <= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[3] = {zexpr, yexpr, xexpr};
   int nterms = 4;
   SCIP_Real offsets[4] = {-0.5321193159399078, -1.1644355645349063, SQRT(17.0909090909091), 3.275663328435811};
   int nnonzeroes[4] = {3, 3, 0, 3};
   SCIP_Real transcoefs[9] = { 0.8403540860560068, -0.8715100046656585, -0.15866733787070855,
                              -1.5723379860224058, -1.2601695810081268, -1.4058990135251308,
                              0.42242364072103755,  0.5895397027601558, -1.0008633075190196};
   int transcoefsidx[9] = {0, 1, 2, 0, 1, 2, 0, 1, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, nterms, 9);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects -x^2 - 5xy - 5xz - 2y^2 - 2yz - 3z^2 - 8 >= 0 as soc expression */
Test(nlhdlrsoc, detectandfree12, .description = "detects complex quadratic constraint")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons,
         (char*) "[expr] <test>: -<x>^2 - 2*<y>^2 - 3*<z>^2 - 5*<x>*<y> - 5*<x>*<z> - 2*<y>*<z> - 10*<x> - 8 >= 0",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* setup expected data */
   SCIP_CONSEXPR_EXPR* vars[3] = {zexpr, yexpr, xexpr};
   int nterms = 4;
   int nnonzeroes[4] = {3, 3, 0, 3};
   SCIP_Real offsets[4] = {-1.1644355645349063, 0.5321193159399078, SQRT(17.0909090909091), 3.275663328435811};
   SCIP_Real transcoefs[9] = {-1.5723379860224058, -1.2601695810081268, -1.4058990135251308,
                              -0.8403540860560068,  0.8715100046656585,  0.15866733787070855,
                              0.42242364072103755,  0.5895397027601558, -1.0008633075190196};
   int transcoefsidx[9] = {0, 1, 2, 0, 1, 2, 0, 1, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, nterms, 9);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects that SQRT(x^2 -4x + 1) <= 2 is not a soc expression */
Test(nlhdlrsoc, detectandfree13, .description = "detects complex quadratic constraint")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons,
         (char*) "[expr] <test>: (<x>^2 - 4*<x> + 1)^0.5 <= 1",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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
   cr_assert_null(nlhdlrexprdata);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects x^2 + y^2 - z*u <= -2 with z not nonnegative but z*u nonnegative */
Test(nlhdlrsoc, detectandfree14, .description = "detects complex quadratic constraint")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons,
         (char*) "[expr] <test>: <x>^2 + <y>^2 - <z> * <u> <= -2",
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

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

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}


/* disaggregates SQRT( 8 + 2*(x + 1)^2 + 3*(y + sin(x) + 2)^2 ) <= -2*(w - 1) */
Test(nlhdlrsoc, disaggregation, .description = "disaggregate soc and check the resulting datastructure")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* normexpr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (8 + 2*(<x> + 1)^2 + 3*(sin(<y>) - 2)^2)^0.5 + 2*(<w> - 1) <= 0", TRUE, TRUE,
            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect_not(infeasible);

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   /* call the separation initialization -> this creates auxvars and creates disaggregation variables and row */
   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeasible) );

   expr = SCIPgetExprConsExpr(scip, cons);
   normexpr = SCIPgetConsExprExprChildren(expr)[1];

   /* find the nlhdlr expr data */
   for( i = 0; i < normexpr->nenfos; ++i )
   {
      if( normexpr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = normexpr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* check disvars */
   cr_expect_not_null(nlhdlrexprdata->disvars[0]);
   cr_expect_not_null(nlhdlrexprdata->disvars[1]);
   cr_expect_not_null(nlhdlrexprdata->disvars[2]);

   cr_expect_eq(SCIProwGetNNonz(nlhdlrexprdata->disrow), 4);

   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[0], 1.0, "expected %f, but got %f\n", 1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[0]);
   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[1], 1.0, "expected %f, but got %f\n", 1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[1]);
   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[2], 1.0, "expected %f, but got %f\n", 1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[2]);
   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[3], -1.0, "expected %f, but got %f\n", -1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[3]);

   cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[0]), nlhdlrexprdata->disvars[0], "expected <%s>, but got <%s>\n", SCIPvarGetName(nlhdlrexprdata->disvars[0]), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[0])));
   cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[1]), nlhdlrexprdata->disvars[1], "expected <%s>, but got <%s>\n", SCIPvarGetName(nlhdlrexprdata->disvars[1]), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[1])));
   cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[2]), nlhdlrexprdata->disvars[2], "expected <%s>, but got <%s>\n", SCIPvarGetName(nlhdlrexprdata->disvars[2]), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[2])));
   cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[3]), SCIPgetConsExprExprAuxVar(nlhdlrexprdata->vars[2]), "expected <%s>, but got <%s>\n", SCIPvarGetName(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->vars[2])), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[3])));

   cr_expect_eq(SCIProwGetLhs(nlhdlrexprdata->disrow), -SCIPinfinity(scip));
   cr_expect_eq(SCIProwGetRhs(nlhdlrexprdata->disrow), 0.0, "expected 0 got %g\n", SCIProwGetRhs(nlhdlrexprdata->disrow));

   /* free row, expr, and cons (freeing of the row doesn't happen automatically because we call createDisaggrRow
    * directly and not through the initlp callback. Thus, the enforce doesn't know that sepainit has been called and so
    * it doesn't call sepaexit in freeEnfoData)
    */
   SCIPreleaseRow(scip, &nlhdlrexprdata->disrow);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* clear sepastorage to get rid of disaggregation row */
   SCIP_CALL( SCIPclearCuts(scip) );
}

/* separates simple norm function from different points */
Test(nlhdlrsoc, separation1, .description = "test separation for simple norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* rootexpr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_SOL* sol;
   SCIP_ROW* cut;
   SCIP_VAR* cutvars[3];
   SCIP_VAR* auxvar;
   SCIP_Real cutvals[3];
   SCIP_Bool infeasible;
   SCIP_Real rhs;
   int i;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*) "exp((<x>^2 + <y>^2 + <z>^2)^0.5)", NULL, &rootexpr) );

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "soc", rootexpr, -SCIPinfinity(scip), 2.0) );

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   /* call the separation initialization -> this creates auxvars and creates disaggregation variables and row */
   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeasible) );

   expr = rootexpr->children[0];

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   auxvar = SCIPgetConsExprExprAuxVar(expr);

   /* create solution */
   SCIPcreateSol(scip, &sol, NULL);
   SCIPsetSolVal(scip, sol, x, 1.0);
   SCIPsetSolVal(scip, sol, y, 2.0);
   SCIPsetSolVal(scip, sol, z, -2.0);
   SCIPsetSolVal(scip, sol, nlhdlrexprdata->disvars[0], 0.0);
   SCIPsetSolVal(scip, sol, nlhdlrexprdata->disvars[1], 1.0);
   SCIPsetSolVal(scip, sol, nlhdlrexprdata->disvars[2], 1.0);
   SCIPsetSolVal(scip, sol, auxvar, 2.0);

   /* check cut w.r.t. x */
   SCIP_CALL( generateCutSolDisagg(scip, expr, cons, sol, nlhdlrexprdata, 0, 0.0, 2.0, &cut) );

   cutvars[0] = nlhdlrexprdata->disvars[0];
   cutvars[1] = x;
   cutvars[2] = auxvar;
   cutvals[0] = -2.0 / SQRT(8.0) - 1.0;
   cutvals[1] = 4.0 / SQRT(8.0);
   cutvals[2] = 2.0 / SQRT(8.0) - 1.0;
   rhs = cutvals[1] + cutvals[2] * 2.0 - SQRT(8.0) + 2.0;

   checkCut(cut, cutvars, cutvals, rhs, 3);
   SCIPreleaseRow(scip, &cut);

   /* check cut w.r.t. y */
   SCIP_CALL( generateCutSolDisagg(scip, expr, cons, sol, nlhdlrexprdata, 1, 0.0, 2.0, &cut) );

   cutvars[0] = auxvar;
   cutvars[1] = y;
   cutvars[2] = nlhdlrexprdata->disvars[1];
   cutvals[0] = 1.0 / SQRT(17.0) - 1.0;
   cutvals[1] = 8.0 / SQRT(17.0);
   cutvals[2] = -1.0 / SQRT(17.0) - 1.0;
   rhs =  cutvals[0] * 2.0 + cutvals[1] * 2.0 + cutvals[2] * 1.0 - SQRT(17.0) + 3.0;

   checkCut(cut, cutvars, cutvals, rhs, 3);
   SCIPreleaseRow(scip, &cut);

   /* check cut w.r.t. z */
   SCIP_CALL( generateCutSolDisagg(scip, expr, cons, sol, nlhdlrexprdata, 2, 0.0, 2.0, &cut) );

   cutvars[0] = auxvar;
   cutvars[1] = z;
   cutvars[2] = nlhdlrexprdata->disvars[2];
   cutvals[0] = 1.0 / SQRT(17.0) - 1.0;
   cutvals[1] = -8.0 / SQRT(17.0);
   cutvals[2] = -1.0 / SQRT(17.0) - 1.0;
   rhs =  cutvals[0] * 2.0 - cutvals[1] * 2.0 + cutvals[2] * 1.0 - SQRT(17.0) + 3.0;

   checkCut(cut, cutvars, cutvals, rhs, 3);
   SCIPreleaseRow(scip, &cut);

   /* free expr and cons */
   SCIPfreeSol(scip, &sol);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &rootexpr) );
}

/* separates simple norm function from different points */
Test(nlhdlrsoc, separation2, .description = "test separation for simple norm expression without disagg")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* rootexpr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_SOL* sol;
   SCIP_ROW* cut;
   SCIP_VAR* cutvars[3];
   SCIP_VAR* auxvar;
   SCIP_Real cutvals[3];
   SCIP_Bool infeasible;
   SCIP_Bool changed;
   SCIP_Real rhs;
   int i;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*) "exp( ((<x> + 0.5)^2 + <y>^2)^0.5 )", NULL, &rootexpr) );

   SCIPsimplifyConsExprExpr(scip, conshdlr, rootexpr, &expr, &changed, &infeasible);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &rootexpr) );
   rootexpr = expr;

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "soc", rootexpr, -SCIPinfinity(scip), 2.0) );

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   /* call the separation initialization -> this creates auxvars and creates disaggregation variables and row */
   SCIP_CALL( initSepa(scip, conshdlr, &cons, 1, &infeasible) );

   /* get norm expression */
   expr = SCIPgetExprConsExpr(scip, cons)->children[0];

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   auxvar = SCIPgetConsExprExprAuxVar(expr);

   /* create solution */
   SCIPcreateSol(scip, &sol, NULL);
   SCIPsetSolVal(scip, sol, x, 1.0);
   SCIPsetSolVal(scip, sol, y, 2.0);
   SCIPsetSolVal(scip, sol, auxvar, 1.0);

   /* check cut */
   SCIP_CALL( generateCutSolSOC(scip, expr, cons, sol, nlhdlrexprdata, 0.0, 1.0, &cut) );

   cutvars[0] = auxvar;
   cutvars[1] = y;
   cutvars[2] = x;
   cutvals[0] = -1.0;
   cutvals[1] = 4.0 / 5;
   cutvals[2] = 3.0 / 5;
   rhs = -3.0 / 10;

   checkCut(cut, cutvars, cutvals, rhs, 3);
   SCIPreleaseRow(scip, &cut);

   /* free expr and cons */
   SCIPfreeSol(scip, &sol);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &rootexpr) );
}

/* separates simple function */
Test(nlhdlrsoc, separation3, .description = "test separation for simple expression without disagg")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_SOL* sol;
   SCIP_ROW* cut;
   SCIP_VAR* cutvars[3];
   SCIP_Real cutvals[3];
   SCIP_Bool infeasible;
   SCIP_Real rhs;
   int i;

   SCIP_CALL( SCIPchgVarLbGlobal(scip, z, 0.0) );

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*) "<x> ^2 - <y> * <z>", NULL, &expr) );

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "soc", expr, -SCIPinfinity(scip), 0.0) );

   /* add locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* call detection method -> this registers the nlhdlr */
   SCIP_CALL( detectNlhdlrs(scip, conshdlr, &cons, 1, &infeasible) );
   cr_assert_not(infeasible);

   /* find the nlhdlr expr data */
   for( i = 0; i < expr->nenfos; ++i )
   {
      if( expr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = expr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* create solution */
   SCIPcreateSol(scip, &sol, NULL);
   SCIPsetSolVal(scip, sol, x, 1.0);
   SCIPsetSolVal(scip, sol, y, 1.0);
   SCIPsetSolVal(scip, sol, z, 0.0);

   /* check cut */
   SCIP_CALL( generateCutSolSOC(scip, expr, cons, sol, nlhdlrexprdata, 0.0, 1.0, &cut) );

   cutvars[0] = x;
   cutvars[1] = y;
   cutvars[2] = z;
   cutvals[0] = 4.0 / SQRT( 5.0 );
   cutvals[1] = -1.0 + 1.0 / SQRT( 5.0 );
   cutvals[2] = -(5.0 + SQRT( 5.0 )) / 5.0;
   rhs = 0.0;

   checkCut(cut, cutvars, cutvals, rhs, 3);
   SCIPreleaseRow(scip, &cut);

   /* free expr and cons */
   SCIPfreeSol(scip, &sol);

   /* remove locks */
   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* disable cons, so it can be deleted */
   SCIP_CALL( consDisableExpr(scip, conshdlr, cons) );

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
