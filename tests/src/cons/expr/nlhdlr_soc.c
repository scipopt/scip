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

/**@file   nlhdlr_soc.c
 * @brief  tests quadratic nonlinear handler methods
 * @author Fabian Wegscheider
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

/* XXX: need the consdata struct because we don't have getNlhdlrs or findNlhdlrs; I don't add those function because I'm unsure
 * we actually need them
 */
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

static SCIP_CONSEXPR_NLHDLR* nlhdlr;
static SCIP_CONSHDLR* conshdlr;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   int h;
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

   /* go to SOLVING stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -0.9, 0.7, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &u, "u", 0.9, 10.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
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
   SCIP_VAR**            vars,
   SCIP_Real*            coefs,
   SCIP_Real*            offsets,
   SCIP_Real*            transcoefs,
   int*                  transcoefsidx,
   int*                  nnonzeroes,
   int                   nvars,
   int                   nterms,
   int                   ntranscoefs,
   SCIP_Real             constant
   )
{
   int i;

   cr_assert_not_null(nlhdlrexprdata->vars);
   cr_assert_not_null(nlhdlrexprdata->coefs);
   cr_assert_not_null(nlhdlrexprdata->offsets);
   cr_assert_not_null(nlhdlrexprdata->nnonzeroes);
   cr_assert_not_null(nlhdlrexprdata->transcoefs);
   cr_assert_not_null(nlhdlrexprdata->transcoefsidx);

   cr_expect_eq(nlhdlrexprdata->nvars, nvars);
   cr_expect_eq(nlhdlrexprdata->nterms, nterms);
   cr_expect_eq(nlhdlrexprdata->ntranscoefs, ntranscoefs);
   cr_expect(SCIPisEQ(scip, nlhdlrexprdata->constant, constant), "expected constant %f, but got %f\n",
      constant, nlhdlrexprdata->constant);

   for( i = 0; i < nvars; ++i )
   {
      cr_assert_not_null(nlhdlrexprdata->vars[i]);
      cr_expect_eq(nlhdlrexprdata->vars[i], vars[i], "expected variable %d to be %s, but got %s\n",
         i+1, SCIPvarGetName(vars[i]), SCIPvarGetName(nlhdlrexprdata->vars[i]));
   }

   for( i = 0; i < nterms; ++i )
   {
      cr_expect(SCIPisEQ(scip, nlhdlrexprdata->coefs[i], coefs[i]), "expected coef %d to be %f, but got %f\n",
         i+1, coefs[i], nlhdlrexprdata->coefs[i]);
   }

   for( i = 0; i < nterms; ++i )
   {
      cr_expect(SCIPisEQ(scip, nlhdlrexprdata->offsets[i], offsets[i]), "expected offset %d to be %f, but got %f\n",
         i+1, offsets[i], nlhdlrexprdata->offsets[i]);
   }

   for( i = 0; i < nterms; ++i )
   {
      cr_expect_eq(nlhdlrexprdata->nnonzeroes[i], nnonzeroes[i], "expected nnonzeroes %d to be %d, but got %d\n",
         i+1, nnonzeroes[i], nlhdlrexprdata->nnonzeroes[i]);
   }

   for( i = 0; i < ntranscoefs; ++i )
   {
      cr_expect(SCIPisEQ(scip, nlhdlrexprdata->transcoefs[i], transcoefs[i]), "expected transcoef %d to be %f, but got %f\n",
         i+1, transcoefs[i], nlhdlrexprdata->transcoefs[i]);
   }

   for( i = 0; i < ntranscoefs; ++i )
   {
      cr_expect_eq(nlhdlrexprdata->transcoefsidx[i], transcoefsidx[i], "expected transcoefsidx %d to be %d, but got %d\n",
         i+1, transcoefsidx[i], nlhdlrexprdata->transcoefsidx[i]);
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
   cr_expect_eq(SCIProwGetRhs(cut), rhs, "expected rhs = %f, but got %f\n", rhs, SCIProwGetRhs(cut));

   for( i = 0; i < nvars; ++i )
   {
      cr_expect_eq(SCIPcolGetVar(SCIProwGetCols(cut)[i]), vars[i], "expected var%d = %s, but got %s\n",
         i+1, SCIPvarGetName(vars[i]), SCIPvarGetName(SCIPcolGetVar(SCIProwGetCols(cut)[i])));
      cr_expect_eq(SCIProwGetVals(cut)[i], vals[i], "expected val%d = %f, but got %f\n", i+1, vals[i],
         SCIProwGetVals(cut)[i]);
   }
}

/* detects ||x|| < 1 as soc expression */
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

   /* setup expected data */
   SCIP_VAR* vars[4] = {x, y, z, SCIPgetConsExprExprAuxVar(expr)};
   SCIP_Real coefs[4] = {1.0, 1.0, 1.0, 1.0};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[4] = {1.0, 1.0, 1.0, 1.0};
   int transcoefsidx[4] = {0, 1, 2, 3};
   int nnonzeroes[4] = {1, 1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 4, 0.0);

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

   /* this also creates the locks */
   SCIP_CALL( SCIPaddCons(scip, cons) );

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
   SCIP_VAR* vars[4] = {x, y, z, SCIPgetConsExprExprAuxVar(normexpr)};
   SCIP_Real coefs[4] = {1.0, 1.0, 1.0, 1.0};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[4] = {1.0, 1.0, 1.0, 1.0};
   int transcoefsidx[4] = {0, 1, 2, 3};
   int nnonzeroes[4] = {1, 1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 4, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects sqrt( 2*(x + 1)^2 + 3*(y + sin(x) + 2)^2 ) as soc expression */
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
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: (8 + 2*(<x> + 1)^2 + 3*(sin(<y>) - 2)^2)^0.5 + 2*(<w> - 1) <= 0",
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
   normexpr = SCIPgetConsExprExprChildren(expr)[1];

   /* find the nlhdlr expr data */
   for( i = 0; i < normexpr->nenfos; ++i )
   {
      if( normexpr->enfos[i]->nlhdlr == nlhdlr )
         nlhdlrexprdata = normexpr->enfos[i]->nlhdlrexprdata;
   }
   cr_assert_not_null(nlhdlrexprdata);

   /* setup expected data */
   SCIP_VAR* sinauxvar = SCIPgetConsExprExprAuxVar(normexpr->children[0]->children[2]);
   SCIP_VAR* vars[3] = {x, sinauxvar, SCIPgetConsExprExprAuxVar(normexpr)};
   SCIP_Real coefs[3] = {2.0, 3.0, 1.0};
   SCIP_Real offsets[3] = {1.0, -2.0, 0.0};
   SCIP_Real transcoefs[3] = {1.0, 1.0, 1.0};
   int transcoefsidx[3] = {0, 1, 2};
   int nnonzeroes[3] = {1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 3, 8.0);

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

   /* setup expected data */
   SCIP_VAR* sinauxvar = SCIPgetConsExprExprAuxVar(expr->children[2]->children[0]);
   SCIP_VAR* vars[3] = {x, sinauxvar, y};
   SCIP_Real coefs[3] = {2.0, 1.0, SQRT(9.0)};
   SCIP_Real offsets[3] = {0.0, 0.0, 0.0};
   SCIP_Real transcoefs[3] = {1.0, 1.0, 1.0};
   int transcoefsidx[3] = {0, 1, 2};
   int nnonzeroes[3] = {1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 3, 0.0);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects 5 - 7cos(x)^2 + y^2 + 2sin(z)^2 <= 4 as soc expression */
Test(nlhdlrsoc, detectandfree5, .description = "detects more complication quadratic expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Bool infeasible;
   SCIP_Bool success;
   int i;

   /* create expression constraint */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: 5 - 7*cos(<x>)^2 + <y>^2 + 2*sin(<z>)^2 <= 4",
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

   /* setup expected data */
   SCIP_VAR* cosauxvar = SCIPgetConsExprExprAuxVar(expr->children[1]->children[0]);
   SCIP_VAR* sinauxvar = SCIPgetConsExprExprAuxVar(expr->children[2]->children[0]);
   SCIP_VAR* vars[3] = {y, sinauxvar, cosauxvar};
   SCIP_Real coefs[3] = {1.0, 2.0, SQRT(7.0)};
   SCIP_Real offsets[3] = {0.0, 0.0, 0.0};
   SCIP_Real transcoefs[3] = {1.0, 1.0, 1.0};
   int transcoefsidx[3] = {0, 1, 2};
   int nnonzeroes[3] = {1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 3, 1);

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
   cr_assert_null(nlhdlrexprdata);

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
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: 1 + 7*<x>^2 - 2*<y>*<u> + 2*<z>^2 <= 0",
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

   /* setup expected data */
   SCIP_VAR* vars[4] = {x, z, y, u};
   SCIP_Real coefs[4] = {14.0, 4.0, 1.0, 1.0};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[6] = {1.0, 1.0, 1.0, -1.0, 1.0, 1.0};
   int transcoefsidx[6] = {0, 1, 2, 3, 2, 3};
   int nnonzeroes[4] = {1, 1, 2, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 6, 2);

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
   cr_assert_null(nlhdlrexprdata);

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
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: -5 + 7*cos(<x>)^2 - <y>^2 - 2*sin(<z>)^2 >= -4",
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

   /* setup expected data */
   SCIP_VAR* cosauxvar = SCIPgetConsExprExprAuxVar(expr->children[1]->children[0]);
   SCIP_VAR* sinauxvar = SCIPgetConsExprExprAuxVar(expr->children[2]->children[0]);
   SCIP_VAR* vars[3] = {y, sinauxvar, cosauxvar};
   SCIP_Real coefs[3] = {1.0, 2.0, SQRT(7.0)};
   SCIP_Real offsets[3] = {0.0, 0.0, 0.0};
   SCIP_Real transcoefs[3] = {1.0, 1.0, 1.0};
   int transcoefsidx[3] = {0, 1, 2};
   int nnonzeroes[3] = {1, 1, 1};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 3, 1);

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
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*) "[expr] <test>: -1 - 7*<x>^2 + 2*<y>*<u> - 2*<z>^2 >= 0",
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

   /* setup expected data */
   SCIP_VAR* vars[4] = {x, z, y, u};
   SCIP_Real coefs[4] = {14.0, 4.0, 1.0, 1.0};
   SCIP_Real offsets[4] = {0.0, 0.0, 0.0, 0.0};
   SCIP_Real transcoefs[6] = {1.0, 1.0, 1.0, -1.0, 1.0, 1.0};
   int transcoefsidx[6] = {0, 1, 2, 3, 2, 3};
   int nnonzeroes[4] = {1, 1, 2, 2};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 4, 4, 6, 2);

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

   /* setup expected data */
   SCIP_VAR* vars[3] = {x, y, z};
   SCIP_Real coefs[3] = {1.221024161223004, 2.45699535030454, -1.2360121980419372};
   SCIP_Real offsets[3] = {0.43579753197260707, 0.4739266455634833, -2.650186894292017};
   SCIP_Real transcoefs[9] = { 0.1299461082832169,  0.7137532837947576, -0.688237065853218,
                               0.5722025535583013,  0.5128905029681603,  0.6399434113001139,
                              -0.3417633267618488, -0.4769691623546202,  0.8097519661250633};
   int transcoefsidx[9] = {0, 1, 2, 0, 1, 2, 2, 1, 0};
   int nnonzeroes[3] = {3, 3, 3};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 9, 17.0909090909091);

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

   /* setup expected data */
   SCIP_VAR* vars[3] = {x, y, z};
   SCIP_Real coefs[3] = {2.45699535030454, 1.221024161223004, -1.2360121980419372};
   SCIP_Real offsets[3] = {0.4739266455634833, 0.43579753197260707, -2.650186894292017};
   SCIP_Real transcoefs[9] = { 0.5722025535583013,  0.5128905029681603,  0.6399434113001139,
                               0.1299461082832169,  0.7137532837947576, -0.688237065853218,
                              -0.3417633267618488, -0.4769691623546202,  0.8097519661250633};
   int transcoefsidx[9] = {0, 1, 2, 0, 1, 2, 2, 1, 0};
   int nnonzeroes[3] = {3, 3, 3};

   /* check nlhdlrexprdata*/
   checkData(nlhdlrexprdata, vars, coefs, offsets, transcoefs, transcoefsidx, nnonzeroes, 3, 3, 9, 17.0909090909091);

   /* free cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* disaggregates sqrt( 2*(x + 1)^2 + 3*(y + sin(x) + 2)^2 ) */
Test(nlhdlrsoc, disaggregation, .description = "detects more complex norm expression")
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

   SCIP_CALL( SCIPaddCons(scip, cons) );

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeasible, NULL, NULL, NULL) );
   cr_expect(!infeasible);

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

   /* check disvars */
   cr_expect(nlhdlrexprdata->disvars[0] != NULL);
   cr_expect(nlhdlrexprdata->disvars[1] != NULL);
   cr_expect(nlhdlrexprdata->disvars[2] != NULL);

   cr_expect(SCIProwGetNNonz(nlhdlrexprdata->disrow) == 4);

   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[0], 1.0, "expected %f, but got %f\n", 1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[0]);
   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[1], 1.0, "expected %f, but got %f\n", 1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[1]);
   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[2], 1.0, "expected %f, but got %f\n", 1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[2]);
   cr_expect_eq(SCIProwGetVals(nlhdlrexprdata->disrow)[3], -1.0, "expected %f, but got %f\n", -1.0, SCIProwGetVals(nlhdlrexprdata->disrow)[3]);

   cr_expect(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[0]) == nlhdlrexprdata->disvars[0]);
   cr_expect(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[1]) == nlhdlrexprdata->disvars[1]);
   cr_expect(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[2]) == nlhdlrexprdata->disvars[2]);
   cr_expect(SCIPcolGetVar(SCIProwGetCols(nlhdlrexprdata->disrow)[3]) == nlhdlrexprdata->vars[2]);

   cr_expect(SCIProwGetLhs(nlhdlrexprdata->disrow) == -SCIPinfinity(scip));
   cr_expect(SCIProwGetRhs(nlhdlrexprdata->disrow) == 0.0);

   /* free expr and cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* separates simple norm function from different points */
Test(nlhdlrsoc, separation1, .description = "test separation for simple norm expression")
{
   SCIP_CONS* cons;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_SOL* sol;
   SCIP_ROW* cut;
   SCIP_VAR* cutvars[3];
   SCIP_VAR* auxvar;
   SCIP_Real cutvals[3];
   SCIP_Bool infeasible;
   int i;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*) "(<x>^2 + <y>^2 + <z>^2)^0.5", NULL, &expr) );

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, "soc", expr, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );

   /* detect */
   SCIP_CALL( SCIPinitlpCons(scip, cons, &infeasible) );
   cr_assert(!infeasible);

   SCIP_CALL( SCIPclearCuts(scip) ); /* we have to clear the separation store */
   SCIP_CALL( SCIPaddCons(scip, cons) );

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
   SCIP_CALL( generateCutSol(scip, expr, conshdlr, sol, nlhdlrexprdata, 0, 0.0, &cut) );

   cr_expect_null(cut);

   /* check cut w.r.t. y */
   SCIP_CALL( generateCutSol(scip, expr, conshdlr, sol, nlhdlrexprdata, 1, 0.0, &cut) );

   cutvars[0] = nlhdlrexprdata->disvars[1];
   cutvars[1] = auxvar;
   cutvars[2] = y;
   cutvals[0] = -1.0 / SQRT(2.0);
   cutvals[1] = -1.0 / (2.0 * SQRT(2.0));
   cutvals[2] = 1.0;

   checkCut(cut, cutvars, cutvals, -4.0 + 2.0 / SQRT(2.0), 3);
   SCIPreleaseRow(scip, &cut);

   /* check cut w.r.t. z */
   SCIP_CALL( generateCutSol(scip, expr, conshdlr, sol, nlhdlrexprdata, 2, 0.0, &cut) );

   cutvars[0] = nlhdlrexprdata->disvars[2];
   cutvars[1] = auxvar;
   cutvars[2] = z;
   cutvals[0] = -1.0 / SQRT(2.0);
   cutvals[1] = -1.0 / (2.0 * SQRT(2.0));
   cutvals[2] = -1.0;

   checkCut(cut, cutvars, cutvals, -4.0 + 2.0 / SQRT(2.0), 3);
   SCIPreleaseRow(scip, &cut);

   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );

   /* free expr and cons */
   SCIPfreeSol(scip, &sol);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}