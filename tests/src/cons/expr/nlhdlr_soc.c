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
#define SCIP_PRIVATE_ROWPREP
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

static SCIP_CONSHDLR* conshdlr;
static SCIP_CONSEXPR_NLHDLR* nlhdlr = NULL;

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

   /* get nlhdlr */
   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), "quadratic") == 0 )
      {
         nlhdlr = conshdlrdata->nlhdlrs[h];
         break;
      }
   cr_assert_not_null(nlhdlr);


   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 1.49, 1.51, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -0.9, 0.7, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   //BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test suite */
TestSuite(nlhdlrsoc, .init = setup, .fini = teardown);

/* detects ||x|| as quadratic expression */
Test(nlhdlrsoc, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_VAR* auxvar
   SCIP_Bool success;
   SCIP_Bool changed = FALSE;

   /* create expression and simplify it: note it fails if not simplified, the order matters! */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"(<x>^2 + <y>^2 + <z>^2)^0.5", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed) );
   cr_expect(!changed);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;

   /* create auxvar */
   SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, &auxvar) );
   cr_assert(auxvar != NULL);

   /* detect */
   success = FALSE;
   SCIP_CALL( detectSOC(scip, expr, auxvar, &nlhdlrexprdata, &success) );
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   /* check nlhdlrexprdata*/
   cr_expect_eq(nlhdlrexprdata->nvars, 4);
   cr_expect_eq(nlhdlrexprdata->ntranscoefs, 4);
   cr_expect_eq(nlhdlrexprdata->nterms, 4);
   cr_expect_eq(nlhdlrexprdata->constant, 0.0);
   cr_expect_eq(nlhdlrexprdata->nnonzeroes[0], 1);
   cr_expect_eq(nlhdlrexprdata->nnonzeroes[1], 1);
   cr_expect_eq(nlhdlrexprdata->nnonzeroes[2], 1);
   cr_expect_eq(nlhdlrexprdata->nnonzeroes[3], 1);
   cr_expect_eq(nlhdlrexprdata->transcoefsidx[0], 0);
   cr_expect_eq(nlhdlrexprdata->transcoefsidx[1], 1);
   cr_expect_eq(nlhdlrexprdata->transcoefsidx[2], 2);
   cr_expect_eq(nlhdlrexprdata->transcoefsidx[3], 3);
   cr_expect_eq(nlhdlrexprdata->transcoefs[0], 1.0);
   cr_expect_eq(nlhdlrexprdata->transcoefs[1], 1.0);
   cr_expect_eq(nlhdlrexprdata->transcoefs[1], 1.0);
   cr_expect_eq(nlhdlrexprdata->transcoefs[2], 1.0);
   cr_expect_eq(nlhdlrexprdata->transcoefs[3], 1.0);
   cr_expect_eq(nlhdlrexprdata->offsets[0], 0.0);
   cr_expect_eq(nlhdlrexprdata->offsets[1], 0.0);
   cr_expect_eq(nlhdlrexprdata->offsets[2], 0.0);
   cr_expect_eq(nlhdlrexprdata->offsets[3], 0.0);
   cr_expect_eq(nlhdlrexprdata->coefs[0], 1.0);
   cr_expect_eq(nlhdlrexprdata->coefs[1], 1.0);
   cr_expect_eq(nlhdlrexprdata->coefs[2], 1.0);
   cr_expect_eq(nlhdlrexprdata->coefs[3], 1.0);
   cr_expect_eq(nlhdlrexprdata->vars[0], x);
   cr_expect_eq(nlhdlrexprdata->vars[1], y);
   cr_expect_eq(nlhdlrexprdata->vars[2], z);
   cr_expect_eq(nlhdlrexprdata->vars[3], auxvar);

   /* register enforcer info in expr and free */
   SCIP_CALL( freeNlhdlrExprData(scip, &nlhdlrexprdata) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}