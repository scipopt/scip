/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_perspective.c
 * @brief  tests perspective nonlinear handler methods
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#define SCIP_PRIVATE_ROWPREP
#include "scip/cons_expr.c"
#include "scip/cons_expr_nlhdlr_perspective.c"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x_1;
static SCIP_VAR* y_1;
static SCIP_VAR* z_1;

static SCIP_CONSHDLR* conshdlr;
static SCIP_CONSEXPR_NLHDLR* nlhdlr = NULL;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   int h;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get perspective handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );

   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* get nlhdlr */
   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
   {
      SCIPinfoMessage(scip, NULL, "\nnlhdlr = %s", SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]));
      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), "perspective") == 0 )
      {
         nlhdlr = conshdlrdata->nlhdlrs[h];
         break;
      }
   }

   cr_assert_not_null(nlhdlr);


   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to PRESOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, TRUE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x_1, "x1", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y_1, "y1", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_1, "z1", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x_1) );
   SCIP_CALL( SCIPaddVar(scip, y_1) );
   SCIP_CALL( SCIPaddVar(scip, z_1) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &y_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_1) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* detects x1^2 + x1 - log(y1) as an on/off expression */
Test(nlhdlrperspective, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONSEXPR_EXPRENFO_METHOD providedexpected;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_Bool changed;
   SCIP_CONS* cons;
   SCIP_CONS* vubcons;
   SCIP_CONS* vlbcons;
   SCIP_SCVARDATA* scv;

   /* create expression and simplify it */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x1>^2 + <x1> - log(<y1>)", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified, &changed) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );

   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );

   /* create upper varbound constraint x1 - 2z <= 0 */
   SCIP_CALL( SCIPcreateConsVarbound(scip, &vubcons, "vub", x_1, z_1, -2.0, -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, vubcons)  );

   /* create lower varbound constraint x1 - z >= 0 */
   SCIP_CALL( SCIPcreateConsVarbound(scip, &vlbcons, "vlb", x_1, z_1, -1.0, 0.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, vlbcons)  );

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);

   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_NONE;
   enforcebelow = FALSE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   providedexpected = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   cr_expect_eq(provided, providedexpected, "expecting provided = %d, got %d\n", providedexpected, provided);
   cr_assert(enforcebelow);
   cr_assert(!enforceabove);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nconvterms, 2, "Expecting 2 convex terms, got %d\n", nlhdlrexprdata->nconvterms);
   cr_expect_eq(nlhdlrexprdata->nonoffterms, 1, "Expecting 1 perspective term, got %d\n", nlhdlrexprdata->nonoffterms);

   scv = SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)x_1);
   cr_expect_not_null(scv, "x should be detected as a semicontinuous variable!\n");
   cr_expect_eq(scv->val0, 0.0, "Expecting val0 to be %g, got %g\n", 0.0, scv->val0);
   cr_expect_eq(scv->bvar, z_1, "Expecting bvar to be %p, got %p\n", z_1, scv->bvar);

   SCIP_CALL( freeAuxVars(scip, conshdlr, &cons, 1) );

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &vubcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &vlbcons) );
}

/* TODO more unittests */










