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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.01, 1.01, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 1.49, 1.51, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0, 1, 0, SCIP_VARTYPE_BINARY) );
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
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* detects x^2 + x as an on/off expression */
Test(nlhdlrperspective, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* simplified;
   SCIP_CONSEXPR_EXPRENFO_METHOD providedexpected;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_CONS* cons;
   SCIP_CONS* vubcons;
   SCIP_CONS* vlbcons;
   SCIP_SCVARDATA* scv;

   /* create expression and simplify it */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>^2 + <x>", NULL, &expr) );
   SCIP_CALL( SCIPsimplifyConsExprExpr(scip, conshdlr, expr, &simplified) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   expr = simplified;

   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );

   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );

   /* create varbound constraint */
   SCIP_CALL( SCIPcreateConsVarbound(scip, &vubcons, "vub", x, z, -2.0, -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, vubcons)  );

   SCIP_CALL( SCIPcreateConsVarbound(scip, &vlbcons, "vlb", x, z, -1.0, 0.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, vlbcons)  );


   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_NONE;
   enforcebelow = FALSE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   providedexpected = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   cr_expect_eq(provided, providedexpected, "expecting %d got %d\n", providedexpected, provided);
   cr_assert(enforcebelow);
   cr_assert(!enforceabove);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nlinterms, 0, "Expecting 0 linear terms, got %d\n", nlhdlrexprdata->nlinterms);
   cr_expect_eq(nlhdlrexprdata->nperspterms, 1, "Expecting 1 perspective term, got %d\n", nlhdlrexprdata->nperspterms);

   scv = SCIPhashmapGetImage(nlhdlrexprdata->scvars, (void*)x);
   cr_expect_not_null(scv, "x should be detected as a semicontinuous variable!\n");
   cr_expect_eq(scv->lb0, 0.0, "Expecting lb0 to be %g, got %g\n", 0.0, scv->lb0);
   cr_expect_eq(scv->ub0, 0.0, "Expecting ub0 to be %g, got %g\n", 0.0, scv->ub0);
   cr_expect_eq(scv->lb1, 1.0, "Expecting lb1 to be %g, got %g\n", 1.0, scv->lb1);
   cr_expect_eq(scv->ub1, 1.01, "Expecting ub1 to be %g, got %g\n", 1.01, scv->ub1);

   SCIP_CALL( freeAuxVars(scip, conshdlr, &cons, 1) );

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
   SCIP_CALL( SCIPreleaseCons(scip, &vubcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &vlbcons) );
}












