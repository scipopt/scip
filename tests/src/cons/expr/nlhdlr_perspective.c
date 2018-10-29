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

/**@file   nlhdlr_quadratic.c
 * @brief  tests quadratic nonlinear handler methods
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_nlhdlr_perspective.c"
#include "scip/cons_expr.h"
#include "scip/struct_cons_expr.h"
#include "scip/cons_varbound.h"
#include "scip/scip.h"

#define SCIP_DEBUG


/*
 * TEST
 */

#include "include/scip_test.h"
#include "scip/type_cons.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;

static SCIP_CONSHDLR* conshdlr;
static SCIP_CONSHDLR* vbdconshdlr;
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
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* get nlhdlr */
   nlhdlr = SCIPfindConsExprNlhdlr(conshdlr, "perspective");
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
   SCIP_CALL( (*nlhdlr->exit)(scip, nlhdlr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   //BMScheckEmptyMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* test detection callback of nonlinear handler */
Test(nlhdlrperspective, detect, .init = setup, .fini = teardown)
{
   SCIP_CONS* cons1;
   SCIP_CONS* cons2;
   SCIP_CONS* vubcons;
   SCIP_CONS* vlbcons;
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;
   SCIP_Bool success;
   SCIP_Real vals[2];
   SCIP_VAR* vars[2];

   /* create expression constraints */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x>^2 + <x> + <z>", NULL, &expr1) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons1, (char*)"nlin1", expr1, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPaddCons(scip, cons1)  );

   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<y>^2 + exp(<y>) + <x>^3 + <x>*<w> + <y> + <x> + exp(<x>) + (<x> + <x>^2 + exp(1/<x>))^3 + 2", NULL, &expr2) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons2, (char*)"nlin2", expr2, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPaddCons(scip, cons2)  );

   /* create varbound constraint */
   vals[0] = 2;
   vals[1] = -4;
   vars[0] = x;
   vars[1] = z;
   SCIP_CALL( SCIPcreateConsVarbound(scip, &vubcons, "vub", x, z, -2.0, -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, vubcons)  );

   vals[0] = 2;
   vals[1] = -2;
   vars[0] = x;
   vars[1] = z;
   SCIP_CALL( SCIPcreateConsVarbound(scip, &vlbcons, "vlb", x, z, -1.0, 0.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, vlbcons)  );

   /* presolve */
   SCIP_CALL( SCIPpresolve(scip) );

   /* release */
   SCIP_CALL( SCIPreleaseCons(scip, &cons1) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons2) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr1) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseCons(scip, &vubcons) );
   SCIP_CALL( SCIPreleaseCons(scip, &vlbcons) );
}
