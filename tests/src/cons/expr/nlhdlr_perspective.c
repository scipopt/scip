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
static SCIP_VAR* z_2;
static SCIP_VAR* z_3;

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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x_1, "x1", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y_1, "y1", 0.07, 0.09, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_1, "z1", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_2, "z2", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_3, "z3", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x_1) );
   SCIP_CALL( SCIPaddVar(scip, y_1) );
   SCIP_CALL( SCIPaddVar(scip, z_1) );
   SCIP_CALL( SCIPaddVar(scip, z_2) );
   SCIP_CALL( SCIPaddVar(scip, z_3) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &y_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_3) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* tests the detection of semicontinuous variables */
Test(nlhdlrperspective, varissc, .init = setup, .fini = teardown)
{
   SCIP_HASHMAP* scvars;
   SCIP_HASHMAPENTRY* entry;
   SCIP_SCVARDATA* scvdata;
   SCIP_Bool result, infeas;
   int nbndchgs, c;

   SCIPinfoMessage(scip, NULL, "\nTesting varIsSemicontinuous");

   /* allocate memory */
   SCIPhashmapCreate(&scvars, SCIPblkmem(scip), 1);

   /* add bound information to the vars */
   SCIPinfoMessage(scip, NULL, "\nvar z1 has bounds %f, %f", SCIPvarGetLbGlobal(z_1), SCIPvarGetUbGlobal(z_1));
   /* z1 <= x1 <= 3*z1 */
   SCIPaddVarVlb(scip, x_1, z_1, 1.0, 0.0, &infeas, &nbndchgs);
   SCIPaddVarVub(scip, x_1, z_1, 3.0, 0.0, &infeas, &nbndchgs);
   SCIPinfoMessage(scip, NULL, "\nvar x1 has bounds %f, %f", SCIPvarGetLbGlobal(x_1), SCIPvarGetUbGlobal(x_1));
   /* x1 <= 2*z2 */
   SCIPaddVarVub(scip, x_1, z_2, 2.0, 0.0, &infeas, &nbndchgs);
   SCIPinfoMessage(scip, NULL, "\nvar x1 has bounds %f, %f", SCIPvarGetLbGlobal(x_1), SCIPvarGetUbGlobal(x_1));
   /* 2z3 <= x1 <= 4*z3 */
   SCIPaddVarVlb(scip, x_1, z_3, 2.0, 0.0, &infeas, &nbndchgs);
   SCIPaddVarVub(scip, x_1, z_3, 4.0, 0.0, &infeas, &nbndchgs);
   SCIPinfoMessage(scip, NULL, "\nvar x1 has bounds %f, %f", SCIPvarGetLbGlobal(x_1), SCIPvarGetUbGlobal(x_1));

   /* check if the var is semicontinuous */
   SCIP_CALL( varIsSemicontinuous(scip, x_1, scvars, &result) );

   /* check result */
   cr_expect_eq(SCIPhashmapGetNElements(scvars), 1, "Expected 1 semicontinuous variable, got %d", SCIPhashmapGetNElements(scvars));
   scvdata = SCIPhashmapGetImage(scvars, (void*)x_1);
   cr_expect_eq(scvdata->nbnds, 3, "Expected 3 on/off bounds for variable x1, got %d", scvdata->nbnds);

   cr_expect_eq(scvdata->bvars[0], (void*)z_1, "bvars[0] expected to be z1, got %s", SCIPvarGetName(scvdata->bvars[0]));
   cr_expect_eq(scvdata->vals0[0], 0.0, "vals0[0] expected to be 0.0, got %f", scvdata->vals0[0]);
   cr_expect_eq(scvdata->bvars[1], (void*)z_2, "bvars[1] expected to be z2, got %s", SCIPvarGetName(scvdata->bvars[1]));
   cr_expect_eq(scvdata->vals0[1], 0.0, "vals0[1] expected to be 0.0, got %f", scvdata->vals0[1]);
   cr_expect_eq(scvdata->bvars[2], (void*)z_3, "bvars[2] expected to be z3, got %s", SCIPvarGetName(scvdata->bvars[2]));
   cr_expect_eq(scvdata->vals0[2], 0.0, "vals0[2] expected to be 0.0, got %f", scvdata->vals0[2]);

   /* free memory */
   for( c = 0; c < SCIPhashmapGetNEntries(scvars); ++c )
   {
      entry = SCIPhashmapGetEntry(scvars, c);
      if( entry != NULL )
      {
         scvdata = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);
         SCIPfreeBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize);
         SCIPfreeBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize);
         SCIPfreeBlockMemory(scip, &scvdata);
      }
   }
   SCIPhashmapFree(&scvars);
}

/* detects x1^2 + x1 - log(y1) as an on/off expression */
#if 0
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
#endif

/* TODO more unittests */










