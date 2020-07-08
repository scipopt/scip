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
#include "scip/cons_expr.c"
#include "scip/cons_expr_nlhdlr_perspective.c"
#include "scip/cons_expr_nlhdlr_convex.h"


/*
 * TEST
 */

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x_1;
static SCIP_VAR* x_2;
static SCIP_VAR* x_3;
static SCIP_VAR* y_1;
static SCIP_VAR* y_2;
static SCIP_VAR* y_3;
static SCIP_VAR* z_1;
static SCIP_VAR* z_2;
static SCIP_VAR* z_3;

static SCIP_CONSHDLR* conshdlr;
static SCIP_CONSEXPR_NLHDLR* nlhdlr = NULL;
static SCIP_CONSEXPR_NLHDLR* nlhdlr_conv = NULL;

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   int h;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers and nonlinear handlers; get perspective handler and conshdlr */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* get perspective and convex nlhdlrs */
   for( h = 0; h < conshdlrdata->nnlhdlrs; ++h )
   {
      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), "perspective") == 0 )
         nlhdlr = conshdlrdata->nlhdlrs[h];

      if( strcmp(SCIPgetConsExprNlhdlrName(conshdlrdata->nlhdlrs[h]), "convex") == 0 )
         nlhdlr_conv = conshdlrdata->nlhdlrs[h];

      if( nlhdlr != NULL && nlhdlr_conv != NULL )
         break;
   }

   cr_assert_not_null(nlhdlr);
   cr_assert_not_null(nlhdlr_conv);


   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* go to SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x_1, "x1", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x_2, "x2", -1.0, 5.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x_3, "x3", -1.0, 5.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y_1, "y1", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y_2, "y2", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y_3, "y3", 0.0, 4.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_1, "z1", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_2, "z2", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z_3, "z3", 0, 1, 0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, x_1) );
   SCIP_CALL( SCIPaddVar(scip, x_2) );
   SCIP_CALL( SCIPaddVar(scip, x_3) );
   SCIP_CALL( SCIPaddVar(scip, y_1) );
   SCIP_CALL( SCIPaddVar(scip, y_2) );
   SCIP_CALL( SCIPaddVar(scip, y_3) );
   SCIP_CALL( SCIPaddVar(scip, z_1) );
   SCIP_CALL( SCIPaddVar(scip, z_2) );
   SCIP_CALL( SCIPaddVar(scip, z_3) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x_2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x_3) );
   SCIP_CALL( SCIPreleaseVar(scip, &y_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &y_2) );
   SCIP_CALL( SCIPreleaseVar(scip, &y_3) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z_3) );
   SCIP_CALL( SCIPfree(&scip) );

   BMSdisplayMemory();
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

static
void checkCut(SCIP_ROW* cut, SCIP_VAR** vars, SCIP_Real* vals, int nvars, SCIP_Real lhs, SCIP_Real rhs)
{
   SCIP_VAR* var;
   SCIP_Real coef;
   SCIP_Bool found;
   int i;
   int j;

   cr_assert(cut != NULL);
   cr_expect_eq(SCIProwGetNNonz(cut), nvars, "\nExpected %d vars, got %d", nvars, SCIProwGetNNonz(cut));
   cr_expect(SCIPisEQ(scip, SCIProwGetLhs(cut), lhs));
   cr_expect(SCIPisEQ(scip, SCIProwGetRhs(cut), rhs));

   for( i = 0; i < SCIProwGetNNonz(cut); ++i )
   {
      var = SCIPcolGetVar(SCIProwGetCols(cut)[i]);
      coef = SCIProwGetVals(cut)[i];
      found = FALSE;

      for( j = 0; j < nvars; ++j )
      {
         if( var == vars[j] )
         {
            cr_expect(SCIPisEQ(scip, coef, vals[j]));
            found = TRUE;
         }
      }

      cr_expect(found, "variable %s must be in the cut", SCIPvarGetName(vars[j]));
   }
}

/* tests the detection of semicontinuous variables */
Test(nlhdlrperspective, varissc, .init = setup, .fini = teardown)
{
   SCIP_HASHMAP* scvars;
   SCIP_HASHMAPENTRY* entry;
   SCVARDATA* scvdata;
   SCIP_Bool result;
   SCIP_Bool infeas;
   int nbndchgs;
   int c;

   /* allocate memory */
   SCIP_CALL( SCIPhashmapCreate(&scvars, SCIPblkmem(scip), 1) );

   /* add bound information to the vars */
   /* z1 <= x1 <= 3*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_1, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_1, 3.0, 0.0, &infeas, &nbndchgs) );
   /* x1 <= 2*z2 */
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_2, 2.0, 0.0, &infeas, &nbndchgs) );
   /* -z3 + 1 <= x1 <= 4*z3 + 1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_3, -1.0, 1.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_3, 4.0, 1.0, &infeas, &nbndchgs) );

   /* check if the var is semicontinuous */
   SCIP_CALL( varIsSemicontinuous(scip, x_1, scvars, &result) );

   /* check result */
   cr_expect_eq(SCIPhashmapGetNElements(scvars), 1, "Expected 1 semicontinuous variable, got %d", SCIPhashmapGetNElements(scvars));
   scvdata = (SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)x_1);
   cr_expect_eq(scvdata->nbnds, 3, "Expected 3 on/off bounds for variable x1, got %d", scvdata->nbnds);

   cr_expect_eq(scvdata->bvars[0], (void*)z_1, "bvars[0] expected to be z1, got %s", SCIPvarGetName(scvdata->bvars[0]));
   cr_expect_eq(scvdata->vals0[0], 0.0, "vals0[0] expected to be 0.0, got %f", scvdata->vals0[0]);
   cr_expect_eq(scvdata->bvars[1], (void*)z_2, "bvars[1] expected to be z2, got %s", SCIPvarGetName(scvdata->bvars[1]));
   cr_expect_eq(scvdata->vals0[1], 0.0, "vals0[1] expected to be 0.0, got %f", scvdata->vals0[1]);
   cr_expect_eq(scvdata->bvars[2], (void*)z_3, "bvars[2] expected to be z3, got %s", SCIPvarGetName(scvdata->bvars[2]));
   cr_expect_eq(scvdata->vals0[2], 1.0, "vals0[2] expected to be 1.0, got %f", scvdata->vals0[2]);

   /* free memory */
   for( c = 0; c < SCIPhashmapGetNEntries(scvars); ++c )
   {
      entry = SCIPhashmapGetEntry(scvars, c);
      if( entry != NULL )
      {
         scvdata = (SCVARDATA*) SCIPhashmapEntryGetImage(entry);
         SCIPfreeBlockMemoryArray(scip, &scvdata->ubs1, scvdata->bndssize);
         SCIPfreeBlockMemoryArray(scip, &scvdata->lbs1, scvdata->bndssize);
         SCIPfreeBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize);
         SCIPfreeBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize);
         SCIPfreeBlockMemory(scip, &scvdata);
      }
   }
   SCIPhashmapFree(&scvars);
}

/* detects x1^2 + x1y1 + x2^2 as an on/off expression with 2 indicator variables */
Test(nlhdlrperspective, detectandfree1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPRENFO_METHOD providedexpected;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x1>^2 + <x1>*<y1> + <x2>^2", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );

   /* add implied variable bounds */
   /* -3z1 + 3 <= x1 <= 3*z1 + 3 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_1, -3.0, 3.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_1, 3.0, 3.0, &infeas, &nbndchgs) );
   /* z1 <= y1 <= 3*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, y_1, z_1, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, y_1, z_1, 3.0, 0.0, &infeas, &nbndchgs) );
   /* -z1 <= x2 <= 5*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_1, -1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_1, 5.0, 0.0, &infeas, &nbndchgs) );

   /* add bounds with z_2 */
   /* z2 <= x1 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_2, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_2, 3.0, 0.0, &infeas, &nbndchgs) );
   /* -z2 <= y1 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, y_1, z_2, -1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, y_1, z_2, 3.0, 0.0, &infeas, &nbndchgs) );
   /* -z2 <= x2 <= 5*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_2, -1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_2, 5.0, 0.0, &infeas, &nbndchgs) );

   /* for a sum, detect wants another (non default) nlhdlr to have detected; pretend that convex has detected */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr->enfos, 2) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &expr->enfos[0]) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &expr->enfos[1]) );
   expr->nenfos = 2;
   expr->enfos[0]->nlhdlr = nlhdlr_conv;
   expr->enfos[1]->nlhdlr = nlhdlr;

   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   enforcebelow = TRUE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   providedexpected = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   cr_expect_eq(provided, providedexpected, "expecting provided = %d, got %d\n", providedexpected, provided);
   cr_assert(enforcebelow);
   cr_assert(!enforceabove);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator variables, got %d\n", nlhdlrexprdata->nindicators);

   /* check the 'off' values */
   cr_assert_eq(nlhdlrexprdata->indicators[0], z_1, "Expecting the first indicator to be z_1, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[0]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[0], 9.0, "Expecting off value = 9.0, got %f\n", nlhdlrexprdata->exprvals0[0]);

   cr_assert_eq(nlhdlrexprdata->indicators[1], z_2, "Expecting the second indicator to be z_2, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[1]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[1], 0.0, "Expecting off value = 0.0, got %f\n", nlhdlrexprdata->exprvals0[1]);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* detects x1^2 as an on/off expression */
Test(nlhdlrperspective, detectandfree2, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPRENFO_METHOD providedexpected;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x1>^2", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );
   SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, NULL) );

   /* z1 <= x1 <= 3*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_1, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_1, 3.0, 0.0, &infeas, &nbndchgs) );
   /* z2 <= x1 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_2, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_2, 3.0, 0.0, &infeas, &nbndchgs) );

   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   enforcebelow = TRUE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   providedexpected = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   cr_expect_eq(provided, providedexpected, "expecting provided = %d, got %d\n", providedexpected, provided);
   cr_assert(enforcebelow);
   cr_assert(!enforceabove);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator variables, got %d\n", nlhdlrexprdata->nindicators);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* detects log(x1+x2) as an on/off expression */
Test(nlhdlrperspective, detectandfree3, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPRENFO_METHOD providedexpected;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"log(<x1>+<x2>+1.0)", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );

   /* z1 <= x1 <= 3*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_1, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_1, 3.0, 0.0, &infeas, &nbndchgs) );
   /* z2 <= x1 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_2, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_2, 3.0, 0.0, &infeas, &nbndchgs) );
   /* z1 <= x2 <= 3*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_1, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_1, 3.0, 0.0, &infeas, &nbndchgs) );
   /* z2 <= x2 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_2, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_2, 3.0, 0.0, &infeas, &nbndchgs) );
   /* z3 <= x2 <= 3*z3 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_3, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_3, 3.0, 0.0, &infeas, &nbndchgs) );

   /* detect */
   provided = SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   enforcebelow = FALSE;
   enforceabove = TRUE;
   success = FALSE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove, &success, &nlhdlrexprdata) );

   providedexpected = SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
   cr_expect_eq(provided, providedexpected, "expecting provided = %d, got %d\n", providedexpected, provided);
   cr_assert(enforceabove);
   cr_assert(!enforcebelow);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator vars, got %d\n", nlhdlrexprdata->nindicators);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/* separates x1^2 + x1y1 + x2^2 for 2 indicator variables */
Test(nlhdlrperspective, sepa1, .init = setup, .fini = teardown)
{
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata_conv = NULL;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPRENFO_METHOD providedexpected;
   SCIP_CONSEXPR_EXPRENFO_METHOD provided;
   SCIP_Bool enforcebelow;
   SCIP_Bool enforceabove;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;
   SCIP_SOL* sol;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_VAR* auxvar;
   SCIP_PTRARRAY* rowpreps;
   SCIP_RESULT result;

   /* skip when no ipopt */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, (char*)"<x1>^2 + <x1>*<x2> + <x2>^2", NULL, &expr) );
   SCIP_CALL( SCIPcreateConsExprBasic(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPcomputeConsExprExprCurvature(scip, expr) );
   SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, expr, &auxvar) );

   /* add implied variable bounds */
   /* -3z1 + 3 <= x1 <= 3*z1 + 3 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_1, -3.0, 3.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_1, 3.0, 3.0, &infeas, &nbndchgs) );
   /* -z1 <= x2 <= 5*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_1, -1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_1, 5.0, 0.0, &infeas, &nbndchgs) );

   /* add bounds with z_2 */
   /* z2 <= x1 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_2, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_2, 3.0, 0.0, &infeas, &nbndchgs) );
   /* -z2 <= x2 <= 5*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_2, z_2, -1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_2, z_2, 5.0, 0.0, &infeas, &nbndchgs) );

   /* detect by convex handler */
   provided = SCIP_CONSEXPR_EXPRENFO_NONE;
   enforcebelow = FALSE;
   enforceabove = FALSE;
   success = FALSE;
   SCIP_CALL( nlhdlr_conv->detect(scip, conshdlr, nlhdlr_conv, expr, FALSE, &provided, &enforcebelow, &enforceabove,
         &success, &nlhdlrexprdata_conv) );
   providedexpected = SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
   cr_expect_eq(provided, providedexpected, "expecting %d got %d\n", providedexpected, provided);
   cr_assert(enforcebelow);
   cr_assert(!enforceabove);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata_conv);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr->enfos, 2) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &expr->enfos[0]) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &expr->enfos[1]) );
   expr->nenfos = 2;
   expr->enfos[0]->nlhdlr = nlhdlr_conv;
   expr->enfos[0]->nlhdlrexprdata = nlhdlrexprdata_conv;
   expr->enfos[1]->nlhdlr = nlhdlr;

   /* detect by perspective handler */
   success = FALSE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &provided, &enforcebelow, &enforceabove,
         &success, &nlhdlrexprdata) );

   cr_expect_eq(provided, providedexpected, "expecting provided = %d, got %d\n", providedexpected, provided);
   cr_assert(enforcebelow);
   cr_assert(!enforceabove);
   cr_assert(success);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator variables, got %d\n", nlhdlrexprdata->nindicators);

   /* check the 'off' values */
   cr_assert_eq(nlhdlrexprdata->indicators[0], z_1, "Expecting the first indicator to be z_1, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[0]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[0], 9.0, "Expecting off value = 9.0, got %f\n", nlhdlrexprdata->exprvals0[0]);

   cr_assert_eq(nlhdlrexprdata->indicators[1], z_2, "Expecting the second indicator to be z_2, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[1]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[1], 0.0, "Expecting off value = 0.0, got %f\n", nlhdlrexprdata->exprvals0[1]);

   /* separate */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x_1, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x_2, 4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z_1, 0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z_2, 0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 10) );

   SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr, expr, nlhdlrexprdata, &(expr->enfos[1]->auxvalue), sol) );
   cr_expect_eq(expr->enfos[1]->auxvalue, 16.0);

   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );

   SCIP_CALL( nlhdlrEnfoPerspective(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, expr->enfos[1]->auxvalue,
         FALSE, FALSE, FALSE, FALSE, &result) );
   cr_expect_eq(result, SCIP_SEPARATED, "Expected result = %d, got %d", SCIP_SEPARATED, result);
   cr_assert(SCIPgetNCuts(scip) == 2);

   /* check the cuts */
   cutvars = (SCIP_VAR*[4]) {z_1, x_2, x_1, auxvar};
   cutvals = (SCIP_Real[4]) {-13.0, 8.0, 4.0, -1.0};
   checkCut(SCIPgetCuts(scip)[0], cutvars, cutvals, 4, -SCIPinfinity(scip), 3.0);

   cutvars = (SCIP_VAR*[4]) {x_2, x_1, auxvar, z_2};
   cutvals = (SCIP_Real[4]) {8.0, 4.0, -1.0, -16.0};
   checkCut(SCIPgetCuts(scip)[1], cutvars, cutvals, 4, -SCIPinfinity(scip), 0.0);

   /* free memory */
   SCIP_CALL( SCIPclearCuts(scip) );
   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps) );
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) ); /* convex nlhdlr is freed when expr is released */
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}
