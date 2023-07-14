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

/**@file   nlhdlr_perspective.c
 * @brief  tests perspective nonlinear handler methods
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"
#include "scip/nlhdlr_perspective.c"


/*
 * TEST
 */

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
static SCIP_NLHDLR* nlhdlr = NULL;
static SCIP_NLHDLR* nlhdlr_conv = NULL;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   cr_assert_not_null(conshdlr);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   cr_assert_not_null(conshdlrdata);

   /* get perspective and convex nlhdlrs */
   nlhdlr = SCIPfindNlhdlrNonlinear(conshdlr, NLHDLR_NAME);
   nlhdlr_conv = SCIPfindNlhdlrNonlinear(conshdlr, "convex");

   cr_assert_not_null(nlhdlr);
   cr_assert_not_null(nlhdlr_conv);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* TODO test also with adjrefpoint enabled (needs adjustment in expected cuts in sepa test) */
   SCIP_CALL( SCIPsetBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/adjrefpoint", FALSE) );

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

      cr_expect(found, "variable %s must be in the cut", SCIPvarGetName(var));
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
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x1>^2 + <x1>*<y1> + <x2>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPcomputeExprCurvature(scip, expr) );
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );

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
   ownerdata = SCIPexprGetOwnerData(expr);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ownerdata->enfos, 2) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[0]) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[1]) );
   ownerdata->nenfos = 2;
   ownerdata->enfos[0]->nlhdlr = nlhdlr_conv;
   ownerdata->enfos[0]->nlhdlrparticipation = SCIP_NLHDLR_METHOD_SEPABELOW;
   ownerdata->enfos[1]->nlhdlr = nlhdlr;

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_SEPABELOW, "expecting sepabelow, got %d\n", participating);
   cr_assert_eq(enforcing, SCIP_NLHDLR_METHOD_NONE);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator variables, got %d\n", nlhdlrexprdata->nindicators);

   /* compute and check the 'off' values */
   SCIP_CALL( nlhdlrInitSepaPerspective(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, TRUE, TRUE, &infeas) );
   cr_assert_eq(nlhdlrexprdata->indicators[0], z_1, "Expecting the first indicator to be z_1, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[0]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[0], 9.0, "Expecting off value = 9.0, got %f\n", nlhdlrexprdata->exprvals0[0]);

   cr_assert_eq(nlhdlrexprdata->indicators[1], z_2, "Expecting the second indicator to be z_2, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[1]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[1], 0.0, "Expecting off value = 0.0, got %f\n", nlhdlrexprdata->exprvals0[1]);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIPfreeBlockMemory(scip, &ownerdata->enfos[1]);
   SCIPfreeBlockMemory(scip, &ownerdata->enfos[0]);
   SCIPfreeBlockMemoryArray(scip, &ownerdata->enfos, 2);
   ownerdata->nenfos = 0;

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects x1^2 as an on/off expression */
Test(nlhdlrperspective, detectandfree2, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x1>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPcomputeExprCurvature(scip, expr) );
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );

   /* z1 <= x1 <= 3*z1 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_1, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_1, 3.0, 0.0, &infeas, &nbndchgs) );
   /* z2 <= x1 <= 3*z2 */
   SCIP_CALL( SCIPaddVarVlb(scip, x_1, z_2, 1.0, 0.0, &infeas, &nbndchgs) );
   SCIP_CALL( SCIPaddVarVub(scip, x_1, z_2, 3.0, 0.0, &infeas, &nbndchgs) );

   /* detect wants another nlhdlr to have detected; pretend that convex has detected */
   ownerdata = SCIPexprGetOwnerData(expr);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ownerdata->enfos, 2) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[0]) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[1]) );
   ownerdata->nenfos = 2;
   ownerdata->enfos[0]->nlhdlr = nlhdlr_conv;
   ownerdata->enfos[0]->nlhdlrparticipation = SCIP_NLHDLR_METHOD_SEPABELOW;
   ownerdata->enfos[1]->nlhdlr = nlhdlr;

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_SEPABELOW, "expecting participating = sepabelow, got %d\n",
         participating);
   cr_assert_eq(enforcing, SCIP_NLHDLR_METHOD_NONE);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator variables, got %d\n", nlhdlrexprdata->nindicators);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIPfreeBlockMemory(scip, &ownerdata->enfos[1]);
   SCIPfreeBlockMemory(scip, &ownerdata->enfos[0]);
   SCIPfreeBlockMemoryArray(scip, &ownerdata->enfos, 2);
   ownerdata->nenfos = 0;

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* detects log(x1+x2+1) as an on/off expression */
Test(nlhdlrperspective, detectandfree3, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_EXPR* expr;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"log(<x1>+<x2>+1.0)", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPcomputeExprCurvature(scip, expr) );
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );

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

   /* detect wants another nlhdlr to have detected; pretend that convex has detected */
   ownerdata = SCIPexprGetOwnerData(expr);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ownerdata->enfos, 2) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[0]) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[1]) );
   ownerdata->nenfos = 2;
   ownerdata->enfos[0]->nlhdlr = nlhdlr_conv;
   ownerdata->enfos[0]->nlhdlrparticipation = SCIP_NLHDLR_METHOD_SEPAABOVE;
   ownerdata->enfos[1]->nlhdlr = nlhdlr;

   /* detect */
   enforcing = SCIP_NLHDLR_METHOD_NONE;
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_SEPAABOVE, "expecting sepaabove, got %d\n", participating);
   cr_assert_eq(enforcing, SCIP_NLHDLR_METHOD_NONE);
   cr_assert_not_null(nlhdlrexprdata);

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator vars, got %d\n", nlhdlrexprdata->nindicators);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, &nlhdlrexprdata);

   SCIPfreeBlockMemory(scip, &ownerdata->enfos[1]);
   SCIPfreeBlockMemory(scip, &ownerdata->enfos[0]);
   SCIPfreeBlockMemoryArray(scip, &ownerdata->enfos, 2);
   ownerdata->nenfos = 0;

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* separates x1^2 + x1y1 + x2^2 for 2 indicator variables */
Test(nlhdlrperspective, sepa, .init = setup, .fini = teardown)
{
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata = NULL;
   SCIP_NLHDLREXPRDATA* nlhdlrexprdata_conv = NULL;
   SCIP_EXPR* expr;
   SCIP_NLHDLR_METHOD enforcing;
   SCIP_NLHDLR_METHOD participating;
   SCIP_NLHDLR_METHOD enforcing_conv;
   SCIP_NLHDLR_METHOD participating_conv;
   SCIP_Bool infeas;
   SCIP_CONS* cons;
   int nbndchgs;
   SCIP_SOL* sol;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_VAR* auxvar;
   SCIP_PTRARRAY* rowpreps;
   SCIP_RESULT result;
   SCIP_EXPR_OWNERDATA* ownerdata;

   /* skip when no ipopt */
   if( ! SCIPisIpoptAvailableIpopt() )
      return;

   /* create expression and constraint */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)"<x1>^2 + <x1>*<x2> + <x2>^2", NULL, NULL, NULL) );
   SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, (char*)"nlin", expr, -SCIPinfinity(scip), 0)  );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPaddConsLocks(scip, cons, 1, 0) );
   SCIP_CALL( SCIPcomputeExprCurvature(scip, expr) );
   SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, expr, TRUE, FALSE, FALSE, FALSE) );

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
   enforcing_conv = SCIP_NLHDLR_METHOD_NONE;
   participating_conv = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( SCIPnlhdlrDetect(scip, conshdlr, nlhdlr_conv, expr, cons, &enforcing_conv, &participating_conv,
         &nlhdlrexprdata_conv) );
   cr_expect_eq(participating_conv, SCIP_NLHDLR_METHOD_SEPABELOW, "expecting participating_conv = sepabelow, got %d\n",
         participating);
   cr_assert_eq(enforcing_conv, SCIP_NLHDLR_METHOD_SEPABELOW);
   cr_assert_not_null(nlhdlrexprdata_conv);

   /* although convex has already detected, it doesn't fill the ownerdata, so we do it here */
   ownerdata = SCIPexprGetOwnerData(expr);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ownerdata->enfos, 2) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[0]) );
   SCIP_CALL( SCIPallocBlockMemory(scip, &ownerdata->enfos[1]) );
   ownerdata->nenfos = 2;
   ownerdata->enfos[0]->nlhdlr = nlhdlr_conv;
   ownerdata->enfos[0]->nlhdlrexprdata = nlhdlrexprdata_conv;
   ownerdata->enfos[0]->nlhdlrparticipation = participating_conv;
   ownerdata->enfos[1]->nlhdlr = nlhdlr;

   /* detect by perspective handler */
   participating = SCIP_NLHDLR_METHOD_NONE;
   SCIP_CALL( nlhdlrDetectPerspective(scip, conshdlr, nlhdlr, expr, cons, &enforcing, &participating, &nlhdlrexprdata) );
   cr_expect_eq(participating, SCIP_NLHDLR_METHOD_SEPABELOW, "expecting sepabelow, got %d\n", participating);
   cr_assert_not_null(nlhdlrexprdata);
   ownerdata->enfos[1]->nlhdlrexprdata = nlhdlrexprdata;

   cr_expect_eq(nlhdlrexprdata->nindicators, 2, "Expecting 2 indicator variables, got %d\n", nlhdlrexprdata->nindicators);

   /* compute and check the 'off' values */
   SCIP_CALL( nlhdlrInitSepaPerspective(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, TRUE, TRUE, &infeas) );
   cr_assert_eq(nlhdlrexprdata->indicators[0], z_1, "Expecting the first indicator to be z_1, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[0]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[0], 9.0, "Expecting off value = 9.0, got %f\n", nlhdlrexprdata->exprvals0[0]);

   cr_assert_eq(nlhdlrexprdata->indicators[1], z_2, "Expecting the second indicator to be z_2, got %s\n", SCIPvarGetName(nlhdlrexprdata->indicators[1]));
   cr_assert_eq(nlhdlrexprdata->exprvals0[1], 0.0, "Expecting off value = 0.0, got %f\n", nlhdlrexprdata->exprvals0[1]);

   /* make sure there is an auxvar; since expr is not part of a constraint, we cannot lean on cons_nonlinear to do that for us TODO but it is part of a constraint? */
   SCIP_CALL( createAuxVar(scip, expr) );
   auxvar = SCIPgetExprAuxVarNonlinear(expr);

   /* initsepa of nlhdlr convex, so its estimate can be called later */
   SCIP_CALL( SCIPnlhdlrInitsepa(scip, conshdlr, cons, nlhdlr_conv, expr, nlhdlrexprdata_conv, FALSE, TRUE, &infeas) );
   /* clear cuts from initsepa, so we can do a cutcount later */
   SCIP_CALL( SCIPclearCuts(scip) );

   /* separate */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x_1, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x_2, 4.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z_1, 0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z_2, 0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, 10) );

   SCIP_CALL( SCIPnlhdlrEvalaux(scip, nlhdlr, expr, nlhdlrexprdata, &(ownerdata->enfos[1]->auxvalue), sol) );
   cr_expect_eq(ownerdata->enfos[1]->auxvalue, 16.0);

   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps) );

   SCIP_CALL( nlhdlrEnfoPerspective(scip, conshdlr, cons, nlhdlr, expr, nlhdlrexprdata, sol, ownerdata->enfos[1]->auxvalue,
         FALSE, FALSE, FALSE, FALSE, &result) );
   cr_expect_eq(result, SCIP_SEPARATED, "Expected enfo result = %d, got %d", SCIP_SEPARATED, result);
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

   SCIP_CALL( SCIPaddConsLocks(scip, cons, -1, 0) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
