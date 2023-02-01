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

/**@file   separation.c
 * @brief  tests rlt cut separation
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "scip/cons_nonlinear.h"
#include "scip/sepastore.h"
#include "scip/scip.h"
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* b1;
static SCIP_VAR* b2;

/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_VAR* x1o;
   SCIP_VAR* x2o;
   SCIP_VAR* x3o;
   SCIP_VAR* x4o;
   SCIP_VAR* b1o;
   SCIP_VAR* b2o;

   SCIP_CALL( SCIPcreate(&scip) );

   /* includes expression handlers as well as nonlinear constraint handler and rlt separator */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* get nonlinear conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   assert(conshdlr != NULL);

   /* get rlt separator */
   sepa = SCIPfindSepa(scip, "rlt");
   assert(sepa != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1o, "x1", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2o, "x2", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3o, "x3", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4o, "x4", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b1o, "b1", 0, 1, 1, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b2o, "b2", 0, 1, 1, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1o) );
   SCIP_CALL( SCIPaddVar(scip, x2o) );
   SCIP_CALL( SCIPaddVar(scip, x3o) );
   SCIP_CALL( SCIPaddVar(scip, x4o) );
   SCIP_CALL( SCIPaddVar(scip, b1o) );
   SCIP_CALL( SCIPaddVar(scip, b2o) );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

   /* get SCIP into SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, x1o, &x1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x2o, &x2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x3o, &x3) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x4o, &x4) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b1o, &b1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b2o, &b2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4o) );
   SCIP_CALL( SCIPreleaseVar(scip, &b1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &b2o) );
   cr_assert(x1 != NULL);
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

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
   cr_expect_eq(SCIProwGetNNonz(cut), nvars, "\nExpected %d nonz, got %d", nvars, SCIProwGetNNonz(cut));
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

      if( !found )
         cr_expect(FALSE, "found an unknown variable");
   }
}

Test(separation, sepadata, .init = setup, .fini = teardown, .description = "test creation and freeing of separator data")
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeasible;
   SCIP_SEPADATA* sepadata;
   ADJACENTVARDATA* adjvardata;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   cr_assert(sepadata->conshdlr != NULL);
   sepadata->maxusedvars = DEFAULT_MAXUSEDVARS;

   /* create a cons with some bilinear expressions */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[nonlinear] <test>: <t_x1>*<t_x2> + <t_x1>*<t_x3> + <t_x4>*<t_x2> + <t_x4>^2 <= 1",
                 TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   SCIP_CALL( SCIPaddCons(scip, cons) ); /* adds locks */

   /* creates auxvars and creates disaggregation variables and row */
   SCIP_CALL( SCIPconstructLP(scip, &infeasible) );
   cr_assert_not(infeasible);

   SCIP_CALL( SCIPcollectBilinTermsNonlinear(scip, conshdlr, &cons, 1) );

   SCIP_CALL( createSepaData(scip, sepadata) );

   cr_expect_eq(sepadata->nbilinvars, 4, "\nExpected 4 bilinear vars, got %d", sepadata->nbilinvars);

   cr_expect_eq(sepadata->varssorted[0], x4, "\nExpected varssorted[0] to be x4, got %s", SCIPvarGetName(sepadata->varssorted[0]));
   cr_expect_eq(sepadata->varssorted[1], x1, "\nExpected varssorted[1] to be x1, got %s", SCIPvarGetName(sepadata->varssorted[1]));
   cr_expect_eq(sepadata->varssorted[2], x2, "\nExpected varssorted[2] to be x2, got %s", SCIPvarGetName(sepadata->varssorted[2]));
   cr_expect_eq(sepadata->varssorted[3], x3, "\nExpected varssorted[3] to be x3, got %s", SCIPvarGetName(sepadata->varssorted[3]));

   adjvardata = (ADJACENTVARDATA*) SCIPhashmapGetImage(sepadata->bilinvardatamap, (void*)(size_t) SCIPvarGetIndex(x1));
   cr_assert(adjvardata != NULL);
   cr_expect_eq(adjvardata->nadjacentvars, 2, "\nExpected 2 bilinear vars for x1, got %d", adjvardata->nadjacentvars);
   cr_expect_eq(adjvardata->adjacentvars[0], x3, "\nBilinear var 0 for x1 should be x3, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[0]));
   cr_expect_eq(adjvardata->adjacentvars[1], x2, "\nBilinear var 1 for x1 should be x2, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[1]));

   adjvardata = (ADJACENTVARDATA*) SCIPhashmapGetImage(sepadata->bilinvardatamap, (void*)(size_t) SCIPvarGetIndex(x2));
   cr_assert(adjvardata != NULL);
   cr_expect_eq(adjvardata->nadjacentvars, 2, "\nExpected 2 bilinear vars for x2, got %d", adjvardata->nadjacentvars);
   cr_expect_eq(adjvardata->adjacentvars[0], x4, "\nBilinear var 0 for x2 should be x4, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[0]));
   cr_expect_eq(adjvardata->adjacentvars[1], x1, "\nBilinear var 1 for x2 should be x1, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[1]));

   adjvardata = (ADJACENTVARDATA*) SCIPhashmapGetImage(sepadata->bilinvardatamap, (void*)(size_t) SCIPvarGetIndex(x3));
   cr_assert(adjvardata != NULL);
   cr_expect_eq(adjvardata->nadjacentvars, 1, "\nExpected 1 bilinear vars for x3, got %d", adjvardata->nadjacentvars);
   cr_expect_eq(adjvardata->adjacentvars[0], x1, "\nBilinear var 0 for x3 should be x1, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[0]));

   adjvardata = (ADJACENTVARDATA*) SCIPhashmapGetImage(sepadata->bilinvardatamap, (void*)(size_t) SCIPvarGetIndex(x4));
   cr_assert(adjvardata != NULL);
   cr_expect_eq(adjvardata->nadjacentvars, 2, "\nExpected 2 bilinear vars for x4, got %d", adjvardata->nadjacentvars);
   cr_expect_eq(adjvardata->adjacentvars[0], x4, "\nBilinear var 0 for x4 should be x4, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[0]));
   cr_expect_eq(adjvardata->adjacentvars[1], x2, "\nBilinear var 1 for x4 should be x2, got %s",
         SCIPvarGetName(adjvardata->adjacentvars[1]));

   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIPfreeBuffer(scip, &sepadata);

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* frees the disaggregation row */
   SCIP_CALL( SCIPclearCuts(scip) );
}

Test(separation, projection, .init = setup, .fini = teardown, .description = "test projection of problem")
{
   SCIP_ROW** rows;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   RLT_SIMPLEROW* projrows;
   SCIP_Bool allcst;

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 3) );

   /* create test row1: -10 <= 4x1 - 7x2 + x3 <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x1, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x2, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x3, 1.0) );
   cr_assert(SCIProwGetNNonz(rows[0]) == 3);

   /* specify solution (only x3 is not at bound) */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   vars[0] = x1; vals[0] = 5.0;
   vars[1] = x2; vals[1] = -6.0;
   vars[2] = x3; vals[2] = 2.0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 3, vars, vals) );

   SCIP_CALL( createProjRows(scip, rows, 1, sol, &projrows, TRUE, &allcst) );

   /* check results */

   /* the projected cut should be: -72 <= x3 <= -57 */
   cr_assert_eq(projrows[0].nnonz, 1, "\nExpected 1 non-zero in the projected row, got %d", projrows[0].nnonz);
   cr_assert_eq(projrows[0].coefs[0], 1.0, "\nExpected coef 0 in projected row 0 to be 1.0, got %f", projrows[0].coefs[0]);
   cr_assert_eq(projrows[0].vars[0], x3, "\nExpected var 0 in projected row 0 to be x3, got %s", SCIPvarGetName(projrows[0].vars[0]));
   cr_assert_eq(projrows[0].cst, 0.0, "\nExpected the const in projected row to be 0.0, got %f", projrows[0].cst);
   cr_assert_eq(projrows[0].lhs, -72.0, "\nExpected the lhs in projected row to be -72.0, got %f", projrows[0].lhs);
   cr_assert_eq(projrows[0].rhs, -57.0, "\nExpected the rhs in projected row to be -57.0, got %f", projrows[0].rhs);

   /* free memory */
   freeProjRows(scip, &projrows, 1);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseRow(scip, &rows[0]) );
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}

Test(separation, compute_projcut, .init = setup, .fini = teardown, .description = "test projected cut computation")
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   RLT_SIMPLEROW projrow;
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* cut;
   SCIP_ROW* row;
   SCIP_Bool success;
   SCIP_Real cut_val;

   /* create row: -10 <= x1 + 2x2 - x3 <= 20 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, "test_row", -10.0, 20.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x1, 1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x2, 2.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row, x3, -1.0) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 3) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 3) );

   /* specify solution (none of the variables is at bound) */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   vars[0] = x1; vals[0] = 0.0;
   vars[1] = x2; vals[1] = -5.0;
   vars[2] = x3; vals[2] = 2.0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 3, vars, vals) );

   /* fill in sepadata */
   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   cr_assert(sepadata->conshdlr != NULL);
   sepadata->maxusedvars = 4;

   /* create projected LP with row -10 <= x1 + 2x2 - x3 <= 20 */
   SCIP_CALL( createProjRow(scip, &(projrow), row, sol, FALSE) );

   /* compute a cut with x1, lb and lhs */
   SCIP_CALL( computeRltCut(scip, sepa, sepadata, &cut, NULL, &projrow, sol, NULL, NULL, x1, &success, TRUE, TRUE,
         FALSE, FALSE, TRUE) );
   assert(success);

   /* the cut should be -8 <= 8x1 */
   cut_val = 8.0;
   checkCut(cut, &x1, &cut_val, 1, -8.0, SCIPinfinity(scip));

   /* free memory */
   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   freeProjRow(scip, &projrow);
   SCIPfreeBuffer(scip, &sepadata);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
}

Test(separation, compute_clique_cuts, .init = setup, .fini = teardown, .description = "test cut computation when cliques are present")
{
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* cut;
   SCIP_ROW** rows;
   SCIP_Bool success;
   SCIP_Bool infeasible;
   SCIP_VAR* clique_vars[2];
   SCIP_Bool clique_vals[2];
   int nbdchgs;
   SCIP_VAR** cut_vars;
   SCIP_Real* cut_vals;

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, 2) );

   /* create test row1: -10 <= b2 <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], b2, 1.0) );

   scip->stat->nnz = 1;

   /* specify solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   vars[0] = b1;
   vals[0] = 1;
   vars[1] = b2;
   vals[1] = 0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 2, vars, vals) );

   /* fill in sepadata */
   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = conshdlr;
   cr_assert(sepadata->conshdlr != NULL);

   /*add a clique (1-b1) + (1-b2) <= 1*/
   clique_vars[0] = b1;
   assert(clique_vars[0] != NULL);
   clique_vars[1] = b2;
   assert(clique_vars[1] != NULL);
   clique_vals[0] = 0;
   clique_vals[1] = 0;
   SCIP_CALL( SCIPaddClique(scip, clique_vars, clique_vals, 2, FALSE, &infeasible, &nbdchgs) );

   /* compute a cut with b1, lb and lhs */
   SCIP_CALL( computeRltCut(scip, sepa, sepadata, &cut, rows[0], NULL, sol, NULL, NULL, b1, &success, TRUE, TRUE,
         FALSE, FALSE, FALSE) );
   assert(success);

   /* the cut should be 1 <= 11b1 + b2 */
   cut_vars = (SCIP_VAR*[4]) {b2, b1};
   cut_vals = (SCIP_Real[4]) {1.0, 11.0};
   checkCut(cut, cut_vars, cut_vals, 2, 1.0, SCIPinfinity(scip));

   /* free memory */
   if( cut != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &cut) );
   }
   SCIP_CALL( SCIPreleaseRow(scip, &rows[0]) );
   SCIPfreeBuffer(scip, &sepadata);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}
