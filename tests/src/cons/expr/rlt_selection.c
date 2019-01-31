/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   rlt.c
 * @brief  tests rlt cut selection
 * @author Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/var.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* y12;
static SCIP_VAR* auxvar;
static SCIP_SOL* sol;

/* creates scip, problem, includes expression constraint handler, creates  and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* include rlt separator */
   SCIP_CALL( SCIPincludeSepaRlt(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* get rlt separator */
   sepa = SCIPfindSepa(scip, "rlt");
   assert(sepa != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3, "x3", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4, "x4", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y12, "y12", 2.0, 4.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, x3) );
   SCIP_CALL( SCIPaddVar(scip, x4) );
   SCIP_CALL( SCIPaddVar(scip, y12) );

   /* get SCIP into SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseVar(scip, &y12) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(rlt_selection, sepadata, .init = setup, .fini = teardown, .description = "test creation and freeing of separator data")
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_SEPADATA* sepadata;
   SCIP_CONSEXPR_EXPR* expr;
   int c;

   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(sepadata->conshdlr != NULL);

   /* create a cons some bilinear expressions */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <x1>*<x2> + <x1>*<x3> + <x4>*<x2> <= 1", TRUE, TRUE,
                 TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   expr = SCIPgetExprConsExpr(scip, cons);
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], NULL) );
   }

   SCIP_CALL( createSepaData(scip, sepadata) );

   cr_expect_eq(sepadata->nbilinvars, 4, "\nExpected 4 bilinear vars, got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nvarbilinvars[0], 2, "\nExpected 2 bilinear vars for entry [0], got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nvarbilinvars[1], 1, "\nExpected 1 bilinear vars for entry [1], got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nvarbilinvars[2], 0, "\nExpected 0 bilinear vars for entry [2], got %d", sepadata->nbilinvars);
   cr_expect_eq(sepadata->nvarbilinvars[3], 0, "\nExpected 0 bilinear vars for entry [3], got %d", sepadata->nbilinvars);
   cr_assert(sepadata->varbilinvars[0][0] == x2);
   cr_assert(sepadata->varbilinvars[0][1] == x3);
   cr_assert(sepadata->varbilinvars[1][0] == x4);

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIPreleaseCons(scip, &cons);
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_selection, projection, .init = setup, .fini = teardown, .description = "test projection of problem")
{
   SCIP_ROW** rows;
   SCIP_ROW** proj_rows;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;

   SCIPallocBufferArray(scip, &rows, 1);
   SCIPallocBufferArray(scip, &vars, 3);
   SCIPallocBufferArray(scip, &vals, 3);
   SCIPallocBufferArray(scip, &proj_rows, 1);

   /* create test row1: -10 <= 4x1 - 7x2 + x3 <= 5 and add it to rows */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x1, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x2, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x3, 1.0) );

   /* specify solution and local bounds */
   SCIPcreateSol(scip, &sol, NULL);
   vars[0] = x1; vals[0] = 5.0;
   vars[1] = x2; vals[1] = -6.0;
   vars[2] = x3; vals[2] = 2.0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 3, vars, vals) );
   cr_assert(SCIProwGetNNonz(rows[0]) == 3);

   createProjectedProb(scip, sepa, rows, 1, sol, proj_rows);

   /* check results */
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPprintRow(scip, rows[0], NULL);
   SCIPprintRow(scip, proj_rows[0], NULL);
   cr_assert_eq(SCIProwGetNNonz(proj_rows[0]), 1, "\nExpected 1 non-zero in the projected row, got %d", SCIProwGetNNonz(proj_rows[0]));

   /* free memory */
   SCIPreleaseRow(scip, &proj_rows[0]);
   SCIPfreeSol(scip, &sol);
   SCIPreleaseRow(scip, &rows[0]);
   SCIPfreeBufferArray(scip, &proj_rows);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}
