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

#define SCIP_DEBUG

#include <scip/sepastore.h>
#include <scip/lp.h>
#include "scip/scip.h"
#include "scip/var.h"
#include "scip/struct_lp.h"
#include "scip/struct_scip.h"
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1o;
static SCIP_VAR* x2o;
static SCIP_VAR* x3o;
static SCIP_VAR* x4o;
static SCIP_VAR* y12o;
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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1o, "x1", -1.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2o, "x2", -6.0, -3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3o, "x3", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4o, "x4", 1.0, 3.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y12o, "y12", 2.0, 4.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x1o) );
   SCIP_CALL( SCIPaddVar(scip, x2o) );
   SCIP_CALL( SCIPaddVar(scip, x3o) );
   SCIP_CALL( SCIPaddVar(scip, x4o) );
   SCIP_CALL( SCIPaddVar(scip, y12o) );

   /* get SCIP into SOLVING stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* get transformed vars and release vars */
   SCIP_CALL( SCIPgetTransformedVar(scip, x1o, &x1) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x2o, &x2) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x3o, &x3) );
   SCIP_CALL( SCIPgetTransformedVar(scip, x4o, &x4) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y12o, &y12) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3o) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4o) );
   SCIP_CALL( SCIPreleaseVar(scip, &y12o) );
   cr_assert(x1 != NULL);

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

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

   /* create a cons with some bilinear expressions */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <t_x1>*<t_x2> + <t_x1>*<t_x3> + <t_x4>*<t_x2> <= 1", TRUE, TRUE,
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
   cr_expect_eq(sepadata->varbilinvars[0][0], x2, "\nExpected varbilinvars[0][0] to be %s, got %s", SCIPvarGetName(x2), SCIPvarGetName(sepadata->varbilinvars[0][0]));
   cr_expect_eq(sepadata->varbilinvars[0][1], x3, "\nExpected varbilinvars[0][1] to be %s, got %s", SCIPvarGetName(x3), SCIPvarGetName(sepadata->varbilinvars[0][1]));
   cr_expect_eq(sepadata->varbilinvars[1][0], x4, "\nExpected varbilinvars[1][0] to be %s, got %s", SCIPvarGetName(x4), SCIPvarGetName(sepadata->varbilinvars[1][0]));

   SCIP_CALL( freeSepaData(scip, sepadata) );

   SCIPreleaseCons(scip, &cons);
   SCIPfreeBuffer(scip, &sepadata);
}

Test(rlt_selection, projection, .init = setup, .fini = teardown, .description = "test projection of problem")
{
   SCIP_ROW** rows;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   PROJLP* projlp;

   SCIPallocBufferArray(scip, &rows, 1);
   SCIPallocBufferArray(scip, &vars, 3);
   SCIPallocBufferArray(scip, &vals, 3);

   /* create test row1: -10 <= 4x1 - 7x2 + x3 <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x1, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x2, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x3, 1.0) );

   /* specify solution (only x3 is not at bound) */
   SCIPcreateSol(scip, &sol, NULL);
   vars[0] = x1; vals[0] = 5.0;
   vars[1] = x2; vals[1] = -6.0;
   vars[2] = x3; vals[2] = 2.0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 3, vars, vals) );
   cr_assert(SCIProwGetNNonz(rows[0]) == 3);

   createProjLP(scip, rows, 1, sol, &projlp);
   printProjLP(scip, projlp, 1, NULL);

   /* check results */
   cr_assert_eq(projlp->nProjNonz[0], 1, "\nExpected 1 non-zero in the projected row, got %d", projlp->nProjNonz[0]);
   cr_assert_eq(projlp->coefs[0][0], 1.0, "\nExpected coef 0 in projected row 0 to be 1.0, got %f", projlp->coefs[0][0]);
   cr_assert_eq(projlp->vars[0][0], x3, "\nExpected var 0 in projected row 0 to be x3, got %s", SCIPvarGetName(projlp->vars[0][0]));

   /* free memory */
   freeProjLP(scip, &projlp, 1);
   SCIPfreeSol(scip, &sol);
   SCIPreleaseRow(scip, &rows[0]);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}

Test(rlt_selection, execlp, .init = setup, .fini = teardown, .description = "test cut separation")
{
   SCIP_RESULT result;
   SCIP_ROW** rows;
   SCIP_SOL* sol;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_SEPADATA* sepadata;
   SCIP_Bool success, infeasible;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONS* cons;
   int c, ncuts;
   SCIP_LP* lp;
   PROJLP* projlp;

   SCIPallocBufferArray(scip, &rows, 1);
   SCIPallocBufferArray(scip, &vars, 6);
   SCIPallocBufferArray(scip, &vals, 6);

   SCIPlpCreate(&lp, scip->set, scip->messagehdlr, scip->stat, "lp");

   /* create a row: -10 <= 4x1 - 7x2 + x3 <= 5 */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &rows[0], "test_row", -10.0, 5.0, FALSE, TRUE, FALSE) );
   SCIPaddRow(scip, rows[0], FALSE, &infeasible);
   SCIPlpAddRow(lp, SCIPblkmem(scip), scip->set, scip->eventqueue, scip->eventfilter, rows[0], 0);
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x1, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x2, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, rows[0], x3, 1.0) );
   SCIPaddRow(scip, rows[0], FALSE, &infeasible);

   /* create a cons with some bilinear expressions */
   SCIP_CALL( SCIPparseCons(scip, &cons, (char*)"[expr] <test>: <t_x1>*<t_x2> + <t_x1>*<t_x3> + <t_x4>*<t_x2> <= 1", TRUE, TRUE,
                            TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   expr = SCIPgetExprConsExpr(scip, cons);
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], NULL) );
   }

   /* specify solution and local bounds */
   SCIPcreateSol(scip, &sol, NULL);
   vars[0] = x1; vals[0] = 5.0;
   vars[1] = x2; vals[1] = -6.0;
   vars[2] = x3; vals[2] = 2.0;
   vars[3] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[0]); vals[3] = -20;
   vars[4] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[1]); vals[4] = 5;
   vars[5] = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(expr)[2]); vals[5] = 0;
   SCIP_CALL( SCIPsetSolVals(scip, sol, 6, vars, vals) );
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIPevalConsExprExpr(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[c], sol, 0);
   }

   /* fill in sepadata */
   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   cr_assert(sepadata->conshdlr != NULL);
   SCIP_CALL( createSepaData(scip, sepadata) );
   sepadata->maxusedvars = 4;

   /* project the problem */
   createProjLP(scip, rows, 1, sol, &projlp);

   /* generate cuts */
   separateRltCuts(scip, sepa, sepadata, sol, rows, 1, TRUE, &ncuts, &result);

   /* TODO check results */
   SCIPinfoMessage(scip, NULL, "\n");

   /* free memory */
   freeProjLP(scip, &projlp, 1);
   SCIPclearCuts(scip);
   SCIP_CALL( freeSepaData(scip, sepadata) );
   SCIPfreeBuffer(scip, &sepadata);
   SCIPfreeSol(scip, &sol);
   SCIPreleaseCons(scip, &cons);
   SCIPreleaseRow(scip, &rows[0]);
   SCIPlpFree(&lp, SCIPblkmem(scip), scip->set, scip->eventqueue, scip->eventfilter);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &rows);
}
