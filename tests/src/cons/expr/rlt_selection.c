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
#include "scip/sepa_rlt.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SEPA* sepa;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
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
   SCIP_CALL( SCIPcreateVarBasic(scip, &y12, "y12", 2.0, 4.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, x3) );
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
   SCIP_CALL( SCIPreleaseVar(scip, &x3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(rlt_selection, projection, .init = setup, .fini = teardown, .description = "test projection of problem")
{
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* row1;
   SCIP_ROW** rows;
   SCIP_ROW** proj_rows;

   SCIPallocBufferArray(scip, &rows, 1);

   /* create test row1: -10 <= 4x1 - 7x2 + x3 <= 5 and add it to rows */
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row1, "test_row", -10.0, 5.0, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, x1, 4.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, x2, -7.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, row1, x3, 1.0) );
   rows[0] = row1;

   /* allocate and fill in sepadata */
   SCIP_CALL( SCIPallocBuffer(scip, &sepadata) );
   sepadata->conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(sepadata->conshdlr != NULL);

   /* add x1 to varssorted */
   SCIP_CALL( SCIPallocBufferArray(scip, &sepadata->varssorted, 1) );
   sepadata->varssorted[0] = x1;
   sepadata->nbilinvars = 1;

   /* add information about bilinear terms to sepadata */
   sepadata->nbilinterms = 1;

   /* TODO fill in sepadata */

   /* TODO specify solution */

   SCIPallocBufferArray(scip, &proj_rows, 1);

   createProjectedProb(scip, sepa, sepadata, rows, 1, NULL, proj_rows);

   /* TODO check results */

   /* free memory */
   SCIPfreeBufferArray(scip, &proj_rows);
   SCIPfreeBufferArray(scip, &sepadata->varssorted);
   SCIPfreeBuffer(scip, &sepadata);
   SCIPreleaseRow(scip, &row1);
   SCIPfreeBufferArray(scip, &rows);
}
