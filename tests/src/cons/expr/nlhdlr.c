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

/**@file   nlhdlr.c
 * @brief  tests basic nonlinear handler methods
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_expr.h"

#include "include/scip_test.h"

static SCIP* testscip;
static SCIP_VAR* x;
static SCIP_VAR* y;

struct SCIP_ConsExpr_NlHdlrData
{
   int dummy;
};

static SCIP_CONSEXPR_NLHDLRDATA testnlhdlrdata = { .dummy = 42 };

/* creates scip, problem, includes expression constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   const char* input1 = "[expr] <test>: (<x>-0)^2 + (<y>-0)^2 <= 1.5;";
   const char* input2 = "[expr] <test>: (<x>-1)^2 + (<y>-1)^2 <= 1.0;";

   /* create scip with all plugins */
   SCIP_CALL( SCIPcreate(&testscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(testscip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(testscip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(testscip, &x, "x", 0.0, 5.0, -1.5, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(testscip, &y, "y", 0.0, 5.0, -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(testscip, x) );
   SCIP_CALL( SCIPaddVar(testscip, y) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(testscip, &consexpr, input1,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(testscip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(testscip, &consexpr) );

   success = FALSE;
   SCIP_CALL( SCIPparseCons(testscip, &consexpr, input2,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(testscip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(testscip, &consexpr) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(testscip, &x) );
   SCIP_CALL( SCIPreleaseVar(testscip, &y) );
   SCIP_CALL( SCIPfree(&testscip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(freeHdlrData)
{
   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(nlhdlrdata != NULL);
   cr_assert(*nlhdlrdata == &testnlhdlrdata, "nlhdlrdata is not the one that was passed in by SCIPincludeConsExprNlHdlrBasic");

   *nlhdlrdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(copyHdlr)
{
   SCIP_CONSEXPR_NLHDLR* targetnlhdlr;

   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourceconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   cr_assert(strcmp(SCIPgetConsExprNlHdlrName(sourcenlhdlr), "testhdlr") == 0, "source nlhdlr is not testhdlr");

   SCIP_CALL( SCIPincludeConsExprNlHdlrBasic(targetscip, targetconsexprhdlr, &targetnlhdlr,
      SCIPgetConsExprNlHdlrName(sourcenlhdlr), SCIPgetConsExprNlHdlrDesc(sourcenlhdlr), SCIPgetConsExprNlHdlrPriority(sourcenlhdlr), SCIPgetConsExprNlHdlrData(sourcenlhdlr)) );
   SCIPsetConsExprNlHdlrFreeHdlrData(targetscip, targetnlhdlr, freeHdlrData);
   SCIPsetConsExprNlHdlrCopyHdlr(testscip, targetnlhdlr, copyHdlr);

   return SCIP_OKAY;
}


Test(conshdlr, nlhdlr, .init = setup, .fini = teardown,
   .description = "test basic functionality of nonlinear handler of the cons_expr constraint handler."
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(testscip, "expr");
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPincludeConsExprNlHdlrBasic(testscip, conshdlr, &nlhdlr, "testhdlr", "tests nonlinear handler functionality", 0, &testnlhdlrdata) );

   SCIPsetConsExprNlHdlrFreeHdlrData(testscip, nlhdlr, freeHdlrData);
   SCIPsetConsExprNlHdlrCopyHdlr(testscip, nlhdlr, copyHdlr);

   SCIP_CALL( SCIPsetIntParam(testscip, "display/verblevel", SCIP_VERBLEVEL_NONE) );
   SCIP_CALL( SCIPsolve(testscip) );
   /* SCIP_CALL( SCIPprintBestSol(testscip, NULL, TRUE) ); */
}
