/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   create.c
 * @brief  tests creation of non-linear rows with expression constraints
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/struct_nlp.h>
#include "scip/cons_expr.h"
#include "scip/cons_expr.c"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_CONS* consexpr;
static SCIP_CONSDATA* consdata;
static SCIP_NLROW* nlrow;
static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* x3;
static SCIP_VAR* x4;
static SCIP_VAR* x5;
static SCIP_Bool success;
static const char* input;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x3, "x3", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x4, "x4", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x5, "x5", 0.0, 1.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, x3) );
   SCIP_CALL( SCIPaddVar(scip, x4) );
   SCIP_CALL( SCIPaddVar(scip, x5) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x3) );
   SCIP_CALL( SCIPreleaseVar(scip, &x4) );
   SCIP_CALL( SCIPreleaseVar(scip, &x5) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* TEST SUITE */
TestSuite(test_create_nlrow, .init = setup, .fini = teardown);

Test(test_create_nlrow, noquad, .init = setup, .fini = teardown)
{
   input = "[expr] <test>: 2*<x1> + 3.2*<x2> + 0.5*<x3>^3 - 4*<x4> + <x5> + 10 <= 2";

   /* create constraint from input string */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add constraint to SCIP and release it */
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   consdata = SCIPconsGetData(consexpr);
   cr_assert(consdata != NULL);

   /* goto presolved stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   SCIP_CALL( createNlRow(scip, consexpr) );

   nlrow = consdata->nlrow;
   cr_assert(nlrow != NULL);

   /* check linear part */
   cr_assert_eq(nlrow->constant, 10);
   cr_assert_eq(nlrow->nlinvars, 4);
   cr_assert_eq(nlrow->linvars[0], x1);
   cr_assert_eq(nlrow->linvars[1], x2);
   cr_assert_eq(nlrow->linvars[2], x4);
   cr_assert_eq(nlrow->linvars[3], x5);

   /* check quadratic part */
   cr_assert_eq(nlrow ->nquadelems, 0);
   cr_assert_eq(nlrow ->nquadvars, 0);
   cr_assert(nlrow ->quadelems == NULL);
   cr_assert(nlrow ->quadvars == NULL);
   cr_assert(nlrow ->quadvarshash == NULL);

   /* check non-quadratic part */
   cr_assert(nlrow->exprtree != NULL);
   cr_assert_eq(SCIPexprtreeGetNVars(nlrow->exprtree), 1);
   cr_assert_eq(SCIPexprtreeGetVars(nlrow->exprtree)[0], x3);

   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

Test(test_create_nlrow, nolin, .init = setup, .fini = teardown)
{
   input = "[expr] <test>: 2*<x1>^2 + 3.2*<x1>*<x2> + 0.5*<x3>^3 - 4*<x4>*<x5> <= 2";

   /* create constraint from input string */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add constraint to SCIP and release it */
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   consdata = SCIPconsGetData(consexpr);
   cr_assert(consdata != NULL);

   /* goto presolved stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   SCIP_CALL( createNlRow(scip, consexpr) );

   nlrow = consdata->nlrow;
   cr_assert(nlrow != NULL);

   /* check linear part */
   cr_assert_eq(nlrow->constant, 0);
   cr_assert_eq(nlrow->nlinvars, 0);

   /* check quadratic part */
   cr_assert_eq(nlrow ->nquadelems, 3);
   cr_assert_eq(nlrow ->nquadvars, 4);
   cr_assert_eq(nlrow ->quadelems[0].coef, 2.0);
   cr_assert_eq(nlrow ->quadelems[0].idx1, 0);
   cr_assert_eq(nlrow ->quadelems[0].idx2, 0);
   cr_assert_eq(nlrow ->quadelems[1].coef, 3.2);
   cr_assert_eq(nlrow ->quadelems[1].idx2, 0);
   cr_assert_eq(nlrow ->quadelems[1].idx2, 1);
   cr_assert_eq(nlrow ->quadelems[2].coef, -4.0);
   cr_assert_eq(nlrow ->quadelems[2].idx2, 3);
   cr_assert_eq(nlrow ->quadelems[2].idx2, 4);
   cr_assert_eq(nlrow ->quadvars[0], x1);
   cr_assert_eq(nlrow ->quadvars[1], x2);
   cr_assert_eq(nlrow ->quadvars[2], x4);
   cr_assert_eq(nlrow ->quadvars[3], x5);

   /* check non-quadratic part */
   cr_assert(nlrow->exprtree != NULL);
   cr_assert_eq(SCIPexprtreeGetNVars(nlrow->exprtree), 1);
   cr_assert_eq(SCIPexprtreeGetVars(nlrow->exprtree)[0], x3);

   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

Test(test_create_nlrow, complex, .init = setup, .fini = teardown)
{
   input = "[expr] <test>: 2*<x1>^2 + <x1> + 3.2*<x1>*<x2> + <x2>^2 +exp(<x2>) + 0.5*<x3>^3 - 4*<x4>*<x5> + 5*<x4> - 10*<x5>^2 + <x1>*<x3>*<x5> - 1 <= 2";

   /* create constraint from input string */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add constraint to SCIP and release it */
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   consdata = SCIPconsGetData(consexpr);
   cr_assert(consdata != NULL);

   cr_assert(consdata->expr != NULL);

   /* goto presolved stage */
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVED, FALSE) );

   SCIP_CALL( createNlRow(scip, consexpr) );

   nlrow = consdata->nlrow;
   cr_assert(nlrow != NULL);

   /* check linear part */
   cr_assert_eq(nlrow->constant, 1);
   cr_assert_eq(nlrow->nlinvars, 2);
   cr_assert_eq(nlrow->linvars[0], x1);
   cr_assert_eq(nlrow->linvars[1], x4);

   /* check quadratic part */
   cr_assert_eq(nlrow ->nquadelems, 5);
   cr_assert_eq(nlrow ->nquadvars, 4);
   cr_assert_eq(nlrow ->quadelems[0].coef, 2.0);
   cr_assert_eq(nlrow ->quadelems[0].idx1, 0);
   cr_assert_eq(nlrow ->quadelems[0].idx2, 0);
   cr_assert_eq(nlrow ->quadelems[1].coef, 3.2);
   cr_assert_eq(nlrow ->quadelems[1].idx2, 0);
   cr_assert_eq(nlrow ->quadelems[1].idx2, 1);
   cr_assert_eq(nlrow ->quadelems[2].coef, 1);
   cr_assert_eq(nlrow ->quadelems[2].idx2, 1);
   cr_assert_eq(nlrow ->quadelems[2].idx2, 1);
   cr_assert_eq(nlrow ->quadelems[3].coef, -4.0);
   cr_assert_eq(nlrow ->quadelems[3].idx2, 3);
   cr_assert_eq(nlrow ->quadelems[3].idx2, 4);
   cr_assert_eq(nlrow ->quadelems[4].coef, 5);
   cr_assert_eq(nlrow ->quadelems[4].idx2, 4);
   cr_assert_eq(nlrow ->quadelems[4].idx2, 4);
   cr_assert_eq(nlrow ->quadvars[0], x1);
   cr_assert_eq(nlrow ->quadvars[1], x2);
   cr_assert_eq(nlrow ->quadvars[2], x4);
   cr_assert_eq(nlrow ->quadvars[3], x5);

   /* check non-quadratic part */
   cr_assert(nlrow->exprtree != NULL);
   cr_assert_eq(SCIPexprtreeGetNVars(nlrow->exprtree), 4);
   cr_assert_eq(SCIPexprtreeGetVars(nlrow->exprtree)[0], x1);
   cr_assert_eq(SCIPexprtreeGetVars(nlrow->exprtree)[1], x2);
   cr_assert_eq(SCIPexprtreeGetVars(nlrow->exprtree)[2], x3);
   cr_assert_eq(SCIPexprtreeGetVars(nlrow->exprtree)[3], x5);

   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}