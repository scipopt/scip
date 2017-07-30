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

/**@file   copy.c
 * @brief  tests copy of expressions when SCIP copies itself
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

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

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(copy, copy, .init = setup, .fini = teardown)
{
   SCIP* subscip;
   SCIP_CONS* consexpr;
   const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 == 2";
   SCIP_Bool success;
   SCIP_Bool valid;

   /* create constraint from input string */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   /* add constraint to SCIP and release it */
   SCIP_CALL( SCIPaddCons(scip, consexpr) );
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );

   /* goto transformed stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_TRANSFORMED, FALSE) );

   /* copy SCIP */
   valid = FALSE;
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPcopy(scip, subscip, NULL, NULL, "copytest_", TRUE, FALSE, FALSE, &valid) );

   /* check copying was valid and some sanity checks */
   cr_assert(valid);
   cr_assert_eq(SCIPgetNConss(subscip), 1);
   cr_assert_neq(scip, subscip);

   /* check that copied constraint is the same as original (i.e. the transformed one!) */
   cr_redirect_stdout();
   SCIP_CALL( SCIPprintCons(subscip, SCIPgetConss(subscip)[0], NULL) );
   SCIPinfoMessage(subscip, NULL, "\n");

   fflush(stdout);

   cr_assert_stdout_eq_str("  [expr] <test>: 1.1*<t_x>*<t_y>*(<t_z>)^-1+3.2*(<t_x>)^2*(<t_y>)^-5*<t_z>+0.5*(<t_z>)^3 == 2\n");

   /* release the copy of SCIP */
   SCIP_CALL( SCIPfree(&subscip) );
}
