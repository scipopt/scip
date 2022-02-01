/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   parse.c
 * @brief  tests nonlinear constraint parsing
 *
 * See also ../../expr/parse.c for expression parsing.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
}

static
void teardown(void)
{

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

TestSuite(parse, .init = setup, .fini = teardown);

Test(parse, constraint_with_spaces)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   const char* input = "[nonlinear] <test1>: <x>[C] / <y>[I] *(5) >= 1;";

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   cr_expect_eq(SCIPgetLhsNonlinear(cons), 1.0);
   cr_expect_eq(SCIPgetRhsNonlinear(cons), SCIPinfinity(scip));
   /* TODO there should be some test that the expression was parsed ok, too */

   /* print constraint */
   SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* release constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(parse, constraint_with_sides)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   const char* input = "[nonlinear] <test2>: 1 <= <x>[C] / <y>[I] *(5) - <x> <= 2;";

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);

   cr_expect_eq(SCIPgetLhsNonlinear(cons), 1.0);
   cr_expect_eq(SCIPgetRhsNonlinear(cons), 2.0);
   /* TODO there should be some test that the expression was parsed ok, too */

   /* print constraint */
   SCIPinfoMessage(scip, NULL, "printing constraint %s after parsing from string:", input);
   SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* release constraint */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
