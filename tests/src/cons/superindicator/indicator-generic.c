/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   indicator-generic.c
 * @brief  unit test for the generic indicator constraint interface
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;

/* creates scip, problem */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(indicatorzero, .init = setup, .fini = teardown);

Test(indicatorzero, feasiblezero, .init = setup, .fini = teardown,
   .description = "check upgrade from quadratic constraint with linear binary variables to SOC"
   )
{
   SCIP_VAR* x1;
   SCIP_VAR* x2;
   SCIP_VAR* z1;
   SCIP_VAR* z2;
   SCIP_CONS* lc1;
   SCIP_CONS* lc2;
   SCIP_CONS* ic1;
   SCIP_CONS* ic2;
   /*
      Unit test based on MathOptInterface.jl
      minimize
        obj: -2 * x1 - 3 * x2
      Subject to
        lc1: x1 + x2 <= 10
        ic2: z1 == 0 => x2 <= 8
        ic3: z2 == 1 => x2 + x1/5 <= 9
        lc2: z1 - z2 <= 0
      Binary
        z1 z2
      End
   */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", 0.0, 10.0, -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", 0.0, 10.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z1, "z1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z2, "z2", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, z1) );
   SCIP_CALL( SCIPaddVar(scip, z2) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc1, "lc1", 0, NULL, NULL, -SCIPinfinity(scip), 10.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x2, 1.0) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc2, "lc2", 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z2, 2.0) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(
       scip,
       &ic1,
       "ic1",
       z1,
       0, NULL, NULL,
       8.0,
       FALSE,
       TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(
       scip,
       ic1,
       x2,
       1.0) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(
       scip,
       &ic2,
       "ic2",
       z2,
       0, NULL, NULL,
       9.0,
       TRUE,
       TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(
       scip,
       ic2,
       x2,
       1.0) );
    SCIP_CALL( SCIPaddVarIndicator(
       scip,
       ic2,
       x1,
       1.0 / 5.0) );

   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_not_null(SCIPgetBestSol(scip));

   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z2) );

}
