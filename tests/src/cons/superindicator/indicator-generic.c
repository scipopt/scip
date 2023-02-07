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

/**@file   indicator-generic.c
 * @brief  unit test for the generic indicator constraint interface
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

#include <math.h>
#include <limits.h>
#include <float.h>

#undef SCIP_CALL
#define SCIP_CALL(x)   do                                                                 \
                       {                                                                  \
                          SCIP_RETCODE _restat_;                                          \
                          if( (_restat_ = (x)) != SCIP_OKAY )                             \
                          {                                                               \
                             cr_assert(FALSE, "Error <%d> in function call\n", _restat_); \
                          }                                                               \
                       }                                                                  \
                       while( FALSE )

/* hardcoded eps equality */
#undef EPSEQ
#define EPSEQ(x, y) (fabs((x) - (y)) <= 0.000005)

static SCIP* scip;

static SCIP_VAR* x1;
static SCIP_VAR* x2;
static SCIP_VAR* z1;
static SCIP_VAR* z2;
static SCIP_CONS* lc1;
static SCIP_CONS* lc2;
static SCIP_CONS* ic1;
static SCIP_CONS* ic2;
static SCIP_SOL* sol;
static SCIP_Real xval;
static SCIP_Real zval;

/* creates scip, problem */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "indicatorProblem") );
}

static
void teardown(void)
{

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

TestSuite(indicator, .init = setup, .fini = teardown);

Test(indicator, feasiblezero1, .init = setup, .fini = teardown,
   .description = "MIP with two indicator constraints"
   )
{
   /*
      Unit test inspired by MathOptInterface.jl
      minimize
        obj: -2 * x1 - 3 * x2
      Subject to
        lc1: x1 + x2 <= 10
        lc2: z1 + z2 >= 1
        ic1: z1 == 1 => x2 <= 8
        ic2: z2 == 1 => x2 + x1/5 <= 9
      Binary
        z1 z2
      Solution:
        x: 1.25, 8.75
        z: 0, 1
        value: -28.75
   */

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z1, "z1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z2, "z2", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, z1) );
   SCIP_CALL( SCIPaddVar(scip, z2) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc1, "lc1", 0, NULL, NULL, -SCIPinfinity(scip), 10.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc1) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc2, "lc2", 0, NULL, NULL, 1.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc2) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic1, "ic1", z1, 0, NULL, NULL, 8.0,
       TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic1, x2, 1.0) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic2, "ic2", z2, 0, NULL, NULL, 9.0,
       TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x2, 1.0) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x1, 1.0 / 5.0) );

   SCIP_CALL( SCIPaddCons(scip, ic1) );
   SCIP_CALL( SCIPaddCons(scip, ic2) );

   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );

   SCIP_CALL( SCIPsolve(scip) );

   sol = SCIPgetBestSol( scip );
   cr_assert_not_null(sol);

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);

   xval = SCIPgetSolVal(scip, sol, x1);
   cr_assert_eq(xval, 1.25);
   xval = SCIPgetSolVal(scip, sol, x2);
   cr_assert_eq(xval, 8.75);

   zval = SCIPgetSolVal(scip, sol, z1);
   cr_assert_eq(zval, 0.0);
   zval = SCIPgetSolVal(scip, sol, z2);
   cr_assert_eq(zval, 1.0);

   SCIP_CALL( SCIPreleaseCons(scip, &ic1) );
   SCIP_CALL( SCIPreleaseCons(scip, &ic2) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc1) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z2) );
}

Test(indicator, feasibleone1, .init = setup, .fini = teardown,
   .description = "MIP with two indicator constraints, one negated variable"
   )
{
   /*
      Unit test inspired by MathOptInterface.jl
      minimize
        obj: -2 * x1 - 3 * x2
      Subject to
        lc1: x1 + x2 <= 10
        lc2: z1 + z2 >= 1
        ic1: z1 == 1 => x2 <= 8
        ic2: z2 == 1 => x2 + x1/5 <= 9
      Binary
        z1 z2 z1neg
      Solution:
        x: 1.25, 8.75
        z: 0, 1, 0
        value: -28.75
   */

   SCIP_VAR* z1neg;
   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z1neg, "z1neg", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z2, "z2", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, z1neg) );
   SCIP_CALL( SCIPaddVar(scip, z2) );

   SCIP_CALL( SCIPgetNegatedVar(scip, z1neg, &z1) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc1, "lc1", 0, NULL, NULL, -SCIPinfinity(scip), 10.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc1) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc2, "lc2", 0, NULL, NULL, 1.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc2) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic1, "ic1", z1, 0, NULL, NULL, 8.0,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarIndicator(scip, ic1, x2, 1.0) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic2, "ic2", z2, 0, NULL, NULL, 9.0,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x2, 1.0) );
   SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x1, 1.0 / 5.0) );

   SCIP_CALL( SCIPaddCons(scip, ic1) );
   SCIP_CALL( SCIPaddCons(scip, ic2) );

   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );

   SCIP_CALL( SCIPsolve(scip) );

   sol = SCIPgetBestSol( scip );
   cr_assert_not_null(sol);

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);

   xval = SCIPgetSolVal(scip, sol, x1);
   cr_assert(EPSEQ(xval, 1.249999));
   xval = SCIPgetSolVal(scip, sol, x2);
   cr_assert(EPSEQ(xval, 8.75));

   zval = SCIPgetSolVal(scip, sol, z1);
   cr_assert(EPSEQ(zval, 0.0));
   zval = SCIPgetSolVal(scip, sol, z2);
   cr_assert(EPSEQ(zval, 1.0));

   SCIP_CALL( SCIPreleaseCons(scip, &ic1) );
   SCIP_CALL( SCIPreleaseCons(scip, &ic2) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc1) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z1neg) );
   SCIP_CALL( SCIPreleaseVar(scip, &z2) );
}

Test(indicator, feasibleone2, .init = setup, .fini = teardown,
   .description = "MIP with two indicator constraints, one negated variable, actual var activated on 0"
   )
{
   /*
      Unit test inspired by MathOptInterface.jl
      minimize
        obj: -2 * x1 - 3 * x2
      Subject to
        lc1: x1 + x2 <= 10
        lc2: z1 - z2 <= 0
        ic1: z1 == 0 => x2 <= 8
        ic2: z2 == 1 => x2 + x1 / 5 <= 9
      Binary
        z1 z2
      Solution:
        x: 1.25, 8.75
        z: 1, 1
        value: -28.75
   */

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z1, "z1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z2, "z2", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, z1) );
   SCIP_CALL( SCIPaddVar(scip, z2) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc1, "lc1", 0, NULL, NULL, -SCIPinfinity(scip), 10.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc1) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc2, "lc2", 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z2, -1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc2) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic1, "ic1", z1, 0, NULL, NULL, 8.0, FALSE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarIndicator(scip, ic1, x2, 1.0) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic2, "ic2", z2, 0, NULL, NULL, 9.0, TRUE,
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x2, 1.0) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x1, 1.0 / 5.0) );

   SCIP_CALL( SCIPaddCons(scip, ic1) );
   SCIP_CALL( SCIPaddCons(scip, ic2) );

   SCIP_CALL( SCIPsolve(scip) );

   sol = SCIPgetBestSol( scip );
   cr_assert_not_null(sol);

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);

   xval = SCIPgetSolVal(scip, sol, x1);
   cr_assert(EPSEQ(xval, 1.25));
   xval = SCIPgetSolVal(scip, sol, x2);
   cr_assert(EPSEQ(xval, 8.75));

   zval = SCIPgetSolVal(scip, sol, z1);
   cr_assert(EPSEQ(zval, 1.0));

   SCIP_CALL( SCIPreleaseCons(scip, &ic1) );
   SCIP_CALL( SCIPreleaseCons(scip, &ic2) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc1) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z2) );
}

Test(indicator, flippedindicator1, .init = setup, .fini = teardown,
   .description = "Indicator constraint with flipped inequality"
   )
{
   /*
      Unit test inspired by MathOptInterface.jl
      minimize
        obj: -2 * x1 - 3 * x2
      Subject to
        lc1: x1 + x2 <= 10
        lc2: z1 + z2 >= 1
        ic1: z1 == 1 => -1.0 * x2 >= -8
        ic2: z2 == 1 => x2 + x1/5 <= 9
      Binary
        z1 z2
      Solution:
        x: 1.25, 8.75
        z: 0, 1
        value: -28.75
   */

   SCIP_CALL( SCIPcreateVarBasic(scip, &x1, "x1", -SCIPinfinity(scip), SCIPinfinity(scip), -2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &x2, "x2", -SCIPinfinity(scip), SCIPinfinity(scip), -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z1, "z1", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z2, "z2", 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

   SCIP_CALL( SCIPaddVar(scip, x1) );
   SCIP_CALL( SCIPaddVar(scip, x2) );
   SCIP_CALL( SCIPaddVar(scip, z1) );
   SCIP_CALL( SCIPaddVar(scip, z2) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc1, "lc1", 0, NULL, NULL, -SCIPinfinity(scip), 10.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc1, x2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc1) );

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &lc2, "lc2", 0, NULL, NULL, 1.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z1, 1.0) );
   SCIP_CALL( SCIPaddCoefLinear(scip, lc2, z2, 1.0) );

   SCIP_CALL( SCIPaddCons(scip, lc2) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic1, "ic1", z1, 0, NULL, NULL, -8.0,
       TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic1, x2, -1.0) );

   SCIP_CALL( SCIPcreateConsIndicatorGeneric(scip, &ic2, "ic2", z2, 0, NULL, NULL, 9.0,
       TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x2, 1.0) );
    SCIP_CALL( SCIPaddVarIndicator(scip, ic2, x1, 1.0 / 5.0) );

   SCIP_CALL( SCIPaddCons(scip, ic1) );
   SCIP_CALL( SCIPaddCons(scip, ic2) );

   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) );

   SCIP_CALL( SCIPsolve(scip) );

   sol = SCIPgetBestSol( scip );
   cr_assert_not_null(sol);

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);

   xval = SCIPgetSolVal(scip, sol, x1);
   cr_assert_eq(xval, 1.25);
   xval = SCIPgetSolVal(scip, sol, x2);
   cr_assert_eq(xval, 8.75);

   zval = SCIPgetSolVal(scip, sol, z1);
   cr_assert_eq(zval, 0.0);
   zval = SCIPgetSolVal(scip, sol, z2);
   cr_assert_eq(zval, 1.0);

   SCIP_CALL( SCIPreleaseCons(scip, &ic1) );
   SCIP_CALL( SCIPreleaseCons(scip, &ic2) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc1) );
   SCIP_CALL( SCIPreleaseCons(scip, &lc2) );
   SCIP_CALL( SCIPreleaseVar(scip, &x1) );
   SCIP_CALL( SCIPreleaseVar(scip, &x2) );
   SCIP_CALL( SCIPreleaseVar(scip, &z1) );
   SCIP_CALL( SCIPreleaseVar(scip, &z2) );
}
