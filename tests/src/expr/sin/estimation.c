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

/**@file   estimation.c
 * @brief  tests separation and estimation of sin()
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_trig.c"
#include "../estimation.h"

static SCIP_Real coefs[SCIP_EXPR_MAXINITESTIMATES];
static SCIP_Real* coefsp[SCIP_EXPR_MAXINITESTIMATES];
static SCIP_Real constants[SCIP_EXPR_MAXINITESTIMATES];
static int nreturned;

static
void sine_setup(void)
{
   int i;

   setup();

   for( i = 0; i < SCIP_EXPR_MAXINITESTIMATES; ++i )
      coefsp[i] = &coefs[i];
}

/* note: for the tests, computeInitialCutsTrig() does:
 * 1. secant -> done if success
 * 2. left mid tangent (left secant) -> if fails tries left tangent
 * 3. right mid tangent (right secant) -> if fails tries right tangent
 * when underestimating.
 * when overestimating 2. and 3. are swapped
 */

TestSuite(separation, .init = sine_setup, .fini = teardown);

/* tests for interval [-1,5] */
Test(separation, sine_x,
   .description = "test separation for a sine expression in large range"
   )
{
   SCIP_EXPR* expr;
   SCIP_Real newtonpoint;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real lincoef;
   SCIP_Real linconst;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreateExprSin(scip, &expr, xexpr, NULL, NULL) );
   childlb = SCIPvarGetLbLocal(x);
   childub = SCIPvarGetUbLocal(x);

   /* sin(x) between -1 and 5 looks like
    *
    *   |                     .... ...
    *   |                   ..        ..
    *   |                 ..            ..
    *   |                /                \
    *   |               /                  \
    *   |              /                    \
    *   |             /                      \
    *   |            /                        \
    *   | ----------/--------------------------\-----------------
    *   |          /                            \
    *   |         /                              ..
    *   |        /                                 \
    *   |      ..                                   \
    *   |     /                                      \
    *   |    /                                        \
    *   |   /                                          ..
    *   | ..                                             \
    *   |                                                 ....../
    *     -1                     2                          5
    *
    * when overestimating we should get a right secant and a left secant, in that order
    * when underestimating we should get a left secant and a right tangent, in that order
    */

   /*
    * test initial overestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, FALSE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2, "expected %d, got %d\n", 2, nreturned);

   /* check right secant */
   newtonpoint = 2.2544608804;
   EXPECTFEQ( coefs[0], cos(newtonpoint) );
   EXPECTFEQ( constants[0], sin(5) - 5 * cos(newtonpoint) );

   /* check lmidtangent */
   newtonpoint = 0.4936608602;
   EXPECTFEQ( coefs[1], cos(newtonpoint) );
   EXPECTFEQ( constants[1], cos(newtonpoint) + sin(-1) );


   /*
    * test initial underestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, TRUE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2, "expected %d, got %d\n", 2, nreturned);

   /* check left secant */
   newtonpoint = 4.6845658560;
   EXPECTFEQ( coefs[0], cos(newtonpoint) );
   EXPECTFEQ( constants[0], cos(newtonpoint) + sin(-1) );

   /* check right tangent */
   EXPECTFEQ( coefs[1], cos(5) );
   EXPECTFEQ( constants[1],  sin(5) - 5 * cos(5) );


   /*
    * test overestimation
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, 1.5, childlb, childub, FALSE);
   cr_expect(success);
   EXPECTFEQ(linconst, -1.5 * cos(1.5) + sin(1.5));
   EXPECTFEQ(lincoef, cos(1.5));

   /*
    * test underestimation
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, 4.8, childlb, childub, TRUE);
   cr_expect(success);
   EXPECTFEQ(linconst, -4.8 * cos(4.8) + sin(4.8));
   EXPECTFEQ(lincoef, cos(4.8));

   /*
    * test point where solution tangent is not feasible
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, 4.0, childlb, childub, TRUE);
   cr_expect(success);

   /* check left tangent */
   newtonpoint = 4.6845658560;
   cr_expect_float_eq(linconst, cos(newtonpoint) + sin(-1), 1e-10);
   cr_expect_float_eq(lincoef, cos(newtonpoint), 1e-10);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* tests for interval [-6,-3] */
Test(separation, sine_y,
   .description = "test separation for a sine expression in mid size range"
)
{
   SCIP_EXPR* expr;
   SCIP_Real newtonpoint;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real lincoef;
   SCIP_Real linconst;
   SCIP_Bool success;


   SCIP_CALL( SCIPcreateExprSin(scip, &expr, yexpr, NULL, NULL) );
   childlb = SCIPvarGetLbLocal(y);
   childub = SCIPvarGetUbLocal(y);

   /* sin(x) between -6 and -3 looks like
    *
    *   |                   ...... ......
    *   |                ...             ..
    *   |              ..                  ..
    *   |            ..                      ..
    *   |          ..                          ..
    *   |         /                              \
    *   |        /                                ..
    *   |      ..                                   \
    *   |     /                                      \
    *   |    /                                        \
    *   |   /                                          ..
    *   |  .                                             \
    *   |                                                 \
    *   |                                                  \
    *   |                                                   \
    *   | ---------------------------------------------------\---
    *   |                                                     \
    *     -6                     -4.5                       -3
    *
    * when overestimating we should get a right secant (sine is convex at -3) and a left tangent, in that order
    * when underestimating we should get a secant
    */

   /*
    * test initial overestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, FALSE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2, "expected %d, got %d\n", 2, nreturned);

   /* check right secant */
   newtonpoint = -3.2123712333;
   EXPECTFEQ( coefs[0], cos(newtonpoint) );
   EXPECTFEQ( constants[0], 3 * cos(newtonpoint) + sin(-3) );

   /* check left tangent */
   EXPECTFEQ( coefs[1], cos(-6) );
   EXPECTFEQ( constants[1], sin(-6) + 6 * cos(-6) );


   /*
    * test initial underestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, TRUE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 1, "expected %d, got %d\n", 1, nreturned);

   /* check secant */
   EXPECTFEQ( coefs[0], (sin(-3) - sin(-6)) / 3.0 );
   EXPECTFEQ( constants[0], 2.0 * sin(-3) - sin(-6) );

   /*
    * test overestimation
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, -4.0, childlb, childub, FALSE);
   cr_expect(success);
   EXPECTFEQ(linconst, 4 * cos(-4) + sin(-4));
   EXPECTFEQ(lincoef, cos(-4));

   /*
    * test underestimation (not possible in this case)
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, -3.1, childlb, childub, TRUE);
   cr_expect(success);
   EXPECTFEQ(linconst, -sin(-6) + 2.0 * sin(-3));
   EXPECTFEQ(lincoef, (sin(-3) - sin(-6)) / 3.0);

   /*
    * test point where solution tangent is not underestimating
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, -3.2, childlb, childub, FALSE);
   cr_expect(success);

   /* check rmidtangent */
   newtonpoint = -3.2123712333;
   cr_expect_float_eq(linconst, 3 * cos(newtonpoint) + sin(-3), 1e-11);
   cr_expect_float_eq(lincoef, cos(newtonpoint), 1e-11);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

/* tests for interval [1,3] */
Test(separation, sine_z,
   .description = "test separation for a sine expression in short range"
)
{
   SCIP_EXPR* expr;
   SCIP_Real childlb;
   SCIP_Real childub;
   SCIP_Real lincoef;
   SCIP_Real linconst;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreateExprSin(scip, &expr, zexpr, NULL, NULL) );
   childlb = SCIPvarGetLbLocal(z);
   childub = SCIPvarGetUbLocal(z);

   /* sin(x) between 1 and 3 looks like
    *   |         ........ ........
    *   |     ....                 ...
    *   |   ..                        ...
    *   | ..                             ..
    *   |                                  ..
    *   |                                    ..
    *   |                                      ..
    *   |                                        ..
    *   |                                          \
    *   |                                           ..
    *   |                                             \
    *   |                                              ..
    *   |                                                \
    *   |                                                 \
    *   |                                                  ..
    *   |                                                    \
    *   |                                                     \
    *   |                                                      ..
    *   |
    *   | -------------------------------------------------------
    *     1                      2                          3
    *
    * when overestimating we should get a right tangent and a left tangent, in that order
    * when underestimating we should get a secant
    */

   /*
    * test initial overestimation
    */
   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, FALSE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2, "expected %d, got %d\n", 2, nreturned);

   /* check right tangent */
   EXPECTFEQ( coefs[0], cos(3) );
   EXPECTFEQ( constants[0], sin(3) - 3 * cos(3) );

   /* check left tangent */
   EXPECTFEQ( coefs[1], cos(1) );
   EXPECTFEQ( constants[1], sin(1) - cos(1) );

   /*
    * test initial underestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, TRUE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 1, "expected %d, got %d\n", 1, nreturned);

   /* check secant */
   EXPECTFEQ( coefs[0], 0.5 * (sin(3) - sin(1)) );
   EXPECTFEQ( constants[0], 1.5 * sin(1) - 0.5 * sin(3) );

   /*
    * test overestimation
    */

   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, 2.0, childlb, childub, FALSE);
   cr_expect(success);
   EXPECTFEQ(linconst, -2 * cos(2) + sin(2));
   EXPECTFEQ(lincoef, cos(2));

   /* note: solution underestimation doesn't make sense in [1,3] */

   /*
    * test point where solution tangent is not underestimating
    */
   success = computeEstimatorsTrig(scip, expr, &lincoef, &linconst, 2.0, childlb, childub, TRUE);
   cr_expect(success);

   /* check secant */
   EXPECTFEQ(linconst, -0.5 * sin(3) + 1.5 * sin(1));
   EXPECTFEQ(lincoef, 0.5 * (sin(3) - sin(1)));

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}


/* tests for interval [-pi,pi] */
Test(separation, sine_v,
   .description = "test separation for a sine expression in large range"
   )
{
   SCIP_EXPR* expr;
   SCIP_Real childlb;
   SCIP_Real childub;

   SCIP_CALL( SCIPcreateExprSin(scip, &expr, vexpr, NULL, NULL) );
   childlb = SCIPvarGetLbLocal(v);
   childub = SCIPvarGetUbLocal(v);

   /* sin(x) between -pi and pi looks like
    *   |                                       ... ....
    *   |                                     ..        ..
    *   |                                   ..            \
    *   |                                  /               \
    *   |                                 /                 \
    *   |                                /                   ..
    *   |                               /                      \
    *   |                              /                        \
    *   | .---------------------------/--------------------------
    *   |                            /
    *   |  .                        /
    *   |   \                      /
    *   |    ..                   /
    *   |      \                 /
    *   |       \               /
    *   |        \            ..
    *   |         ..        ..
    *   |           ........
    *     -3.14159               0                          3.14159
    *
    * when overestimating we should get a right tangent and a left secant
    * when underestimating we should get left tangent and a right secant
    */

   /*
    * test initial overestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, FALSE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2, "expected %d, got %d\n", 2, nreturned);

   /* check right tangent */
   EXPECTFEQ( coefs[0], -1.0 );
   EXPECTFEQ( constants[0], M_PI );

   /* check left secant */
   EXPECTFEQ( coefs[1], 0.217233628211222 );
   EXPECTFEQ( constants[1], 0.682459570501030 );

   /*
    * test initial underestimation
    */

   SCIP_CALL( computeInitialCutsTrig(scip, expr, childlb, childub, TRUE, coefsp, constants, &nreturned) );
   cr_expect_eq(nreturned, 2, "expected %d, got %d\n", 2, nreturned);

   /* check left tangent */
   EXPECTFEQ( coefs[0], -1.0 );
   EXPECTFEQ( constants[0], -M_PI );

   /* check right secant */
   EXPECTFEQ( coefs[1], 0.217233628211222 );
   EXPECTFEQ( constants[1], -0.682459570501030 );


   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
