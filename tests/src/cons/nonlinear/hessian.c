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

/**@file   hessian.c
 * @brief  tests computation of Hessian direction for expression in constraint
 *
 * This is similar to expr/hessian.c, but canonicalizes expressions and uses
 * cons_nonlinear's function to access Hessian direction
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"

#define EXPECTFEQ(a,b) cr_expect_float_eq(a, b, 1e-6, "%s = %g != %g (dif %g)", #a, a, b, ABS(a-b))


static SCIP* scip;
static SCIP_SOL* sol;
static SCIP_SOL* dir;
static SCIP_VAR* x;
static SCIP_VAR* y;

static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", -1) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_PRESOLVING, FALSE) );
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPcreateSol(scip, &dir, NULL) );

}

static
void teardown(void)
{
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPfreeSol(scip, &dir) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   assert(BMSgetMemoryUsed() == 0);
}

Test(hess, hessian1, .init = setup, .fini = teardown)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_Real expected;
   SCIP_Real xv = 2.3;
   SCIP_Real yv = -4.0;
   SCIP_VAR* tx;
   SCIP_VAR* ty;
   int dummy;
   SCIP_EXPR* expr;
   const char* input = "[nonlinear] <test>: sin(<t_x>[C]^2 * <t_y>[C]) + <t_x>[C]^2 <= 2;";

   SCIP_CALL( SCIPgetTransformedVar(scip, x, &tx) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y, &ty) );

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);

   SCIP_CALL( canonicalizeConstraints(scip, SCIPfindConshdlr(scip, "nonlinear"), &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeas, &dummy, &dummy, &dummy) );
   assert(!infeas);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, tx, xv) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, ty, yv) );

   /* set direction values */
   SCIP_CALL( SCIPsetSolVal(scip, dir, tx, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, ty, 1.0) );

   expr = SCIPgetExprNonlinear(cons);

   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   expected = 2.0 + cos(xv * xv * yv) * (2 * xv + 2 * yv) - sin(xv * xv * yv) * (4 * xv * xv * yv * yv + 2 * xv * xv * xv * yv);
   EXPECTFEQ(SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, tx), expected);

   expected = 2.0 * xv * cos(xv * xv * yv) - sin(xv * xv * yv) * ( xv * xv * xv * xv + 2 * xv * xv * xv * yv);
   EXPECTFEQ(SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, ty), expected);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(hess, hessian2, .init = setup, .fini = teardown)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_Real xv = 2.3;
   SCIP_Real yv = -4.0;
   SCIP_Real fxx;
   SCIP_Real fxy;
   SCIP_Real fyy;
   SCIP_VAR* tx;
   SCIP_VAR* ty;
   int dummy;
   SCIP_EXPR* expr;
   const char* input = "[nonlinear] <test>: (3.0 +  2.8 * <t_x>[C])^3 * <t_y>[C] + 1.0 <= 2;";

   SCIP_CALL( SCIPgetTransformedVar(scip, x, &tx) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y, &ty) );

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);

   SCIP_CALL( canonicalizeConstraints(scip, SCIPfindConshdlr(scip, "nonlinear"), &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeas, &dummy, &dummy, &dummy) );
   assert(!infeas);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, yv) );

   /* get first column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 0.0) );

   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   fxx = SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, tx);
   fxy = SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, ty);

   /* get second column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 1.0) );

   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   fyy = SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, ty);

   EXPECTFEQ(SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, tx), fxy);

   EXPECTFEQ(fxx, 6*(3+xv*2.8)*2.8*2.8*yv);
   EXPECTFEQ(fxy, 3*(3+xv*2.8)*(3+xv*2.8)*2.8);
   EXPECTFEQ(fyy, 0.0);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

Test(hess, hessian3, .init = setup, .fini = teardown)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Bool infeas;
   SCIP_Real expected;
   SCIP_Real xv = 2.3;
   SCIP_Real fxx;
   SCIP_VAR* tx;
   SCIP_VAR* ty;
   int dummy;
   SCIP_EXPR* expr;
   const char* input = "[nonlinear] <test>:  <t_x>[C]^3 <= 2.0;";

   SCIP_CALL( SCIPgetTransformedVar(scip, x, &tx) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y, &ty) );

   SCIP_CALL( SCIPparseCons(scip, &cons, input, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);

   SCIP_CALL( canonicalizeConstraints(scip, SCIPfindConshdlr(scip, "nonlinear"), &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeas, &dummy, &dummy, &dummy) );
   assert(!infeas);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, tx, xv) );

   /* get first column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, tx, 1.0) );

   expr = SCIPgetExprNonlinear(cons);
   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   fxx = SCIPgetExprPartialDiffGradientDirNonlinear(scip, expr, tx);
   expected = 3 * 2 * xv;

   EXPECTFEQ(fxx, expected);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
