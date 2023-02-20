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
 * @brief  tests computation of Hessian direction
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

#define EXPECTFEQ(a,b) cr_expect_float_eq(a, b, 1e-6, "%s = %g != %g (dif %g)", #a, a, b, ABS(a-b))

static SCIP* scip;
static SCIP_SOL* sol;
static SCIP_SOL* dir;
static SCIP_VAR* x;
static SCIP_VAR* y;

/* gives bardot expr that belongs to var */
static
SCIP_Real getPartialDiffGradientDir(
   SCIP_EXPR* expr,
   SCIP_VAR*           var
   )
{
   SCIP_EXPRITER* it;
   SCIP_Real deriv = 0.0;

   SCIP_CALL_ABORT( SCIPcreateExpriter(scip, &it) );

   /* we need to sum-up the derivative value from all expr that represent variable var; but only visit each expr once */
   for( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      if( SCIPisExprVar(scip, expr) && SCIPgetVarExprVar(expr) == var )
         deriv += SCIPexprGetBardot(expr);

   SCIPfreeExpriter(&it);

   return deriv;
}


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
   SCIP_Real expected;
   SCIP_EXPR* expr;
   const char* input = "sin(<x>[C]^2 * <y>[C]) + <x>[C]^2";
   SCIP_Real xv = 2.3;
   SCIP_Real yv = -4.0;

   SCIP_CALL( SCIPparseExpr(scip, &expr, input, NULL, NULL, NULL) );

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, yv) );

   /* set direction values */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 1.0) );

   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   expected = 2.0 + cos(xv * xv * yv) * (2 * xv + 2 * yv) - sin(xv * xv * yv) * (4 * xv * xv * yv * yv + 2 * xv * xv * xv * yv);
   EXPECTFEQ(getPartialDiffGradientDir(expr, x), expected);

   expected = 2.0 * xv * cos(xv * xv * yv) - sin(xv * xv * yv) * ( xv * xv * xv * xv + 2 * xv * xv * xv * yv);
   EXPECTFEQ(getPartialDiffGradientDir(expr, y), expected);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(hess, hessian2, .init = setup, .fini = teardown)
{
   SCIP_Real xv = 2.3;
   SCIP_Real yv = -4.0;
   SCIP_Real fxx;
   SCIP_Real fxy;
   SCIP_Real fyy;
   SCIP_EXPR* expr;
   const char* input = "(3.0 + 2.8 * <x>[C])^3 * <y>[C] + 1.0";

   SCIP_CALL( SCIPparseExpr(scip, &expr, input, NULL, NULL, NULL) );

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, yv) );

   /* get first column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 0.0) );

   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   fxx = getPartialDiffGradientDir(expr, x);
   fxy = getPartialDiffGradientDir(expr, y);

   /* get second column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 1.0) );

   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   fyy = getPartialDiffGradientDir(expr, y);

   EXPECTFEQ(getPartialDiffGradientDir(expr, x), fxy);

   EXPECTFEQ(fxx, 6*(3+xv*2.8)*2.8*2.8*yv);
   EXPECTFEQ(fxy, 3*(3+xv*2.8)*(3+xv*2.8)*2.8);
   EXPECTFEQ(fyy, 0.0);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(hess, hessian3, .init = setup, .fini = teardown)
{
   SCIP_Real expected;
   SCIP_Real xv = 2.3;
   SCIP_Real fxx;
   SCIP_EXPR* expr;
   const char* input = "<x>[C]^3";

   SCIP_CALL( SCIPparseExpr(scip, &expr, input, NULL, NULL, NULL) );

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );

   /* get first column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );

   SCIP_CALL( SCIPevalExprHessianDir(scip, expr, sol, 0, dir) );

   fxx = getPartialDiffGradientDir(expr, x);
   expected = 3 * 2 * xv;

   EXPECTFEQ(fxx, expected);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}
