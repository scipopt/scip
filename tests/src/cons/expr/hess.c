/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   hess.c
 * @brief  tests computation of hessian
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.c"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_pow.h"
#include "include/scip_test.h"

#define EXPECTFEQ(a,b) cr_expect_float_eq(a, b, 1e-6, "%s = %g != %g (dif %g)", #a, a, b, ABS(a-b))


static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_SOL* sol;
static SCIP_SOL* dir;
static SCIP_VAR* x;
static SCIP_VAR* y;

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


   SCIP_CALL( SCIPgetTransformedVar(scip, x, &tx) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y, &ty) );

   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "[expr] <test>: sin(<t_x>[C]^2 * <t_y>[C]) + <t_x>[C]^2 <= 2;";

   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeas, &dummy, &dummy, &dummy) );
   assert(!infeas);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, yv) );

   /* set direction values */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 1.0) );

   expr = SCIPgetExprConsExpr(scip, cons);
   SCIP_CALL( SCIPcomputeConsExprHessianDir(scip, conshdlr, cons, sol, dir, 0) );

   expected = 2.0 + cos(xv * xv * yv) * (2 * xv + 2 * yv) - sin(xv * xv * yv) * (4 * xv * xv * yv * yv + 2 * xv * xv * xv * yv);
   EXPECTFEQ(SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, tx), expected);

   expected =  2 * xv * cos(xv * xv * yv) - sin(xv * xv * yv) * ( xv * xv * xv * xv + 2 * xv * xv * xv * yv);
   EXPECTFEQ(SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, ty), expected);


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

   SCIP_CALL( SCIPgetTransformedVar(scip, x, &tx) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y, &ty) );

   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "[expr] <test>: (3.0 +  2.8 * <t_x>[C])^3 * <t_y>[C] + 1.0 <= 2;";

   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeas, &dummy, &dummy, &dummy) );
   assert(!infeas);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, yv) );

   expr = SCIPgetExprConsExpr(scip, cons);

   /* get first column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 0.0) );

   SCIP_CALL( SCIPcomputeConsExprHessianDir(scip, conshdlr, cons, sol, dir, 0) );

   fxx = SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, tx);
   fxy = SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, ty);

   /* get second column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, dir, y, 1.0) );

   SCIP_CALL( SCIPcomputeConsExprHessianDir(scip, conshdlr, cons, sol, dir, 0) );

   fyy = SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, ty);

   EXPECTFEQ(SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, tx), fxy);

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

   SCIP_CALL( SCIPgetTransformedVar(scip, x, &tx) );
   SCIP_CALL( SCIPgetTransformedVar(scip, y, &ty) );

   SCIP_CONSEXPR_EXPR* expr;
   const char* input = "[expr] <test>:  <t_x>[C]^3 <= 2.0;";

   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   assert(success);

   SCIP_CALL( canonicalizeConstraints(scip, conshdlr, &cons, 1, SCIP_PRESOLTIMING_ALWAYS, &infeas, &dummy, &dummy, &dummy) );
   assert(!infeas);

   /* set solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, xv) );

   expr = SCIPgetExprConsExpr(scip, cons);

   /* get first column of hessian */
   SCIP_CALL( SCIPsetSolVal(scip, dir, x, 1.0) );

   SCIP_CALL( SCIPcomputeConsExprHessianDir(scip, conshdlr, cons, sol, dir, 0) );

   fxx = SCIPgetConsExprExprPartialDiffGradientDir(scip, conshdlr, expr, tx);
   expected = 3 * 2 * xv;

   EXPECTFEQ(fxx, expected);

   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}
