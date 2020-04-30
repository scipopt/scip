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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   conscheck.c
 * @brief  tests check
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_expr.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_SOL* sol;

/* creates scip, problem, includes expression constraint handler, creates  and adds variables */
static
void setup(void)
{
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include cons_expr: this adds the operator handlers */
   SCIP_CALL( SCIPincludeConshdlrExpr(scip) );

   /* get expr conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 5.0, 2.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 5.0, 0.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* release solution */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(conshdlr, conscheck, .init = setup, .fini = teardown,
   .description = "test feasibility check of the cons_expr constraint handler."
   )
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );


   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "an infeasible solution has been accepted");

   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success, "a feasible solution has been declined");

   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "undefined solution has been accepted");

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

static
SCIP_Real getActivity(void)
{
   SCIP_Real xval;
   SCIP_Real yval;
   SCIP_Real zval;

   xval = SCIPgetSolVal(scip, sol, x);
   yval = SCIPgetSolVal(scip, sol, y);
   zval = SCIPgetSolVal(scip, sol, z);

   if( zval == 0.0 )
      return SCIP_INVALID;

   return 1.1*xval*yval/zval + 3.2*xval*xval*pow(yval,-5)*zval + 0.5*pow(zval,3);
}

static
SCIP_Real getGradNorm(void)
{
   SCIP_Real xval;
   SCIP_Real yval;
   SCIP_Real zval;

   xval = SCIPgetSolVal(scip, sol, x);
   yval = SCIPgetSolVal(scip, sol, y);
   zval = SCIPgetSolVal(scip, sol, z);

   if( zval == 0.0 )
      return SCIP_INVALID;

   /* Gradient:
    *  dx:  1.1*y/z + 3.2*2*x*<y>^(-5)*<z>
    *  dy:  1.1*x/z + 3.2*(-5)*x^2*y^(-6)*z
    *  dz:  -1.1*x*y/z^2 + 3.2*x^2*y^(-5) + 1.5*z^2
    */
   return sqrt(
      SQR(1.1*yval/zval + 6.4*xval*pow(yval,-5)*zval) +
      SQR(1.1*xval/zval - 5*3.2*xval*xval*pow(yval,-6)*zval) +
      SQR(-1.1*xval*yval/(zval*zval) + 3.2*xval*xval*pow(yval,-5) + 1.5*zval*zval) );
}

Test(conshdlr, consrelviol_n, .init = setup, .fini = teardown,
   .description = "test relative constraint violation of the cons_expr constraint handler for violscale=n."
   )
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* no scaling, relative should be equal to absolute violation */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/expr/violscale", 'n') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "an infeasible solution has been accepted");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_float_eq(absviol, activity-2.0, SCIPepsilon(scip), "activity: %g, rhs: 2.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_eq(relviol, absviol);


   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success, "a feasible solution has been declined");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, 0.0);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_eq(relviol, absviol);


   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "undefined solution has been accepted");

   activity = getActivity();
   cr_expect_eq(activity, SCIP_INVALID);
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, SCIPinfinity(scip));

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_eq(relviol, absviol);

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

Test(conshdlr, consrelviol_a, .init = setup, .fini = teardown,
   .description = "test relative constraint violation of the cons_expr constraint handler for violscale=a."
   )
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* scale by max(activity, rhs) */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/expr/violscale", 'a') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "an infeasible solution has been accepted");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_float_eq(absviol, activity-2.0, SCIPepsilon(scip), "activity: %g, rhs: 2.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_float_eq(relviol, absviol / MAX(activity, 2.0), SCIPepsilon(scip));


   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success, "a feasible solution has been declined");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, 0.0);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_float_eq(relviol, absviol / MAX(activity, 2.0), SCIPepsilon(scip));


   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "undefined solution has been accepted");

   activity = getActivity();
   cr_expect_eq(activity, SCIP_INVALID);
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, SCIPinfinity(scip));

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_eq(relviol, absviol);

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

Test(conshdlr, consrelviol_g, .init = setup, .fini = teardown,
   .description = "test relative constraint violation of the cons_expr constraint handler for violscale=g."
   )
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real gradnorm;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[expr] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* scale by norm of gradient */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/expr/violscale", 'g') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "an infeasible solution has been accepted");

   activity = getActivity();
   gradnorm = getGradNorm();
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_float_eq(absviol, activity-2.0, SCIPepsilon(scip), "activity: %g, rhs: 2.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_float_eq(relviol, absviol / MAX(gradnorm, 1.0), SCIPepsilon(scip));


   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect(success, "a feasible solution has been declined");

   activity = getActivity();
   gradnorm = getGradNorm();
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, 0.0);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_float_eq(relviol, absviol / MAX(gradnorm, 1.0), SCIPepsilon(scip));


   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "undefined solution has been accepted");

   activity = getActivity();
   cr_expect_eq(activity, SCIP_INVALID);
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, SCIPinfinity(scip));

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_eq(relviol, absviol);

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}

Test(conshdlr, consrelviol_g2, .init = setup, .fini = teardown,
   .description = "test relative constraint violation of the cons_expr constraint handler for violscale=g."
   )
{
   SCIP_CONS* consexpr;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real gradnorm;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[expr] <test>: <x>^0.5*<y>^0.5 >= 1;";
   /* dx: 0.5*sqrt(y/x)  dy: 0.5*sqrt(x/y) */

   /* scale by norm of gradient */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/expr/violscale", 'g') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &consexpr, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   cr_assert(success);
   SCIP_CALL( SCIPaddCons(scip, consexpr) );

   /* create an infeasible solution that can be evaluated and scaled */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.5) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "an infeasible solution has been accepted");

   activity = sqrt(0.25);
   gradnorm = sqrt(SQR(0.5)+SQR(0.5));
   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_float_eq(absviol, 1.0-activity, SCIPepsilon(scip), "activity: %g, lhs: 1.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   cr_expect_float_eq(relviol, absviol / MAX(gradnorm, 1.0), SCIPepsilon(scip));


   /* create an infeasible solution that can be evaluated but gradient doesn't exist */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.5) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   cr_expect_not(success, "an infeasible solution has been accepted");

   SCIP_CALL( SCIPgetAbsViolationConsExpr(scip, consexpr, sol, &absviol) );
   cr_expect_eq(absviol, 1.0);

   SCIP_CALL( SCIPgetRelViolationConsExpr(scip, consexpr, sol, &relviol) );
   /* if gradient cannot be evaluated, then no scaling should happen */
   cr_expect_eq(relviol, absviol);


   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &consexpr) );
}
