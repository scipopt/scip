/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   check.c
 * @brief  tests check callback
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;
static SCIP_SOL* sol;

/* creates scip, problem, includes nonlinear constraint handler, creates  and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 5.0, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 5.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
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

   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/** @brief test feasibility check of the nonlinear constraint handler. */
void test_conshdlr_check(void)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   const char* input = "[nonlinear] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   TEST_ASSERT(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );


   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "an infeasible solution has been accepted");

   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, FALSE, FALSE, FALSE, &success) );
   /* soft assert, here and in the tests below: the violations queried afterwards
    * remain meaningful if the solution was declined, and they say why it was
    */
   SOFT_ASSERT(success, "a feasible solution has been declined");

   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, TRUE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "undefined solution has been accepted");

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
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

/** @brief test relative constraint violation of the nonlinear constraint handler for violscale=n. */
void test_conshdlr_relviol_n(void)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[nonlinear] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* no scaling, relative should be equal to absolute violation */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/nonlinear/violscale", 'n') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   TEST_ASSERT(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "an infeasible solution has been accepted");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(absviol, activity-2.0, SCIPepsilon(scip), "activity: %g, rhs: 2.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_EQUAL(relviol, absviol);


   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT(success, "a feasible solution has been declined");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, 0.0);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_EQUAL(relviol, absviol);


   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "undefined solution has been accepted");

   activity = getActivity();
   SOFT_ASSERT_EQUAL(activity, SCIP_INVALID);
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, SCIPinfinity(scip));

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_EQUAL(relviol, absviol);

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/** @brief test relative constraint violation of the nonlinear constraint handler for violscale=a. */
void test_conshdlr_relviol_a(void)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[nonlinear] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* scale by max(activity, rhs) */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/nonlinear/violscale", 'a') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   TEST_ASSERT(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "an infeasible solution has been accepted");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(absviol, activity-2.0, SCIPepsilon(scip), "activity: %g, rhs: 2.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(relviol, absviol / MAX(activity, 2.0), SCIPepsilon(scip));


   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT(success, "a feasible solution has been declined");

   activity = getActivity();
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, 0.0);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(relviol, absviol / MAX(activity, 2.0), SCIPepsilon(scip));


   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "undefined solution has been accepted");

   activity = getActivity();
   SOFT_ASSERT_EQUAL(activity, SCIP_INVALID);
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, SCIPinfinity(scip));

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_EQUAL(relviol, absviol);

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/** @brief test relative constraint violation of the nonlinear constraint handler for violscale=g. */
void test_conshdlr_relviol_g(void)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real gradnorm;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[nonlinear] <test>: 1.1*<x>*<y>/<z> + 3.2*<x>^2*<y>^(-5)*<z> + 0.5*<z>^3 <= 2;";

   /* scale by norm of gradient */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/nonlinear/violscale", 'g') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   TEST_ASSERT(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* create an infeasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 2) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 3) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "an infeasible solution has been accepted");

   activity = getActivity();
   gradnorm = getGradNorm();
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(absviol, activity-2.0, SCIPepsilon(scip), "activity: %g, rhs: 2.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(relviol, absviol / MAX(gradnorm, 1.0), SCIPepsilon(scip));


   /* create a feasible solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 1) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT(success, "a feasible solution has been declined");

   activity = getActivity();
   gradnorm = getGradNorm();
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, 0.0);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(relviol, absviol / MAX(gradnorm, 1.0), SCIPepsilon(scip));


   /* create an undefined solution */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, z, 0) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "undefined solution has been accepted");

   activity = getActivity();
   SOFT_ASSERT_EQUAL(activity, SCIP_INVALID);
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, SCIPinfinity(scip));

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_EQUAL(relviol, absviol);

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/** @brief test relative constraint violation of the nonlinear constraint handler for violscale=g. */
void test_conshdlr_relviol_g2(void)
{
   SCIP_CONS* cons;
   SCIP_Bool success;
   SCIP_Real activity;
   SCIP_Real gradnorm;
   SCIP_Real absviol;
   SCIP_Real relviol;
   const char* input = "[nonlinear] <test>: <x>^0.5*<y>^0.5 >= 1;";
   /* dx: 0.5*sqrt(y/x)  dy: 0.5*sqrt(x/y) */

   /* scale by norm of gradient */
   SCIP_CALL( SCIPsetCharParam(scip, "constraints/nonlinear/violscale", 'g') );

   /* parse constraint */
   success = FALSE;
   SCIP_CALL( SCIPparseCons(scip, &cons, input,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   TEST_ASSERT(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );

   /* create an infeasible solution that can be evaluated and scaled */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.5) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.5) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "an infeasible solution has been accepted");

   activity = sqrt(0.25);
   gradnorm = sqrt(SQR(0.5)+SQR(0.5));
   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(absviol, 1.0-activity, SCIPepsilon(scip), "activity: %g, lhs: 1.0, but absviol: %g", activity, absviol);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   SOFT_ASSERT_DOUBLE_WITHIN(relviol, absviol / MAX(gradnorm, 1.0), SCIPepsilon(scip));


   /* create an infeasible solution that can be evaluated but gradient doesn't exist */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 0.5) );
   SCIP_CALL( SCIPcheckSol(scip, sol, FALSE, TRUE, FALSE, FALSE, FALSE, &success) );
   SOFT_ASSERT_NOT(success, "an infeasible solution has been accepted");

   SCIP_CALL( SCIPgetAbsViolationNonlinear(scip, cons, sol, &absviol) );
   SOFT_ASSERT_EQUAL(absviol, 1.0);

   SCIP_CALL( SCIPgetRelViolationNonlinear(scip, cons, sol, &relviol) );
   /* if gradient cannot be evaluated, then no scaling should happen */
   SOFT_ASSERT_EQUAL(relviol, absviol);


   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
}

/* test violation in expression in nonlinear constraint */
/** @brief Tests expression violation. */
void test_conshdlr_exprviol(void)
{
   SCIP_EXPR* mainexpr;
   SCIP_VAR* auxvar;
   SCIP_Real viol;
   SCIP_Bool violunder;
   SCIP_Bool violover;
   SCIP_Bool success;
   SCIP_Bool cutoff;
   SCIP_CONS* cons;

   SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE);
   SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE);

   /* dualfix will fix y to its lower bound, but 0 isn't in the domain of y^1, so propagation would declare the problem to be infeasible
    * (this indicates that the locks in cons_nonlinear miss to take domain-restrictions into account)
    */
   SCIP_CALL( SCIPchgVarLbGlobal(scip, y, 0.1) );

   SCIP_CALL( SCIPparseCons(scip, &cons, "[nonlinear] <test>: 0.5 * (<x>^2*<y>^(-1)*5^(-4))^2 * (2*<x> + 1)^(-1) >= 0", TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
   TEST_ASSERT(success);
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPfreeSol(scip, &sol) );  /* we need a sol in transformed space */

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* initialize solution values */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, 4.0) );

   TEST_ASSERT_EQUAL(SCIPgetNConss(scip), 1);
   mainexpr = SCIPgetExprNonlinear(SCIPgetConss(scip)[0]);
   SCIP_CALL( SCIPdismantleExpr(scip, NULL, mainexpr) );

   /* construct LP to get auxvar */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   SOFT_ASSERT_NOT(cutoff);

   auxvar = SCIPgetExprAuxVarNonlinear(mainexpr);
   TEST_ASSERT_NOT_NULL(auxvar);

   SCIP_CALL( SCIPevalExpr(scip, mainexpr, sol, 1) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, SCIPexprGetEvalValue(mainexpr) + 5.0) );
   SCIP_CALL( SCIPgetExprAbsOrigViolationNonlinear(scip, mainexpr, sol, 1, &viol, &violunder, &violover) );
   SOFT_ASSERT(!violunder);
   SOFT_ASSERT(violover);
   SOFT_ASSERT_DOUBLE_WITHIN(viol, 5.0, 1e-12, "got violation %g, but expected 5", viol);

   /* because we don't have a positive lock, there should be no violation "on the other side" */
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, SCIPexprGetEvalValue(mainexpr) - 5.0) );
   SCIP_CALL( SCIPgetExprAbsOrigViolationNonlinear(scip, mainexpr, sol, 1, &viol, &violunder, &violover) );
   SOFT_ASSERT(!violunder);
   SOFT_ASSERT(!violover);
   SOFT_ASSERT_EQUAL(viol, 0.0, "got violation %g, but none was expected", viol);

   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, SCIPexprGetEvalValue(mainexpr) + 5.0) );
   SCIP_CALL( SCIPgetExprAbsAuxViolationNonlinear(scip, mainexpr, SCIPexprGetEvalValue(mainexpr)-3.0, sol, &viol, &violunder, &violover) );
   SOFT_ASSERT(!violunder);
   SOFT_ASSERT(violover);
   SOFT_ASSERT_DOUBLE_WITHIN(viol, 8.0, 1e-12, "got violation %g, but expected 8", viol);

   SCIP_CALL( SCIPgetExprRelAuxViolationNonlinear(scip, mainexpr, SCIPexprGetEvalValue(mainexpr)-3.0, sol, &viol, &violunder, &violover) );
   SOFT_ASSERT(!violunder);
   SOFT_ASSERT(violover);
   SOFT_ASSERT_DOUBLE_WITHIN(viol, 8.0 / REALABS(SCIPexprGetEvalValue(mainexpr)-3.0), 1e-12, "got violation %g, but expected %g", viol, viol / REALABS(SCIPexprGetEvalValue(mainexpr)-3.0));

   /* because we don't have a positive lock, there should be no violation "on the other side" */
   SCIP_CALL( SCIPsetSolVal(scip, sol, auxvar, SCIPexprGetEvalValue(mainexpr) - 5.0) );
   SCIP_CALL( SCIPgetExprAbsAuxViolationNonlinear(scip, mainexpr, SCIPexprGetEvalValue(mainexpr)-3.0, sol, &viol, &violunder, &violover) );
   SOFT_ASSERT(!violunder);
   SOFT_ASSERT(!violover);
   SOFT_ASSERT_EQUAL(viol, 0.0, "got violation %g, but none was expected", viol);
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_conshdlr_check);
   RUN_TEST(test_conshdlr_relviol_n);
   RUN_TEST(test_conshdlr_relviol_a);
   RUN_TEST(test_conshdlr_relviol_g);
   RUN_TEST(test_conshdlr_relviol_g2);
   RUN_TEST(test_conshdlr_exprviol);
   return UNITY_END();
}
