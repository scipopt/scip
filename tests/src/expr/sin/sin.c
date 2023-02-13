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

/**@file   sin.c
 * @brief  tests expression handler functions of sine expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/expr.h"
#include "scip/expr_trig.h"
#include "scip/expr_var.h"
#include "scip/expr_value.h"
#include "scip/scipdefplugins.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_EXPR* sinexpr;
static SCIP_EXPR* xexpr;
static SCIP_EXPR* yexpr;
static SCIP_RANDNUMGEN* rndgen;


/* creates scip, problem, includes nonlinear constraint handler, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   /* includes expr handlers */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create variable and sine expressions */
   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &sinexpr, xexpr, NULL, NULL) );

   /* create solution */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &rndgen, 1, TRUE) );
}

/* releases variables, frees scip */
static
void teardown(void)
{
   /* free random number generator */
   SCIPfreeRandom(scip, &rndgen);

   /* free allocated memory */
   SCIP_CALL( SCIPfreeSol(scip, &sol) );

   SCIP_CALL( SCIPreleaseExpr(scip, &sinexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/* test suite */
TestSuite(sin, .init = setup, .fini = teardown);

/*
 * TESTS
 */

Test(sin, creation, .description = "Tests the expression creation.")
{
   SCIP_EXPR* expr;

   /* create sine expression */
   SCIP_CALL( SCIPcreateExprSin(scip, &expr, xexpr, NULL, NULL) );

   cr_assert(expr != NULL);
   cr_expect(SCIPexprGetNChildren(expr) == 1);
   cr_expect(SCIPexprGetChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(sin, parse, .description = "Tests the expression parsing.")
{
   SCIP_EXPR* expr;
   const char* input = "sin(<x>[C])";

   /* create sine expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL) );

   cr_assert(expr != NULL);
   cr_expect(SCIPexprGetNChildren(expr) == 1);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(expr)[0]));
   cr_expect(SCIPgetVarExprVar(SCIPexprGetChildren(expr)[0]) == x);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(sin, eval, .description = "Tests the expression evaluation.")
{
   int i;
   SCIP_Real randnum;
   SCIP_Real testvalues[5] = {0.0, M_PI, 2.5*M_PI, -0.5*M_PI, -2.0*M_PI};
   SCIP_Real results[5] = {0.0, 0.0, 1.0, -1.0, 0.0};

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalExpr(scip, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(sinexpr), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i )
   {
      randnum = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );
      SCIP_CALL( SCIPevalExpr(scip, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(sinexpr), sin(randnum)));
   }
}

Test(sin, inteval, .description = "Tests the expression interval evaluation.")
{
   int i;
   SCIP_INTERVAL interval;
   SCIP_Real rndlb[5];
   SCIP_Real rndub[5];
   SCIP_Real rndreslb[5];
   SCIP_Real rndresub[5];

   /* pick 5 special cases with well known results */
   SCIP_Real detlb[5] = {0.0, -0.5*M_PI, -M_PI, 2.0*M_PI, 0.5*M_PI};
   SCIP_Real detub[5] = {2.0*M_PI, 0.5*M_PI, -0.5*M_PI, 3.0*M_PI, 1.5*M_PI};
   SCIP_Real detreslb[5] = {-1.0, -1.0, -1.0, 0.0, -1.0};
   SCIP_Real detresub[5] = {1.0, 1.0, 0.0, 1.0, 1.0};

   /* create 5 random cases within specific bounds that have non-trivial results */
   rndlb[0]    = SCIPrandomGetReal(rndgen, 0.0, 0.5*M_PI);
   rndub[0]    = SCIPrandomGetReal(rndgen, 0.5*M_PI, M_PI);
   rndreslb[0] = MIN(sin(rndlb[0]), sin(rndub[0]));
   rndresub[0] = 1.0;

   rndlb[1]    = SCIPrandomGetReal(rndgen, 0.0, 0.5*M_PI);
   rndub[1]    = SCIPrandomGetReal(rndgen, M_PI, 1.5*M_PI);
   rndreslb[1] = sin(rndub[1]);
   rndresub[1] = 1.0;

   rndlb[2]    = SCIPrandomGetReal(rndgen, 0.5*M_PI, M_PI);
   rndub[2]    = SCIPrandomGetReal(rndgen, M_PI, 1.5*M_PI);
   rndreslb[2] = sin(rndub[2]);
   rndresub[2] = sin(rndlb[2]);

   rndlb[3]    = SCIPrandomGetReal(rndgen, M_PI, 1.5*M_PI);
   rndub[3]    = SCIPrandomGetReal(rndgen, 1.5*M_PI, 2.0*M_PI);
   rndreslb[3] = -1.0;
   rndresub[3] = MAX(sin(rndlb[3]), sin(rndub[3]));

   rndlb[4]    = SCIPrandomGetReal(rndgen, 0.0, 0.5*M_PI);
   rndub[4]    = SCIPrandomGetReal(rndgen, 1.5*M_PI, 2.0*M_PI);
   rndreslb[4] = -1.0;
   rndresub[4] = 1.0;

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, detlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, detub[i]) );
      SCIP_CALL( SCIPevalExpr(scip, sinexpr, sol, 0) );
      SCIP_CALL( SCIPevalExprActivity(scip, sinexpr) );
      interval = SCIPexprGetActivity(sinexpr);

      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), detreslb[i]));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), detresub[i]));
   }

   /* random part */
   for ( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, rndlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, rndub[i]) );
      SCIP_CALL( SCIPevalExpr(scip, sinexpr, sol, 0) );
      SCIP_CALL( SCIPevalExprActivity(scip, sinexpr) );
      interval = SCIPexprGetActivity(sinexpr);

      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetInf(interval), rndreslb[i]));
      cr_expect(SCIPisFeasEQ(scip, SCIPintervalGetSup(interval), rndresub[i]));
   }
}

Test(sin, derivative, .description = "Tests the expression derivation.")
{
   int i;
   SCIP_Real randnum;
   SCIP_Real testvalues[5] = {0.0, M_PI, 2.0*M_PI, 2.5*M_PI, -0.5*M_PI};
   SCIP_Real results[5] = {1.0, -1.0, 1.0, 0.0, 0.0};

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalExprGradient(scip, sinexpr, sol, 0) );

      cr_expect(SCIPisEQ(scip, SCIPexprGetDerivative(xexpr), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i)
   {
      randnum = SCIPrandomGetReal(rndgen, -10.0, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );
      SCIP_CALL( SCIPevalExprGradient(scip, sinexpr, sol, 0) );

      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetDerivative(xexpr), cos(randnum)));
   }
}

Test(sin, hash, .description = "Tests the expression hash.")
{
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_EXPR* expr3;
   unsigned int hashkey1;
   unsigned int hashkey2;
   unsigned int hashkey3;

   SCIP_CALL( SCIPcreateExprSin(scip, &expr1, xexpr, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr2, xexpr, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr3, yexpr, NULL, NULL) );

   SCIP_CALL( SCIPhashExpr(scip, expr1, &hashkey1) );
   SCIP_CALL( SCIPhashExpr(scip, expr2, &hashkey2) );
   SCIP_CALL( SCIPhashExpr(scip, expr3, &hashkey3) );

   cr_expect(hashkey1 != 0);
   cr_expect(hashkey2 != 0);
   cr_expect(hashkey3 != 0);
   cr_expect(hashkey1 == hashkey2);
   cr_expect(hashkey1 != hashkey3);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );
}

Test(sin, simplify, .description = "Tests the expression simplification.")
{
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_EXPR* expr3;
   SCIP_Bool changed;
   SCIP_Bool infeasible;

   /* expr1 = <5.0>, expr2 = sin(<5.0>), expr3 is buffer for simplification */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr1, 5.0, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprSin(scip, &expr2, expr1, NULL, NULL) );
   SCIP_CALL( SCIPsimplifyExpr(scip, expr2, &expr3, &changed, &infeasible, NULL, NULL) );
   SCIP_CALL( SCIPevalExpr(scip, expr2, sol, 0) );

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPisExprValue(scip, expr3));
   cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr2), sin(5.0)));

   SCIP_CALL( SCIPreleaseExpr(scip, &expr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );
}
