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

/**@file   entropy.c
 * @brief  tests expression handler functions of entropy expressions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/expr_var.h"
#include "scip/expr_log.h"
#include "scip/expr_value.h"
#include "scip/expr_entropy.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_SOL* sol;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_EXPR* entropyexpr;
static SCIP_EXPR* prodexpr;          /* xlogx as product expression */
static SCIP_EXPR* negprodexpr;       /* -xlogx as product expression */
static SCIP_EXPR* xexpr;
static SCIP_EXPR* yexpr;
static SCIP_RANDNUMGEN* rndgen;

/* creates scip, problem, includes expression handlers, creates and adds variables */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );

   /* create variable and entropy expressions */
   SCIP_CALL( SCIPcreateExprVar(scip, &xexpr, x, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprEntropy(scip, &entropyexpr, xexpr, NULL, NULL) );
   SCIP_CALL( SCIPparseExpr(scip, &prodexpr, "<x>[C]*log(<x>[C])", NULL, NULL, NULL) );
   SCIP_CALL( SCIPparseExpr(scip, &negprodexpr, "-<x>[C]*log(<x>[C])", NULL, NULL, NULL) );

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
   SCIP_CALL( SCIPreleaseExpr(scip, &negprodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &prodexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &entropyexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
   SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

/* test suite */
TestSuite(entropy, .init = setup, .fini = teardown);

/*
 * TESTS
 */

Test(entropy, creation, .description = "Tests the expression creation.")
{
   SCIP_EXPR* expr;

   /* create entropy expression */
   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr, xexpr, NULL, NULL) );

   cr_assert(expr != NULL);
   cr_expect(SCIPexprGetNChildren(expr) == 1);
   cr_expect(SCIPexprGetChildren(expr)[0] == xexpr);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(entropy, parse, .description = "Tests the expression parsing.")
{
   SCIP_EXPR* expr;
   const char* input = "entropy(<x>[C])";

   /* create entropy expression */
   SCIP_CALL( SCIPparseExpr(scip, &expr, (char*)input, NULL, NULL, NULL) );

   cr_assert(expr != NULL);
   cr_expect(SCIPexprGetNChildren(expr) == 1);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(expr)[0]));
   cr_expect(SCIPgetVarExprVar(SCIPexprGetChildren(expr)[0]) == x);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}

Test(entropy, eval, .description = "Tests the expression evaluation.")
{
   SCIP_Real testvalues[4] = {-1.0, 0.0, exp(-1.0), 1.0};
   SCIP_Real results[4] = {SCIP_INVALID, 0.0, exp(-1.0), 0.0};
   SCIP_Real randnum;
   int i;

   /* deterministic part */
   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalExpr(scip, entropyexpr, sol, 0) );

      cr_expect(SCIPisEQ(scip, SCIPexprGetEvalValue(entropyexpr), results[i]));
   }

   /* random part */
   for( i = 0; i < 100; ++i )
   {
      randnum = SCIPrandomGetReal(rndgen, 0.0, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );

      SCIP_CALL( SCIPevalExpr(scip, entropyexpr, sol, 0) );
      cr_expect(SCIPisEQ(scip, SCIPexprGetEvalValue(entropyexpr), -randnum * log(randnum)));

      SCIP_CALL( SCIPevalExpr(scip, negprodexpr, sol, 0) );
      cr_expect(SCIPisEQ(scip, SCIPexprGetEvalValue(negprodexpr), -randnum * log(randnum)));
   }
}

Test(entropy, inteval, .description = "Tests the expression interval evaluation.")
{
   SCIP_INTERVAL intervalEntropy;
   SCIP_Real rndlb[4];
   SCIP_Real rndub[4];
   SCIP_Real rndreslb[4];
   SCIP_Real rndresub[4];
   int i;

   /* pick 5 special cases with well known results */
   SCIP_Real detlb[4] = {0.0, 0.0, exp(-1.0), 1.0};
   SCIP_Real detub[4] = {exp(-1.0), 1.0, 1.0, 2.0};
   SCIP_Real detreslb[4] = {0.0, 0.0, 0.0, -2.0 * log(2.0)};
   SCIP_Real detresub[4] = {exp(-1.0), exp(-1.0), exp(-1.0), 0.0};

   /* create 5 random cases within specific bounds that have non-trivial results */
   rndlb[0]    = SCIPrandomGetReal(rndgen, 0.0, exp(-1.0));
   rndub[0]    = SCIPrandomGetReal(rndgen, exp(-1.0), 1.0);
   rndreslb[0] = MIN(-rndlb[0] * log(rndlb[0]), -rndub[0] * log(rndub[0]));
   rndresub[0] = exp(-1.0);

   rndlb[1]    = SCIPrandomGetReal(rndgen, 0.0, exp(-1.0));
   rndub[1]    = SCIPrandomGetReal(rndgen, 1.0, 10.0);
   rndreslb[1] = -rndub[1] * log(rndub[1]);
   rndresub[1] = exp(-1.0);

   rndlb[2]    = SCIPrandomGetReal(rndgen, exp(-1.0), 1.0);
   rndub[2]    = SCIPrandomGetReal(rndgen, 1.0, 10.0);
   rndreslb[2] = -rndub[2] * log(rndub[2]);
   rndresub[2] = -rndlb[2] * log(rndlb[2]);

   rndlb[3]    = SCIPrandomGetReal(rndgen, 1.0, 10.0);
   rndub[3]    = SCIPrandomGetReal(rndgen, rndlb[3], 10.0);
   rndreslb[3] = -rndub[3] * log(rndub[3]);
   rndresub[3] = -rndlb[3] * log(rndlb[3]);

   /* deterministic part */
   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, detlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, detub[i]) );
      SCIP_CALL( SCIPevalExprActivity(scip, entropyexpr) );
      intervalEntropy = SCIPexprGetActivity(entropyexpr);

      cr_expect(SCIPisEQ(scip, intervalEntropy.inf, detreslb[i]));
      cr_expect(SCIPisEQ(scip, intervalEntropy.sup, detresub[i]));
   }

   /* random part */
   for ( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPchgVarLb(scip, x, rndlb[i]) );
      SCIP_CALL( SCIPchgVarUb(scip, x, rndub[i]) );
      SCIP_CALL( SCIPevalExprActivity(scip, entropyexpr) );
      intervalEntropy = SCIPexprGetActivity(entropyexpr);

      cr_expect(SCIPisEQ(scip, intervalEntropy.inf, rndreslb[i]));
      cr_expect(SCIPisEQ(scip, intervalEntropy.sup, rndresub[i]));
   }
}

Test(entropy, derivative, .description = "Tests the expression derivation.")
{
   SCIP_Real testvalues[5] = {-1.0, 0.0, exp(-1.0), 1.0, 2.0};
   SCIP_Real results[5] = {SCIP_INVALID, SCIP_INVALID, 0.0, -1.0, -1.0 - log(2.0)};
   SCIP_Real randnum;
   int i;

   /* deterministic part */
   for( i = 0; i < 5; ++i )
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, testvalues[i]) );
      SCIP_CALL( SCIPevalExprGradient(scip, entropyexpr, sol, 0) );

      /* iff cannot be differentiated, then entropyexpr->derivative is SCIP_INVALID */
      cr_expect((SCIPexprGetDerivative(entropyexpr) == SCIP_INVALID) == (results[i] == SCIP_INVALID));
      /* if entropyexpr cannot be differentiated, then the equality may not hold since
       * xexpr->derivative may not be SCIP_INVALID */
      cr_expect(SCIPisEQ(scip, SCIPexprGetDerivative(xexpr), results[i]) || (results[i] == SCIP_INVALID));
   }

   /* random part */
   for( i = 0; i < 100; ++i)
   {
      randnum = SCIPrandomGetReal(rndgen, 1e-12, 10.0);
      SCIP_CALL( SCIPsetSolVal(scip, sol, x, randnum) );

      SCIP_CALL( SCIPevalExprGradient(scip, entropyexpr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetDerivative(xexpr), -1.0 - log(randnum)));

      SCIP_CALL( SCIPevalExprGradient(scip, negprodexpr, sol, 0) );
      cr_expect(SCIPisFeasEQ(scip, SCIPexprGetDerivative(xexpr), -1.0 - log(randnum)));
   }
}

Test(entropy, hash, .description = "Tests the expression hash.")
{
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_EXPR* expr3;
   unsigned int hashkey1;
   unsigned int hashkey2;
   unsigned int hashkey3;

   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr1, xexpr, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr2, xexpr, NULL, NULL) );
   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr3, yexpr, NULL, NULL) );

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

Test(entropy, simplify, .description = "Tests the expression simplification.")
{
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;
   SCIP_EXPR* expr3;
   SCIP_Bool changed = FALSE;
   SCIP_Bool infeasible;

   /* expr1 = <5.0>, expr2 = entropy(<5.0>), expr3 is buffer for simplification */
   SCIP_CALL( SCIPcreateExprValue(scip, &expr1, 5.0, NULL, NULL));
   SCIP_CALL( SCIPcreateExprEntropy(scip, &expr2, expr1, NULL, NULL));
   SCIP_CALL( SCIPsimplifyExpr(scip, expr2, &expr3, &changed, &infeasible, NULL, NULL));
   SCIP_CALL( SCIPevalExpr(scip, expr2, sol, 0) );

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPexprGetHdlr(expr3) == SCIPgetExprhdlrValue(scip));
   cr_expect(SCIPisFeasEQ(scip, SCIPexprGetEvalValue(expr2), -5.0 * log(5.0)));

   SCIP_CALL( SCIPreleaseExpr(scip, &expr3) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );

   /* test if product of x and log(x) is transformed to sum of entropy expression
    * expr1 is buffer for simplification and expr2 is used to store children
    */
   changed = FALSE;
   SCIP_CALL( SCIPsimplifyExpr(scip, prodexpr, &expr1, &changed, &infeasible, NULL, NULL));

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPexprGetHdlr(expr1) == SCIPgetExprhdlrSum(scip));
   cr_expect(SCIPexprGetNChildren(expr1) == 1);
   cr_expect(SCIPgetCoefsExprSum(expr1)[0] == -1.0);

   expr2 = SCIPexprGetChildren(expr1)[0];
   cr_expect(SCIPexprGetHdlr(expr2) == SCIPfindExprhdlr(scip, "entropy"));
   cr_expect(SCIPexprGetNChildren(expr2) == 1);
   cr_expect(SCIPisExprVar(scip, SCIPexprGetChildren(expr2)[0]));
   cr_expect(SCIPgetVarExprVar(SCIPexprGetChildren(expr2)[0]) == x);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );


   /* test if product of -x and log(x) is transformed to entropy expression
    * expr1 is buffer for simplification
    */
   changed = FALSE;
   SCIP_CALL( SCIPsimplifyExpr(scip, negprodexpr, &expr1, &changed, &infeasible, NULL, NULL));

   cr_expect(changed);
   cr_assert_not(infeasible);
   cr_expect(SCIPexprGetHdlr(expr1) == SCIPfindExprhdlr(scip, "entropy"));
   cr_expect(SCIPexprGetNChildren(expr1) == 1);
   cr_expect(SCIPgetVarExprVar(SCIPexprGetChildren(expr1)[0]) == x);

   SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );
}
