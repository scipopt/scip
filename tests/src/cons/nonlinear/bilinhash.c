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

/**@file   bilinhash.c
 * @brief  tests functionalities to access all bilinear terms
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

static SCIP* scip;
static SCIP_CONSHDLR* conshdlr;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* z;

static
void setup(void)
{
   /* create SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* store nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "nonlinear");
   TEST_ASSERT(conshdlr != NULL);

   /* create problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", 0.0, 2.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   SCIP_CALL( SCIPsetPresolving(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, TRUE) );
}

static
void teardown(void)
{
   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &z) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &x) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory leak!!");
}

/* define test test suite */
/* tests the creating and release of the hash table using non-API methods from cons_nonlinear.c */
void test_bilinhash_createInsert(void)
{
   SCIP_CONSNONLINEAR_BILINTERM* bilinterms;

   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* inserts two bilinear terms into the hash table */
   SCIP_CALL( SCIPinsertBilinearTermExistingNonlinear(scip, conshdlr, x, y, NULL, 0, 0) );
   SCIP_CALL( SCIPinsertBilinearTermExistingNonlinear(scip, conshdlr, y, z, NULL, 0, 0) );

   TEST_ASSERT_EQUAL(SCIPgetNBilinTermsNonlinear(conshdlr), 2);

   bilinterms = SCIPgetBilinTermsNonlinear(conshdlr);

   SOFT_ASSERT_EQUAL(bilinterms[0].x, x);
   SOFT_ASSERT_EQUAL(bilinterms[0].y, y);
   SOFT_ASSERT_EQUAL(bilinterms[1].x, y);
   SOFT_ASSERT_EQUAL(bilinterms[1].y, z);
}

/* tests API methods for a simple problem containing two nonlinear constraints */
void test_bilinhash_api_methods(void)
{
   const char* inputs[2] = {"[nonlinear] <c1>: (<x>[C])^2 + <x>[C] * <y>[C] <= 4;",
      "[nonlinear] <c2>: abs(<y>[C] * <z>[C] + <x>[C] * <y>[C]) - (log(<x>[C] + <z>[C]))^2 <= 1;"};
   SCIP_CONSNONLINEAR_BILINTERM* bilinterms;
   SCIP_VAR* tx;
   SCIP_VAR* ty;
   SCIP_VAR* tz;
   SCIP_Bool cutoff;
   int i;

   /* create, add, and release nonlinear constraints */
   for( i = 0; i < 2; ++i )
   {
      SCIP_CONS* cons;
      SCIP_Bool success;

      SCIP_CALL( SCIPparseCons(scip, &cons, inputs[i],
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, &success) );
      TEST_ASSERT(success);
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* go to the solving stage */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   /* collect transformed variables */
   tx = SCIPvarGetTransVar(x);
   ty = SCIPvarGetTransVar(y);
   tz = SCIPvarGetTransVar(z);

   /* collect all bilinear terms by getting CONSINITLP called */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   SOFT_ASSERT_NOT(cutoff);

   /*
    * because auxiliary variables are present, there are four bilinear terms: xx, xy, yz, log()^2
    */
   SOFT_ASSERT_EQUAL(SCIPgetNBilinTermsNonlinear(conshdlr), 4);

   bilinterms = SCIPgetBilinTermsNonlinear(conshdlr);
   TEST_ASSERT(bilinterms != NULL);
   SOFT_ASSERT_EQUAL(bilinterms[0].x, tx);
   SOFT_ASSERT_EQUAL(bilinterms[0].y, tx);
   SOFT_ASSERT_NOT_NULL(bilinterms[0].aux.var);
   SOFT_ASSERT_EQUAL(bilinterms[1].x, tx);
   SOFT_ASSERT_EQUAL(bilinterms[1].y, ty);
   SOFT_ASSERT_NOT_NULL(bilinterms[1].aux.var);
   SOFT_ASSERT_EQUAL(bilinterms[2].x, ty);
   SOFT_ASSERT_EQUAL(bilinterms[2].y, tz);
   SOFT_ASSERT_NOT_NULL(bilinterms[2].aux.var);
   SOFT_ASSERT_EQUAL(bilinterms[3].x, bilinterms[3].y);

   /* xx exists */
   SOFT_ASSERT_NOT_NULL(SCIPgetBilinTermNonlinear(conshdlr, tx, tx));

   /* xy exists */
   SOFT_ASSERT_NOT_NULL(SCIPgetBilinTermNonlinear(conshdlr, tx, ty));

   /* yx = xy exists */
   SOFT_ASSERT_NOT_NULL(SCIPgetBilinTermNonlinear(conshdlr, ty, tx));

   /* yz exists */
   SOFT_ASSERT_NOT_NULL(SCIPgetBilinTermNonlinear(conshdlr, ty, tz));

   /* xz does not exist */
   SOFT_ASSERT_NULL(SCIPgetBilinTermNonlinear(conshdlr, tx, tz));

   /* zz does not exist */
   SOFT_ASSERT_NULL(SCIPgetBilinTermNonlinear(conshdlr, tz, tz));
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_bilinhash_createInsert);
   RUN_TEST(test_bilinhash_api_methods);
   return UNITY_END();
}
