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

/**@file   numerics.c
 * @brief  unit test for numeric parameter consistency checks
 * @author João Dionísio
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* UNIT TEST */

#include "include/scip_test.h"

/* GLOBAL VARIABLES */
static SCIP* scip;

/* TEST SUITES */
static
void setup(void)
{
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );
   TEST_ASSERT_NOT_NULL(scip);
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/* TESTS */

/** @brief default parameters should pass and allow transforming */
void test_numerics_defaultParamsTransform(void)
{
   SCIP_CALL( SCIPtransformProb(scip) );
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_TRANSFORMED);
}

/** @brief inconsistent settings can be set in PROBLEM stage (deferred check) */
void test_numerics_inconsistentSetInProblemStage(void)
{
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** @brief epsilon > feastol rejected at transform, stays in PROBLEM */
void test_numerics_epsilonExceedsFeastolRejectsTransform(void)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );

   retcode = SCIPtransformProb(scip);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** @brief epsilon > sumepsilon rejected at transform */
void test_numerics_epsilonExceedsSumepsilonRejectsTransform(void)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-5) );

   retcode = SCIPtransformProb(scip);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** @brief sumepsilon > feastol rejected at transform */
void test_numerics_sumepsilonExceedsFeastolRejectsTransform(void)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-4) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-5) );

   retcode = SCIPtransformProb(scip);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** @brief consistent params in any order should allow transforming */
void test_numerics_consistentParamsAnyOrder(void)
{
   /* epsilon set before feastol -- would fail with per-param callbacks */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-3) );

   SCIP_CALL( SCIPtransformProb(scip) );
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_TRANSFORMED);
}

/** @brief after transform, setting epsilon too high should be rejected and reverted */
void test_numerics_rejectInvalidChangeAfterTransform(void)
{
   SCIP_RETCODE retcode;
   SCIP_Real oldepsilon;

   SCIP_CALL( SCIPtransformProb(scip) );

   oldepsilon = SCIPepsilon(scip);

   /* try to set epsilon larger than feastol */
   retcode = SCIPsetRealParam(scip, "numerics/epsilon", 1e-4);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);

   /* value should be reverted */
   TEST_ASSERT_EQUAL(SCIPepsilon(scip), oldepsilon); /*lint !e777*/
}

/** @brief after transform, lowering feastol below sumepsilon should be rejected and reverted */
void test_numerics_rejectFeastolBelowSumepsilonAfterTransform(void)
{
   SCIP_RETCODE retcode;
   SCIP_Real oldfeastol;

   SCIP_CALL( SCIPtransformProb(scip) );

   oldfeastol = SCIPfeastol(scip);

   /* try to set feastol below default sumepsilon (1e-6) */
   retcode = SCIPsetRealParam(scip, "numerics/feastol", 1e-8);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);

   /* value should be reverted */
   TEST_ASSERT_EQUAL(SCIPfeastol(scip), oldfeastol); /*lint !e777*/
}

/** @brief after transform, lowering dualfeastol below epsilon should be rejected and reverted */
void test_numerics_rejectDualfeastolBelowEpsilonAfterTransform(void)
{
   SCIP_RETCODE retcode;
   SCIP_Real olddualfeastol;

   SCIP_CALL( SCIPtransformProb(scip) );

   olddualfeastol = SCIPdualfeastol(scip);

   /* try to set dualfeastol below default epsilon (1e-9) */
   retcode = SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-11);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);

   /* value should be reverted */
   TEST_ASSERT_EQUAL(SCIPdualfeastol(scip), olddualfeastol); /*lint !e777*/
}

/** @brief after transform, consistent change should be accepted */
void test_numerics_acceptValidChangeAfterTransform(void)
{
   SCIP_CALL( SCIPtransformProb(scip) );

   /* epsilon = 1e-12 is below all defaults */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-12) );
   TEST_ASSERT_EQUAL(SCIPepsilon(scip), (SCIP_Real)1e-12); /*lint !e777*/
}

/** @brief user can fix settings after rejected transform and retry */
void test_numerics_fixAndRetryTransform(void)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );

   /* first attempt should fail */
   retcode = SCIPtransformProb(scip);
   TEST_ASSERT_EQUAL(retcode, SCIP_PARAMETERWRONGVAL);
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);

   /* fix: raise feastol, sumepsilon, dualfeastol above epsilon */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-3) );

   /* retry should succeed */
   SCIP_CALL( SCIPtransformProb(scip) );
   TEST_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_TRANSFORMED);
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_numerics_acceptValidChangeAfterTransform);
   RUN_TEST(test_numerics_consistentParamsAnyOrder);
   RUN_TEST(test_numerics_defaultParamsTransform);
   RUN_TEST(test_numerics_epsilonExceedsFeastolRejectsTransform);
   RUN_TEST(test_numerics_epsilonExceedsSumepsilonRejectsTransform);
   RUN_TEST(test_numerics_fixAndRetryTransform);
   RUN_TEST(test_numerics_inconsistentSetInProblemStage);
   RUN_TEST(test_numerics_rejectDualfeastolBelowEpsilonAfterTransform);
   RUN_TEST(test_numerics_rejectFeastolBelowSumepsilonAfterTransform);
   RUN_TEST(test_numerics_rejectInvalidChangeAfterTransform);
   RUN_TEST(test_numerics_sumepsilonExceedsFeastolRejectsTransform);
   return UNITY_END();
}
