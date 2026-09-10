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
   cr_assert_not_null(scip);
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(numerics, .init = setup, .fini = teardown);

/* TESTS */

/** default parameters should pass and allow transforming */
Test(numerics, defaultParamsTransform)
{
   SCIP_CALL( SCIPtransformProb(scip) );
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_TRANSFORMED);
}

/** inconsistent settings can be set in PROBLEM stage (deferred check) */
Test(numerics, inconsistentSetInProblemStage)
{
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** epsilon > feastol rejected at transform, stays in PROBLEM */
Test(numerics, epsilonExceedsFeastolRejectsTransform)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );

   retcode = SCIPtransformProb(scip);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** epsilon > sumepsilon rejected at transform */
Test(numerics, epsilonExceedsSumepsilonRejectsTransform)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-5) );

   retcode = SCIPtransformProb(scip);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** sumepsilon > feastol rejected at transform */
Test(numerics, sumepsilonExceedsFeastolRejectsTransform)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-4) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-5) );

   retcode = SCIPtransformProb(scip);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);
}

/** consistent params in any order should allow transforming */
Test(numerics, consistentParamsAnyOrder)
{
   /* epsilon set before feastol -- would fail with per-param callbacks */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-3) );

   SCIP_CALL( SCIPtransformProb(scip) );
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_TRANSFORMED);
}

/** after transform, setting epsilon too high should be rejected and reverted */
Test(numerics, rejectInvalidChangeAfterTransform)
{
   SCIP_RETCODE retcode;
   SCIP_Real oldepsilon;

   SCIP_CALL( SCIPtransformProb(scip) );

   oldepsilon = SCIPepsilon(scip);

   /* try to set epsilon larger than feastol */
   retcode = SCIPsetRealParam(scip, "numerics/epsilon", 1e-4);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);

   /* value should be reverted */
   cr_assert_eq(SCIPepsilon(scip), oldepsilon); /*lint !e777*/
}

/** after transform, lowering feastol below sumepsilon should be rejected and reverted */
Test(numerics, rejectFeastolBelowSumepsilonAfterTransform)
{
   SCIP_RETCODE retcode;
   SCIP_Real oldfeastol;

   SCIP_CALL( SCIPtransformProb(scip) );

   oldfeastol = SCIPfeastol(scip);

   /* try to set feastol below default sumepsilon (1e-6) */
   retcode = SCIPsetRealParam(scip, "numerics/feastol", 1e-8);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);

   /* value should be reverted */
   cr_assert_eq(SCIPfeastol(scip), oldfeastol); /*lint !e777*/
}

/** after transform, lowering dualfeastol below epsilon should be rejected and reverted */
Test(numerics, rejectDualfeastolBelowEpsilonAfterTransform)
{
   SCIP_RETCODE retcode;
   SCIP_Real olddualfeastol;

   SCIP_CALL( SCIPtransformProb(scip) );

   olddualfeastol = SCIPdualfeastol(scip);

   /* try to set dualfeastol below default epsilon (1e-9) */
   retcode = SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-11);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);

   /* value should be reverted */
   cr_assert_eq(SCIPdualfeastol(scip), olddualfeastol); /*lint !e777*/
}

/** after transform, consistent change should be accepted */
Test(numerics, acceptValidChangeAfterTransform)
{
   SCIP_CALL( SCIPtransformProb(scip) );

   /* epsilon = 1e-12 is below all defaults */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-12) );
   cr_assert_eq(SCIPepsilon(scip), (SCIP_Real)1e-12); /*lint !e777*/
}

/** user can fix settings after rejected transform and retry */
Test(numerics, fixAndRetryTransform)
{
   SCIP_RETCODE retcode;

   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-4) );

   /* first attempt should fail */
   retcode = SCIPtransformProb(scip);
   cr_assert_eq(retcode, SCIP_PARAMETERWRONGVAL);
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_PROBLEM);

   /* fix: raise feastol, sumepsilon, dualfeastol above epsilon */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-3) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/dualfeastol", 1e-3) );

   /* retry should succeed */
   SCIP_CALL( SCIPtransformProb(scip) );
   cr_assert_eq(SCIPgetStage(scip), SCIP_STAGE_TRANSFORMED);
}
