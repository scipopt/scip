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

/**@file   verification.c
 * @brief  unit test for testing circular pattern verification
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "probdata_rpa.h"
#include "pricer_rpa.h"
#include "pattern.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_PATTERN* pattern;
static SCIP_PROBDATA* probdata;

/** setup of test run */
static
void setup(void)
{
   SCIP_Real rexts[4] = {2.414213562373095, 1.0, 0.6, 0.5};
   SCIP_Real rints[4] = {2.414213562373095, 1.0, 0.5, 0.0};
   int demands[4] = {1, 5, 5, 5};

   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include ringpacking pricer  */
   SCIP_CALL( SCIPincludePricerRpa(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create problem data */
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 4, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   TEST_ASSERT(probdata != NULL);

   /* creates a circular pattern */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, 1) );
}

/** deinitialization method */
static
void teardown(void)
{
   /* release circular pattern */
   SCIPpatternRelease(scip, &pattern);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(0, BMSgetMemoryUsed(), "There is a memory leak!!");
}

void setUp(void) { setup(); }
void tearDown(void) { teardown(); }

/** verifies empty circular pattern with NLP */
void test_verification_nlp_empty(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing a single element with NLP */
void test_verification_nlp_single(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing two elements with NLP */
void test_verification_nlp_two(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 3 and element of type 2 -> not packable */
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPpatternAddElement(pattern, 2, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_NO);

   /* remove element of type 2 and add another element of type 3 -> packable */
   SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);
   SCIPpatternRemoveLastElements(pattern, 1);
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing two elements with NLP */
void test_verification_nlp_complex(void)
{
   SCIP_PATTERN* p;
   int i;

   /* pattern with inner radius of 1+sqrt(2) */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &p, 0) );

   /* add four circles with radius 1 */
   for( i = 0; i < 4; ++i )
      SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);

   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, p, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_YES);

   /* add one more element -> not packable anymore */
   SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);
   SCIPpatternSetPackableStatus(p, SCIP_PACKABLE_UNKNOWN);
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, p, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_NO);

   /* release pattern */
   SCIPpatternRelease(scip, &p);
}

/** verifies empty circular pattern with heuristic */
void test_verification_heur_empty(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing a single element with heuristic */
void test_verification_heur_single(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verify circular pattern containing two elements with heuristic */
void test_verification_heur_two(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 3 and element of type 2 -> not packable */
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPpatternAddElement(pattern, 2, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* remove element of type 3 and add another element of type 2 -> packable */
   SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);
   SCIPpatternRemoveLastElements(pattern, 1);
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing four and five elements with heuristic */
void test_verification_heuer_four(void)
{
   SCIP_PATTERN* p;
   int i;

   /* pattern with inner radius of 1+sqrt(2) */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &p, 0) );

   /* add four circles with radius 1 */
   for( i = 0; i < 4; ++i )
      SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);

   /* easy enough for any heuristic */
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, p, SCIPinfinity(scip), 1) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_YES);

   /* add one more element -> not packable anymore (heuristic should only detect UNKNOWN) */
   SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);
   SCIPpatternSetPackableStatus(p, SCIP_PACKABLE_UNKNOWN);
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, p, SCIPinfinity(scip), 10) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_UNKNOWN);

   /* release pattern */
   SCIPpatternRelease(scip, &p);
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_verification_nlp_empty);
   RUN_TEST(test_verification_nlp_single);
   RUN_TEST(test_verification_nlp_two);
   RUN_TEST(test_verification_nlp_complex);
   RUN_TEST(test_verification_heur_empty);
   RUN_TEST(test_verification_heur_single);
   RUN_TEST(test_verification_heur_two);
   RUN_TEST(test_verification_heuer_four);
   return UNITY_END();
}
