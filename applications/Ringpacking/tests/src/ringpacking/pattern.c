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

/**@file   pattern.c
 * @brief  unit test for testing pattern interface functions
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
static SCIP_PROBDATA* probdata;
static SCIP_PATTERN* rpattern;
static SCIP_PATTERN* cpattern;

/** setup of test run */
static
void setup(void)
{
   SCIP_Real rexts[3] = {1.0, 0.6, 0.5};
   SCIP_Real rints[3] = {1.0, 0.5, 0.0};
   int demands[3] = {100, 100, 100};

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
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 3, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   TEST_ASSERT(probdata != NULL);

   /* creates circular and rectangular pattern */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &cpattern, 1) );
   SCIP_CALL( SCIPpatternCreateRectangular(scip, &rpattern) );
}

/** deinitialization method */
static
void teardown(void)
{
   /* release patterns */
   SCIPpatternRelease(scip, &rpattern);
   SCIPpatternRelease(scip, &cpattern);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(0, BMSgetMemoryUsed(), "There is a memory leak!!");
}

void setUp(void) { setup(); }
void tearDown(void) { teardown(); }

/* checks the pattern */
void test_pattern_patterntype(void)
{
   SOFT_ASSERT(SCIPpatternGetPatternType(cpattern) == SCIP_PATTERNTYPE_CIRCULAR);
   SOFT_ASSERT(SCIPpatternGetPatternType(rpattern) == SCIP_PATTERNTYPE_RECTANGULAR);
}

/* checks the type of a circular pattern */
void test_pattern_type(void)
{
   SOFT_ASSERT(SCIPpatternGetCircleType(cpattern) == 1);
}

/* checks the position of an element */
void test_pattern_position(void)
{
   SCIP_CALL( SCIPpatternAddElement(cpattern, 0, -1.0, 1.0) );
   SOFT_ASSERT(SCIPpatternGetElementPosX(cpattern, 0) == -1.0);
   SOFT_ASSERT(SCIPpatternGetElementPosY(cpattern, 0) == 1.0);

   SCIP_CALL( SCIPpatternAddElement(cpattern, 0, -2.0, 2.0) );
   SOFT_ASSERT(SCIPpatternGetElementPosX(cpattern, 1) == -2.0);
   SOFT_ASSERT(SCIPpatternGetElementPosY(cpattern, 1) == 2.0);
}

/* checks the packable status */
void test_pattern_packable(void)
{
   SOFT_ASSERT(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_UNKNOWN);

   SCIPpatternSetPackableStatus(rpattern, SCIP_PACKABLE_YES);
   SOFT_ASSERT(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_YES);

   /* adding an element does not change packable status */
   SCIP_CALL( SCIPpatternAddElement(rpattern, 0, 0.0, 0.0) );
   SOFT_ASSERT(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_YES);

   /* removing an element does not change packable status */
   SCIPpatternRemoveLastElements(rpattern, 1);
   SOFT_ASSERT(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_YES);
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_pattern_patterntype);
   RUN_TEST(test_pattern_type);
   RUN_TEST(test_pattern_position);
   RUN_TEST(test_pattern_packable);
   return UNITY_END();
}
