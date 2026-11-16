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

/**@file   empty.c
 * @brief  unit test for counting solutions of the trivial problem
 * @author Dominik Kamp
 */

#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;

/* TEST SUITE */

/** create SCIP instance */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "empty") );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIPfree(&scip);
}


/* TESTS  */
/** @brief test checking that the empty solution is counted */
void test_empty_counter(void)
{
   SCIP_CALL( SCIPsetParamsCountsols(scip) );
   SCIP_CALL( SCIPcount(scip) );
   SCIP_Bool valid;
   SCIP_Longint count = SCIPgetNCountedSols(scip, &valid);
   int nsols = SCIPgetNSols(scip);
   TEST_ASSERT(valid);
   TEST_ASSERT_EQUAL(count, 1);
   TEST_ASSERT_EQUAL(nsols, 0);
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_empty_counter);
   return UNITY_END();
}
