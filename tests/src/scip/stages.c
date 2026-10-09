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

/**@file   stages.c
 * @brief  unit test for checking setters on scip.c
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"


/** GLOBAL VARIABLES **/
static SCIP* scip;

/** TEST SUITES **/
static
void setup(void)
{

   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

void setUp(void)
{
   SOFT_ASSERT_RESET();
   setup();
}

void tearDown(void)
{
   teardown();
   SOFT_ASSERT_CHECK();
}

/* test that we get to each of the following stages:
 *  SCIP_STAGE_TRANSFORMED
 *  SCIP_STAGE_PRESOLVING
 *  SCIP_STAGE_PRESOLVED
 *  SCIP_STAGE_SOLVING
 *  SCIP_STAGE_SOLVED
 */

/* helper method */
static
void gotoStage(SCIP_STAGE stage)
{
   SCIP_CALL( TESTscipSetStage(scip, stage, FALSE) );
   SOFT_ASSERT_EQUAL(SCIPgetStage(scip), stage, "got stage %d, expected %d", SCIPgetStage(scip), stage);
}

void test_stages_transformed(void)
{
   gotoStage(SCIP_STAGE_TRANSFORMED);
}

void test_stages_presolving(void)
{
   gotoStage(SCIP_STAGE_PRESOLVING);
}

void test_stages_presolved(void)
{
   gotoStage(SCIP_STAGE_PRESOLVED);
}

void test_stages_solving(void)
{
   gotoStage(SCIP_STAGE_SOLVING);
}

void test_stages_solved(void)
{
   gotoStage(SCIP_STAGE_SOLVED);
}

void test_stages_solving_with_nlp(void)
{
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );
   SOFT_ASSERT_EQUAL(SCIPgetStage(scip), SCIP_STAGE_SOLVING, "got stage %d, expected %d", SCIPgetStage(scip), SCIP_STAGE_SOLVING);

   /* check that NLP is created */
   SOFT_ASSERT(SCIPisNLPConstructed(scip), "NLP is not constructed");
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_stages_transformed);
   RUN_TEST(test_stages_presolving);
   RUN_TEST(test_stages_presolved);
   RUN_TEST(test_stages_solving);
   RUN_TEST(test_stages_solved);
   RUN_TEST(test_stages_solving_with_nlp);
   return UNITY_END();
}
