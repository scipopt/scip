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

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(stages, .init = setup, .fini = teardown);


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
   cr_expect_eq(SCIPgetStage(scip), stage, "got stage %d, expected %d", SCIPgetStage(scip), stage);
}

Test(stages, transformed)
{
   gotoStage(SCIP_STAGE_TRANSFORMED);
}

Test(stages, presolving)
{
   gotoStage(SCIP_STAGE_PRESOLVING);
}

Test(stages, presolved)
{
   gotoStage(SCIP_STAGE_PRESOLVED);
}

Test(stages, solving)
{
   gotoStage(SCIP_STAGE_SOLVING);
}

Test(stages, solved)
{
   gotoStage(SCIP_STAGE_SOLVED);
}

Test(stages, solving_with_nlp)
{
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, TRUE) );
   cr_expect_eq(SCIPgetStage(scip), SCIP_STAGE_SOLVING, "got stage %d, expected %d", SCIPgetStage(scip), SCIP_STAGE_SOLVING);

   /* check that NLP is created */
   cr_expect(SCIPisNLPConstructed(scip), "NLP is not constructed");
}
