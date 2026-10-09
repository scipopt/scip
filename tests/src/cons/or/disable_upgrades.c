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

/**@file   disable_upgrades.c
 * @brief  tests for or-constraint methods with upgrades to and-constraints disabled
 * @author Gregor Hendel
 * @author Helena Müller
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_or.h"
#include "scip/pub_cons.h"
#include <stdio.h>

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;


/* TEST SUITE */

/** the setup creates the necessary source SCIP, sets parameter */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPsetIntParam(scip, "constraints/or/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "branching/inference/priority", 1000000) );

   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_OFF, FALSE) );
}

/** free all allocated memory */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

}



/* TESTS  */
void test_disable_upgrades_create_and_free(void)
{
   /* calls setup and teardown */
}

/** @brief disable upgrades of or-constraints to and-constraints to keep or-constraints during solution process */
void test_disable_upgrades_disable_upgrades_or(void)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   int nconss;
   int i;

   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Or/or_constraint.cip", "cip") );

   conshdlr = SCIPfindConshdlr(scip, "or");

   nconss = SCIPgetNOrigConss(scip);
   conss = SCIPgetOrigConss(scip);

   TEST_ASSERT_EQUAL(nconss, 8);

   /* set or-constraints to be modifiable */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPconsGetHdlr(conss[i]) == conshdlr )
         SCIPsetConsModifiable(scip, conss[i], TRUE);
   }

   SCIP_CALL( SCIPchgVarUbGlobal(scip, SCIPgetVars(scip)[3], 0.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, SCIPgetVars(scip)[4], 0.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, SCIPgetVars(scip)[2], 0.0) );
   SCIP_CALL( SCIPchgVarUbGlobal(scip, SCIPgetVars(scip)[1], 0.0) );

   SCIP_CALL( SCIPsetIntParam(scip, "constraints/or/maxprerounds", 1) );

   SCIP_CALL( SCIPsolve(scip) );
}

/** @brief test that CIP write method works for or constraints */
void test_disable_upgrades_write_problem(void)
{
   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Or/or_constraint.cip", "cip") );

   SCIP_CALL( SCIPwriteOrigProblem(scip, NULL, "cip", FALSE) );
}

/** @brief test copying of or-constraints */
void test_disable_upgrades_copy_problem(void)
{
   SCIP* targetscip;
   SCIP_Bool valid;
   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Or/or_constraint.cip", "cip") );

   SCIP_CALL( SCIPcreate(&targetscip) );

   SCIP_CALL( SCIPcopy(scip, targetscip, NULL, NULL, "copy_of_prob", TRUE, TRUE, FALSE, FALSE, FALSE, &valid) );

   TEST_ASSERT(valid);

   SCIP_CALL( SCIPsolve(targetscip) );

   SCIP_CALL( SCIPfree(&targetscip) );
}

/** @brief test SCIPcreateConsBasicOr */
void test_disable_upgrades_test_basic_creation(void)
{
   SCIP_CONS* newcons;
   SCIP_VAR** origvars;

   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Or/or_constraint.cip", "cip") );

   TEST_ASSERT(SCIPgetNVars(scip) == 24);
   TEST_ASSERT(SCIPgetNBinVars(scip) == 24);

   origvars = SCIPgetOrigVars(scip);

   SCIP_CALL( SCIPcreateConsBasicOr(scip, &newcons, "new_or_constraint", origvars[3], 6, &origvars[4]) );

   SCIP_CALL( SCIPaddCons(scip, newcons) );

   SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPprintBestSol(scip, NULL, TRUE) );

}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_disable_upgrades_create_and_free);
   RUN_TEST(test_disable_upgrades_disable_upgrades_or);
   RUN_TEST(test_disable_upgrades_write_problem);
   RUN_TEST(test_disable_upgrades_copy_problem);
   RUN_TEST(test_disable_upgrades_test_basic_creation);
   return UNITY_END();
}
