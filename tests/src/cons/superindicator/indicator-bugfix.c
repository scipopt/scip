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

/**@file   indicator-bugfix.c
 * @brief  unit test for bugfix on indicator constraint presolving
 * @author Mathieu Besançon
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/scipdefplugins.h"

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;

/* TEST SUITE */

/** create SCIP instance */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "indicator-bugcheck") );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(bugfixxor, .init = setup, .fini = teardown);

/* TESTS  */
Test(bugfixxor, demofix, .description = "test checking that a problem with indicator constraint is not unbounded")
{
   SCIP_CALL( SCIPreadProb(scip, "src/cons/superindicator/indicator-bugfix-instance.cip", NULL) );
   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_geq(SCIPgetNSols(scip), 1);
   SCIP_SOL* sol = SCIPgetBestSol(scip);
   cr_assert(sol != NULL);
   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);
   cr_assert(EPSEQ(SCIPgetSolOrigObj(scip, sol), 28.75, SCIPfeastol(scip)));
}
