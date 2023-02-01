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

/**@file   bugfix-fixing.c
 * @brief  unit test for bugfix on XOR constraint and bounds
 * @author Mathieu Besançon
 */

#include "scip/scip.h"
#include "include/scip_test.h"
#include "scip/cons_xor.h"
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
   SCIP_CALL( SCIPcreateProbBasic(scip, "xor-bugcheck") );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(bugfixxor, .init = setup, .fini = teardown);

/* TESTS  */
Test(bugfixxor, demofix, .description = "test checking that the optimal solution of a problem with XOR constraints is not cut off")
{
   SCIP_CALL( SCIPreadProb(scip, "src/cons/xor/Demo5.cip", NULL) );
   SCIP_CALL( SCIPsolve(scip) );
   int nsols = SCIPgetNSols(scip);
   SCIP_SOL** sols = SCIPgetSols(scip);
   cr_assert_geq(nsols, 1);
   cr_assert_float_eq(SCIPgetSolOrigObj(scip, sols[0]), 2.0, 1e-5);
}
