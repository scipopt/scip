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

/**@file   varboundnumerics.c
 * @brief  unit test for problem with difficult numerics
 * @author Gioni Mexi
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
   SCIP_CALL( SCIPcreateProbBasic(scip, "varboundnumerics") );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(problem, .init = setup, .fini = teardown);

/* TESTS  */
Test(problem, tolerance, .description = "test whether propagating a var bound constraint does not falsely remove it")
{
   SCIP_CALL( SCIPreadProb(scip, "src/bugs/varboundnumerics.cip", NULL) );
   SCIP_CALL( SCIPsolve(scip) );
   int nvars = SCIPgetNOrigVars(scip);
   int nconss = SCIPgetNOrigConss(scip);
   cr_assert_eq(nvars, 19);
   cr_assert_eq(nconss, 7);
   SCIP_SOL* sol = SCIPgetSols(scip)[0];
   cr_assert_float_eq(SCIPgetSolOrigObj(scip, sol), 57138036.3824479, 1e-5);
}
