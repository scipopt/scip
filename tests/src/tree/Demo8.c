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

/**@file   Demo8.c
 * @brief  unit test for root node repropagation during tree path switch
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
   SCIP_CALL( SCIPcreateProbBasic(scip, "Demo8") );
}

/** free SCIP instance */
static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(Demo8, .init = setup, .fini = teardown);

/* TESTS  */
Test(Demo8, repropagation, .description = "test checking that no infeasible solutions are accepted during root node repropagation")
{
   SCIP_CALL( SCIPreadParams(scip, "src/tree/Demo8.set") );
   SCIP_CALL( SCIPreadProb(scip, "src/tree/Demo8.cip", NULL) );
   SCIP_CALL( SCIPsolve(scip) );
   int nvars = SCIPgetNOrigVars(scip);
   int nsols = SCIPgetNSols(scip);
   cr_assert_eq(nvars, 38);
   cr_assert_geq(nsols, 1);
   SCIP_VAR** vars = SCIPgetOrigVars(scip);
   SCIP_SOL* sol = SCIPgetSols(scip)[0];
   int i;
   for(i = 0; i < 4    ; ++i) cr_assert_float_eq(SCIPgetSolVal(scip, sol, var[i]), 1.0, 1e-5);
   for(i = 4; i < nvars; ++i) cr_assert_float_eq(SCIPgetSolVal(scip, sol, var[i]), 0.0, 1e-5);
   cr_assert_float_eq(SCIPgetSolOrigObj(scip, sol), 3.0, 1e-5);
}
