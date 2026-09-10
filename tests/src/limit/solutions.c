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

/**@file   solutions.c
 * @brief  unit test for continued solve after solution limit
 * @author Dominik Kamp
 */

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
static SCIP* scip = NULL;

/* TEST SUITE */
static
void setup(void)
{
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );
   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(solutions, .init = setup, .fini = teardown);

Test(solutions, continued, .description = "tests continued solve after solution limit")
{
   cr_assert_not_null(scip);
   SCIP_CALL( SCIPsetIntParam(scip, "limits/solutions", 3) );
   SCIP_CALL( SCIPreadProb(scip, "../check/instances/MIP/MANN_a9.clq.lp", NULL) );
   SCIP_CALL( SCIPsolve(scip) );
   cr_assert(!SCIPisInfinity(scip, -SCIPgetPrimalbound(scip)),
         "Primal bound should be finite after solution limit");
   cr_assert(SCIPisFeasLE(scip, SCIPgetPrimalbound(scip), 16.0),
         "Primal bound is %f but should be <= 16.0", SCIPgetPrimalbound(scip));
   cr_assert(SCIPisFeasGE(scip, SCIPgetDualbound(scip), 16.0),
         "Dual bound is %f but should be >= 16.0", SCIPgetDualbound(scip));

   SCIP_CALL( SCIPsetIntParam(scip, "limits/solutions", -1) );
   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL,
         "SCIP terminated with status %d but should have terminated with status %d", SCIPgetStatus(scip), SCIP_STATUS_OPTIMAL);
   cr_assert(SCIPisFeasEQ(scip, SCIPgetPrimalbound(scip), 16.0),
         "Primal bound is %f but should be 16.0", SCIPgetPrimalbound(scip));
   cr_assert(SCIPisFeasEQ(scip, SCIPgetDualbound(scip), 16.0),
         "Dual bound is %f but should be 16.0", SCIPgetDualbound(scip));
}
