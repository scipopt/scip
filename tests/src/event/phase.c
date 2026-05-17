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

/**@file   phase.c
 * @brief  unit tests for the solving phase API
 * @author João Dionísio
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "include/scip_test.h"
#include "scip/scipdefplugins.h"
#include "scip/event_solvingphase.h"

/* known optimum of misc03 */
#define OPTIMAL_VALUE 3360.0

/** GLOBAL VARIABLES **/
static SCIP* scip_test = NULL;

/** TEST SUITES **/
static
void setup(void)
{
   char filename[SCIP_MAXSTRLEN];

   cr_assert(scip_test == NULL);

   SCIP_CALL( SCIPcreate(&scip_test) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip_test) );

   /* test on misc03 with optimum 3360, for which nodes are explored after an optimal solution is found */
   TESTsetTestfilename(filename, __FILE__, "../../../check/instances/MIP/misc03.mps");
   SCIP_CALL( SCIPreadProb(scip_test, filename, NULL) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip_test) );

   cr_assert(scip_test == NULL);
   cr_assert(BMSgetMemoryUsed() == 0, "There is a memory leak!");
}

TestSuite(solvingphase, .init = setup, .fini = teardown);

/* TESTS */

Test(solvingphase, disabled, .description = "solving phase getters report uninitialized state when the event handler is disabled")
{
   SCIP_CALL( SCIPsolve(scip_test) );

   cr_assert_eq(SCIPgetStatus(scip_test), SCIP_STATUS_OPTIMAL, "expected optimal status");

   cr_expect_eq(SCIPgetSolvingPhase(scip_test), SCIP_SOLVINGPHASE_UNINITIALIZED,
      "phase should remain uninitialized while the event handler is disabled");
   cr_expect_eq(SCIPgetSolvingPhaseFlags(scip_test), SCIP_SOLVINGPHASEFLAG_NONE,
      "no transition flags should be set while the event handler is disabled");
}

Test(solvingphase, incumbent, .description = "enabled mode tracks the solving phase and populates the transition flag bit-field")
{
   SCIP_SOLVINGPHASEFLAG flags;

   /* with transition method 'o' and provided optimal value, final phase is PROOF after OPTIMAL state is reached */
   SCIP_CALL( SCIPsetBoolParam(scip_test, "solvingphases/enabled", TRUE) );
   SCIP_CALL( SCIPsetCharParam(scip_test, "solvingphases/transitionmethod", 'o') );
   SCIP_CALL( SCIPsetRealParam(scip_test, "solvingphases/optimalvalue", OPTIMAL_VALUE) );

   SCIP_CALL( SCIPsolve(scip_test) );

   cr_assert_eq(SCIPgetStatus(scip_test), SCIP_STATUS_OPTIMAL, "expected optimal status");
   cr_assert(SCIPisEQ(scip_test, SCIPgetPrimalbound(scip_test), OPTIMAL_VALUE),
      "primal bound %.15g does not match expected optimum %.15g",
      SCIPgetPrimalbound(scip_test), OPTIMAL_VALUE);

   cr_expect_eq(SCIPgetSolvingPhase(scip_test), SCIP_SOLVINGPHASE_PROOF,
      "expected SCIP_SOLVINGPHASE_PROOF, got %d", (int)SCIPgetSolvingPhase(scip_test));

   flags = SCIPgetSolvingPhaseFlags(scip_test);
   cr_expect((flags & SCIP_SOLVINGPHASEFLAG_OPTIMAL) != 0,
      "expected SCIP_SOLVINGPHASEFLAG_OPTIMAL bit to be set under enabled mode, got flags = %#04x", (unsigned)flags);
}

Test(solvingphase, reachesproof, .description = "solving phase reaches SCIP_SOLVINGPHASE_PROOF when the optimal value is provided")
{
   SCIP_SOLVINGPHASEFLAG flags;

   /* testmode alone is enough to track both the solving phase and the reached transition criteria */
   SCIP_CALL( SCIPsetBoolParam(scip_test, "solvingphases/testmode", TRUE) );
   SCIP_CALL( SCIPsetCharParam(scip_test, "solvingphases/transitionmethod", 'o') );
   SCIP_CALL( SCIPsetRealParam(scip_test, "solvingphases/optimalvalue", OPTIMAL_VALUE) );

   SCIP_CALL( SCIPsolve(scip_test) );

   cr_assert_eq(SCIPgetStatus(scip_test), SCIP_STATUS_OPTIMAL, "expected optimal status");
   cr_assert(SCIPisEQ(scip_test, SCIPgetPrimalbound(scip_test), OPTIMAL_VALUE),
      "primal bound %.15g does not match expected optimum %.15g",
      SCIPgetPrimalbound(scip_test), OPTIMAL_VALUE);

   cr_expect_eq(SCIPgetSolvingPhase(scip_test), SCIP_SOLVINGPHASE_PROOF,
      "expected SCIP_SOLVINGPHASE_PROOF, got %d", (int)SCIPgetSolvingPhase(scip_test));

   flags = SCIPgetSolvingPhaseFlags(scip_test);
   cr_expect((flags & SCIP_SOLVINGPHASEFLAG_OPTIMAL) != 0,
      "expected SCIP_SOLVINGPHASEFLAG_OPTIMAL bit to be set, got flags = %#04x", (unsigned)flags);
}
