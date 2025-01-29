/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   primal.c
 * @brief  unit tests for primal limit
 * @author Leon Eifler
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

TestSuite(primal, .init = setup, .fini = teardown);

Test(primal, primallimit_min, .description = "tests primal limit for minimization")
{
   const SCIP_Real target = 6040;
   SCIP_Real primallimit;

   cr_assert_not_null(scip);
   SCIP_CALL( SCIPsetRealParam(scip, "limits/primal", target) );
   /* turn on aggressive heuristics to get more limit candidates */
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
   /* turn off separators to oppose optimality before the limit */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPreadProb(scip, "../check/instances/MIP/rgn.mps", NULL) );
   SCIP_CALL( SCIPsolve(scip) );

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_PRIMALLIMIT, "SCIP terminated with status %d but should have terminated with status %d", SCIPgetStatus(scip), SCIP_STATUS_PRIMALLIMIT);
   primallimit = SCIPgetPrimalbound(scip);
   cr_assert_leq(primallimit, target, "Primal bound is %f but should be at most %f", primallimit, target);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/primal", &primallimit) );
   cr_assert_eq(primallimit, target, "Primal limit is %f but should be %f", primallimit, target);
}

Test(primal, primallimit_max, .description = "tests primal limit for maximization")
{
   const SCIP_Real target = 14;
   SCIP_Real primallimit;

   cr_assert_not_null(scip);
   SCIP_CALL( SCIPsetRealParam(scip, "limits/primal", target) );
   /* turn on aggressive heuristics to get more limit candidates */
   SCIP_CALL( SCIPsetHeuristics(scip, SCIP_PARAMSETTING_AGGRESSIVE, TRUE) );
   /* turn off separators to oppose optimality before the limit */
   SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
   SCIP_CALL( SCIPreadProb(scip, "../check/instances/Symmetry/packorb_1-FullIns_3.cip", NULL) );
   SCIP_CALL( SCIPsolve(scip) );

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_PRIMALLIMIT, "SCIP terminated with status %d but should have terminated with status %d", SCIPgetStatus(scip), SCIP_STATUS_PRIMALLIMIT);
   primallimit = SCIPgetPrimalbound(scip);
   cr_assert_geq(primallimit, target, "Primal bound is %f but should be at least %f", primallimit, target);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/primal", &primallimit) );
   cr_assert_eq(primallimit, target, "Primal limit is %f but should be %f", primallimit, target);
}
