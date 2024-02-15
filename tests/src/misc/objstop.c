/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   objstop.c
 * @brief  unit tests for objective stop
 * @author Leon Eifler
 */

#include<stdio.h>

#include "scip/pub_misc.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "include/scip_test.h"

/* helper methods */

/** GLOBAL VARIABLES **/
static SCIP* scip;

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

TestSuite(objstop, .init = setup, .fini = teardown);

Test(objstop, objstop_min, .description = "tests objective stop for minimization")
{
   SCIP_Real objstop;
   cr_assert_not_null(scip);

   SCIPreadProb(scip, "../check/instances/MIP/rgn.mps", NULL);

   SCIPsetRealParam(scip, "limits/objectivestop", 352);
   SCIPgetRealParam(scip, "limits/objectivestop", &objstop);
   /* at the moment, the first solution found with objective value at most 352 is an optimal solution
    * we turn off separation to avoid that the dual bound is already at the optimal value when this happens, since SCIP would then return with status SCIP_STATUS_OPTIMAL
    */
   SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE);
   cr_assert_eq(objstop, 352, "Objective stop should be 352.0, but is %f", objstop);

   SCIPsolve(scip);

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_USERINTERRUPT, "SCIP should have stopped at objective value. Stopped with %d", SCIPgetStatus(scip));

}

Test(objstop, objstop_max, .description = "tests objective stop for maximization")
{
   SCIP_Real objstop;
   cr_assert_not_null(scip);

   SCIPreadProb(scip, "../check/instances/Symmetry/packorb_1-FullIns_3.cip", NULL);

   SCIPsetRealParam(scip, "limits/objectivestop", 24);
   SCIPgetRealParam(scip, "limits/objectivestop", &objstop);
   cr_assert_eq(objstop, 24, "Objective stop should be 24.0, but is %f", objstop);

   SCIPsolve(scip);

   cr_assert_eq(SCIPgetStatus(scip), SCIP_STATUS_USERINTERRUPT, "SCIP should have stopped at objective value. Stopped with %d", SCIPgetStatus(scip));

}
