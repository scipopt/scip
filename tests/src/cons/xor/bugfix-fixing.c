/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
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
