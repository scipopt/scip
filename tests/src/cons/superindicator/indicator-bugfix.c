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
