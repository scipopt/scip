/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probingobj.c
 * @brief  unit test for changing objective in probing
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include <signal.h>

#include "include/scip_test.h"

#define EPS 1e-6

static SCIP* scip = NULL;
static SCIP_VAR* x;
static SCIP_VAR* y;

/* methods to test: */


static
void setup(void)
{
   SCIP_VAR* a;
   SCIP_VAR* b;
   SCIP_CALL( SCIPcreate(&scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create and add variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &a, "x", 0.0, 1.0, 2.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &b, "y", 0.0, 1.0, -3.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, a) );
   SCIP_CALL( SCIPaddVar(scip, b) );

   /* change SCIP's stage to be able to create nlrows and rows */
   SCIP_CALL( TESTscipSetStage(scip, SCIP_STAGE_SOLVING, FALSE) );

   SCIP_CALL( SCIPgetTransformedVar(scip, a, &x) );
   SCIP_CALL( SCIPgetTransformedVar(scip, b, &y) );

   SCIP_CALL( SCIPreleaseVar(scip, &a) );
   SCIP_CALL( SCIPreleaseVar(scip, &b) );
}

static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}


/* TEST SUITE */
TestSuite(probingobj, .init = setup, .fini = teardown);

Test(probingobj, create_before_probing, .description="create a solution, start probing, change objective, end probing, and test that the objective value stays the same all the time")
{
   SCIP_SOL* sol;

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), 2.0, "expected 2.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPchgVarObjProbing(scip, x, 100.0) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), 2.0, "expected 2.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPendProbing(scip) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), 2.0, "expected 2.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

Test(probingobj, create_during_probing, .description="start probing, create a solution, change objective, end probing, and test that the objective value stays the same all the time")
{
   SCIP_SOL* sol;

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPchgVarObjProbing(scip, x, 100.0) );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 1.0) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), 2.0, "expected 2.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPendProbing(scip) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), 2.0, "expected 2.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

Test(probingobj, unlink_lpsol_during_probing, .description="start probing, solve LP, link solution to the LP, unlink solution, change objective, end probing, and test that the objective value stays the same all the time")
{
   SCIP_SOL* sol;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;

   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   cr_assert_not(lperror);
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL) );

   SCIP_CALL( SCIPunlinkSol(scip, sol) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPchgVarObjProbing(scip, y, 100.0) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPendProbing(scip) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

Test(probingobj, unlink_lpsol_after_objchange, .description="start probing, solve LP, link solution to the LP, change objective, unlink solution, end probing, and test that the objective value stays the same all the time")
{
   SCIP_SOL* sol;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;

   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   cr_assert_not(lperror);
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPchgVarObjProbing(scip, y, 100.0) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPunlinkSol(scip, sol) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPendProbing(scip) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}

/* @todo: this should actually be in another file, because it just tests that you cannot link a solution to an unsolved LP */
Test(probingobj, link_lpsol_after_objchange, .description="start probing, solve LP, change objective, link solution to the LP, unlink solution, end probing, and test that the objective value stays the same all the time", .signal = SIGABRT)
{
   SCIP_SOL* sol;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;

   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   cr_assert_not(lperror);
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPchgVarObjProbing(scip, y, 100.0) );

   SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL) );

#ifdef NDEBUG
   abort(); /* return SIGABORT in opt mode so it passes */
#endif

}

Test(probingobj, solve_lp_after_objchange, .description="start probing, change objective, solve LP, link solution to the LP, unlink solution, end probing, and test that the objective value stays the same all the time")
{
   SCIP_SOL* sol;
   SCIP_Bool cutoff;
   SCIP_Bool lperror;

   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPstartProbing(scip) );

   SCIP_CALL( SCIPnewProbingNode(scip) );

   SCIP_CALL( SCIPchgVarObjProbing(scip, y, -100.0) );

   SCIP_CALL( SCIPsolveProbingLP(scip, -1, &lperror, &cutoff) );
   cr_assert_not(lperror);
   cr_assert_not(cutoff);

   SCIP_CALL( SCIPcreateLPSol(scip, &sol, NULL) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPunlinkSol(scip, sol) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPendProbing(scip) );

   cr_expect_eq(SCIPgetSolOrigObj(scip, sol), -3.0, "expected -3.0, got %g\n", SCIPgetSolOrigObj(scip, sol));

   SCIP_CALL( SCIPfreeSol(scip, &sol) );
}
