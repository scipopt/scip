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

/**@file   presol.c
 * @brief  unit test for checking presolvers
 * @author Franziska Schloesser
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_linear.h"


/* UNIT TEST */
/* presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecUnittest)
{
   SCIP_VAR** vars;
   SCIP_Bool infeasible;
   SCIP_Bool fixed;

   assert(scip != NULL);

   vars = SCIPgetVars(scip);


   SCIP_CALL( SCIPfixVar(scip, vars[0], 2.0, &infeasible, &fixed) );
   SCIPdebugMessage("FIXED: %d INFEASIBLE: %d\n" , fixed, infeasible);

   if( fixed )
   {
      (*nfixedvars)++;
   }


   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#include "include/scip_test.h"

/* GLOBAL VARIABLES */
static SCIP_PRESOL* presol;
static SCIP* scip;

/** create bounded problem */
static
void initProb(void)
{
   SCIP_VAR* xvar;
   SCIP_VAR* yvar;
   SCIP_CONS* cons;
   SCIP_Real vals[2];
   SCIP_VAR* vars[2];
   char name[SCIP_MAXSTRLEN];

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &xvar, "x", 0, 3, -1.0, SCIP_VARTYPE_INTEGER) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &yvar, "y", 0, 3, -1.0, SCIP_VARTYPE_INTEGER) );

   SCIP_CALL( SCIPaddVar(scip, xvar) );
   SCIP_CALL( SCIPaddVar(scip, yvar) );

   /* add a constraint x + y <= 2 */
   vals[0] = 1.0;
   vals[1] = 1.0;
   vars[0] = xvar;
   vars[1] = yvar;
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "1. Cons");
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 2, vars, vals, 0.0, 2.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   SCIP_CALL( SCIPreleaseVar(scip, &xvar) );
   SCIP_CALL( SCIPreleaseVar(scip, &yvar) );
}


/* TEST SUITES */
/** setup of test run */
static
void setup(void)
{
   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );
   cr_assert_not_null(scip);

   /* include unittest presolver */
   presol = NULL;
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, "unittest", "presolver to test", 20010001, 20, SCIP_PRESOLTIMING_ALWAYS, presolExecUnittest, NULL) );
   SCIP_CALL( SCIPincludeConshdlrLinear(scip) );
   SCIP_CALL( SCIPincludeNodeselDfs(scip) );
   cr_assert_not_null(presol);

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
   initProb();
}

/** deinitialization method */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(presol, .init = setup, .fini = teardown);


/* TESTS */
Test(presol, checkPresolGetName)
{
   cr_assert_str_eq("unittest", SCIPpresolGetName(presol));
}

Test(presol, CheckPresolGetDesc)
{
   cr_assert_str_eq(SCIPpresolGetDesc(presol), "presolver to test");
}

Test(presol, checkPresolGetPriority)
{
   cr_assert_eq(SCIPpresolGetPriority(presol), 20010001 );
}

Test(presol, checkPresolIsInitialized)
{
   cr_assert_eq(SCIPpresolIsInitialized(presol), FALSE, "was expecting the presolver not to be initialized");

   SCIP_CALL( SCIPpresolve(scip) );
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   cr_assert_eq(SCIPpresolIsInitialized(presol), TRUE, "was expecting the presolver to be initialized");
}
