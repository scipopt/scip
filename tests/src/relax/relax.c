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

/**@file   relax.c
 * @brief  unit test for checking relaxator
 * @author Franziska Schloesser
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* UNIT TEST */

/** relaxator data */
struct SCIP_RelaxData
{
   int ncalls;
};

/**! [SnippetRelaxFreeUnittest] */
static
SCIP_DECL_RELAXFREE(relaxFreeUnittest)
{  /*lint --e{715}*/
   /* call destructor of relaxation handler */

   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   SCIPfreeMemory(scip, &relaxdata);
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}
/**! [SnippetRelaxFreeUnittest] */


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecUnittest)
{
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);
   relaxdata->ncalls++;

   return SCIP_OKAY;
}

#include "include/scip_test.h"

/* GLOBAL VARIABLES */
static SCIP* scip;
static SCIP_RELAX* relax;

/** get the number of calls of the relaxator */
static
int SCIPgetNcallsUnittest(void)
{
   SCIP_RELAXDATA* relaxdata;

   relax = SCIPfindRelax(scip, "unittest");
   assert(relax != NULL);
   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   return relaxdata->ncalls;
}

/** creates the unittest relaxator and includes it in SCIP */
static
void SCIPincludeRelaxUnittest(void)
{
   SCIP_RELAXDATA* relaxdata;

   /* create unittest relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
   relaxdata->ncalls = 0;

   relax = NULL;

   /* use SCIPincludeRelaxBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, "unittest", "relaxator for unittest", 101, 2, relaxExecUnittest, relaxdata) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeUnittest) );

   assert(relax != NULL);
}


/* TEST SUITES */
/** setup of test run */
static
void setup(void)
{
   scip = NULL;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include binpacking reader */
   //SCIP_CALL( SCIPincludeRelaxUnittest(scip) );
   SCIPincludeRelaxUnittest();

   /* create a problem and disable the presolver */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );
}

/** deinitialization method */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}



TestSuite(relax, .init = setup, .fini = teardown);

/* TESTS */

Test(relax, relaxCheckName)
{
   cr_assert_str_eq(SCIPrelaxGetName(relax), "unittest");
}

Test(relax, find)
{
   cr_assert_eq(relax, SCIPfindRelax(scip, "unittest"));
}

Test(relax, relaxCheckDesc)
{
   cr_assert_str_eq(SCIPrelaxGetDesc(relax), "relaxator for unittest");
}

Test(relax, relaxCheckPriority)
{
   cr_assert_eq(SCIPrelaxGetPriority(relax), 101);
}

Test(relax, relaxCheckFreq)
{
   cr_assert_eq(SCIPrelaxGetFreq(relax), 2);
}

/*@todo how to check this? */
Test(relax, relaxCheckSetupTime)
{
   cr_assert_geq(SCIPrelaxGetSetupTime(relax), 0.0);
}

/*@todo how to check this? */
Test(relax, relaxCheckTime)
{
   cr_assert_geq(SCIPrelaxGetTime(relax), 0.0);
}

/*@todo how to check this? */
Test(relax, relaxCheckNCalls)
{
   cr_assert_eq( SCIPrelaxGetNCalls(relax), SCIPgetNcallsUnittest() );
}

Test(relax, relaxCheckInitialized)
{
   cr_assert_eq( SCIPrelaxIsInitialized(relax), FALSE );
   SCIP_CALL( SCIPsolve(scip) );
   cr_assert_eq( SCIPrelaxIsInitialized(relax), TRUE );
}
