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

/**@file   enum.c
 * @brief  unit test for testing circular pattern enumeration
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "pricer_rpa.h"
#include "probdata_rpa.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_PROBDATA* probdata;

/** setup of test run */
static
void setup(void)
{
   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include ringpacking pricer  */
   SCIP_CALL( SCIPincludePricerRpa(scip) );
}

/** deinitialization method */
static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/* test suite */
TestSuite(enumerate, .init = setup, .fini = teardown);

/** single empty pattern */
Test(enumerate, empty)
{
   SCIP_Real rexts[1] = {1.0};
   SCIP_Real rints[1] = {0.0};
   int demands[1] = {100};
   SCIP_PATTERN** patterns;
   int npatterns;

   /* create problem data */
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 1, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   cr_assert(probdata != NULL);

   /* compute circular patterns */
   SCIP_CALL( SCIPprobdataEnumeratePatterns(scip, probdata, SCIPinfinity(scip), SCIPinfinity(scip), SCIPinfinity(scip),
      SCIP_LONGINT_MAX, INT_MAX) );

   /* get circular pattern information */
   SCIPprobdataGetCInfos(probdata, &patterns, NULL, &npatterns);
   cr_assert(patterns != NULL);
   cr_expect(npatterns == 1);
   cr_expect(SCIPpatternGetCircleType(patterns[0]) == 0);
   cr_expect(SCIPpatternGetNElemens(patterns[0]) == 0);
}

/** two ring types */
Test(enumerate, two)
{
   SCIP_Real rexts[2] = {1.2, 0.5};
   SCIP_Real rints[2] = {1.0, 0.0};
   int demands[2] = {100, 100};
   SCIP_PATTERN** patterns;
   int npatterns;

   /* create problem data */
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 2, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   cr_assert(probdata != NULL);

   /* compute circular patterns */
   SCIP_CALL( SCIPprobdataEnumeratePatterns(scip, probdata, SCIPinfinity(scip), SCIPinfinity(scip), SCIPinfinity(scip),
      SCIP_LONGINT_MAX, INT_MAX) );

   /* get circular pattern information */
   SCIPprobdataGetCInfos(probdata, &patterns, NULL, &npatterns);
   cr_assert(patterns != NULL);
   cr_expect(npatterns == 2);
   cr_expect(SCIPpatternGetCircleType(patterns[0]) == 0);
   cr_expect(SCIPpatternGetNElemens(patterns[0]) == 2);
   cr_expect(SCIPpatternGetElementType(patterns[0], 0) == 1);
   cr_expect(SCIPpatternGetElementType(patterns[0], 1) == 1);
   cr_expect(SCIPpatternGetCircleType(patterns[1]) == 1);
   cr_expect(SCIPpatternGetNElemens(patterns[1]) == 0);
}

/** three ring types */
Test(enumerate, three)
{
   SCIP_Real rexts[3] = {3.0, 1.0, 0.45};
   SCIP_Real rints[3] = {2.414213562373095, 0.9, 0.0};
   int demands[3] = {100, 3, 3};
   SCIP_PATTERN** patterns;
   int npatterns;

   /* create problem data */
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 3, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   cr_assert(probdata != NULL);

   /* compute circular patterns */
   SCIP_CALL( SCIPprobdataEnumeratePatterns(scip, probdata, SCIPinfinity(scip), SCIPinfinity(scip), SCIPinfinity(scip),
      SCIP_LONGINT_MAX, 10) );

   /* get circular pattern information */
   SCIPprobdataGetCInfos(probdata, &patterns, NULL, &npatterns);
   cr_assert(patterns != NULL);
   cr_expect(npatterns == 3);
   cr_expect(SCIPpatternGetCircleType(patterns[0]) == 0);
   cr_expect(SCIPpatternGetNElemens(patterns[0]) == 6);
   cr_expect(SCIPpatternGetCircleType(patterns[1]) == 1);
   cr_expect(SCIPpatternGetNElemens(patterns[1]) == 2);
   cr_expect(SCIPpatternGetCircleType(patterns[2]) == 2);
   cr_expect(SCIPpatternGetNElemens(patterns[2]) == 0);
}
