/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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
