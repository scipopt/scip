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

/**@file   verification.c
 * @brief  unit test for testing circular pattern verification
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "probdata_rpa.h"
#include "pricer_rpa.h"
#include "pattern.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_PATTERN* pattern;
static SCIP_PROBDATA* probdata;

/** setup of test run */
static
void setup(void)
{
   SCIP_Real rexts[4] = {2.414213562373095, 1.0, 0.6, 0.5};
   SCIP_Real rints[4] = {2.414213562373095, 1.0, 0.5, 0.0};
   int demands[4] = {1, 5, 5, 5};

   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include ringpacking pricer  */
   SCIP_CALL( SCIPincludePricerRpa(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create problem data */
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 4, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   cr_assert(probdata != NULL);

   /* creates a circular pattern */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, 1) );
}

/** deinitialization method */
static
void teardown(void)
{
   /* release circular pattern */
   SCIPpatternRelease(scip, &pattern);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/* test suite */
TestSuite(verification, .init = setup, .fini = teardown);

/** verifies empty circular pattern with NLP */
Test(verification, nlp_empty)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing a single element with NLP */
Test(verification, nlp_single)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing two elements with NLP */
Test(verification, nlp_two)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 3 and element of type 2 -> not packable */
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPpatternAddElement(pattern, 2, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_NO);

   /* remove element of type 2 and add another element of type 3 -> packable */
   SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);
   SCIPpatternRemoveLastElements(pattern, 1);
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing two elements with NLP */
Test(verification, nlp_complex)
{
   SCIP_PATTERN* p;
   int i;

   /* pattern with inner radius of 1+sqrt(2) */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &p, 0) );

   /* add four circles with radius 1 */
   for( i = 0; i < 4; ++i )
      SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);

   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, p, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_YES);

   /* add one more element -> not packable anymore */
   SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);
   SCIPpatternSetPackableStatus(p, SCIP_PACKABLE_UNKNOWN);
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, p, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_NO);

   /* release pattern */
   SCIPpatternRelease(scip, &p);
}

/** verifies empty circular pattern with heuristic */
Test(verification, heur_empty)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   cr_expect(SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing a single element with heuristic */
Test(verification, heur_single)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verify circular pattern containing two elements with heuristic */
Test(verification, heur_two)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 3 and element of type 2 -> not packable */
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPpatternAddElement(pattern, 2, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* remove element of type 3 and add another element of type 2 -> packable */
   SCIPpatternSetPackableStatus(pattern, SCIP_PACKABLE_UNKNOWN);
   SCIPpatternRemoveLastElements(pattern, 1);
   SCIP_CALL( SCIPpatternAddElement(pattern, 3, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verifies circular pattern containing four and five elements with heuristic */
Test(verification, heuer_four)
{
   SCIP_PATTERN* p;
   int i;

   /* pattern with inner radius of 1+sqrt(2) */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &p, 0) );

   /* add four circles with radius 1 */
   for( i = 0; i < 4; ++i )
      SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);

   /* easy enough for any heuristic */
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, p, SCIPinfinity(scip), 1) );
   cr_expect(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_YES);

   /* add one more element -> not packable anymore (heuristic should only detect UNKNOWN) */
   SCIPpatternAddElement(p, 1, SCIP_INVALID, SCIP_INVALID);
   SCIPpatternSetPackableStatus(p, SCIP_PACKABLE_UNKNOWN);
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, p, SCIPinfinity(scip), 10) );
   cr_expect(SCIPpatternGetPackableStatus(p) == SCIP_PACKABLE_UNKNOWN);

   /* release pattern */
   SCIPpatternRelease(scip, &p);
}
