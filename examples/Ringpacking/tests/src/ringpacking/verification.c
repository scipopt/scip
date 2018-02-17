/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
   SCIP_Real rexts[3] = {1.0, 0.5, 0.6};
   SCIP_Real rints[3] = {1.0, 0.0, 0.5};
   int demands[3] = {100, 100, 100};

   /* initialize SCIP */
   scip = NULL;
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include ringpacking pricer  */
   SCIP_CALL( SCIPincludePricerRingpacking(scip) );

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create problem data */
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 3, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   cr_assert(probdata != NULL);

   /* creates a circular pattern */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &pattern, 0) );
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

/** verify empty circular pattern with NLP */
Test(verification, nlp_empty)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIP_PACKABLE_YES);
}

/** verify circular pattern containing a single element with NLP */
Test(verification, nlp_single)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verify circular pattern containing two elements with NLP */
Test(verification, nlp_two)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 and element of type 2 -> not packable */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPpatternAddElement(pattern, 2, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_NO);

   /* remove element of type 2 and add another element of type 1 -> packable */
   SCIPpatternRemoveLastElement(pattern);
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternNLP(scip, probdata, pattern, SCIPinfinity(scip), SCIP_LONGINT_MAX) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verify empty circular pattern with heuristic */
Test(verification, heur_empty)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   cr_expect(SCIP_PACKABLE_YES);
}

/** verify circular pattern containing a single element with heuristic */
Test(verification, heur_single)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES);
}

/** verify circular pattern containing a single element with heuristic */
Test(verification, heur_two)
{
   cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_UNKNOWN);

   /* add element of type 1 and element of type 2 -> not packable */
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPpatternAddElement(pattern, 2, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );

   /* TODO */
   /* cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_NO); */

   /* remove element of type 2 and add another element of type 1 -> packable */
   SCIPpatternRemoveLastElement(pattern);
   SCIP_CALL( SCIPpatternAddElement(pattern, 1, SCIP_INVALID, SCIP_INVALID) );
   SCIP_CALL( SCIPverifyCircularPatternHeuristic(scip, probdata, pattern, SCIPinfinity(scip), 1) );

   /* TODO */
   /* cr_expect(SCIPpatternGetPackableStatus(pattern) == SCIP_PACKABLE_YES); */
}
