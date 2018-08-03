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

/**@file   pattern.c
 * @brief  unit test for testing pattern interface functions
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
static SCIP_PROBDATA* probdata;
static SCIP_PATTERN* rpattern;
static SCIP_PATTERN* cpattern;

/** setup of test run */
static
void setup(void)
{
   SCIP_Real rexts[3] = {1.0, 0.6, 0.5};
   SCIP_Real rints[3] = {1.0, 0.5, 0.0};
   int demands[3] = {100, 100, 100};

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
   SCIP_CALL( SCIPprobdataCreate(scip, "unit test", demands, rints, rexts, 3, 100.0, 100.0) );
   probdata = SCIPgetProbData(scip);
   cr_assert(probdata != NULL);

   /* creates circular and rectangular pattern */
   SCIP_CALL( SCIPpatternCreateCircular(scip, &cpattern, 1) );
   SCIP_CALL( SCIPpatternCreateRectangular(scip, &rpattern) );
}

/** deinitialization method */
static
void teardown(void)
{
   /* release patterns */
   SCIPpatternRelease(scip, &rpattern);
   SCIPpatternRelease(scip, &cpattern);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_null(scip);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

/* test suite */
TestSuite(pattern, .init = setup, .fini = teardown);

/* checks the pattern */
Test(pattern, patterntype)
{
   cr_expect(SCIPpatternGetPatternType(cpattern) == SCIP_PATTERNTYPE_CIRCULAR);
   cr_expect(SCIPpatternGetPatternType(rpattern) == SCIP_PATTERNTYPE_RECTANGULAR);
}

/* checks the type of a circular pattern */
Test(pattern, type)
{
   cr_expect(SCIPpatternGetCircleType(cpattern) == 1);
}

/* checks the position of an element */
Test(pattern, position)
{
   SCIP_CALL( SCIPpatternAddElement(cpattern, 0, -1.0, 1.0) );
   cr_expect(SCIPpatternGetElementPosX(cpattern, 0) == -1.0);
   cr_expect(SCIPpatternGetElementPosY(cpattern, 0) == 1.0);

   SCIP_CALL( SCIPpatternAddElement(cpattern, 0, -2.0, 2.0) );
   cr_expect(SCIPpatternGetElementPosX(cpattern, 1) == -2.0);
   cr_expect(SCIPpatternGetElementPosY(cpattern, 1) == 2.0);
}

/* checks the packable status */
Test(pattern, packable)
{
   cr_expect(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_UNKNOWN);

   SCIPpatternSetPackableStatus(rpattern, SCIP_PACKABLE_YES);
   cr_expect(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_YES);

   /* adding an element does not change packable status */
   SCIP_CALL( SCIPpatternAddElement(rpattern, 0, 0.0, 0.0) );
   cr_expect(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_YES);

   /* removing an element does not change packable status */
   SCIPpatternRemoveLastElements(rpattern, 1);
   cr_expect(SCIPpatternGetPackableStatus(rpattern) == SCIP_PACKABLE_YES);
}
