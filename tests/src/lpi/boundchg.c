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

/**@file   boundchg.c
 * @brief  unit test for checking bound changes
 * @author Marc Pfetsch
 *
 * We perform two tests:
 * - We change the bounds by a very small value and check whether this has an effect on the LP solver interface.
 * - We fix a variable to infinity.
 *
 * In both cases it is unclear what happens. Also LP-solvers might react differently. CPLEX simply ignores small bound changes.
 *
 * For fixing variables to infinity the following variants happen:
 * - SoPlex allows fixing of variables to infinity.
 * - Gurobi returns an error.
 * - CPLEX simply ignores such changes.
 *
 * These tests can be used for debugging or checking the behavior of LP-solvers.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <lpi/lpi.h>

#include "include/scip_test.h"

#define EPS 1e-6

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;

/* TEST SUITE */
static
void setup(void)
{
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 1.0;

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

}
static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

TestSuite(boundchg, .init = setup, .fini = teardown);

/** TESTS **/
Test(boundchg, simple_bound_test)
{
   SCIP_Real lb = 1.0;
   SCIP_Real ub = 2.0;
   SCIP_Real lbnew;
   SCIP_Real ubnew;
   int ind = 0;

   /* change bounds to some value */
   SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );

   /* get bounds and compare */
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 0, &lbnew, &ubnew) );

   cr_assert_float_eq(lb, lbnew, EPS, "Violation of lower bounds: %g != %g\n", lb, lbnew);
   cr_assert_float_eq(ub, ubnew, EPS, "Violation of upper bounds: %g != %g\n", ub, ubnew);
}

Test(boundchg, change_bound_by_small_value)
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lbnew;
   SCIP_Real ubnew;
   int ind = 0;

   /* change bound to small value */
   lb = 1e-11;
   ub = 1.0 - 1e-11;
   SCIP_CALL( SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub) );

   /* get bounds and compare */
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, 0, &lbnew, &ubnew) );

   cr_assert_float_eq(lb, lbnew, EPS, "Violation of lower bounds: %g != %g\n", lb, lbnew);
   cr_assert_float_eq(ub, ubnew, EPS, "Violation of upper bounds: %g != %g\n", ub, ubnew);
}

Test(boundchg, fix_to_infinity)
{
   SCIP_RETCODE retcode;
   SCIP_Real lb;
   SCIP_Real ub;
   int ind = 0;

   /* try to fix variables to infinity */
   lb = SCIPlpiInfinity(lpi);
   ub = SCIPlpiInfinity(lpi);

   /* calling should return an LPERROR */
   retcode = SCIPlpiChgBounds(lpi, 1, &ind, &lb, &ub);

   cr_assert(retcode == SCIP_LPERROR, "Fixing variables to infinity does not return an error.");
}
