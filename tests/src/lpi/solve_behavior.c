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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solve_behavior.c
 * @brief  unit tests for testing the behaviour of the LP solver when calling solving methods and
 *    getting the solving state afterwards.  *
 *    Methods tested are:
 *       SCIPlpiSolveBarier(), SCIPlpiWasSolved()
 *       SCIPlpiHas{Simplex,Barrier}Solve()
 *
 * @author Franziska Schloesser
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <lpi/lpi.h>

#include <signal.h>
#include "include/scip_test.h"

#define EPS 1e-6

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;

/** write ncols, nrows and objsen into variables to check later */
static
SCIP_Bool initProb(int* ncols, int* nrows, int* nnonz, SCIP_OBJSEN* objsen)
{
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   SCIP_Real obj[100] = { 1.0, 1.0 };
   SCIP_Real  lb[100] = { 0.0, 0.0 };
   SCIP_Real  ub[100] = { SCIPlpiInfinity(lpi),  SCIPlpiInfinity(lpi) };
   SCIP_Real lhs[100] = {-SCIPlpiInfinity(lpi), -SCIPlpiInfinity(lpi) };
   SCIP_Real rhs[100] = { 1.0, 1.0 };
   SCIP_Real val[100] = { 1.0, 1.0 };
   int beg[100] = { 0, 1 };
   int ind[100] = { 0, 1 };

   /* optimal - optimal
    * (P):  max x+y
    * x     <= 1 (constr)
    *     y <= 1 (constr)
    *
    *  0 <= x (bound)
    *  0 <= y (bound)
    *
    * (D):  min x+y
    * 1 <= x    (constr)
    * 1 <=    y (constr)
    *
    * 0 <= x (bound)
    * 0 <= y (bound)
    j* */
   *ncols = 2;
   *nrows = 2;
   *nnonz = 2;
   *objsen = SCIP_OBJSEN_MAXIMIZE;

   SCIP_CALL( SCIPlpiChgObjsen(lpi, *objsen) );
   SCIP_CALL( SCIPlpiAddCols(lpi, *ncols, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddRows(lpi, *nrows, lhs, rhs, NULL, *nnonz, beg, ind, val) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   cr_assert( SCIPlpiWasSolved(lpi) );

   return TRUE;
}

/* TEST SUITE */
static
void setup(void)
{
   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
}

static
void teardown(void)
{
   cr_assert( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(solve_behavior, .init = setup, .fini = teardown);

/** Tests */

/** SCIPlpiSolveBarrier should always solve the problem
 *  if the lp solver does not have a barrier method it should default to dualsimplex. */
Test(solve_behavior, testbarriersolve)
{
   int nrows, ncols, nnonz;
   SCIP_OBJSEN sense;
   initProb(&ncols, &nrows, &nnonz, &sense);

   SCIP_CALL( SCIPlpiSolveBarrier(lpi, true) );

   /* Problem is optimal - optimal, should be solved */
   cr_assert( SCIPlpiWasSolved(lpi) );

   SCIP_CALL( SCIPlpiSolveBarrier(lpi, false) );

   /* Problem is optimal - optimal, should be solved */
   cr_assert( SCIPlpiWasSolved(lpi) );
}

/** Test if the two method giving information about availability of solve methods do not crash. */
Test(solve_behavior, testhassolve)
{
   /* try calling all methods at least once */
   SCIPlpiHasSimplexSolve(lpi);
   SCIPlpiHasBarrierSolve(lpi);
}
