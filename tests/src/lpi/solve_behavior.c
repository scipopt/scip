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

/** initialize a 3x3 maximization problem, write nrows and ncols in variables for later use. */
static
void initProb(int* ncols, int* nrows)
{
   int beg = 0;
   int inds[2];
   SCIP_Real vals[2];
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_Real lhs;
   SCIP_Real rhs;
   *ncols = 3;
   *nrows = 3;

   /* clean up after old problem */
   if ( lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&lpi) );
   }

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* use the following LP (same as in bases.c):
    * max 1 x1 + 1 x2 + 1 x3
    *       -8 <= -x1 -          x3 <= -1
    *       -7 <= -x1 -   x2        <= -1
    *              x1 + 2 x2        <= 12
    *              x1,    x2,    x3 >= 0
    */
   /* add columns */
   lb = 0.0;
   ub = SCIPlpiInfinity(lpi);
   obj = 1.0;

   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add rows */
   lhs = -8.0;
   rhs = -1.0;
   inds[0] = 0;
   inds[1] = 2;
   vals[0] = -1.0;
   vals[1] = -1.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals) );

   lhs = -7.0;
   rhs = -1.0;
   inds[0] = 0;
   inds[1] = 1;
   vals[0] = -1.0;
   vals[1] = -1.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals) );

   lhs = -SCIPlpiInfinity(lpi);
   rhs = 12.0;
   inds[0] = 0;
   inds[1] = 1;
   vals[0] = 1.0;
   vals[1] = 2.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 2, &beg, inds, vals) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
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
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(solve_behavior, .init = setup, .fini = teardown);

/** Tests */

/** SCIPlpiSolveBarrier should always solve the problem
 *  if the lp solver does not have a barrier method it should default to dualsimplex. */
Test(solve_behavior, testbarriersolve)
{
   SCIP_Real objval;
   int i;
   bool crossover = true;
   int nrows, ncols;
   SCIP_Real exp_primsol[3] = { 0.0, 6.0, 8.0};
   SCIP_Real exp_dualsol[3] = {-1.0, 0.0, 0.5};
   SCIP_Real* primsol;
   SCIP_Real* dualsol;

   /* do this twice, with and without crossover */
   for (i = 0; i < 2; ++i)
   {
      /* initialize */
      initProb(&nrows, &ncols);

      /* solve problem */
      SCIP_CALL( SCIPlpiSolveBarrier(lpi, crossover) );

      /* very short check (subset of complex test1 in bases.c) */
      cr_assert( SCIPlpiWasSolved(lpi) );
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      cr_assert_float_eq(objval, 14.0, EPS);

      /* allocate storage for solution */
      BMSallocMemoryArray(&primsol, ncols);
      BMSallocMemoryArray(&dualsol, nrows);

      /* get solution */
      SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primsol, dualsol, NULL, NULL) );

      for (i = 0; i < ncols; ++i)
      {
         cr_assert_float_eq(primsol[i], exp_primsol[i], EPS, "Violation of primal solution %d: %g != %g\n", i, primsol[i], exp_primsol[i]);
      }

      for (i = 0; i < nrows; ++i)
      {
         cr_assert_float_eq(dualsol[i], exp_dualsol[i], EPS, "Violation of dual solution %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
      }
      /* free up memory */
      BMSfreeMemoryArray(&primsol);
      BMSfreeMemoryArray(&dualsol);

      /* prepare for second run */
      crossover = !crossover;
   }
}

/** Test if the two method giving information about availability of solve methods do not crash. */
Test(solve_behavior, testhassolve)
{
   /* try calling all methods at least once */
   SCIPlpiHasPrimalSolve();
   SCIPlpiHasDualSolve();
   SCIPlpiHasBarrierSolve();
}
