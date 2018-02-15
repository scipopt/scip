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
 *       SCIPlpiSolveDual(), SCIPlpiWasSolved()
 *       SCIPlpiHas{Primal,Dual,Barrier}Solve()
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
SCIP_Bool initProb(int pos, int* ncols, int* nrows, int* nnonz, SCIP_OBJSEN* objsen)
{
   /* this is necessary because the theories don't setup and teardown after each call, but once before and after */
   if ( lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&lpi) );
   }
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   SCIP_Real obj[100] = { 1.0, 1.0 };
   SCIP_Real  lb[100] = { 0.0, 0.0 };
   SCIP_Real  ub[100] = { SCIPlpiInfinity(lpi),  SCIPlpiInfinity(lpi) };
   SCIP_Real lhs[100] = {-SCIPlpiInfinity(lpi), -SCIPlpiInfinity(lpi) };
   SCIP_Real rhs[100] = { 1.0, 1.0 };
   SCIP_Real val[100] = { 1.0, 1.0 };
   int beg[100] = { 0, 1 };
   int ind[100] = { 0, 1 };

   switch ( pos )
   {

   case 0:
      /* maximization problems, ncols is 2, *nrows is 2 */
      /* unbounded - infeasible
       * (P):  max x+y
       * -x    <= 1 (constr)
       *    -y <= 1 (constr)
       *
       *  0 <= x (bound)
       *  0 <= y (bound)
       *
       * (D):  min x+y
       * 1 <= -x   (constr)
       * 1 <=   -y (constr)
       *
       * 0 <= x (bound)
       * 0 <= y (bound)
       * */
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      val[0] = -1.0;
      val[1] = -1.0;
      break;

   case 1:
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
      break;

   case 2:
      /* infeasible - infeasible
       * (P):  max x+y
       * -x    <= -1 (constr)
       *     y <= -1 (constr)
       *
       *  0 <= x (bound)
       *  0 <= y (bound)
       *
       * (D):  min -x-y
       * 1 <= -x    (constr)
       * 1 <=     y (constr)
       *
       * 0 <= x (bound)
       * 0 <= y (bound)
       */
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      rhs[0] = -1.0;
      rhs[1] = -1.0;
      val[0] = -1.0;
      break;

   case 3:
      /* minimization problems (duals of the above) */
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MINIMIZE;
      rhs[0] = SCIPlpiInfinity(lpi);
      rhs[1] = SCIPlpiInfinity(lpi);
      lhs[0] = 1.0;
      lhs[1] = 1.0;
      val[0] = -1.0;
      val[1] = -1.0;
      break;

   case 4:
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MINIMIZE;
      rhs[0] = SCIPlpiInfinity(lpi);
      rhs[1] = SCIPlpiInfinity(lpi);
      lhs[0] = 1.0;
      lhs[1] = 1.0;
      break;

   case 5:
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MINIMIZE;
      rhs[0] = SCIPlpiInfinity(lpi);
      rhs[1] = SCIPlpiInfinity(lpi);
      lhs[0] = 1.0;
      lhs[1] = 1.0;
      obj[0] = -1.0;
      obj[1] = -1.0;
      val[0] = -1.0;
      break;

   default:
      return FALSE;
   }

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

/** Test */

/** Test adding zero coeffs in rows, expecting an assert in debug mode */
/* This test should fail with an assert from the LPI, which causes SIGABRT to be issued. Thus, this test should pass. */
Test(solve_behavior, testbarriersolve, .signal = SIGABRT)
{
   int nrows, ncols, nnonz;
   SCIP_OBJSEN sense;
   cr_assume( initProb(4, &ncols, &nrows, &nnonz, &sense) );
   cr_assume_eq( 3, nrows );
   cr_assume_eq( 3, ncols );

   SCIP_CALL( SCIPlpiSolveBarrier(lpi, true) );
   SCIP_CALL( SCIPlpiSolveBarrier(lpi, false) );

   /* Problem is optimal - optimal, should be solved */
   cr_assert( SCIPlpiWasSolved(lpi) );

   /* Does it make sense not to have a primal solve as an lp solver? */
   cr_assert( SCIPlpiHasPrimalSolve() );
   /* try calling all methods at least once */
   SCIPlpiHasDualSolve();

   if( !SCIPlpiHasBarrierSolve() )
   {
      /* make sure problem was solved with dual */
      if ( !SCIPlpiHasDualSolve() )
      {
         /* make sure problem was solved with primal */
      }
   }
}
