/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solve.c
 * @brief  unit test for checking lpi solution
 * @author Marc Pfetsch
 *
 * We perform tests with solving several examples. These are inspired by the unit tests of OSI in COIN-OR.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <lpi/lpi.h>

#include "include/scip_test.h"

#define EPS 1e-6

/** expected feasibility status for primal or dual problem */
enum SCIPfeasStatus
{
   SCIPfeas      = 0,    /**< the problem is feasible */
   SCIPunbounded = 1,    /**< the problem is unbounded */
   SCIPinfeas    = 2     /**< the problem is infeasible */
};
typedef enum SCIPfeasStatus SCIPFEASSTATUS;

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;

/* TEST SUITE */
static
void setup(void)
{
   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) );
#endif
}

static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(solve, .init = setup, .fini = teardown);

/* local functions */

/** perform basic test for the given problem */
static
SCIP_RETCODE performTest(
   SCIP_Bool             solveprimal,        /**< use primal simplex */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val,                /**< values of constraint matrix entries */
   SCIPFEASSTATUS        exp_primalfeas,     /**< expected primal feasibility status */
   SCIPFEASSTATUS        exp_dualfeas,       /**< expected primal feasibility status */
   const SCIP_Real*      exp_primsol,        /**< expected primal optimal solution or primal ray if primal is unbounded or NULL */
   const SCIP_Real*      exp_dualsol,        /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
   const SCIP_Real*      exp_activity,       /**< expected activity of optimal solution or NULL */
   const SCIP_Real*      exp_redcost         /**< expected reduced cost of optimal solution or NULL */
   )
{
   /* solution data */
   SCIP_Real objval;
   SCIP_Real* primsol;
   SCIP_Real* dualsol;
   SCIP_Real* activity;
   SCIP_Real* redcost;

   /* auxiliary data */
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   int ntmprows;
   int ntmpcols;
   int i;
   int j;

   /* load problem */
   SCIP_CALL( SCIPlpiLoadColLP(lpi, objsen, 2, obj, lb, ub, NULL, 2, lhs, rhs, NULL, 4, beg, ind, val) );

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &ntmprows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ntmpcols) );
   cr_assert( nrows == ntmprows );
   cr_assert( ncols == ntmpcols );

   cr_assert( ! SCIPlpiWasSolved(lpi) );

   /* solve problem */
   if ( solveprimal )
   {
      SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   }
   else
   {
      SCIP_CALL( SCIPlpiSolveDual(lpi) );
   }

   /* check status */
   cr_assert( SCIPlpiWasSolved(lpi) );
   cr_assert( ! SCIPlpiIsObjlimExc(lpi) );
   cr_assert( ! SCIPlpiIsIterlimExc(lpi) );
   cr_assert( ! SCIPlpiIsTimelimExc(lpi) );

   /* check feasibility status */
   SCIP_CALL( SCIPlpiGetSolFeasibility(lpi, &primalfeasible, &dualfeasible) );

   /* CPLEX has a strange definition of stability such that we cannot use the test here. */
   /* cr_assert( SCIPlpiIsStable(lpi) ); */

   /* if we are feasible, we should be optimal */
   if ( exp_primalfeas == SCIPfeas && exp_dualfeas == SCIPfeas )
   {
      cr_assert( SCIPlpiIsOptimal(lpi) );
   }

   /* check more primal statuses */
   switch ( exp_primalfeas )
   {
   case SCIPfeas:
      cr_assert( primalfeasible );
      cr_assert( ! SCIPlpiExistsPrimalRay(lpi) );
      cr_assert( ! SCIPlpiHasPrimalRay(lpi) );
      cr_assert( ! SCIPlpiIsPrimalUnbounded(lpi) );
      cr_assert( ! SCIPlpiIsPrimalInfeasible(lpi) );
      cr_assert( SCIPlpiIsPrimalFeasible(lpi) );
      break;

   case SCIPunbounded:
      /* Because of SoPlex, cannot always determine feasibility status here, even if we want to apply the primal
       * simplex. In any case, the results of primalfeasible and SCIPlpiIsPrimalFeasible(lpi) should coincide. */
      cr_assert( primalfeasible == SCIPlpiIsPrimalFeasible(lpi) );

      /* It seems that we cannot guarantee that the primal is shown to be unbounded. */
      /* cr_assert( SCIPlpiIsPrimalUnbounded(lpi) ); */

      cr_assert( SCIPlpiExistsPrimalRay(lpi) );
      cr_assert( ! SCIPlpiIsPrimalInfeasible(lpi) );
      break;

   case SCIPinfeas:
      cr_assert( ! primalfeasible );
      /* It seems that we cannot always prove that primal is infeasible. */
      /* cr_assert( SCIPlpiIsPrimalInfeasible(lpi) ); */

      /* It seems that we cannot always prove that primal is not unbounded. */
      /* cr_assert( ! SCIPlpiIsPrimalUnbounded(lpi) ); */
      cr_assert( ! SCIPlpiIsPrimalFeasible(lpi) );
      break;

   default:
      abort();
   }

   /* check more dual statuses */
   switch ( exp_dualfeas )
   {
   case SCIPfeas:
      cr_assert( dualfeasible );
      cr_assert( ! SCIPlpiExistsDualRay(lpi) );
      cr_assert( ! SCIPlpiHasDualRay(lpi) );
      cr_assert( ! SCIPlpiIsDualUnbounded(lpi) );
      cr_assert( ! SCIPlpiIsDualInfeasible(lpi) );
      cr_assert( SCIPlpiIsDualFeasible(lpi) );
      break;

   case SCIPunbounded:
      /* Because of SoPlex, cannot always determine feasibility status here, even if we want to apply the dual
       * simplex. In any case, the results of dualfeasible and SCIPlpiIsDualFeasible(lpi) should coincide. */
      cr_assert( dualfeasible == SCIPlpiIsDualFeasible(lpi) );

      /* It seems that we cannot guarantee that the dual is shown to be unbounded. */
      /* cr_assert( SCIPlpiIsDualUnbounded(lpi) ); */

      cr_assert( SCIPlpiExistsDualRay(lpi) );
      cr_assert( ! SCIPlpiIsDualInfeasible(lpi) );
      break;

   case SCIPinfeas:
      cr_assert( ! dualfeasible );
      cr_assert( ! SCIPlpiIsDualUnbounded(lpi) );
      /* It seems that we cannot always prove that dual is infeasible. */
      /* cr_assert( SCIPlpiIsDualInfeasible(lpi) ); */
      cr_assert( ! SCIPlpiIsDualFeasible(lpi) );
      break;

   default:
      abort();
   }

   /* allocate storage for solution */
   BMSallocMemoryArray(&primsol, ncols);
   BMSallocMemoryArray(&dualsol, nrows);
   BMSallocMemoryArray(&activity, nrows);
   BMSallocMemoryArray(&redcost, ncols);

   /* get solution */
   SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost) );

   /* check solution */
   if ( exp_primalfeas == SCIPfeas )
   {
      assert( exp_primsol != NULL && exp_redcost != NULL );
      for (j = 0; j < ncols; ++j)
      {
         cr_assert_float_eq(primsol[j], exp_primsol[j], EPS, "Violation of primal solution %d: %g != %g\n", j, primsol[j], exp_primsol[j]);
         cr_assert_float_eq(redcost[j], exp_redcost[j], EPS, "Violation of reduced cost of solution %d: %g != %g\n", j, redcost[j], exp_redcost[j]);
      }
   }
   else if ( exp_primalfeas == SCIPunbounded )
   {
      assert( exp_primsol != NULL );
      if ( SCIPlpiHasPrimalRay(lpi) )
      {
         SCIP_Real scalingfactor = 1.0;

         SCIP_CALL( SCIPlpiGetPrimalRay(lpi, primsol) );

         /* loop until scaling factor can be determined */
         for (j = 0; j < ncols; ++j)
         {
            if ( REALABS(exp_primsol[j]) < EPS )
               cr_assert_float_eq(primsol[j], exp_primsol[j], EPS, "Violation of primal ray %d: %g != %g\n", j, primsol[j], exp_primsol[j]);
            else
            {
               scalingfactor = primsol[j]/exp_primsol[j];
               break;
            }
         }

         /* again loop over ray */
         for (j = 0; j < ncols; ++j)
         {
            cr_assert_float_eq(primsol[j], scalingfactor * exp_primsol[j], EPS, "Violation of primal ray %d: %g != %g\n", j, primsol[j], scalingfactor * exp_primsol[j]);
         }
      }
   }

   if ( exp_dualfeas == SCIPfeas )
   {
      assert( exp_dualsol != NULL && exp_activity != NULL );
      for (i = 0; i < nrows; ++i)
      {
         cr_assert_float_eq(dualsol[i], exp_dualsol[i], EPS, "Violation of dual solution %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
         cr_assert_float_eq(activity[i], exp_activity[i], EPS, "Violation of activity of solution %d: %g != %g\n", i, activity[i], exp_activity[i]);
      }
   }
   else if ( exp_dualfeas == SCIPunbounded )
   {
      assert( exp_dualsol != NULL );
      if ( SCIPlpiHasDualRay(lpi) )
      {
         SCIP_Real scalingfactor = 1.0;

         SCIP_CALL( SCIPlpiGetDualfarkas(lpi, dualsol) );

         /* loop until scaling factor can be determined */
         for (i = 0; i < nrows; ++i)
         {
            if ( REALABS(exp_dualsol[i]) < EPS )
               cr_assert_float_eq(dualsol[i], exp_dualsol[i], EPS, "Violation of dual ray %d: %g != %g\n", i, dualsol[i], exp_dualsol[i]);
            else
            {
               scalingfactor = dualsol[i]/exp_dualsol[i];
               break;
            }
         }

         /* again loop over ray */
         for (i = 0; i < nrows; ++i)
         {
            cr_assert_float_eq(dualsol[i], scalingfactor * exp_dualsol[i], EPS, "Violation of dual ray %d: %g != %g\n", i, dualsol[i], scalingfactor * exp_dualsol[i]);
         }
      }
   }

   BMSfreeMemoryArray(&primsol);
   BMSfreeMemoryArray(&dualsol);
   BMSfreeMemoryArray(&activity);
   BMSfreeMemoryArray(&redcost);

   return SCIP_OKAY;
}


/** TESTS **/

/* Test 1:
 *
 * max 3 x1 +   x2
 *     2 x1 +   x2 <= 10
 *       x1 + 3 x2 <= 15
 *       x1,    x2 >= 0
 *
 * with primal optimal solution (5, 0), dual optimal solution (1.5, 0), activity (10, 5), and redcost (0, -0.5).
 *
 * Then use objective (1, 1) with primal optimal solution (3,4), dual optimal solution (0.4, 0.2), activity (10, 15), and redcost (0, 0).
 */
Test(solve, test1)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {3, 1};
   SCIP_Real lb[2] = {0, 0};
   SCIP_Real rhs[2] = {10, 15};
   int beg[2] = {0, 2};
   int ind[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real ub[2];
   SCIP_Real lhs[2];

   /* expected solutions */
   SCIP_Real exp_primsol[2] = {5, 0};
   SCIP_Real exp_dualsol[2] = {1.5, 0};
   SCIP_Real exp_activity[2] = {10, 5};
   SCIP_Real exp_redcost[2] = {0, -0.5};

   /* fill variable data */
   ub[0] = SCIPlpiInfinity(lpi);
   ub[1] = SCIPlpiInfinity(lpi);
   lhs[0] = -SCIPlpiInfinity(lpi);
   lhs[1] = -SCIPlpiInfinity(lpi);

   /* check problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   /* check problem with dual simplex */
   SCIP_CALL( performTest(FALSE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   /* change objective and expected solution */
   obj[0] = 1.0;
   exp_primsol[0] = 3;
   exp_primsol[1] = 4;
   exp_dualsol[0] = 0.4;
   exp_dualsol[1] = 0.2;
   exp_activity[0] = 10;
   exp_activity[1] = 15;
   exp_redcost[1] = 0;

   /* check changed problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );
}


/* Test 2:
 *
 * max 3 x1 +   x2
 *     2 x1 +   x2 <= 10
 *       x1 + 3 x2 <= 15
 *       x1, x2 free
 *
 * which is unbounded (the only difference to Test 1 is that the variables are free).
 *
 * Then use objective (1, 1) with primal optimal solution (3,4), dual optimal solution (0.4, 0.2), activity (10, 15), and redcost (0, 0).
 */
Test(solve, test2)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {3, 1};
   SCIP_Real rhs[2] = {10, 15};
   int beg[2] = {0, 2};
   int ind[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real lb[2];
   SCIP_Real ub[2];
   SCIP_Real lhs[2];

   /* expected ray for first run */
   SCIP_Real exp_primray[2] = {0.5, -1};

   /* expected solutions for second run */
   SCIP_Real exp_primsol[2] = {3, 4};
   SCIP_Real exp_dualsol[2] = {0.4, 0.2};
   SCIP_Real exp_activity[2] = {10, 15};
   SCIP_Real exp_redcost[2] = {0, 0};

   /* fill variable data */
   lb[0] = -SCIPlpiInfinity(lpi);
   lb[1] = -SCIPlpiInfinity(lpi);
   ub[0] = SCIPlpiInfinity(lpi);
   ub[1] = SCIPlpiInfinity(lpi);
   lhs[0] = -SCIPlpiInfinity(lpi);
   lhs[1] = -SCIPlpiInfinity(lpi);

   /* check problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPunbounded, SCIPinfeas, exp_primray, NULL, NULL, NULL) );

   /* check problem with dual simplex */
   SCIP_CALL( performTest(FALSE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPunbounded, SCIPinfeas, exp_primray, NULL, NULL, NULL) );

   /* change objective */
   obj[0] = 1.0;
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );
}


/* Test 3:
 *
 * min 10 y1 + 15 y2
 *      2 y1 +   y2 == 3
 *        y1 + 3 y2 == 1
 *        y1,    y2 >= 0
 *
 * which is dual unbounded (this is the dual of the problem in Test 2).
 *
 * Then use rhs (1, 1) with primal optimal solution (0.4,0.2), dual optimal solution (3, 4), activity (0, 0), and redcost (0, 0).
 */
Test(solve, test3)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {10, 15};
   SCIP_Real lb[2] = {0, 0};
   SCIP_Real lhs[2] = {3, 1};
   SCIP_Real rhs[2] = {3, 1};
   int beg[2] = {0, 2};
   int ind[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real ub[2];

   /* expected ray */
   SCIP_Real exp_dualray[2] = {0.5, -1};

   /* expected solutions for second run */
   SCIP_Real exp_primsol[2] = {0.4, 0.2};
   SCIP_Real exp_dualsol[2] = {3, 4};
   SCIP_Real exp_activity[2] = {1, 1};
   SCIP_Real exp_redcost[2] = {0, 0};

   /* fill variable data */
   ub[0] = SCIPlpiInfinity(lpi);
   ub[1] = SCIPlpiInfinity(lpi);

   /* check problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPinfeas, SCIPunbounded, NULL, exp_dualray, NULL, NULL) );

   /* check problem with dual simplex */
   SCIP_CALL( performTest(FALSE, SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPinfeas, SCIPunbounded, NULL, exp_dualray, NULL, NULL) );

   /* change objective */
   lhs[0] = 1.0;
   rhs[0] = 1.0;
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );
}


/* Test 4:
 *
 * max x1 + x2
 *     x1 - x2 <= 0
 *   - x1 + x2 <= -1
 *     x1,  x2 free
 *
 * which primal and dual infeasible.
 */
Test(solve, test4)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {1, 1};
   SCIP_Real rhs[2] = {0, -1};
   int beg[2] = {0, 2};
   int ind[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {1, -1, -1, 1};

   /* data to be filled */
   SCIP_Real lhs[2];
   SCIP_Real lb[2];
   SCIP_Real ub[2];

   /* fill variable data */
   lb[0] = -SCIPlpiInfinity(lpi);
   lb[1] = -SCIPlpiInfinity(lpi);
   ub[0] = SCIPlpiInfinity(lpi);
   ub[1] = SCIPlpiInfinity(lpi);
   lhs[0] = -SCIPlpiInfinity(lpi);
   lhs[1] = -SCIPlpiInfinity(lpi);

   /* check problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPinfeas, SCIPinfeas, NULL, NULL, NULL, NULL) );

   /* check problem with dual simplex */
   SCIP_CALL( performTest(FALSE, SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPinfeas, SCIPinfeas, NULL, NULL, NULL, NULL) );
}


/* Test 5: Test objective limit
 *
 * Use second problem from Test 1 and set objective limit.
 *
 * This is quite weak test. For instance SoPlex directly finishes with the optimal solution.
 */
Test(solve, test5)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {1, 1};
   SCIP_Real lb[2] = {0, 0};
   SCIP_Real rhs[2] = {10, 15};
   int beg[2] = {0, 2};
   int ind[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};
   SCIP_Real objval;
   int cstat[2] = {SCIP_BASESTAT_LOWER, SCIP_BASESTAT_LOWER};
   int rstat[2] = {SCIP_BASESTAT_BASIC, SCIP_BASESTAT_BASIC};
   SCIP_Real exp_objval = 5.0;

   /* data to be filled */
   SCIP_Real ub[2];
   SCIP_Real lhs[2];

   /* fill variable data */
   ub[0] = SCIPlpiInfinity(lpi);
   ub[1] = SCIPlpiInfinity(lpi);
   lhs[0] = -SCIPlpiInfinity(lpi);
   lhs[1] = -SCIPlpiInfinity(lpi);

   /* load problem */
   SCIP_CALL( SCIPlpiLoadColLP(lpi, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, NULL, 2, lhs, rhs, NULL, 4, beg, ind, val) );

   /* set objective limit */
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FROMSCRATCH, 1) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_PRESOLVING, 0) );
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, 0.0) );
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_LOBJLIM, 0.0) );

   /* set basis */
   SCIP_CALL( SCIPlpiSetBase(lpi, cstat, rstat) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* check status */
   cr_assert( SCIPlpiWasSolved(lpi) );
   cr_assert( SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsOptimal(lpi) );
   cr_assert( ! SCIPlpiIsIterlimExc(lpi) );
   cr_assert( ! SCIPlpiIsTimelimExc(lpi) );

   /* the objective should be equal to the objective limit */
   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   cr_assert_geq(objval, exp_objval, "Objective value not equal to objective limit: %g != %g\n", objval, exp_objval);
}
