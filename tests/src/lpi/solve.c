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

/* global variables */
static SCIP_LPI* lpi = NULL;

/* macro for parameters */
#define SCIP_CALL_PARAM(x) /*lint -e527 */ do                                                   \
{                                                                                               \
   SCIP_RETCODE _restat_;                                                                       \
   if ( (_restat_ = (x)) != SCIP_OKAY && (_restat_ != SCIP_PARAMETERUNKNOWN) )                  \
   {                                                                                            \
      SCIPerrorMessage("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_);  \
      abort();                                                                                  \
   }                                                                                            \
}                                                                                               \
while ( FALSE )

/** setup of test suite */
static
void setup(void)
{
   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

#ifdef SCIP_DEBUG
   /* turn on output */
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) );
#endif
}

/** deinitialization method of test */
static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(solve, .init = setup, .fini = teardown);

/* local functions */

/** solve problem */
static
SCIP_RETCODE solveTest(
   SCIP_Bool             solveprimal,        /**< use primal simplex */
   int                   ncols,              /**< number of columns */
   int                   nrows,              /**< number of rows */
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

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &ntmprows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ntmpcols) );
   cr_assert( nrows == ntmprows );
   cr_assert( ncols == ntmpcols );

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

      /* primal ray should exist if the primal simplex ran */
      cr_assert( ! solveprimal || SCIPlpiExistsPrimalRay(lpi) );
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

      /* dual ray should exist if the dual simplex ran */
      cr_assert( solveprimal || SCIPlpiExistsDualRay(lpi) );
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

   /* check solution */
   if ( exp_primalfeas == SCIPfeas )
   {
      /* get solution */
      SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost) );

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
      /* get solution */
      SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost) );

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
         SCIP_Real* lhs;
         SCIP_Real* rhs;

         /* get lhs/rhs for check of dual ray */
         BMSallocMemoryArray(&lhs, nrows);
         BMSallocMemoryArray(&rhs, nrows);
         SCIP_CALL( SCIPlpiGetSides(lpi, 0, nrows-1, lhs, rhs) );

         /* get dual ray */
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
            cr_assert( ! SCIPlpiIsInfinity(lpi, -lhs[i]) || dualsol[i] <= -EPS );
            cr_assert( ! SCIPlpiIsInfinity(lpi, rhs[i]) || dualsol[i] >= EPS );
         }

         BMSfreeMemoryArray(&rhs);
         BMSfreeMemoryArray(&lhs);
      }
   }

   BMSfreeMemoryArray(&primsol);
   BMSfreeMemoryArray(&dualsol);
   BMSfreeMemoryArray(&activity);
   BMSfreeMemoryArray(&redcost);

   return SCIP_OKAY;
}

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
   /* load problem */
   SCIP_CALL( SCIPlpiLoadColLP(lpi, objsen, 2, obj, lb, ub, NULL, 2, lhs, rhs, NULL, 4, beg, ind, val) );
   cr_assert( ! SCIPlpiWasSolved(lpi) );

   /* solve problem */
   SCIP_CALL( solveTest(solveprimal, ncols, nrows, exp_primalfeas, exp_dualfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   return SCIP_OKAY;
}

/** check whether data in LP solver aggrees with original data */
static
SCIP_RETCODE checkData(
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
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   SCIP_OBJSEN lpiobjsen;
   SCIP_Real* lpival;
   SCIP_Real* lpilb;
   SCIP_Real* lpiub;
   SCIP_Real* lpiobj;
   SCIP_Real* lpilhs;
   SCIP_Real* lpirhs;
   int* lpibeg;
   int* lpiind;
   int lpincols;
   int lpinrows;
   int lpinnonz;
   int lpinnonz2;
   int i;
   int j;

   /* check number of rows and columns */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &lpinrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &lpincols) );
   cr_assert( lpinrows == nrows );
   cr_assert( lpincols == ncols );

   /* check objective sense */
   SCIP_CALL( SCIPlpiGetObjsen(lpi, &lpiobjsen) ) ;
   cr_assert( objsen == lpiobjsen );

   /* get number of nonzeros in matrix */
   SCIP_CALL( SCIPlpiGetNNonz(lpi, &lpinnonz) );
   cr_assert( lpinnonz == nnonz );

   /* allocate storage for data */
   BMSallocMemoryArray(&lpilb, ncols);
   BMSallocMemoryArray(&lpiub, ncols);
   BMSallocMemoryArray(&lpibeg, ncols);
   BMSallocMemoryArray(&lpiind, lpinnonz);
   BMSallocMemoryArray(&lpival, lpinnonz);
   BMSallocMemoryArray(&lpiobj, ncols);

   /* get matrix data */
   SCIP_CALL( SCIPlpiGetCols(lpi, 0, ncols-1, lpilb, lpiub, &lpinnonz2, lpibeg, lpiind, lpival) );
   SCIP_CALL( SCIPlpiGetObj(lpi, 0, ncols-1, lpiobj) );

   /* compare data */
   for (j = 0; j < ncols; ++j)
   {
      cr_assert_float_eq(lpilb[j], lb[j], EPS, "Violation of lower bound %d: %g != %g\n", j, lpilb[j], lb[j]);
      cr_assert_float_eq(lpiub[j], ub[j], EPS, "Violation of upper bound %d: %g != %g\n", j, lpiub[j], ub[j]);

      cr_assert_float_eq(lpiobj[j], obj[j], EPS, "Violation of objective coefficient %d: %g != %g\n", j, lpiobj[j], obj[j]);

      cr_assert( lpibeg[j] == beg[j] );
   }

   /* compare matrix */
   for (j = 0; j < nnonz; ++j)
   {
      cr_assert( lpiind[j] == ind[j] );
      cr_assert_float_eq(lpival[j], val[j], EPS, "Violation of matrix entry (%d, %d): %g != %g\n", ind[j], j, lpival[j], val[j]);
   }

   BMSfreeMemoryArray(&lpiobj);
   BMSfreeMemoryArray(&lpival);
   BMSfreeMemoryArray(&lpiind);
   BMSfreeMemoryArray(&lpibeg);
   BMSfreeMemoryArray(&lpiub);
   BMSfreeMemoryArray(&lpilb);

   /* compare lhs/rhs */
   BMSallocMemoryArray(&lpilhs, nrows);
   BMSallocMemoryArray(&lpirhs, nrows);

   SCIP_CALL( SCIPlpiGetSides(lpi, 0, nrows-1, lpilhs, lpirhs) );

   for (i = 0; i < nrows; ++i)
   {
      cr_assert_float_eq(lpilhs[i], lhs[i], EPS, "Violation of lhs %d: %g != %g\n", i, lpilhs[i], lhs[i]);
      cr_assert_float_eq(lpilhs[i], lhs[i], EPS, "Violation of rhs %d: %g != %g\n", i, lpirhs[i], rhs[i]);
   }

   BMSfreeMemoryArray(&lpirhs);
   BMSfreeMemoryArray(&lpilhs);

   return SCIP_OKAY;
}


/** TESTS **/

/** Test 1
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

   /* solve problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* solve problem with dual simplex */
   SCIP_CALL( solveTest(FALSE, 2, 2, SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* change objective */
   obj[0] = 1.0;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, ind, obj) );

   /* change expected solution */
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

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );
}


/** Test 2
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

   /* solve problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPunbounded, SCIPinfeas, exp_primray, NULL, NULL, NULL) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* solve problem with dual simplex */
   SCIP_CALL( solveTest(FALSE, 2, 2, SCIPunbounded, SCIPinfeas, exp_primray, NULL, NULL, NULL) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* change objective */
   obj[0] = 1.0;
   SCIP_CALL( SCIPlpiChgObj(lpi, 1, ind, obj) );

   /* solve with primal simplex */
   SCIP_CALL( solveTest(TRUE, 2, 2, SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );
}


/** Test 3
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

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* check problem with dual simplex */
   SCIP_CALL( solveTest(FALSE, 2, 2, SCIPinfeas, SCIPunbounded, NULL, exp_dualray, NULL, NULL) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* change lhs/rhs */
   lhs[0] = 1.0;
   rhs[0] = 1.0;
   SCIP_CALL( SCIPlpiChgSides(lpi, 1, ind, lhs, rhs) );
   SCIP_CALL( solveTest(TRUE, 2, 2, SCIPfeas, SCIPfeas, exp_primsol, exp_dualsol, exp_activity, exp_redcost) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );
}


/** Test 4
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

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* check problem with dual simplex */
   SCIP_CALL( solveTest(FALSE, 2,  2, SCIPinfeas, SCIPinfeas, NULL, NULL, NULL, NULL) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );
}


/** Test 5: Test objective limit
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
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FROMSCRATCH, 1) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_PRESOLVING, 0) );
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, 0.0) );

   /* set basis */
   SCIP_CALL( SCIPlpiSetBase(lpi, cstat, rstat) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolveDual(lpi) );

   /* check status */
   cr_assert( SCIPlpiWasSolved(lpi) );
   cr_assert( SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsOptimal(lpi) );
   cr_assert( ! SCIPlpiIsIterlimExc(lpi) );
   cr_assert( ! SCIPlpiIsTimelimExc(lpi) );

   /* the objective should be equal to the objective limit */
   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   cr_assert_geq(objval, exp_objval, "Objective value not equal to objective limit: %g != %g\n", objval, exp_objval);

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MAXIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );
}


/** Test 6: More complex example
 *
 * The original problem was the following (arising from the qpkktref unit test), which displays a bug in CPLEX 12.7.0
 * w.r.t. scaling:
 *   Minimize t_objvar
 *   Subject To
 *     KKTBinary1_y:                 - t_dual_y_bin1 + t_dual_y_bin2 + t_dual_y_slackbin1 = 0
 *     KKTlin_lower_1:               - t_x - t_y + t_slack_lhs_lower + t_slack_ub_z       = 0.75
 *     KKTBinary1_x:                 - t_dual_x_bin1 + t_dual_x_bin2 + t_dual_x_slackbin1 = 0
 *     KKTlin_lower_0:               - t_x - t_y + t_slack_ub_z - t_slack_rhs_lower       = 0.25
 *     quadratic_side1_estimation_0: 2.75 t_x - 3.75 t_y + t_objvar + 2.28 t_slack_ub_z  <= 5.0496
 *     quadratic_side0_estimation_0: 1.25 t_x - 0.25 t_y + t_objvar + 2 t_slack_ub_z     >= 2.6875
 *     quadratic_side1_estimation_0: 0.75 t_x - 0.25 t_y + t_objvar + 0.68 t_slack_ub_z  <= 4.2056
 *     quadratic_side0_estimation_0: 2.75 t_x - 0.25 t_y + t_objvar + 3 t_slack_ub_z     >= 4.4375
 *   Bounds
 *     t_x = 1
 *     t_y = 0
 *     -2.562500001 <= t_objvar <= -0.0624999989999999
 *     0 <= t_slack_lhs_lower <= 0.5
 *     t_dual_x_bin1 Free
 *     t_dual_x_bin2 Free
 *     t_dual_x_slackbin1 = 0
 *     t_dual_y_bin1 = 0
 *     t_dual_y_bin2 Free
 *     t_dual_y_slackbin1 Free
 *     1.25 <= t_slack_ub_z <= 1.75
 *     0 <= t_slack_rhs_lower <= 0.5
 *   End
 *
 *  We use the following mapping between variables and indices:
 *  0:  t_x = 1
 *  1:  t_y = 0
 *  2:  t_objvar
 *  3:  t_slack_lhs_lower
 *  4:  t_dual_x_bin1
 *  5:  t_dual_x_bin2
 *  6:  t_dual_x_slackbin1
 *  7:  t_dual_y_bin1
 *  8:  t_dual_y_bin2
 *  9:  t_dual_y_slackbin1
 *  10: t_slack_ub_z
 *  11: t_slack_rhs_lower
 */
Test(solve, test6)
{
   SCIP_Real exp_objval = -2.0625;
   SCIP_Real objval;
   int varind[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

   /* LP data: */
   SCIP_Real obj[12] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   SCIP_Real lb[12] = {1, 0, -2.5625, 0, -1e20, -1e20, 0.0, 0.0, -1e20, -1e20, 1.25, 0};
   SCIP_Real ub[12] = {1, 0, -0.0625, 0.5, 1e20, 1e20, 0.0, 0.0, 1e20, 1e20, 1.75, 0.5};

   SCIP_Real lhs[8] = {0, 0.75, 0, 0.25, -1e20, 2.6875, -1e20, 4.4375};
   SCIP_Real rhs[8] = {0, 0.75, 0, 0.25, 5.0496, 1e20, 4.2056, 1e20};

   /* matrix */
   int beg[12] = {0, 6, 12, 16, 17, 18, 19, 20, 21, 22, 23, 29};
   /*             x0                x1                x2          x3 x4 x5 x6 x7 x8 x9 x10               x11 */
   int ind[30] = {1, 3, 4, 5, 6, 7, 1, 3, 4, 5, 6, 7, 4, 5, 6, 7, 1, 2, 2, 2, 0, 0, 0, 1, 3, 4, 5, 6, 7, 3};
   /*                   x0                              x1                                  x2          x3 x4  x5 x6 x7  x8 x9 x10                     x11 */
   SCIP_Real val[30] = {-1, -1, 2.75, 1.25, 0.75, 2.75, -1, -1, -3.75, -0.25, -0.25, -0.25, 1, 1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 2.28, 2, 0.68, 3, -1.0};
   int j;

   /* possibly convert |1e20| to infinity of LPI */
   for (j = 0; j < 12; ++j)
   {
      if ( lb[j] == -1e20 )
         lb[j] = -SCIPlpiInfinity(lpi);
      if ( ub[j] == 1e20 )
         ub[j] = SCIPlpiInfinity(lpi);
   }
   for (j = 0; j < 8; ++j)
   {
      if ( lhs[j] == -1e20 )
         lhs[j] = -SCIPlpiInfinity(lpi);
      if ( rhs[j] == 1e20 )
         rhs[j] = SCIPlpiInfinity(lpi);
   }

   /* load problem */
   SCIP_CALL( SCIPlpiLoadColLP(lpi, SCIP_OBJSEN_MINIMIZE, 12, obj, lb, ub, NULL, 8, lhs, rhs, NULL, 30, beg, ind, val) );

   /* set some parameters - simulate settings in SCIP */
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FROMSCRATCH, 0) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_SCALING, 1) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_PRESOLVING, 1) );
   SCIP_CALL_PARAM( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_PRICING, 0) );

   SCIP_CALL_PARAM( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, 1e-06) );
   SCIP_CALL_PARAM( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_DUALFEASTOL, 1e-07) );

   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* set objlimit */
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, 4.320412501) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolveDual(lpi) );

   /* check status */
   cr_assert( SCIPlpiWasSolved(lpi) );

   /* the objective should be equal to the objective limit */
   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   cr_assert_geq(objval, exp_objval, "Objective value not equal to objective limit: %g != %g\n", objval, exp_objval);

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 12, obj, lb, ub, 8, lhs, rhs, 30, beg, ind, val) );

   /* change some bounds */
   lb[0] = 1;
   ub[0] = 1;
   lb[1] = 0;
   ub[1] = 0;
   lb[2] = -2.06255;
   ub[2] = -2.0625;
   lb[3] = 0;
   ub[3] = 4.94694e-05;
   lb[6] = 0;
   ub[6] = 0;
   lb[7] = 0;
   ub[7] = 0;
   lb[10] = 1.74995;
   ub[10] = 1.750;
   lb[11] = 0.499951;
   ub[11] = 0.5;
   SCIP_CALL( SCIPlpiChgBounds(lpi, 12, varind, lb, ub) );

   /* set objlimit */
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_OBJLIM, -2.0625) );

   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolveDual(lpi) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 12, obj, lb, ub, 8, lhs, rhs, 30, beg, ind, val) );
}


/** Test 7
 *
 *  min 10 x1 + 15 x2
 *       2 x1 +   x2 >= 3
 *         x1 + 3 x2 <= 1
 *         x1,    x2 >= 0
 *
 *  which is dual unbounded (this is a variant of Test 3 in which the equations have been replaced by inequalities).
 *
 *  The dual is:
 *  max  3 y1 +   y2
 *       2 y1 +   y2 <= 10
 *         y1 + 3 y2 <= 15
 *         y1 >= 0, y2 <= 0
 */
Test(solve, test7)
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

   /* fill data */
   rhs[0] = SCIPlpiInfinity(lpi);
   lhs[1] = -SCIPlpiInfinity(lpi);
   ub[0] = SCIPlpiInfinity(lpi);
   ub[1] = SCIPlpiInfinity(lpi);

   /* check problem with primal simplex */
   SCIP_CALL( performTest(TRUE, SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val,
         SCIPinfeas, SCIPunbounded, NULL, exp_dualray, NULL, NULL) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );

   /* clear basis status */
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* check problem with dual simplex */
   SCIP_CALL( solveTest(FALSE, 2, 2, SCIPinfeas, SCIPunbounded, NULL, exp_dualray, NULL, NULL) );

   /* check that data stored in lpi is still the same */
   SCIP_CALL( checkData(SCIP_OBJSEN_MINIMIZE, 2, obj, lb, ub, 2, lhs, rhs, 4, beg, ind, val) );
}
