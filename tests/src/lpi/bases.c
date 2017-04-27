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

/**@file   bases.c
 * @brief  unit test for checking the settings of slack variables in a basis of the lpi
 * @author Marc Pfetsch
 *
 * The behavior of different LP solvers w.r.t. the slack variables should not differ, if interfaced by LPI.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <scip/scip.h>
#include <lpi/lpi.h>
#include "include/scip_test.h"

#define EPS 1e-6

static SCIP_LPI* lpi;
static SCIP_Real obj;
static SCIP_Real vals[2];
static SCIP_Real lhs;
static SCIP_Real rhs;
static SCIP_Real val;
static SCIP_Real objval;
static int ind;


/*** TEST SUITE SIMPLE ***/
static
void setup_simple(void)
{
   int nrows;
   int ncols;
   int beg = 0;
   SCIP_Real lb;
   SCIP_Real ub;

   ind = 0;
   lpi = NULL;
   obj = 1.0;
   lb = 0.0;
   ub = 3.0;
   lhs = 1.0;
   rhs = 2.0;
   val = 1.0;

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* use the following LP as base:
    *   max x
    *       1 <= x <= 2  (linear constraint)
    *       0 <= x <= 3  (bounds)
    */
   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add one row */
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   cr_assert( nrows == 1 );
   cr_assert( ncols == 1 );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}
TestSuite(simple, .init = setup_simple, .fini = teardown);

/*** TESTS ***/
Test(simple, test1)
{
   int cstat;
   int rstat;

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_UPPER );
}

Test(simple, test2)
{
   int cstat;
   int rstat;

   /* modify LP to:
    *   min x
    *       1 <= x <= 2  (linear constraint)
    *       0 <= x <= 3  (bounds)
    */
   /* change sense */
   SCIP_CALL( SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MINIMIZE) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_LOWER );
}

Test(simple, test3)
{
   int cstat;
   int rstat;

   /* modify LP to:
    *   min x
    *       1 <= x       (linear constraint)
    *       0 <= x <= 3  (bounds)
    */
   /* change sense */
   SCIP_CALL( SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MINIMIZE) );

   /* change row side */
   rhs = SCIPlpiInfinity(lpi);
   SCIP_CALL( SCIPlpiChgSides(lpi, 1, &ind, &lhs, &rhs) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_LOWER );
}

Test(simple, test4)
{

   int cstat;
   int rstat;

   /* modify LP to:
    *   max x
    *       x <= 1       (linear constraint)
    *       0 <= x <= 3  (bounds)
    */

   /* change row sides */
   lhs = -SCIPlpiInfinity(lpi);
   rhs = 1.0;
   SCIP_CALL( SCIPlpiChgSides(lpi, 1, &ind, &lhs, &rhs) );

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   /* get basis */
   SCIP_CALL( SCIPlpiGetBase(lpi, &cstat, &rstat) );

   /* the variable should be basic and the slack variable at the upper bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_UPPER );
}

/*** TEST SUITE COMPLEX ***/
static
void setup_complex(void)
{
   int ncols;
   int beg = 0;
   int inds[2];
   int nrows;
   SCIP_Real lb;
   SCIP_Real ub;

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* use the following LP:
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

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   cr_assert_eq(nrows, 3);
   cr_assert_eq(ncols, 3);
}
TestSuite(complex, .init = setup_complex, .fini = teardown);

/*** TESTS ***/
Test(complex, test1)
{
   int i;
   SCIP_Real binvrow[3];
   SCIP_Real coef[3];
   int cstats[3];
   int nrows;
   int rstats[3];
   int basinds[3];

   /* solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   cr_assert_float_eq(objval, 14, EPS);

   /* the optimal basis should be: {x2, x3, slack for second row} */
   SCIP_CALL( SCIPlpiGetBase(lpi, cstats, rstats) );
   cr_assert(cstats[0] == SCIP_BASESTAT_LOWER);
   cr_assert(cstats[1] == SCIP_BASESTAT_BASIC);
   cr_assert(cstats[2] == SCIP_BASESTAT_BASIC);

   cr_assert(rstats[0] == SCIP_BASESTAT_LOWER);
   cr_assert(rstats[1] == SCIP_BASESTAT_BASIC);
   cr_assert(rstats[2] == SCIP_BASESTAT_UPPER);

   /* get basis indices */
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, basinds) );

   /* search for slack variable in basis */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   for (i = 0; i < nrows; ++i)
   {
      if ( basinds[i] < 0 )
         break;
   }
   cr_assert_lt(i, nrows);

   /* check basis inverse */
   SCIP_CALL( SCIPlpiGetBInvRow(lpi, i, binvrow, NULL, NULL) );

   /* row of basis inverse should be (0, 1, 0.5) */
   cr_expect_float_eq(binvrow[0], 0.0, EPS);
   cr_expect_float_eq(binvrow[1], 1.0, EPS);
   cr_expect_float_eq(binvrow[2], 0.5, EPS);

   /* check basis inverse times nonbasic matrix */
   SCIP_CALL( SCIPlpiGetBInvARow(lpi, i, binvrow, coef, NULL, NULL) );

   /* row of basis inverse times nonbasic matrix should be (-0.5, 0, 0) */
   cr_expect_float_eq(coef[0], -0.5, EPS);
   cr_expect_float_eq(coef[1], 0.0, EPS);
   cr_expect_float_eq(coef[2], 0.0, EPS);
}
