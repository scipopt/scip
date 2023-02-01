/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bases.c
 * @brief  unit test for checking the settings of slack variables in a basis of the lpi
 * @author Marc Pfetsch
 * @author Franziska Schloesser
 * @author Felipe Serrano
 *
 * The behavior of different LP solvers w.r.t. the slack variables should not differ, if interfaced by LPI.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include <scip/scip.h>
#include <lpi/lpi.h>
#include "include/scip_test.h"

#define EPS 1e-6

/* global variable for LPI */
static SCIP_LPI* lpi;


/*** TEST SUITE SIMPLE ***/
static
void setup_simple(void)
{
   int nrows;
   int ncols;
   int beg = 0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 3.0;
   SCIP_Real lhs = 1.0;
   SCIP_Real rhs = 2.0;
   SCIP_Real obj = 1.0;
   SCIP_Real val = 1.0;
   int ind = 0;

   lpi = NULL;

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

#ifdef SCIP_DEBUG
   /* turn on output */
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) );
#endif
}

static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
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

   /* the variable should be basic and the slack variable at the lower bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_LOWER );
}

Test(simple, test3)
{
   int cstat;
   int rstat;
   SCIP_Real lhs = 1.0;
   SCIP_Real rhs;
   int ind = 0;

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

   /* the variable should be basic and the slack variable at the lower bound */
   assert( cstat == SCIP_BASESTAT_BASIC );
   assert( rstat == SCIP_BASESTAT_LOWER );
}

Test(simple, test4)
{
   int cstat;
   int rstat;
   SCIP_Real lhs;
   SCIP_Real rhs = 1.0;
   int ind = 0;

   /* modify LP to:
    *   max x
    *       x <= 1       (linear constraint)
    *       0 <= x <= 3  (bounds)
    */

   /* change row sides */
   lhs = -SCIPlpiInfinity(lpi);
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
   SCIP_Real vals[2];
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_Real lhs;
   SCIP_Real rhs;

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

#ifdef SCIP_DEBUG
   /* turn on output */
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) );
#endif
}
TestSuite(complex, .init = setup_complex, .fini = teardown);

/*** TESTS ***/
Test(complex, test1)
{
   SCIP_Real binvrow[3];
   SCIP_Real binvcol[3];
   SCIP_Real coef[3];
   SCIP_Real coeftwo[3];
   SCIP_Real objval;
   int cstats[3];
   int nrows;
   int rstats[3];
   int basinds[3];
   int inds[3];
   int ninds;
   int idx;
   int entry;
   int i;

   /* expected values for the first column of BInv with corresponding variables */
   int exp_vars[] = {-2, 1, 2};
   float exp_vals[] = {0.0, 0.0, -1.0};

   /* expected values for the first column of BAInv with corresponding variables */
   float exp_avals[] = {-0.5, 0.5, 1.0};

   /* ------------------------------------- */
   /* first solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   cr_assert_float_eq(objval, 14.0, EPS);

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
   /* assert that we found the slack variable in the basis */
   cr_assert_lt(i, nrows);

   /* check basis inverse for the row corresponding to the basic slack variable */
   SCIP_CALL( SCIPlpiGetBInvRow(lpi, i, binvrow, NULL, NULL) );

   /* row of basis inverse should be (0, 1, 0.5) */
   cr_expect_float_eq(binvrow[0], 0.0, EPS, "BInvRow[%d] = %g != %g\n", 0, binvrow[0], 0.0);
   cr_expect_float_eq(binvrow[1], 1.0, EPS, "BInvRow[%d] = %g != %g\n", 1, binvrow[1], 1.0);
   cr_expect_float_eq(binvrow[2], 0.5, EPS, "BInvRow[%d] = %g != %g\n", 2, binvrow[2], 0.5);

   /* check whether sparse version is available and the same */
   SCIP_CALL( SCIPlpiGetBInvRow(lpi, i, coef, inds, &ninds) );
   if ( ninds >= 0 )
   {
      cr_assert( ninds == 2 );
      for (entry = 0; entry < ninds; ++entry)
      {
         idx = inds[entry];
         cr_assert( 0 <= idx && idx < 3 );
         cr_expect_float_eq(coef[idx], binvrow[idx], EPS);
      }
   }

   /* check first column of basis inverse */
   SCIP_CALL( SCIPlpiGetBInvCol(lpi, 0, binvcol, NULL, NULL) );

   /* The columns will be in the same order, however, the rows might be permuted.
    * For each row/entry we check that it corresponds to the value of the corresponding variable.
    * The correspondance variable to row/entry is given by basinds. */
   for( entry = 0; entry < nrows; entry++ )
   {
      /* for the given entry try each variable in exp_vars */
      for( idx = 0; idx < nrows; idx++)
      {
         /* Check that the value is the expected one if the column corresponds to the current variable given in exp_vars. */
         if ( exp_vars[idx] == basinds[entry] )
         {
            cr_expect_float_eq(binvcol[entry], exp_vals[idx], EPS);
         }
      }
   }

   /* check whether number of nonzeros fits */
   SCIP_CALL( SCIPlpiGetBInvCol(lpi, 0, coef, inds, &ninds) );
   if ( ninds >= 0 )
   {
      cr_assert( ninds == 1 );
   }

   /* check basis inverse times nonbasic matrix for row corresponding to the basic slack variable */
   cr_assert_geq(i, 0);
   cr_assert_lt(i, nrows);
   SCIP_CALL( SCIPlpiGetBInvARow(lpi, i, NULL, coef, NULL, NULL) );

   /* row of basis inverse times nonbasic matrix should be (-0.5, 0, 0) */
   cr_expect_float_eq(coef[0], -0.5, EPS, "BInvARow[%d] = %g != %g\n", 0, coef[0], -0.5);
   cr_expect_float_eq(coef[1], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 1, coef[1], 0.0);
   cr_expect_float_eq(coef[2], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 2, coef[2], 0.0);

   /* check nonzeros */
   SCIP_CALL( SCIPlpiGetBInvARow(lpi, i, NULL, coeftwo, inds, &ninds) );
   if ( ninds >= 0 )
   {
      cr_assert( ninds == 1 );
      for (entry = 0; entry < ninds; ++entry)
      {
         idx = inds[entry];
         cr_assert( 0 <= idx && idx < 3 );
         cr_expect_float_eq(coeftwo[idx], coef[idx], EPS);
      }
   }

   /* check first column of basis inverse times nonbasic matrix */
   SCIP_CALL( SCIPlpiGetBInvACol(lpi, 0, coef, NULL, NULL) );

   /* The columns will be in the same order, however, the rows will be permuted.
    * For each row/entry we check that it corresponds to the value of the corresponding variable.
    * The correspondance variable to row/entry is given by basinds. */
   for( entry = 0; entry < nrows; entry++ )
   {
      /* for the given entry try each variable in exp_vars */
      for( idx = 0; idx < nrows; idx++)
      {
         /* Check that the value is the expected one if the column corresponds to the current variable given in exp_vars. */
         if ( exp_vars[idx] == basinds[entry] )
         {
            cr_expect_float_eq(coef[entry], exp_avals[idx], EPS);
         }
      }
   }

   /* check nonzeros */
   SCIP_CALL( SCIPlpiGetBInvACol(lpi, 0, coef, inds, &ninds) );
   if ( ninds >= 0 )
   {
      cr_assert( ninds == 3 );
   }
}

/*** TEST SUITE MORE VARS THAN ROWS ***/
static
void setup_more_vars(void)
{
   int ncols;
   int beg = 0;
   int inds[4];
   int nrows;
   SCIP_Real vals[4];
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   SCIP_Real lhs;
   SCIP_Real rhs;

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );

   /* use the following LP:
    * max 1 x1 + 1 x2 + 1 x3 + x4
    *       +x1          + x3 +   x4 >=   1
    *       +x1 +   x2   - x3 -   x4 >=   1
    *       -x1 - 2 x2        - 3 x4 >= -12
    *        x1,    x2,    x3,    x4 >=   0
    */
   /* add columns */
   lb = 0.0;
   ub = SCIPlpiInfinity(lpi);
   obj = 1.0;

   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add rows */
   lhs = 1.0;
   rhs = SCIPlpiInfinity(lpi);
   inds[0] = 0; inds[1] = 2; inds[2] = 3;
   vals[0] = 1.0; vals[1] = 1.0; vals[2] = 1.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 3, &beg, inds, vals) );

   lhs = 1.0;
   rhs = SCIPlpiInfinity(lpi);
   inds[0] = 0; inds[1] = 1; inds[2] = 2; inds[3] = 3;
   vals[0] = 1.0; vals[1] = 1.0; vals[2] = -1.0; vals[3] = -1.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 4, &beg, inds, vals) );

   lhs = -12;
   rhs = SCIPlpiInfinity(lpi);
   inds[0] = 0; inds[1] = 1; inds[2] = 3;
   vals[0] = -1.0; vals[1] = -2.0; vals[2] = -3.0;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 3, &beg, inds, vals) );

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   cr_assert_eq(nrows, 3);
   cr_assert_eq(ncols, 4);

#ifdef SCIP_DEBUG
   /* turn on output */
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) );
#endif
}
TestSuite(more_vars, .init = setup_more_vars, .fini = teardown);

/*** TESTS ***/
Test(more_vars, test1)
{
   SCIP_Real binvarow[4];
   SCIP_Real objval;
   int cstats[4];
   int rstats[3];
   int basinds[3];
   int basicvarpos;

   SCIP* scip;
   SCIP_CALL( SCIPcreate(&scip) );
   SCIPprintVersion(scip, 0);
   SCIP_CALL( SCIPfree(&scip) );

   /* ------------------------------------- */
   /* first solve problem */
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
   cr_assert_float_eq(objval, 23.0, EPS);

   /* the optimal basis should be: {x1, x3, s1 = slack for first row} */
   SCIP_CALL( SCIPlpiGetBase(lpi, cstats, rstats) );
   cr_assert(cstats[0] == SCIP_BASESTAT_BASIC);
   cr_assert(cstats[1] == SCIP_BASESTAT_LOWER);
   cr_assert(cstats[2] == SCIP_BASESTAT_BASIC);
   cr_assert(cstats[3] == SCIP_BASESTAT_LOWER);

   cr_assert(rstats[0] == SCIP_BASESTAT_BASIC);
   cr_assert(rstats[1] == SCIP_BASESTAT_LOWER);
   cr_assert(rstats[2] == SCIP_BASESTAT_LOWER);

   /* binvarow should be
    * 1.0   2.0  0.0   3.0  <- basic var x1
    * 0.0   1.0  1.0   4.0  <- basic var x3
    * 0.0  -3.0  0.0  -6.0  <- basic var s1
    */

   /* get basis indices */
   SCIP_CALL( SCIPlpiGetBasisInd(lpi, basinds) );

   /* find position of x1 in basis indices; check binvarow of row where x1 is basic */
   for( basicvarpos = 0; basicvarpos < 3; ++basicvarpos )
   {
      if( basinds[basicvarpos] == 0 )
         break;
   }
   cr_assert(basicvarpos < 3); /* assert that we found the variable */

   SCIP_CALL( SCIPlpiGetBInvARow(lpi, basicvarpos, NULL, binvarow, NULL, NULL) );
   cr_expect_float_eq(binvarow[0], 1.0, EPS);
   cr_expect_float_eq(binvarow[1], 2.0, EPS);
   cr_expect_float_eq(binvarow[2], 0.0, EPS);
   cr_expect_float_eq(binvarow[3], 3.0, EPS);

   /* find position of x3 in basis indices; check binvarow of row where x3 is basic */
   for( basicvarpos = 0; basicvarpos < 3; ++basicvarpos )
   {
      if( basinds[basicvarpos] == 2 )
         break;
   }
   cr_assert(basicvarpos < 3); /* assert that we found the variable */

   SCIP_CALL( SCIPlpiGetBInvARow(lpi, basicvarpos, NULL, binvarow, NULL, NULL) );
   cr_expect_float_eq(binvarow[0], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 0, binvarow[0], 0.0);
   cr_expect_float_eq(binvarow[1], 1.0, EPS, "BInvARow[%d] = %g != %g\n", 1, binvarow[1], 0.0);
   cr_expect_float_eq(binvarow[2], 1.0, EPS, "BInvARow[%d] = %g != %g\n", 2, binvarow[2], 0.0);
   cr_expect_float_eq(binvarow[3], 4.0, EPS, "BInvARow[%d] = %g != %g\n", 3, binvarow[3], 0.0);

   /* find position of s1 in basis indices; check binvarow of row where s1 is basic */
   for( basicvarpos = 0; basicvarpos < 3; ++basicvarpos )
   {
      if( basinds[basicvarpos] == -1 )
         break;
   }
   cr_assert(basicvarpos < 3); /* assert that we found the variable */

   SCIP_CALL( SCIPlpiGetBInvARow(lpi, basicvarpos, NULL, binvarow, NULL, NULL) );
   cr_expect_float_eq(binvarow[0], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 0, binvarow[0], 0.0);
   cr_expect_float_eq(binvarow[1], -3.0, EPS, "BInvARow[%d] = %g != %g\n", 1, binvarow[1], -3.0);
   cr_expect_float_eq(binvarow[2], 0.0, EPS, "BInvARow[%d] = %g != %g\n", 2, binvarow[2], 0.0);
   cr_expect_float_eq(binvarow[3], -6.0, EPS, "BInvARow[%d] = %g != %g\n", 3, binvarow[3], -6.0);
}
