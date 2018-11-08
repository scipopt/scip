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

/**@file   matrix.c
 * @brief  unit test for checking lpi matrix changes
 * @author Marc Pfetsch
 *
 * We perform two tests:
 * - We set up a matrix, also using new row entries (this usually fails).
 * - We perform several changes to the matrix.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <lpi/lpi.h>

#include <signal.h>
#include "include/scip_test.h"

#define EPS 1e-6

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;

/* TEST SUITE */
static
void setup(void)
{
   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MAXIMIZE) );
}

static
void teardown(void)
{
   SCIP_CALL( SCIPlpiFree(&lpi) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!");
}

TestSuite(matrix, .init = setup, .fini = teardown);

/** TESTS **/

/* This test should fail with an assert from the LPI, which causes SIGABRT to be issued. Thus, this test should pass. */
Test(matrix, create_matrix, .signal = SIGABRT)
{
   SCIP_Real obj = 0.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 1.0;
   SCIP_Real lhs = 1.0;
   SCIP_Real rhs = 2.0;
   SCIP_Real val = 1.0;
   SCIP_Real matval[2];
   SCIP_Real matlhs[2];
   SCIP_Real matrhs[2];
   int nnonz;
   int nrows;
   int ncols;
   int matbeg[2];
   int matind[2];
   int beg = 0;
   int ind = 0;

   /* this test can only work in debug mode, so we make it pass in opt mode */
#ifdef NDEBUG
   abort(); /* return SIGABORT */
#endif

   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add one row */
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* add one more row using a new variable - often undefined behavior */
   ind = 1;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* add corresponding column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* ------------------------------------------------------------ */

   /* check size */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   cr_assert( nrows == 2 );
   cr_assert( ncols == 2 );

   /* get rows */
   SCIP_CALL( SCIPlpiGetRows(lpi, 0, 1, matlhs, matrhs, &nnonz, matbeg, matind, matval) );
   cr_assert( nnonz == 2 );

   cr_assert_float_eq(matlhs[0], 1.0, EPS, "Violation of lhs: %g != %g\n", matlhs[0], 1.0);
   cr_assert_float_eq(matlhs[1], 1.0, EPS, "Violation of lhs: %g != %g\n", matlhs[1], 1.0);

   cr_assert_float_eq(matrhs[0], 2.0, EPS, "Violation of lhs: %g != %g\n", matrhs[0], 2.0);
   cr_assert_float_eq(matrhs[1], 2.0, EPS, "Violation of lhs: %g != %g\n", matrhs[1], 2.0);

   cr_assert( matbeg[0] == 0 );
   cr_assert( matbeg[1] == 1 );

   cr_assert( matind[0] == 0 );
   cr_assert( matind[1] == 1 );

   cr_assert_float_eq(matval[0], 1.0, EPS, "Violation of matrix entry (1,1): %g != %g\n", matval[0], 1.0);
   cr_assert_float_eq(matval[1], 1.0, EPS, "Violation of matrix entry (2,2): %g != %g\n", matval[1], 1.0);
}

Test(matrix, change_matrix)
{
   SCIP_Real obj = 0.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 1.0;
   SCIP_Real lhs = 1.0;
   SCIP_Real rhs = 2.0;
   SCIP_Real val = 1.0;
   SCIP_Real matval[2];
   SCIP_Real matlhs;
   SCIP_Real matrhs;
   SCIP_Real row1[2];
   SCIP_Real row2[2];
   int nnonz;
   int nrows;
   int ncols;
   int matbeg[2];
   int matind[2];
   int beg = 0;
   int ind = 0;
   int i;

   /* add one column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add one more column */
   SCIP_CALL( SCIPlpiAddCols(lpi, 1, &obj, &lb, &ub, NULL, 0, NULL, NULL, NULL) );

   /* add one row */
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* add one more row using second variable */
   ind = 1;
   val = 2;
   SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, 1, &beg, &ind, &val) );

   /* ------------------------------------------------------------ */

   /* change matrix */
   SCIP_CALL( SCIPlpiChgCoef(lpi, 1, 0, 3.0) );

   /* setup matrix */
   row1[0] = 1.0;
   row1[1] = 0.0;
   row2[0] = 3.0;
   row2[1] = 2.0;

   /* ------------------------------------------------------------ */

   /* check matrix */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncols) );
   cr_assert( nrows == 2 );
   cr_assert( ncols == 2 );

   /* check first row */
   SCIP_CALL( SCIPlpiGetRows(lpi, 0, 0, &matlhs, &matrhs, &nnonz, matbeg, matind, matval) );

   cr_assert_float_eq(matlhs, 1.0, EPS, "Violation of lhs: %g != %g\n", matlhs, 1.0);
   cr_assert_float_eq(matrhs, 2.0, EPS, "Violation of rhs: %g != %g\n", matlhs, 2.0);

   for (i = 0; i < nnonz; ++i)
   {
      cr_assert( 0 <= matind[i] && matind[i] <= 1 );
      cr_assert_float_eq(matval[i], row1[matind[i]], EPS, "Violation of matrix entry (1,%d): %g != %g\n", matind[i], matval[i], row1[matind[i]]);
   }

   /* check second row */
   SCIP_CALL( SCIPlpiGetRows(lpi, 1, 1, &matlhs, &matrhs, &nnonz, matbeg, matind, matval) );

   cr_assert_float_eq(matlhs, 1.0, EPS, "Violation of lhs: %g != %g\n", matlhs, 1.0);
   cr_assert_float_eq(matrhs, 2.0, EPS, "Violation of rhs: %g != %g\n", matlhs, 2.0);

   for (i = 0; i < nnonz; ++i)
   {
      cr_assert( 0 <= matind[i] && matind[i] <= 1 );
      cr_assert_float_eq(matval[i], row2[matind[i]], EPS, "Violation of matrix entry (2,%d): %g != %g\n", matind[i], matval[i], row2[matind[i]]);
   }

   /* check get method for first row */
   for (i = 0; i < 2; ++i)
   {
      SCIP_CALL( SCIPlpiGetCoef(lpi, 0, i, &val) );
      cr_assert_float_eq(val, row1[i], EPS, "Violation of matrix entry (1,%d): %g != %g\n", i, val, row1[i]);
   }
   /* check get method for second row */
   for (i = 0; i < 2; ++i)
   {
      SCIP_CALL( SCIPlpiGetCoef(lpi, 1, i, &val) );
      cr_assert_float_eq(val, row2[i], EPS, "Violation of matrix entry (1,%d): %g != %g\n", i, val, row2[i]);
   }
}
