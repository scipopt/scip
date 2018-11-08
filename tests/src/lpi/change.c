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

/**@file   change.c
 * @brief  unit tests for testing the methods that change and get coefficients in the lpi.
 * @author Franziska Schloesser
 *
 * The tested methods are:
 * SCIPlpiChgCoef, SCIPlpiGetCoef,
 * SCIPlpiChgObj, SCIPlpiGetObj,
 * SCIPlpiChgBounds, SCIPlpiGetBounds,
 * SCIPlpiChgSides, SCIPlpiGetSides,
 * SCIPlpiChgObjsen, SCIPlpiGetObjsen,
 * SCIPlpiGetNCols, SCIPlpiGetNRows, SCIPlpiGetNNonz,
 * SCIPlpiGetCols, SCIPlpiGetRows,
 * SCIPlpiWriteState, SCIPlpiClearState,
 * SCIPlpiWriteLP, SCIPlpiReadLP, SCIPlpiClear
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <scip/message_default.h>
#include <lpi/lpi.h>

#include <signal.h>
#include "include/scip_test.h"

#define EPS 1e-6

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;
static SCIP_MESSAGEHDLR* messagehdlr = NULL;

/** write ncols, nrows and objsen into variables to check later */
static
SCIP_Bool initProb(int pos, int* ncols, int* nrows, int* nnonz, SCIP_OBJSEN* objsen)
{
   SCIP_Real obj[100] = { 1.0, 1.0 };
   SCIP_Real  lb[100] = { 0.0, 0.0 };
   SCIP_Real  ub[100] = { SCIPlpiInfinity(lpi),  SCIPlpiInfinity(lpi) };
   SCIP_Real lhs[100] = {-SCIPlpiInfinity(lpi), -SCIPlpiInfinity(lpi) };
   SCIP_Real rhs[100] = { 1.0, 1.0 };
   SCIP_Real val[100] = { 1.0, 1.0 };
   int beg[100] = { 0, 1 };
   int ind[100] = { 0, 1 };

   /* this is necessary because the theories don't setup and teardown after each call, but once before and after */
   if ( lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&lpi) );
   }

   /* create message handler if necessary */
   if ( messagehdlr == NULL )
   {
      /* create default message handler */
      SCIP_CALL( SCIPcreateMessagehdlrDefault(&messagehdlr, TRUE, NULL, FALSE) );
   }

   /* This name is necessary because if CPLEX reads a problem from a file its problemname will be the filename. */
   SCIP_CALL( SCIPlpiCreate(&lpi, messagehdlr, "lpi_change_test_problem.lp", SCIP_OBJSEN_MAXIMIZE) );

   /* maximization problems, ncols is 1, nrows is 1*/
   switch ( pos )
   {
   case 0:
      /* unbounded - infeasible
       * (P):  max x
       * -x <= 1 (constr)
       *  0 <= x (bound)
       *
       * (D):  min y
       * 1 <= -y (constr)
       * 0 <= y (bound)
       * */
      *ncols = 1;
      *nrows = 1;
      *nnonz = 1;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      val[0] = -1.0;
      break;

   case 1:
      /* optimal - optimal
       * (P):  max x
       *  x <= 0 (constr)
       *  0 <= x (bound)
       *
       * (D):  min 0
       * 1 <= y (constr)
       * 0 <= y (bound)
       * */
      *ncols = 1;
      *nrows = 1;
      *nnonz = 1;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      rhs[0] =  0.0;
      break;

   case 2:
      /* minimization problems (duals of the above) */
      *ncols = 1;
      *nrows = 1;
      *nnonz = 1;
      *objsen = SCIP_OBJSEN_MINIMIZE;
      rhs[0] = SCIPlpiInfinity(lpi);
      lhs[0] = 1;
      val[0] = -1.0;
      break;

   case 3:
      *ncols = 1;
      *nrows = 1;
      *nnonz = 1;
      *objsen = SCIP_OBJSEN_MINIMIZE;
      rhs[0] = SCIPlpiInfinity(lpi);
      lhs[0] = 1;
      obj[0] = 0.0;
      break;

   case 4:
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

   case 5:
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

   case 6:
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

   case 7:
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

   case 8:
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MINIMIZE;
      rhs[0] = SCIPlpiInfinity(lpi);
      rhs[1] = SCIPlpiInfinity(lpi);
      lhs[0] = 1.0;
      lhs[1] = 1.0;
      break;

   case 9:
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
   /* create default message handler */
   assert( messagehdlr == NULL );
   SCIP_CALL( SCIPcreateMessagehdlrDefault(&messagehdlr, TRUE, NULL, FALSE) );

   /* create LPI */
   SCIP_CALL( SCIPlpiCreate(&lpi, messagehdlr, "prob", SCIP_OBJSEN_MAXIMIZE) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
}

static
void teardown(void)
{
   cr_assert( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiFree(&lpi) );

   if ( messagehdlr != NULL )
   {
      SCIP_CALL( SCIPmessagehdlrRelease(&messagehdlr) );
   }

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(change, .init = setup, .fini = teardown);

/* TESTS **/

/* We treat -2 and 2 as -infinity and infinity resp., as SCIPlpiInfinity is not accepted by the TheoryDataPoints */
static
SCIP_Real substituteInfinity(SCIP_Real inf)
{
   if( inf == 2 )
      return SCIPlpiInfinity(lpi);
   else if( inf == -2 )
      return -SCIPlpiInfinity(lpi);
   else
      return inf;
}

/** Test SCIPlpiChgCoef */
static
void checkChgCoef(int row, int col, SCIP_Real newval)
{
   SCIP_Real val;

   SCIP_CALL( SCIPlpiChgCoef(lpi, row, col, newval) );
   cr_assert( !SCIPlpiWasSolved(lpi) );

   SCIP_CALL( SCIPlpiGetCoef(lpi, row, col, &val) );
   cr_assert_float_eq(newval, val, EPS);
}

TheoryDataPoints(change, testchgcoef) =
{
   DataPoints(int, 0, 1),
   DataPoints(int, 0, 1),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((int row, int col, SCIP_Real newval, int prob), change, testchgcoef)
{
   int nrows, ncols, nnonz;
   SCIP_OBJSEN sense;

   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );
   cr_assume_lt( row, nrows );
   cr_assume_lt( col, ncols );
   newval = substituteInfinity( newval );
   checkChgCoef( row, col, newval );
}

/** Test SCIPlpiChgObj */
static
void checkChgObj(int dim, int* ind, SCIP_Real* setobj)
{
   SCIP_Real obj[100];

   assert( dim < 100 );
   SCIP_CALL( SCIPlpiChgObj(lpi, dim, ind, setobj) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiGetObj(lpi, 0, dim-1, obj) );

   cr_assert_arr_eq(obj, setobj, dim*sizeof(SCIP_Real));
}

TheoryDataPoints(change, testchgobjectives) =
{
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((SCIP_Real first, SCIP_Real second, int prob), change, testchgobjectives)
{
   int nrows, ncols, nnonz;
   SCIP_OBJSEN sense;
   int ind[2] = { 0, 1 };
   SCIP_Real setobj[2];

   first = substituteInfinity( first );
   second = substituteInfinity( second );
   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );

   setobj[0] = first;
   setobj[1] = second;
   checkChgObj(ncols, ind, setobj);
}

/** Test SCIPlpiChgBounds */
static
void checkChgBounds(int dim, int* ind, SCIP_Real* setlb, SCIP_Real* setub)
{
   SCIP_Real ub[100];
   SCIP_Real lb[100];

   assert( dim < 100 );
   SCIP_CALL( SCIPlpiChgBounds(lpi, dim, ind, setlb, setub) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, dim - 1, lb, ub) );

   cr_assert_arr_eq(ub, setub, dim*sizeof(SCIP_Real));
   cr_assert_arr_eq(lb, setlb, dim*sizeof(SCIP_Real));
}

TheoryDataPoints(change, testchgbounds) =
{
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((SCIP_Real upper1, SCIP_Real upper2, SCIP_Real lower1, SCIP_Real lower2, int prob), change, testchgbounds)
{
   int nrows, ncols, nnonz;
   int ind[2] = {0, 1};
   SCIP_OBJSEN sense;
   SCIP_Real setub[2];
   SCIP_Real setlb[2];

   lower1 = substituteInfinity( lower1 );
   lower2 = substituteInfinity( lower2 );
   upper1 = substituteInfinity( upper1 );
   upper2 = substituteInfinity( upper2 );
   cr_assume_lt(lower1, upper1);
   cr_assume_lt(lower2, upper2);

   setub[0] = upper1;
   setub[1] = upper2;
   setlb[0] = lower1;
   setlb[1] = lower2;

   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );

   checkChgBounds(ncols, ind, setlb, setub);
}

/** Test SCIPlpiChgSides */
static
void checkChgSides(int dim, int* ind, SCIP_Real* setls, SCIP_Real* setrs)
{
   SCIP_Real ls[100];
   SCIP_Real rs[100];

   assert( dim < 100 );
   SCIP_CALL( SCIPlpiChgSides(lpi, dim, ind, setls, setrs) );
   cr_assert( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiGetSides(lpi, 0, dim - 1, ls, rs) );

   cr_assert_arr_eq(ls, setls, dim*sizeof(SCIP_Real));
   cr_assert_arr_eq(rs, setrs, dim*sizeof(SCIP_Real));
}

TheoryDataPoints(change, testchgsides) =
{
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(SCIP_Real, 0, 1, -1, 2, -2),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((SCIP_Real left1, SCIP_Real left2, SCIP_Real right1, SCIP_Real right2, int prob), change, testchgsides)
{
   int nrows, ncols, nnonz;
   SCIP_OBJSEN sense;
   SCIP_Real setrhs[2];
   SCIP_Real setlhs[2];
   int ind[2] = { 0, 1 };

   left1 = substituteInfinity( left1 );
   left2 = substituteInfinity( left2 );
   right1 = substituteInfinity( right1 );
   right2 = substituteInfinity( right2 );

   setrhs[0] = left1;
   setrhs[1] = left2;
   setlhs[0] = right1;
   setlhs[1] = right2;

   cr_assume_lt(right1, left1);
   cr_assume_lt(right2, left2);

   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );

   checkChgSides(nrows, ind, setlhs, setrhs);
}

/** Test SCIPlpiChgObjsen */
TheoryDataPoints(change, testchgobjsen) =
{
   DataPoints(SCIP_OBJSEN, SCIP_OBJSEN_MAXIMIZE, SCIP_OBJSEN_MINIMIZE),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((SCIP_OBJSEN newsense, int prob), change, testchgobjsen)
{
   int nrows, ncols, nnonz;
   SCIP_OBJSEN sense;
   SCIP_OBJSEN probsense;

   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );

   SCIP_CALL( SCIPlpiChgObjsen(lpi, newsense) );
   cr_assert( !SCIPlpiWasSolved(lpi) );

   SCIP_CALL( SCIPlpiGetObjsen(lpi, &probsense) );

   cr_assert_eq( newsense, probsense, "Expected: %d, got %d\n", newsense, probsense );
}

/** Test SCIPlpiScaleCol */
TheoryDataPoints(change, testscalecol) =
{
   DataPoints(SCIP_Real, 1e10, 1e-10, 1, -1, 2, -2),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((SCIP_Real scale, int prob), change, testscalecol)
{
   int nrows, ncols, cols, rows, nnonz;
   SCIP_OBJSEN sense;
   int col;
   int i;

   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );

   SCIP_CALL( SCIPlpiGetNCols(lpi, &cols) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &rows) );
   cr_assert_eq(ncols, cols);
   cr_assert_eq(nrows, rows);

   assert( nrows < 100 );
   for( col = 0; col < ncols; col++ )
   {
      SCIP_Real colbefore[100];

      for( i = 0; i < nrows; i++ )
      {
         SCIP_Real coef;

         SCIP_CALL( SCIPlpiGetCoef(lpi, i, col, &coef) );
         colbefore[i] = coef;
      }

      SCIP_CALL( SCIPlpiScaleCol(lpi, col, scale) );
      for( i = 0; i < nrows; i++ )
      {
         SCIP_Real colafter;

         SCIP_CALL( SCIPlpiGetCoef(lpi, i, col, &colafter) );
         cr_assert_float_eq( colbefore[i] * scale, colafter, EPS,
            "Found values: scale %.20f, colbefore[i] %.20f, colafter %.20f, i %d\n",
            scale, colbefore[i], colafter, i );
      }
   }
}

/** Test SCIPlpiScaleRow */
TheoryDataPoints(change, testscalerow) =
{
   DataPoints(SCIP_Real, 1e10, 1e-10, 1, -1, 2, -2),
   DataPoints(int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
};

Theory((SCIP_Real scale, int prob), change, testscalerow)
{
   int ncols, nrows, rows, cols, nnonz;
   SCIP_OBJSEN sense;
   int row;
   int i;

   cr_assume( initProb(prob, &ncols, &nrows, &nnonz, &sense) );

   SCIP_CALL( SCIPlpiGetNRows(lpi, &rows) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &cols) );
   cr_assert_eq(nrows, rows);
   cr_assert_eq(ncols, cols);

   assert( nrows < 100 );
   for( row = 0; row < nrows; row++ )
   {
      SCIP_Real rowbefore[100];

      for( i = 0; i < ncols; i++ )
      {
         SCIP_CALL( SCIPlpiGetCoef(lpi, row, i, &rowbefore[i]) );
      }
      SCIP_CALL( SCIPlpiScaleRow(lpi, row, scale) );

      for( i = 0; i < ncols; i++ )
      {
         SCIP_Real rowafter;

         SCIP_CALL( SCIPlpiGetCoef(lpi, row, i, &rowafter) );
         cr_assert_float_eq( rowbefore[i] * scale, rowafter, EPS,
            "Found values: scale %.20f, rowbefore[i] %.20f, rowafter %.20f, i %d\n",
            scale, rowbefore[i], rowafter, i );
      }
   }
}

/** Test for SCIPlpiAddRows, SCIPlpiDelRowset, SCIPlpiDelRows */
Test(change, testrowmethods)
{
   /* problem data */
   SCIP_Real obj[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
   SCIP_Real  lb[5] = { -1.0, -SCIPlpiInfinity(lpi), 0.0, -SCIPlpiInfinity(lpi), 0.0 };
   SCIP_Real  ub[5] = { 10.0, SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi), 29.0, 0.0 };
   int ncolsbefore, ncolsafter;
   int nrowsbefore, nrowsafter;
   SCIP_Real lhsvals[6] = { -SCIPlpiInfinity(lpi), -1.0,   -3e-10, 0.0, 1.0,  3e10 };
   SCIP_Real rhsvals[6] = { -1.0,                  -3e-10, 0.0,    1.0, 3e10, SCIPlpiInfinity(lpi) };
   int     nnonzs[6]  = { 1, 10, -1, 6, -1 };
   int    begvals[6]  = { 0, 2, 3, 5, 8, 9 };
   int    indvals[10] = { 0, 1, 3, 2, 1, 1, 2, 4, 0, 3 };
   SCIP_Real vals[10] = { 1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2 };

   int iterations = 5;
   int k[5] = { 1, 6, -1, 4, -2 };
   int nnonzsdiff[5] = {1, 10, -1, 6, -3 };
   int i;
   int j;

   /* create original lp */
   SCIP_CALL( SCIPlpiAddCols(lpi, 5, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );

   for( i = 0; i < iterations; i++ )
   {
      /* setup row values */
      int nrows;
      int nnonzsbefore;
      int nnonzsafter;

      nrows = k[i];

      /* get data before modification */
      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsbefore) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );

      if (nrows < 0)
      {
         SCIP_CALL( SCIPlpiDelRows(lpi, 0, -(1 + nrows)) );
      }
      else
      {  /* nrows >= 0 */
         SCIP_Real lhs[100];
         SCIP_Real rhs[100];
         int beg[100];

         int nnonz = nnonzs[i];
         int ind[100];
         SCIP_Real val[100];

         SCIP_Real newlhs[100];
         SCIP_Real newval[100];
         SCIP_Real newrhs[100];
         int newbeg[100];
         int newind[100];
         int newnnonz;
         int indold;
         int indnew;

         assert( nrows < 100 );
         for( j = 0; j < nrows; j++ )
         {
            lhs[j] = lhsvals[j];
            rhs[j] = rhsvals[j];
            beg[j] = begvals[j];
         }

         assert( nnonz < 100 );
         for( j = 0; j < nnonz; j++ )
         {
            ind[j] = indvals[j];
            val[j] = vals[j];
         }
         SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, NULL, nnonz, beg, ind, val) );

         /* checks */
         SCIP_CALL( SCIPlpiGetRows(lpi, nrowsbefore, nrowsbefore - 1 + nrows, newlhs, newrhs, &newnnonz, newbeg, newind, newval) );
         cr_assert_eq(nnonz, newnnonz, "expecting %d, got %d\n", nnonz, newnnonz);

         cr_assert_arr_eq(beg, newbeg, nrows*sizeof(int));

         beg[nrows] = nnonz;
         newbeg[nrows] = newnnonz;

         /* check each row seperately */
         for( j = 0; j < nrows; j++ )
         {
            cr_assert_float_eq(lhs[j], newlhs[j], 1e-16);
            cr_assert_float_eq(rhs[j], newrhs[j], 1e-16);

            /* We add a row where the indices are not sorted, some lp solvers give them back sorted (e.g. soplex), some others don't (e.g. cplex).
             * Therefore we cannot simply assert the ind and val arrays to be equal, but have to search for and check each value individually. */
            for( indold = beg[j]; indold < beg[j+1]; indold++ )
            {
               int occurrences = 0;

               /* for each value ind associated to the current row search for it in the newind array */
               for( indnew = beg[j]; indnew < beg[j+1]; indnew++ )
               {
                  if( ind[indold] == newind[indnew] )
                  {
                     occurrences = occurrences + 1;
                     cr_assert_float_eq( val[indold], newval[indnew], 1e-16, "expected %g got %g\n", val[indold], newval[indnew] );
                  }
               }
               /* assert that we found only one occurrence in the current row */
               cr_assert_eq(occurrences, 1);
            }
         }
      }

      /* checks */
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );
      cr_assert_eq(nrowsbefore + nrows, nrowsafter);

      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsafter) );
      cr_assert_eq(nnonzsbefore + nnonzsdiff[i], nnonzsafter, "nnonzsbefore %d, nnonzsafter %d, nnonzsdiff[i] %d, in iteration %d\n",
         nnonzsbefore, nnonzsafter, nnonzsdiff[i], i);

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );
      cr_assert_eq(ncolsbefore, ncolsafter);
   }

   /* delete rowsets */
   /* should have 8 rows now */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );
   cr_assert_eq(8, nrowsbefore);
   for( i = 3; i > 0; i-- )
   {
      int rows[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

      for( j = 0; j < i; j++ )
         rows[(2 * j) + 1] = 1;

      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );
      SCIP_CALL( SCIPlpiDelRowset(lpi, rows) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );

      cr_assert_eq(nrowsbefore - i, nrowsafter);
      /* assert that the rows that are left are the ones I intended */
   }
}

/** Test for SCIPlpiAddCols, SCIPlpiDelColset, SCIPlpiDelCols */
Test(change, testcolmethods)
{
   /* problem data */
   SCIP_Real obj[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
   SCIP_Real lhs[5] = { -1.0, -SCIPlpiInfinity(lpi), 0.0, -SCIPlpiInfinity(lpi), 0.0 };
   SCIP_Real rhs[5] = { 10.0, SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi), 29.0, 0.0 };
   int ncolsbefore, ncolsafter;
   int nrowsbefore, nrowsafter;
   SCIP_Real lbvals[6] = { -SCIPlpiInfinity(lpi), -1.0, -3e-10, 0.0, 1.0, 3e10 };
   SCIP_Real ubvals[6] = { -1.0, -3e-10, 0.0, 1.0, 3e10, SCIPlpiInfinity(lpi) };
   SCIP_Real   vals[10] = { 1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2 };
   int  nnonzs[6]  = { 1, 10, -1, 6, -1 };
   int begvals[6]  = { 0, 2, 3, 5, 8, 9 };
   int indvals[10] = { 0, 1, 3, 2, 1, 1, 2, 4, 0, 3 };

   int iterations = 5;
   int k[5] = { 1, 6, -1, 4, -2 };
   int nnonzsdiff[5] = {1, 10, -1, 6, -3 };
   int i;
   int j;

   /* create original lp */
   SCIP_CALL( SCIPlpiAddRows(lpi, 5, lhs, rhs, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );

   for( i = 0; i < iterations; i++ )
   {
      /* setup col values */
      int ncols;
      int nnonzsbefore;
      int nnonzsafter;

      ncols = k[i];

      /* get data before modification */
      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsbefore) );
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );

      if (ncols < 0)
      {
         SCIP_CALL( SCIPlpiDelCols(lpi, 0, -(1 + ncols)) );
      }
      else
      {  /* ncols >= 0 */
         SCIP_Real lb[100];
         SCIP_Real ub[100];
         int beg[100];

         int nnonz = nnonzs[i];
         int ind[100];
         SCIP_Real val[100];

         SCIP_Real newlb[100];
         SCIP_Real newval[100];
         SCIP_Real newub[100];
         int newbeg[100];
         int newind[100];
         int newnnonz;

         assert( ncols < 100 );
         for( j = 0; j < ncols; j++ )
         {
            lb[j] = lbvals[j];
            ub[j] = ubvals[j];
            beg[j] = begvals[j];
         }

         assert( nnonz < 100 );
         for( j = 0; j < nnonz; j++ )
         {
            ind[j] = indvals[j];
            val[j] = vals[j];
         }
         SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, NULL, nnonz, beg, ind, val) );

         /* checks */
         SCIP_CALL( SCIPlpiGetCols(lpi, ncolsbefore, ncolsbefore-1+ncols, newlb, newub, &newnnonz, newbeg, newind, newval) );
         cr_assert_eq(nnonz, newnnonz, "expecting %d, got %d\n", nnonz, newnnonz);

         cr_assert_arr_eq(lb, newlb, ncols*sizeof(SCIP_Real));
         cr_assert_arr_eq(ub, newub, ncols*sizeof(SCIP_Real));
         cr_assert_arr_eq(beg, newbeg, ncols*sizeof(int));
         cr_assert_arr_eq(ind, newind, nnonz*sizeof(int));
         cr_assert_arr_eq(val, newval, nnonz*sizeof(SCIP_Real));
      }

      /* checks */
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );
      cr_assert_eq(nrowsbefore, nrowsafter);

      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsafter) );
      cr_assert_eq(nnonzsbefore+nnonzsdiff[i], nnonzsafter, "nnonzsbefore %d, nnonzsafter %d, nnonzsdiff[i] %d, in iteration %d\n",
         nnonzsbefore, nnonzsafter, nnonzsdiff[i], i);

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );
      cr_assert_eq(ncolsbefore+ncols, ncolsafter);
   }

   /* delete rowsets */
   /* should have 8 rows now */
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );
   cr_assert_eq(8, ncolsbefore);
   for( i = 3; i > 0; i-- )
   {
      int cols[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

      for( j = 0; j < i; j++ )
         cols[(2 * j) + 1] = 1;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );
      SCIP_CALL( SCIPlpiDelColset(lpi, cols) );
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );

      cr_assert_eq(ncolsbefore - i, ncolsafter);
      /* assert that the rows that are left are the ones I intended */
   }
}

/** Test adding zero coeffs cols */
Test(change, testzerosincols, .signal = SIGABRT)
{
   int nrows;
   int ncols;
   int nnonz = 2;
   SCIP_OBJSEN sense;
   SCIP_Real lb[100] = { 0 };
   SCIP_Real ub[100] = { 20 };
   int beg[1]  = { 0 };
   int ind[2] = { 0, 1 };
   SCIP_Real val[2] = { 0, 3 };
   SCIP_Real obj[1] = { 1 };

   /* 2x2 problem */
   cr_assume( initProb(4, &ncols, &nrows, &nnonz, &sense) );
   cr_assume_eq( 2, nrows );
   cr_assume_eq( 2, ncols );

   SCIP_CALL( SCIPlpiAddCols(lpi, 1, obj, lb, ub, NULL, nnonz, beg, ind, val) );

   /* this test can only work in debug mode, so we make it pass in opt mode */
#ifdef NDEBUG
   abort(); /* return SIGABORT */
#endif
}

/** Test adding zero coeffs in rows, expecting an assert in debug mode
 *
 *  This test should fail with an assert from the LPI, which causes SIGABRT to be issued. Thus, this test should pass.
 */
Test(change, testzerosinrows, .signal = SIGABRT)
{
   int nrows;
   int ncols;
   int nnonz = 2;
   SCIP_OBJSEN sense;
   SCIP_Real lhs[100] = { 0 };
   SCIP_Real rhs[100] = { 20 };
   int beg[1]  = { 0 };
   int ind[2] = { 0, 1 };
   SCIP_Real val[2] = { 0, 3 };

   /* 2x2 problem */
   cr_assume( initProb(4, &ncols, &nrows, &nnonz, &sense) );
   cr_assume_eq( 2, nrows );
   cr_assume_eq( 2, ncols );

   SCIP_CALL( SCIPlpiAddRows(lpi, 1, lhs, rhs, NULL, nnonz, beg, ind, val) );

   /* this test can only work in debug mode, so we make it pass in opt mode */
#ifdef NDEBUG
   abort(); /* return SIGABORT */
#endif
}

/** test SCIPlpiWriteState, SCIPlpiClearState */
Test(change, testlpiwritestatemethods)
{
   int nrows, ncols, nnonz;
   int cstat[2];
   int rstat[2];
   SCIP_OBJSEN sense;

   /* 2x2 problem */
   cr_assume( initProb(5, &ncols, &nrows, &nnonz, &sense) );

   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   SCIP_CALL( SCIPlpiWriteState(lpi, "testlpiwriteandreadstate.bas") );
   SCIP_CALL( SCIPlpiGetBase(lpi, cstat, rstat) );
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* Ideally we want to check in the same matter as in testlpiwritereadlpmethods, but this is
    * not possible at the moment because we have no names for columns and rows. */
   SCIP_CALL( SCIPlpiClear(lpi) );
}

/** test SCIPlpiWriteLP, SCIPlpiReadLP, SCIPlpiClear */
Test(change, testlpiwritereadlpmethods)
{
   int nrows, ncols, nnonz;
   SCIP_Real objval;
   SCIP_Real primsol[2];
   SCIP_Real dualsol[2];
   SCIP_Real activity[2];
   SCIP_Real redcost[2];
   SCIP_Real objval2;
   SCIP_Real primsol2[2];
   SCIP_Real dualsol2[2];
   SCIP_Real activity2[2];
   SCIP_Real redcost2[2];
   SCIP_OBJSEN sense;
   FILE* file;
   FILE* file2;

   /* 2x2 problem */
   cr_assume( initProb(5, &ncols, &nrows, &nnonz, &sense) );

   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost) );

   SCIP_CALL( SCIPlpiWriteLP(lpi, "lpi_change_test_problem.lp") );
   SCIP_CALL( SCIPlpiClear(lpi) );

   SCIP_CALL( SCIPlpiReadLP(lpi, "lpi_change_test_problem.lp") );

   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   SCIP_CALL( SCIPlpiGetSol(lpi, &objval2, primsol2, dualsol2, activity2, redcost2) );
   cr_assert_float_eq( objval, objval2, EPS );
   cr_assert_arr_eq( primsol, primsol2, 2*sizeof(SCIP_Real) );
   cr_assert_arr_eq( dualsol, dualsol2, 2*sizeof(SCIP_Real) );
   cr_assert_arr_eq( activity, activity2, 2*sizeof(SCIP_Real) );
   cr_assert_arr_eq( redcost, redcost2, 2*sizeof(SCIP_Real) );

   SCIP_CALL( SCIPlpiWriteLP(lpi, "lpi_change_test_problem2.lp") );
   SCIP_CALL( SCIPlpiClear(lpi) );

   file = fopen("lpi_change_test_problem.lp", "r");
   file2 = fopen("lpi_change_test_problem2.lp", "r");
   cr_assert_file_contents_eq(file, file2);
}
