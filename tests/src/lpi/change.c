/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/** macro to keep control of infinity
 *
 *  Some LPIs use std::numeric_limits<SCIP_Real>::infinity() as infinity value. Comparing two infinty values then yields
 *  nan. This is a workaround. */
#undef TEST_ASSERT_DOUBLE_WITHIN_INF
#undef TEST_ASSERT_DOUBLE_WITHIN_INF_
#define TEST_ASSERT_DOUBLE_WITHIN_INF(Actual, Expected, Epsilon) \
   if ( fabs(Actual) > 1e30 || fabs(Expected) > 1e30 )    \
      TEST_ASSERT( Actual == Expected );                    \
   else                                                   \
      TEST_ASSERT_DOUBLE_WITHIN(Actual, Expected, Epsilon);

#define EPS 1e-6

/* GLOBAL VARIABLES */
static SCIP_LPI* lpi = NULL;
static SCIP_MESSAGEHDLR* messagehdlr = NULL;

/** write ncols, nrows and objsen into variables to check later */
static
SCIP_Bool initProb(int pos, int* ncols, int* nrows, int* nnonz, SCIP_OBJSEN* objsen)
{
   const char* rownames[] = { "row1", "row2" };
   const char* colnames[] = { "x1", "x2" };
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
      *ncols = 1;
      *nrows = 1;
      *nnonz = 1;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      val[0] = -1.0;
      break;

   case 1:
      *ncols = 1;
      *nrows = 1;
      *nnonz = 1;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      rhs[0] =  0.0;
      break;

   case 2:
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
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      val[0] = -1.0;
      val[1] = -1.0;
      break;

   case 5:
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      break;

   case 6:
      *ncols = 2;
      *nrows = 2;
      *nnonz = 2;
      *objsen = SCIP_OBJSEN_MAXIMIZE;
      rhs[0] = -1.0;
      rhs[1] = -1.0;
      val[0] = -1.0;
      break;

   case 7:
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
   SCIP_CALL( SCIPlpiAddCols(lpi, *ncols, obj, lb, ub, (char**) colnames, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddRows(lpi, *nrows, lhs, rhs, (char**) rownames, *nnonz, beg, ind, val) );
   TEST_ASSERT( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   TEST_ASSERT( SCIPlpiWasSolved(lpi) );

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
   TEST_ASSERT( !SCIPlpiWasSolved(lpi) );
}

static
void teardown(void)
{
   /* theory-loop tests call theoryCleanup() which already frees lpi/messagehdlr */
   if( lpi != NULL )
   {
      TEST_ASSERT( !SCIPlpiWasSolved(lpi) );
      SCIP_CALL( SCIPlpiFree(&lpi) );
   }

   if ( messagehdlr != NULL )
   {
      SCIP_CALL( SCIPmessagehdlrRelease(&messagehdlr) );
   }

   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

/** cleanup helper for theory-loop functions */
static
void theoryCleanup(void)
{
   if( lpi != NULL )
   {
      SCIPlpiFree(&lpi);
      lpi = NULL;
   }
   if( messagehdlr != NULL )
   {
      SCIPmessagehdlrRelease(&messagehdlr);
      messagehdlr = NULL;
   }
   TEST_ASSERT_EQUAL(BMSgetMemoryUsed(), 0, "Memory leak!");
}


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
   TEST_ASSERT( !SCIPlpiWasSolved(lpi) );

   SCIP_CALL( SCIPlpiGetCoef(lpi, row, col, &val) );
   TEST_ASSERT_DOUBLE_WITHIN_INF(newval, val, EPS);
}

/** Test SCIPlpiChgCoef: 2 rows x 2 cols x 3 values x 10 problems = 120 combos */
void test_change_testchgcoef(void)
{
   int rowvals[2] = {0, 1};
   int colvals[2] = {0, 1};
   SCIP_Real newvalvals[3] = {0, 1, -1};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int ri, ci, vi, pi;

   for( ri = 0; ri < 2; ri++ )
   {
      for( ci = 0; ci < 2; ci++ )
      {
         for( vi = 0; vi < 3; vi++ )
         {
            for( pi = 0; pi < 10; pi++ )
            {
               int row = rowvals[ri];
               int col = colvals[ci];
               SCIP_Real newval = newvalvals[vi];
               int prob = probvals[pi];
               int nrows, ncols, nnonz;
               SCIP_OBJSEN sense;

               if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
                  continue;
               if( !(row < nrows) )
                  continue;
               if( !(col < ncols) )
                  continue;
               newval = substituteInfinity(newval);
               checkChgCoef(row, col, newval);
            }
         }
      }
   }
   theoryCleanup();
}

/** Test SCIPlpiChgObj */
static
void checkChgObj(int dim, int* ind, SCIP_Real* setobj)
{
   SCIP_Real obj[100];

   assert( dim < 100 );
   SCIP_CALL( SCIPlpiChgObj(lpi, dim, ind, setobj) );
   TEST_ASSERT( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiGetObj(lpi, 0, dim-1, obj) );

   TEST_ASSERT_EQUAL_MEMORY(obj, setobj, dim*sizeof(SCIP_Real));
}

/** Test SCIPlpiChgObj: 5 x 5 x 10 = 250 combos */
void test_change_testchgobjectives(void)
{
   SCIP_Real firstvals[5] = {0, 1, -1, 2, -2};
   SCIP_Real secondvals[5] = {0, 1, -1, 2, -2};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int fi, si, pi;

   for( fi = 0; fi < 5; fi++ )
   {
      for( si = 0; si < 5; si++ )
      {
         for( pi = 0; pi < 10; pi++ )
         {
            SCIP_Real first = firstvals[fi];
            SCIP_Real second = secondvals[si];
            int prob = probvals[pi];
            int nrows, ncols, nnonz;
            SCIP_OBJSEN sense;
            int ind[2] = {0, 1};
            SCIP_Real setobj[2];

            first = substituteInfinity(first);
            second = substituteInfinity(second);
            if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
               continue;

            setobj[0] = first;
            setobj[1] = second;
            checkChgObj(ncols, ind, setobj);
         }
      }
   }
   theoryCleanup();
}

/** Test SCIPlpiChgBounds */
static
void checkChgBounds(int dim, int* ind, SCIP_Real* setlb, SCIP_Real* setub)
{
   SCIP_Real ub[100];
   SCIP_Real lb[100];

   assert( dim < 100 );
   SCIP_CALL( SCIPlpiChgBounds(lpi, dim, ind, setlb, setub) );
   TEST_ASSERT( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiGetBounds(lpi, 0, dim - 1, lb, ub) );

   TEST_ASSERT_EQUAL_MEMORY(ub, setub, dim*sizeof(SCIP_Real));
   TEST_ASSERT_EQUAL_MEMORY(lb, setlb, dim*sizeof(SCIP_Real));
}

/** Test SCIPlpiChgBounds: 5^4 x 10 = 6250 combos (with filters) */
void test_change_testchgbounds(void)
{
   SCIP_Real vals[5] = {0, 1, -1, 2, -2};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int u1i, u2i, l1i, l2i, pi;

   for( u1i = 0; u1i < 5; u1i++ )
   {
      for( u2i = 0; u2i < 5; u2i++ )
      {
         for( l1i = 0; l1i < 5; l1i++ )
         {
            for( l2i = 0; l2i < 5; l2i++ )
            {
               for( pi = 0; pi < 10; pi++ )
               {
                  SCIP_Real upper1 = vals[u1i];
                  SCIP_Real upper2 = vals[u2i];
                  SCIP_Real lower1 = vals[l1i];
                  SCIP_Real lower2 = vals[l2i];
                  int prob = probvals[pi];
                  int nrows, ncols, nnonz;
                  int ind[2] = {0, 1};
                  SCIP_OBJSEN sense;
                  SCIP_Real setub[2];
                  SCIP_Real setlb[2];

                  lower1 = substituteInfinity(lower1);
                  lower2 = substituteInfinity(lower2);
                  upper1 = substituteInfinity(upper1);
                  upper2 = substituteInfinity(upper2);
                  if( !(lower1 < upper1) )
                     continue;
                  if( !(lower2 < upper2) )
                     continue;

                  setub[0] = upper1;
                  setub[1] = upper2;
                  setlb[0] = lower1;
                  setlb[1] = lower2;

                  if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
                     continue;

                  checkChgBounds(ncols, ind, setlb, setub);
               }
            }
         }
      }
   }
   theoryCleanup();
}

/** Test SCIPlpiChgSides */
static
void checkChgSides(int dim, int* ind, SCIP_Real* setls, SCIP_Real* setrs)
{
   SCIP_Real ls[100];
   SCIP_Real rs[100];

   assert( dim < 100 );
   SCIP_CALL( SCIPlpiChgSides(lpi, dim, ind, setls, setrs) );
   TEST_ASSERT( !SCIPlpiWasSolved(lpi) );
   SCIP_CALL( SCIPlpiGetSides(lpi, 0, dim - 1, ls, rs) );

   TEST_ASSERT_EQUAL_MEMORY(ls, setls, dim*sizeof(SCIP_Real));
   TEST_ASSERT_EQUAL_MEMORY(rs, setrs, dim*sizeof(SCIP_Real));
}

/** Test SCIPlpiChgSides: 5^4 x 10 = 6250 combos (with filters) */
void test_change_testchgsides(void)
{
   SCIP_Real vals[5] = {0, 1, -1, 2, -2};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int l1i, l2i, r1i, r2i, pi;

   for( l1i = 0; l1i < 5; l1i++ )
   {
      for( l2i = 0; l2i < 5; l2i++ )
      {
         for( r1i = 0; r1i < 5; r1i++ )
         {
            for( r2i = 0; r2i < 5; r2i++ )
            {
               for( pi = 0; pi < 10; pi++ )
               {
                  SCIP_Real left1 = vals[l1i];
                  SCIP_Real left2 = vals[l2i];
                  SCIP_Real right1 = vals[r1i];
                  SCIP_Real right2 = vals[r2i];
                  int prob = probvals[pi];
                  int nrows, ncols, nnonz;
                  SCIP_OBJSEN sense;
                  SCIP_Real setrhs[2];
                  SCIP_Real setlhs[2];
                  int ind[2] = {0, 1};

                  left1 = substituteInfinity(left1);
                  left2 = substituteInfinity(left2);
                  right1 = substituteInfinity(right1);
                  right2 = substituteInfinity(right2);

                  setrhs[0] = left1;
                  setrhs[1] = left2;
                  setlhs[0] = right1;
                  setlhs[1] = right2;

                  if( !(right1 < left1) )
                     continue;
                  if( !(right2 < left2) )
                     continue;

                  if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
                     continue;

                  checkChgSides(nrows, ind, setlhs, setrhs);
               }
            }
         }
      }
   }
   theoryCleanup();
}

/** Test SCIPlpiChgObjsen: 2 x 10 = 20 combos */
void test_change_testchgobjsen(void)
{
   SCIP_OBJSEN sensvals[2] = {SCIP_OBJSEN_MAXIMIZE, SCIP_OBJSEN_MINIMIZE};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int si, pi;

   for( si = 0; si < 2; si++ )
   {
      for( pi = 0; pi < 10; pi++ )
      {
         SCIP_OBJSEN newsense = sensvals[si];
         int prob = probvals[pi];
         int nrows, ncols, nnonz;
         SCIP_OBJSEN sense;
         SCIP_OBJSEN probsense;

         if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
            continue;

         SCIP_CALL( SCIPlpiChgObjsen(lpi, newsense) );
         TEST_ASSERT( !SCIPlpiWasSolved(lpi) );

         SCIP_CALL( SCIPlpiGetObjsen(lpi, &probsense) );

         TEST_ASSERT_EQUAL(newsense, probsense, "Expected: %d, got %d\n", newsense, probsense);
      }
   }
   theoryCleanup();
}

/** Helper method that is used in the scaling tests */
static
void assertScaledBounds(SCIP_Real boundbefore, SCIP_Real boundafter, SCIP_Real oppositeboundafter, SCIP_Real scale)
{
   if( SCIPlpiIsInfinity(lpi, boundbefore) || SCIPlpiIsInfinity(lpi, -boundbefore) )
   {
      /* if infinite, check for equality */
      if( scale > 0 )
      {
         TEST_ASSERT_EQUAL( boundbefore, boundafter,
         "bound before scaling %.20f, bound after scaling %.20f, scale %.20f\n",
         boundbefore, boundafter, scale );
      }
      else
      {
         TEST_ASSERT_EQUAL( -boundbefore, oppositeboundafter,
         "bound before scaling %.20f, bound after scaling %.20f, scale %.20f\n",
         boundbefore, oppositeboundafter, scale );
      }
   }
   else
   {
      /* if finite, check for correct scaling */
      if( scale > 0 )
      {
         TEST_ASSERT_DOUBLE_WITHIN( boundbefore / scale, boundafter, EPS,
         "bound before scaling %.20f, bound after scaling %.20f, scale %.20f\n",
         boundbefore, boundafter, scale );
      }
      else
      {
         TEST_ASSERT_DOUBLE_WITHIN( boundbefore / scale, oppositeboundafter, EPS,
         "bound before scaling %.20f, bound after scaling %.20f, scale %.20f\n",
         boundbefore, oppositeboundafter, scale );
      }
   }
}

/** Helper method that is used in the scaling tests */
static
void assertScaledSides(SCIP_Real sidebefore, SCIP_Real sideafter, SCIP_Real oppositesideafter, SCIP_Real scale)
{
   if( SCIPlpiIsInfinity(lpi, sidebefore) || SCIPlpiIsInfinity(lpi, -sidebefore) )
   {
      /* if infinite, check for equality */
      if( scale > 0 )
      {
         TEST_ASSERT_EQUAL( sidebefore, sideafter,
         "side before scaling %.20f, side after scaling %.20f, scale %.20f\n",
         sidebefore, sideafter, scale );
      }
      else
      {
         TEST_ASSERT_EQUAL( -sidebefore, oppositesideafter,
         "side before scaling %.20f, side after scaling %.20f, scale %.20f\n",
         sidebefore, oppositesideafter, scale );
      }
   }
   else
   {
      /* if finite, check for correct scaling */
      if( scale > 0 )
      {
         TEST_ASSERT_DOUBLE_WITHIN( sidebefore * scale, sideafter, EPS,
         "side before scaling %.20f, side after scaling %.20f, scale %.20f\n",
         sidebefore, sideafter, scale );
      }
      else
      {
         TEST_ASSERT_DOUBLE_WITHIN( sidebefore * scale, oppositesideafter, EPS,
         "side before scaling %.20f, side after scaling %.20f, scale %.20f\n",
         sidebefore, oppositesideafter, scale );
      }
   }
}

/** Test SCIPlpiScaleCol: 6 x 10 = 60 combos */
void test_change_testscalecol(void)
{
   SCIP_Real scalevals[6] = {1e10, 1e-10, 1, -1, 2, -2};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int sci, pi;

   for( sci = 0; sci < 6; sci++ )
   {
      for( pi = 0; pi < 10; pi++ )
      {
         SCIP_Real scale = scalevals[sci];
         int prob = probvals[pi];
         int nrows, ncols, cols, rows, nnonz;
         SCIP_OBJSEN sense;
         int col;
         int i;

         if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
            continue;

         SCIP_CALL( SCIPlpiGetNCols(lpi, &cols) );
         SCIP_CALL( SCIPlpiGetNRows(lpi, &rows) );
         TEST_ASSERT_EQUAL(ncols, cols);
         TEST_ASSERT_EQUAL(nrows, rows);

         assert( nrows < 100 );
         for( col = 0; col < ncols; col++ )
         {
            SCIP_Real colbefore[100];
            SCIP_Real objbefore;
            SCIP_Real lbbefore;
            SCIP_Real ubbefore;
            SCIP_Real objafter;
            SCIP_Real lbafter;
            SCIP_Real ubafter;

            for( i = 0; i < nrows; i++ )
            {
               SCIP_Real coef;

               SCIP_CALL( SCIPlpiGetCoef(lpi, i, col, &coef) );
               colbefore[i] = coef;
            }

            SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &objbefore) );

            SCIP_CALL( SCIPlpiGetBounds(lpi, col, col, &lbbefore, &ubbefore) );

            SCIP_CALL( SCIPlpiScaleCol(lpi, col, scale) );
            for( i = 0; i < nrows; i++ )
            {
               SCIP_Real colafter;

               SCIP_CALL( SCIPlpiGetCoef(lpi, i, col, &colafter) );
               TEST_ASSERT_DOUBLE_WITHIN( colbefore[i] * scale, colafter, EPS,
                  "Found values: scale %.20f, colbefore[i] %.20f, colafter %.20f, i %d\n",
                  scale, colbefore[i], colafter, i );
            }

            SCIP_CALL( SCIPlpiGetObj(lpi, col, col, &objafter) );

            TEST_ASSERT_DOUBLE_WITHIN( objbefore * scale, objafter, EPS,
                  "Found values: scale %.20f, objbefore %.20f, objafter %.20f",
                  scale, objbefore, objafter );

            SCIP_CALL( SCIPlpiGetBounds(lpi, col, col, &lbafter, &ubafter) );

            assertScaledBounds(lbbefore, lbafter, ubafter, scale);
            assertScaledBounds(ubbefore, ubafter, lbafter, scale);
         }
      }
   }
   theoryCleanup();
}

/** Test SCIPlpiScaleRow: 6 x 10 = 60 combos */
void test_change_testscalerow(void)
{
   SCIP_Real scalevals[6] = {1e10, 1e-10, 1, -1, 2, -2};
   int probvals[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int sci, pi;

   for( sci = 0; sci < 6; sci++ )
   {
      for( pi = 0; pi < 10; pi++ )
      {
         SCIP_Real scale = scalevals[sci];
         int prob = probvals[pi];
         int ncols, nrows, rows, cols, nnonz;
         SCIP_OBJSEN sense;
         int row;
         int i;

         if( !initProb(prob, &ncols, &nrows, &nnonz, &sense) )
            continue;

         SCIP_CALL( SCIPlpiGetNRows(lpi, &rows) );
         SCIP_CALL( SCIPlpiGetNCols(lpi, &cols) );
         TEST_ASSERT_EQUAL(nrows, rows);
         TEST_ASSERT_EQUAL(ncols, cols);

         assert( nrows < 100 );
         for( row = 0; row < nrows; row++ )
         {
            SCIP_Real rowbefore[100];
            SCIP_Real lhsbefore;
            SCIP_Real rhsbefore;
            SCIP_Real lhsafter;
            SCIP_Real rhsafter;

            for( i = 0; i < ncols; i++ )
            {
               SCIP_CALL( SCIPlpiGetCoef(lpi, row, i, &rowbefore[i]) );
            }

            SCIP_CALL( SCIPlpiGetSides(lpi, row, row, &lhsbefore, &rhsbefore) );

            SCIP_CALL( SCIPlpiScaleRow(lpi, row, scale) );

            for( i = 0; i < ncols; i++ )
            {
               SCIP_Real rowafter;

               SCIP_CALL( SCIPlpiGetCoef(lpi, row, i, &rowafter) );
               TEST_ASSERT_DOUBLE_WITHIN( rowbefore[i] * scale, rowafter, EPS,
                  "Found values: scale %.20f, rowbefore[i] %.20f, rowafter %.20f, i %d\n",
                  scale, rowbefore[i], rowafter, i );
            }

            SCIP_CALL( SCIPlpiGetSides(lpi, row, row, &lhsafter, &rhsafter) );

            assertScaledSides(lhsbefore, lhsafter, rhsafter, scale);
            assertScaledSides(rhsbefore, rhsafter, lhsafter, scale);
         }
      }
   }
   theoryCleanup();
}

/** Test for SCIPlpiAddRows, SCIPlpiDelRowset, SCIPlpiDelRows */
void test_change_testrowmethods(void)
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
         TEST_ASSERT_EQUAL(nnonz, newnnonz, "expecting %d, got %d\n", nnonz, newnnonz);

         TEST_ASSERT_EQUAL_MEMORY(beg, newbeg, nrows*sizeof(int));

         beg[nrows] = nnonz;
         newbeg[nrows] = newnnonz;

         /* check each row seperately */
         for( j = 0; j < nrows; j++ )
         {
            TEST_ASSERT_DOUBLE_WITHIN_INF(lhs[j], newlhs[j], 1e-16);
            TEST_ASSERT_DOUBLE_WITHIN_INF(rhs[j], newrhs[j], 1e-16);

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
                     TEST_ASSERT_DOUBLE_WITHIN( val[indold], newval[indnew], 1e-16, "expected %g got %g\n", val[indold], newval[indnew] );
                  }
               }
               /* assert that we found only one occurrence in the current row */
               TEST_ASSERT_EQUAL(occurrences, 1);
            }
         }
      }

      /* checks */
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );
      TEST_ASSERT_EQUAL(nrowsbefore + nrows, nrowsafter);

      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsafter) );
      TEST_ASSERT_EQUAL(nnonzsbefore + nnonzsdiff[i], nnonzsafter, "nnonzsbefore %d, nnonzsafter %d, nnonzsdiff[i] %d, in iteration %d\n",
         nnonzsbefore, nnonzsafter, nnonzsdiff[i], i);

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );
      TEST_ASSERT_EQUAL(ncolsbefore, ncolsafter);
   }

   /* delete rowsets */
   /* should have 8 rows now */
   SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );
   TEST_ASSERT_EQUAL(8, nrowsbefore);
   for( i = 3; i > 0; i-- )
   {
      int rows[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

      for( j = 0; j < i; j++ )
         rows[(2 * j) + 1] = 1;

      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsbefore) );
      SCIP_CALL( SCIPlpiDelRowset(lpi, rows) );
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );

      TEST_ASSERT_EQUAL(nrowsbefore - i, nrowsafter);
      /* assert that the rows that are left are the ones I intended */
   }
}

/** Test for SCIPlpiAddCols, SCIPlpiDelColset, SCIPlpiDelCols */
void test_change_testcolmethods(void)
{
   /* problem data */
   SCIP_Real lhs[5] = { -1.0, -SCIPlpiInfinity(lpi), 0.0, -SCIPlpiInfinity(lpi), 0.0 };
   SCIP_Real rhs[5] = { 10.0, SCIPlpiInfinity(lpi), SCIPlpiInfinity(lpi), 29.0, 0.0 };
   int ncolsbefore, ncolsafter;
   int nrowsbefore, nrowsafter;
   SCIP_Real obj[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
   SCIP_Real lbvals[6] = { -SCIPlpiInfinity(lpi), -1.0, -3e-10, 0.0, 1.0, 3e10 };
   SCIP_Real ubvals[6] = { -1.0, -3e-10, 0.0, 1.0, 3e10, SCIPlpiInfinity(lpi) };
   SCIP_Real vals[10] = { 1.0, 5.0, -1.0, 3e5, 2.0, 1.0, 20, 10, -1.9, 1e-2 };
   int nnonzs[5]  = { 1, 10, -1, 6, -1 };
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

         assert( nnonz > 0 );
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
         TEST_ASSERT_EQUAL(nnonz, newnnonz, "expecting %d, got %d\n", nnonz, newnnonz);

         TEST_ASSERT_EQUAL_MEMORY(lb, newlb, ncols*sizeof(SCIP_Real));
         TEST_ASSERT_EQUAL_MEMORY(ub, newub, ncols*sizeof(SCIP_Real));
         TEST_ASSERT_EQUAL_MEMORY(beg, newbeg, ncols*sizeof(int));
         TEST_ASSERT_EQUAL_MEMORY(ind, newind, nnonz*sizeof(int));
         TEST_ASSERT_EQUAL_MEMORY(val, newval, nnonz*sizeof(SCIP_Real));
      }

      /* checks */
      SCIP_CALL( SCIPlpiGetNRows(lpi, &nrowsafter) );
      TEST_ASSERT_EQUAL(nrowsbefore, nrowsafter);

      SCIP_CALL( SCIPlpiGetNNonz(lpi, &nnonzsafter) );
      TEST_ASSERT_EQUAL(nnonzsbefore+nnonzsdiff[i], nnonzsafter, "nnonzsbefore %d, nnonzsafter %d, nnonzsdiff[i] %d, in iteration %d\n",
         nnonzsbefore, nnonzsafter, nnonzsdiff[i], i);

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );
      TEST_ASSERT_EQUAL(ncolsbefore+ncols, ncolsafter);
   }

   /* delete rowsets */
   /* should have 8 rows now */
   SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );
   TEST_ASSERT_EQUAL(8, ncolsbefore);
   for( i = 3; i > 0; i-- )
   {
      int cols[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

      for( j = 0; j < i; j++ )
         cols[(2 * j) + 1] = 1;

      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsbefore) );
      SCIP_CALL( SCIPlpiDelColset(lpi, cols) );
      SCIP_CALL( SCIPlpiGetNCols(lpi, &ncolsafter) );

      TEST_ASSERT_EQUAL(ncolsbefore - i, ncolsafter);
      /* assert that the rows that are left are the ones I intended */
   }
}

/** Test adding zero coeffs cols */
void test_change_testzerosincols(void)
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
   if( !initProb(4, &ncols, &nrows, &nnonz, &sense) )
      return;
   if( 2 != nrows )
      return;
   if( 2 != ncols )
      return;

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
void test_change_testzerosinrows(void)
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
   if( !initProb(4, &ncols, &nrows, &nnonz, &sense) )
      return;
   if( 2 != nrows )
      return;
   if( 2 != ncols )
      return;

   SCIP_CALL( SCIPlpiAddRows(lpi, 1, lhs, rhs, NULL, nnonz, beg, ind, val) );

   /* this test can only work in debug mode, so we make it pass in opt mode */
#ifdef NDEBUG
   abort(); /* return SIGABORT */
#endif
}

/** test SCIPlpiWriteState, SCIPlpiClearState */
void test_change_testlpiwritestatemethods(void)
{
   int nrows, ncols, nnonz;
   int cstat[2];
   int rstat[2];
   SCIP_OBJSEN sense;
   SCIP_RETCODE retcode;

   /* 2x2 problem */
   if( !initProb(5, &ncols, &nrows, &nnonz, &sense) )
      return;

   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );

   retcode = SCIPlpiWriteState(lpi, "testlpiwriteandreadstate.bas");
   if ( retcode != SCIP_OKAY && retcode != SCIP_NOTIMPLEMENTED )
   {
      SCIP_CALL( retcode );
   }
   SCIP_CALL( SCIPlpiGetBase(lpi, cstat, rstat) );
   SCIP_CALL( SCIPlpiClearState(lpi) );

   /* Ideally we want to check in the same matter as in testlpiwritereadlpmethods, but this is
    * not possible at the moment because we have no names for columns and rows. */
   SCIP_CALL( SCIPlpiClear(lpi) );

   remove("testlpiwriteandreadstate.bas");
}

/** test SCIPlpiWriteLP, SCIPlpiReadLP, SCIPlpiClear */
void test_change_testlpiwritereadlpmethods(void)
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

   /* 2x2 problem */
   if( !initProb(5, &ncols, &nrows, &nnonz, &sense) )
      return;

   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primsol, dualsol, activity, redcost) );

   SCIP_CALL( SCIPlpiWriteLP(lpi, "lpi_change_test_problem.lp") );
   SCIP_CALL( SCIPlpiClear(lpi) );

   SCIP_CALL( SCIPlpiReadLP(lpi, "lpi_change_test_problem.lp") );

   SCIP_CALL( SCIPlpiSolvePrimal(lpi) );
   SCIP_CALL( SCIPlpiGetSol(lpi, &objval2, primsol2, dualsol2, activity2, redcost2) );
   TEST_ASSERT_DOUBLE_WITHIN( objval, objval2, EPS );
   TEST_ASSERT_EQUAL_MEMORY( primsol, primsol2, 2*sizeof(SCIP_Real) );
   TEST_ASSERT_EQUAL_MEMORY( dualsol, dualsol2, 2*sizeof(SCIP_Real) );
   TEST_ASSERT_EQUAL_MEMORY( activity, activity2, 2*sizeof(SCIP_Real) );
   TEST_ASSERT_EQUAL_MEMORY( redcost, redcost2, 2*sizeof(SCIP_Real) );

   SCIP_CALL( SCIPlpiWriteLP(lpi, "lpi_change_test_problem2.lp") );
   SCIP_CALL( SCIPlpiClear(lpi) );

   /* skip comparison of file content if LP solver is Xpress
    * Xpress writes a timestamp to the .lp file, which makes a simple check for
    * file content to be equal not reliable
    */
   if( strncmp(SCIPlpiGetSolverName(), "Xpress", 6) != 0 )
   {
      FILE* file;
      FILE* file2;

      file = fopen("lpi_change_test_problem.lp", "r");
      file2 = fopen("lpi_change_test_problem2.lp", "r");
      TEST_ASSERT_FILES_EQUAL(file, file2);

      fclose(file);
      fclose(file2);
   }

   remove("lpi_change_test_problem.lp");
   remove("lpi_change_test_problem2.lp");
}

void setUp(void) { setup(); }

void tearDown(void) { teardown(); }

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_change_testrowmethods);
   RUN_TEST(test_change_testcolmethods);
   /* testzerosincols and testzerosinrows expect SIGABRT which Unity can't handle */
   RUN_TEST(test_change_testlpiwritestatemethods);
   RUN_TEST(test_change_testlpiwritereadlpmethods);
   RUN_TEST(test_change_testchgcoef);
   RUN_TEST(test_change_testchgobjectives);
   RUN_TEST(test_change_testchgbounds);
   RUN_TEST(test_change_testchgsides);
   RUN_TEST(test_change_testchgobjsen);
   RUN_TEST(test_change_testscalecol);
   RUN_TEST(test_change_testscalerow);
   return UNITY_END();
}
