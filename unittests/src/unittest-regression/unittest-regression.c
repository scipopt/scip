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

/**@file   unittest-regression.c
 * @brief  unit tests for incremental best-fit computation
 * @author Gregor Hendel
 */

#include "scip/pub_misc.h"
#include "scip/scip.h"

/** macro to check the return of tests */
#define CHECK_TEST(x)                            \
   do                                            \
   {                                             \
      SCIP_RETCODE retcode = (x);                \
      if( retcode != SCIP_OKAY )                 \
      {                                          \
         printf("Unit test " #x " failed\n");    \
         SCIPprintError(retcode);                \
         return -1;                              \
      }                                          \
   }                                             \
   while( FALSE )

/** try to create and free a regression */
static
SCIP_RETCODE regressionUnittestCreateAndFree(
   void
   )
{
   SCIP_REGRESSION* regression;

   printf("Testing creation and freeing a regression\n");
   SCIP_CALL( SCIPregressionCreate(&regression) );
   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

static
SCIP_RETCODE testRegressionSlope(
   SCIP_REGRESSION*      regression,
   SCIP_Real             value
   )
{
   /* test slope of regression */
   if( !EPSEQ(value, SCIPregressionGetSlope(regression), 1e-4) )
   {
      printf("Error: Slope <%.1f> and regression slope <%.4f> are not equal within tolerance.\n", value, SCIPregressionGetSlope(regression));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE testRegressionIntercept(
   SCIP_REGRESSION*      regression,
   SCIP_Real             value
   )
{
   if( !EPSEQ(value, SCIPregressionGetIntercept(regression), 1e-4) )
   {
      printf("Error: Y intercept <%.1f> and regression intercept <%.4f> are not equal within tolerance.\n", value, SCIPregressionGetIntercept(regression));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE testRegressionSlopeAndIntercept(
   SCIP_REGRESSION*      regression,
   SCIP_Real             slope,
   SCIP_Real             intercept
   )
{
   SCIP_CALL( testRegressionSlope(regression, slope) );
   SCIP_CALL( testRegressionIntercept(regression, intercept) );

   return SCIP_OKAY;
}

/** try determine the regression and slope of a cloud of points that is actually a line */
static
SCIP_RETCODE regressionUnittestLineProperties(
   void
   )
{
   SCIP_REGRESSION* regression;
   int x;
   SCIP_RETCODE result;

   SCIP_Real slope = .5;
   SCIP_Real yintercept = 2;
   int xmax = 5;

   regression = NULL;
   SCIP_CALL( SCIPregressionCreate(&regression) );
   assert(regression != NULL);

   result = SCIP_OKAY;
   printf("Testing regression on a linear function\n");
   /* create points from an actual line and add them one by one to the regression */
   for( x = 1;  x <= xmax; ++x )
   {
      SCIP_Real y = slope * x + yintercept;
      SCIPregressionAddObservation(regression, (SCIP_Real)x, y);
   }

   /* test slope of regression */
   SCIP_CALL( testRegressionSlope(regression, slope) );

   /* test intercept of regression */
   SCIP_CALL( testRegressionIntercept(regression, yintercept) );


   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

/* if all points are actually equal, there is no best-fit line */
static
SCIP_RETCODE regressionUnittestManyEqualPoints(
   void
   )
{
   SCIP_REGRESSION* regression;
   int i;
   SCIP_Real x = 0;
   SCIP_Real y = 500;

   SCIP_CALL( SCIPregressionCreate(&regression) );

   printf("Testing regression slope and intercept on many equal points\n");
   /* add 10000 times the same point */
   for (i = 0; i < 10000; ++i) {
      SCIPregressionAddObservation(regression, x, y);
   }

   SCIP_CALL( testRegressionSlopeAndIntercept(regression, SCIP_INVALID, SCIP_INVALID) );

   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

/* test for same Y's over many observations. In this case, the best fit line is horizontal, i.e. it has zero slope
 * and a y-intercept equal to the average y-value
 */
static
SCIP_RETCODE regressionUnittestSameY(
   void
   )
{
   SCIP_REGRESSION* regression;
   int i;
   SCIP_Real x = 0;
   SCIP_Real y = 500;

   SCIP_CALL( SCIPregressionCreate(&regression) );

   printf("Testing same y points\n");
   /* add 10000 times a point with a new x coordinate, but the same y */
   for (i = 0; i < 10000; ++i) {

      /* use positive and negative x */
      x = i % 2 == 0 ? i : -1;

      SCIPregressionAddObservation(regression, x, y);
   }

   SCIP_CALL( testRegressionSlopeAndIntercept(regression, 0.0, y) );

   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

/* test for same X values over many points. In this case, a best fit line is the vertical line going through
 * the mean value of the X observations, with infinite slope and no y-intercept
 */
static
SCIP_RETCODE regressionUnittestSameX(
   void
   )
{
   SCIP_REGRESSION* regression;
   int i;
   SCIP_Real x = 42;
   SCIP_Real y;

   SCIP_CALL( SCIPregressionCreate(&regression) );

   printf("Testing same X points\n");
   /* add 10000 times a point with a new x coordinate, but the same y */
   for (i = 0; i < 10000; ++i) {

      /* use positive and negative x */
      y = i % 2 == 0 ? i : -1;

      SCIPregressionAddObservation(regression, x, y);
   }

   SCIP_CALL( testRegressionSlopeAndIntercept(regression, SCIP_INVALID, SCIP_INVALID) );

   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

/* test that no slope and intercept SCIP_INVALID for zero or one observation */
static
SCIP_RETCODE regressionUnittestZeroOrOneObservation(
   void
   )
{
   SCIP_REGRESSION* regression;

   SCIP_CALL( SCIPregressionCreate(&regression) );

   printf("Testing regression slope and intercept on zero and one observation\n");
   do
   {
      SCIP_CALL( testRegressionSlopeAndIntercept(regression, SCIP_INVALID, SCIP_INVALID) );
      SCIPregressionAddObservation(regression, 14, 1986);
   }
   while ( SCIPregressionGetNObservations(regression) <= 1);

   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

/* test that observations on a straight line from which you remove some observations still fits
 * the remaining points
 */
static
SCIP_RETCODE regressionUnittestRemovingPoints(
   void
   )
{
   SCIP_REGRESSION* regression;
   int i;
   SCIP_Real slope = 0.4;
   SCIP_Real intercept = -83;

   SCIP_CALL( SCIPregressionCreate(&regression) );

   printf("Testing removal of points\n");
   for( i = 0; i < 20; ++i )
      SCIPregressionAddObservation(regression, i, slope * i + intercept);

   /* remove every second previously added point */
   for( i = 19; i >= 1; i -= 2 )
      SCIPregressionRemoveObservation(regression, i, slope * i + intercept);

   SCIP_CALL( testRegressionSlopeAndIntercept(regression, slope, intercept) );
   SCIPregressionFree(&regression);

   return SCIP_OKAY;
}

int main(
   int argc,
   char **argv
   )
{
   CHECK_TEST( regressionUnittestCreateAndFree() );
   CHECK_TEST( regressionUnittestZeroOrOneObservation() );
   CHECK_TEST( regressionUnittestLineProperties() );
   CHECK_TEST( regressionUnittestManyEqualPoints() );
   CHECK_TEST( regressionUnittestSameX() );
   CHECK_TEST( regressionUnittestSameY() );
   CHECK_TEST( regressionUnittestRemovingPoints() );
   printf("Unit test for regression passed\n");

   /* for automatic testing output the following */
   printf("SCIP Status        : all tests passed\n");
   printf("Ignore the following:\n");
   printf("  solving          : 0.00\n");
   printf("  nodes (total)    : 0\n");
   printf("  Primal Bound     : 0.0\n");
   printf("  Dual Bound       : 0.0\n");

   return 0;
}
