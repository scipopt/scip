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

/**@file   regression.c
 * @brief  unit tests for incremental best-fit computation
 * @author Gregor Hendel
 */

#include "scip/pub_misc.h"
#include "scip/scip.h"

#include "include/scip_test.h"

#define EPS 1e-04

/* helper methods */
static
void testRegressionSlope(
   SCIP_REGRESSION*      regression,
   SCIP_Real             value
   )
{
   SCIP_Real regression_slope = SCIPregressionGetSlope(regression);

   /* test slope of regression */
   cr_assert_float_eq(value, regression_slope, EPS,
         "Slope <%.1f> and regression slope <%.4f> are not equal within tolerance.\n", value, regression_slope);
}

static
void testRegressionIntercept(
   SCIP_REGRESSION*      regression,
   SCIP_Real             value
   )
{
   SCIP_Real regression_intercept = SCIPregressionGetIntercept(regression);

   /* test intercept of regression */
   cr_assert_float_eq(value, regression_intercept, EPS,
         "Y intercept <%.1f> and regression intercept <%.4f> are not equal within tolerance.\n", value, regression_intercept);
}

static
void testRegressionSlopeAndIntercept(
   SCIP_REGRESSION*      regression,
   SCIP_Real             slope,
   SCIP_Real             intercept
   )
{
   testRegressionSlope(regression, slope);
   testRegressionIntercept(regression, intercept);
}

/** GLOBAL VARIABLES **/
static SCIP_REGRESSION* regression;

/* TEST SUITE */
static
void setup(void)
{
   regression = NULL;
   SCIP_CALL( SCIPregressionCreate(&regression) );
   cr_assert_not_null(regression);
}

static
void teardown(void)
{
   SCIPregressionFree(&regression);

   cr_assert_null(regression);
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is are memory leak!!");
}

TestSuite(reg, .init = setup, .fini = teardown);

/* TESTS  */
Test(reg, create_and_free)
{
   /* calls setup and teardown */
}

Test(reg, line_properties, .description = "determine regression of a cloud of points that is actually a line")
{
   int x;
   SCIP_Real slope = .5;
   SCIP_Real yintercept = 2;
   int xmax = 5;

   /* create points from an actual line and add them one by one to the regression */
   for( x = 1;  x <= xmax; ++x )
   {
      SCIP_Real y = slope * x + yintercept;
      SCIPregressionAddObservation(regression, (SCIP_Real)x, y);
   }

   testRegressionSlopeAndIntercept(regression, slope, yintercept);
}

Test(reg, all_point_equal, .description = "tests that there is no best-fit line when all points are equal")
{
   int i;
   SCIP_Real x = 0;
   SCIP_Real y = 500;

   /* add 10000 times the same point */
   for (i = 0; i < 10000; ++i)
      SCIPregressionAddObservation(regression, x, y);

   testRegressionSlopeAndIntercept(regression, SCIP_INVALID, SCIP_INVALID);
}

/* test for same Y's over many observations. In this case, the best fit line is horizontal, i.e. it has zero slope
 * and a y-intercept equal to the average y-value
 */
Test(reg, horizontal_line, .description = "tests properties of regression on a horizontal line")
{
   int i;
   SCIP_Real x = 0;
   SCIP_Real y = 500;

   /* add 10000 times a point with a new x coordinate, but the same y */
   for (i = 0; i < 10000; ++i) {

      /* use positive and negative x */
      x = i % 2 == 0 ? i : -1;

      SCIPregressionAddObservation(regression, x, y);
   }

   testRegressionSlopeAndIntercept(regression, 0.0, y);
}

/* test for same X values over many points. In this case, a best fit line is the vertical line going through
 * the mean value of the X observations, with infinite slope and no y-intercept
 */
Test(reg, vertical_line, .description = "tests properties of regression on a horizontal line")
{
   int i;
   SCIP_Real x = 42;
   SCIP_Real y;

   /* add 10000 times a point with a new y coordinate, but the same x */
   for (i = 0; i < 10000; ++i) {

      /* use positive and negative y */
      y = i % 2 == 0 ? i : -1;

      SCIPregressionAddObservation(regression, x, y);
   }

   testRegressionSlopeAndIntercept(regression, SCIP_INVALID, SCIP_INVALID);
}

Test(reg, zero_one_obs, .description = "tests that slope and intercept are SCIP_INVALID for zero or one observation")
{
   do
   {
      testRegressionSlopeAndIntercept(regression, SCIP_INVALID, SCIP_INVALID);
      SCIPregressionAddObservation(regression, 14, 1986);
   }
   while ( SCIPregressionGetNObservations(regression) <= 1);
}

/* test that observations on a straight line from which you remove some observations still fits
 * the remaining points
 */
Test(reg, removing_points, .description = "tests removal of regression observations")
{
   int i;
   SCIP_Real slope = 0.4;
   SCIP_Real intercept = -83;

   for( i = 0; i < 20; ++i )
      SCIPregressionAddObservation(regression, i, slope * i + intercept);

   testRegressionSlopeAndIntercept(regression, slope, intercept);

   /* remove every second previously added point */
   for( i = 19; i >= 1; i -= 2 )
      SCIPregressionRemoveObservation(regression, i, slope * i + intercept);

   testRegressionSlopeAndIntercept(regression, slope, intercept);
}
