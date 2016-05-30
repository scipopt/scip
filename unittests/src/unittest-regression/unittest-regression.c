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
   if( !EPSEQ(slope, SCIPregressionGetSlope(regression), 1e-4) )
   {
      printf("Error: Slope <%.1f> and regression slope <%.4f> are not equal within tolerance.\n", slope, SCIPregressionGetSlope(regression));
      result = SCIP_ERROR;
   }

   /* test intercept of regression */
   else if( !EPSEQ(yintercept, SCIPregressionGetIntercept(regression), 1e-4) )
   {
      printf("Error: Y intercept <%.1f> and regression intercept <%.4f> are not equal within tolerance.\n", yintercept, SCIPregressionGetIntercept(regression));
      result = SCIP_ERROR;
   }


   SCIPregressionFree(&regression);
   return result;
}


int main(
   int argc,
   char **argv
   )
{
   /* todo unit test for 0,1 observations such that it correctly returns SCIP_INVALID */
   /* todo check for many points but all the same */
   /* todo check for slope 0 */
   /* todo check for infinite slope */
   /* todo add if removing points produces correct results */
   CHECK_TEST( regressionUnittestCreateAndFree() );
   CHECK_TEST( regressionUnittestLineProperties() );

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
