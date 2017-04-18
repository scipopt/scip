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

/**@file   dbldblarith.c
 * @brief  unit tests for double double arithmetic
 * @author Robert Lion Gottwald;
 */

#include "scip/dbldblarith.h"

#include "include/scip_test.h"

/** GLOBAL VARIABLES **/

/* corresponds to exact rational number 5223617350827361/151115727451828646838272 */
volatile double a = 0.000000034567;

/* corresponds to exact rational number 1320494817241769/562949953421312 */
volatile double b = 2.34567;

/* corresponds to exact rational number 6628035890302157/536870912 */
volatile double c = 12345678.9;

/* TEST SUITE */
TestSuite(dbldblarith);

Test(dbldblarith, product, .description = "tests product with double double arithmetic")
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = (a * b) * c;

   SCIPdbldblProd(dbldblres, dbldblreserr, a, b);
   SCIPdbldblProd21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, c);

   cr_assert_neq(result, 1.0010219031129229441, "double precision should not be exact");
   cr_assert_eq((dbldblres + dbldblreserr), 1.0010219031129229441, "double double arithmetic did not give exact double precision result");
   cr_assert_lt(fabs(result - 1.0010219031129229441), 1e-14, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 1.0010219031129229441), fabs((dbldblres + dbldblreserr) - 1.0010219031129229441));
}


Test(dbldblarith, sum, .description = "tests sum with double double arithmetic")
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = ((c + b) + a) + a;

   SCIPdbldblSum(dbldblres, dbldblreserr, c, b);
   SCIPdbldblSum21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, a);
   SCIPdbldblSum21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, a);

   cr_assert_neq(result, 1.2345681245670069507e7, "double precision should not be exact");
   cr_assert_eq((dbldblres + dbldblreserr), 1.2345681245670069507e7, "double double arithmetic did not give exact double precision result");
   cr_assert_lt(fabs(result - 1.2345681245670069507e7), 1e-7, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 1.2345681245670069507e7), fabs((dbldblres + dbldblreserr) - 1.2345681245670069507e7));
}


Test(dbldblarith, division, .description = "tests division with double double arithmetic")
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = ((((b / a) / c) / a) / c);

   SCIPdbldblDiv(dbldblres, dbldblreserr, b, a);
   SCIPdbldblDiv21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, c);
   SCIPdbldblDiv21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, a);
   SCIPdbldblDiv21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, c);

   cr_assert_neq(result, 12.879932287432127447, "double precision should not be exact");
   cr_assert_eq((dbldblres + dbldblreserr), 12.879932287432127447, "double double arithmetic did not give exact double precision result");
   cr_assert_lt(fabs(result - 12.879932287432127447), 1e-14, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 12.879932287432127447), fabs((dbldblres + dbldblreserr) - 12.879932287432127447));
}


Test(dbldblarith, square, .description = "tests squaring with double double arithmetic")
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = a * a;
   result = result * result;
   /* scale with 2^60 and square again */
   result *= 1152921504606846976.0;
   result = result * result;

   SCIPdbldblSquare(dbldblres, dbldblreserr, a);
   SCIPdbldblSquare2(dbldblres, dbldblreserr, dbldblres, dbldblreserr);

   /* scale with 2^60 and square again */
   dbldblres *= 1152921504606846976.0;
   dbldblreserr *= 1152921504606846976.0;
   SCIPdbldblSquare2(dbldblres, dbldblreserr, dbldblres, dbldblreserr);

   cr_assert_neq(result, 2.7095239662690573102e-24, "double precision should not be exact");
   cr_assert_eq((dbldblres + dbldblreserr), 2.7095239662690573102e-24, "double double arithmetic did not give exact double precision result");
   cr_assert_lt(fabs(result - 2.7095239662690573102e-24), 1e-39, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 2.7095239662690573102e-24), fabs((dbldblres + dbldblreserr) - 2.7095239662690573102e-24));
}


Test(dbldblarith, sqrt, .description = "tests sqrt with double double arithmetic")
{
   int i;
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = sqrt(c);

   for( i = 0; i < 50; ++i )
      result = sqrt(result);

   SCIPdbldblSqrt(dbldblres, dbldblreserr, c);

   for( i = 0; i < 50; ++i )
      SCIPdbldblSqrt2(dbldblres, dbldblreserr, dbldblres, dbldblreserr);

   cr_assert_neq(result, 1.0000000000000072514, "double precision should not be exact");
   cr_assert_eq((dbldblres + dbldblreserr), 1.0000000000000072514, "double double arithmetic did not give exact double precision result");
   cr_assert_lt(fabs(result - 1.0000000000000072514), 1e-14, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 1.0000000000000072514), fabs((dbldblres + dbldblreserr) - 1.0000000000000072514));
}
