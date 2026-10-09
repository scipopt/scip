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

/**@file   dbldblarith.c
 * @brief  unit tests for double double arithmetic
 * @author Leona Gottwald;
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

/* example for a number where floor is failing in double precision (see issue #1873) */
volatile double xhi = -2.881951920481567e+17;
volatile double xlo = 6.0078125;

/* TEST SUITE */
void setUp(void)
{
}

void tearDown(void)
{
}

/* TESTS  */
/** @brief tests product with double double arithmetic */
void test_dbldblarith_product(void)
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = (a * b) * c;

   SCIPdbldblProd(dbldblres, dbldblreserr, a, b);
   SCIPdbldblProd21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, c);

   TEST_ASSERT(result != 1.0010219031129229441, "double precision should not be exact");
   TEST_ASSERT_EQUAL_DOUBLE(1.0010219031129229441, (dbldblres + dbldblreserr), "double double arithmetic did not give exact double precision result");
   TEST_ASSERT(fabs(result - 1.0010219031129229441) < 1e-14, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 1.0010219031129229441), fabs((dbldblres + dbldblreserr) - 1.0010219031129229441));
}


/** @brief tests sum with double double arithmetic */
void test_dbldblarith_sum(void)
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = ((c + b) + a) + a;

   SCIPdbldblSum(dbldblres, dbldblreserr, c, b);
   SCIPdbldblSum21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, a);
   SCIPdbldblSum21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, a);

   TEST_ASSERT(result != 1.2345681245670069507e7, "double precision should not be exact");
   TEST_ASSERT_EQUAL_DOUBLE(1.2345681245670069507e7, (dbldblres + dbldblreserr), "double double arithmetic did not give exact double precision result");
   TEST_ASSERT(fabs(result - 1.2345681245670069507e7) < 1e-7, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 1.2345681245670069507e7), fabs((dbldblres + dbldblreserr) - 1.2345681245670069507e7));
}


/** @brief tests division with double double arithmetic */
void test_dbldblarith_division(void)
{
   volatile double result;
   volatile double dbldblres;
   volatile double dbldblreserr;

   result = ((((b / a) / c) / a) / c);

   SCIPdbldblDiv(dbldblres, dbldblreserr, b, a);
   SCIPdbldblDiv21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, c);
   SCIPdbldblDiv21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, a);
   SCIPdbldblDiv21(dbldblres, dbldblreserr, dbldblres, dbldblreserr, c);

   TEST_ASSERT(result != 12.879932287432127447, "double precision should not be exact");
   TEST_ASSERT_EQUAL_DOUBLE(12.879932287432127447, (dbldblres + dbldblreserr), "double double arithmetic did not give exact double precision result");
   TEST_ASSERT(fabs(result - 12.879932287432127447) < 1e-14, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 12.879932287432127447), fabs((dbldblres + dbldblreserr) - 12.879932287432127447));
}


/** @brief tests squaring with double double arithmetic */
void test_dbldblarith_square(void)
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

   TEST_ASSERT(result != 2.7095239662690573102e-24, "double precision should not be exact");
   TEST_ASSERT_EQUAL_DOUBLE(2.7095239662690573102e-24, (dbldblres + dbldblreserr), "double double arithmetic did not give exact double precision result");
   TEST_ASSERT(fabs(result - 2.7095239662690573102e-24) < 1e-39, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 2.7095239662690573102e-24), fabs((dbldblres + dbldblreserr) - 2.7095239662690573102e-24));
}


/** @brief tests sqrt with double double arithmetic */
void test_dbldblarith_sqrt(void)
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

   TEST_ASSERT(result != 1.0000000000000072514, "double precision should not be exact");
   TEST_ASSERT_EQUAL_DOUBLE(1.0000000000000072514, (dbldblres + dbldblreserr), "double double arithmetic did not give exact double precision result");
   TEST_ASSERT(fabs(result - 1.0000000000000072514) < 1e-14, "double precision has to large error");

   printf("error with double: %.18g, error with dbldbl: %.18g\n", fabs(result - 1.0000000000000072514), fabs((dbldblres + dbldblreserr) - 1.0000000000000072514));
}


/** @brief tests floor/ceil with double double arithmetic */
void test_dbldblarith_floor_ceil(void)
{
   double tmphi;
   double tmplo;
   volatile double resultfloor;
   volatile double resultceil;
   volatile double dbldblfloorres;
   volatile double dbldblfloorreserr;
   volatile double dbldblceilres;
   volatile double dbldblceilreserr;

   resultfloor = floor(xhi + xlo);
   resultceil = ceil(-xhi - xlo);

   TEST_ASSERT_EQUAL(resultfloor, (-resultceil), "floor(x) should be equal to -ceil(-x)");

   SCIPdbldblFloor2(dbldblfloorres, dbldblfloorreserr, xhi, xlo);
   SCIPdbldblCeil2(dbldblceilres, dbldblceilreserr, -xhi, -xlo);

   printf("ceil(- (%.16g + %.16g)) = %.16g + %.16g  floor(%.16g + %.16g) = %.16g + %.16g\n", xhi, xlo, dbldblceilres, dbldblceilreserr, xhi, xlo, dbldblfloorres, dbldblfloorreserr);

   TEST_ASSERT_EQUAL(dbldblfloorres, -dbldblceilres, "floor(x) should be equal to -ceil(-x)");
   TEST_ASSERT_EQUAL(dbldblfloorreserr, -dbldblceilreserr, "floor(x) should be equal to -ceil(-x)");

   SCIPdbldblSum22(tmphi, tmplo, xhi, xlo, -dbldblfloorres, -dbldblfloorreserr);
   TEST_ASSERT_EQUAL((tmphi + tmplo), 0.0078125, "double double arithmetic floor should give the correct result 0.0078125");

   SCIPdbldblSum21(tmphi, tmplo, xhi, xlo, -resultfloor);
   TEST_ASSERT_EQUAL((tmphi + tmplo),  6.0078125, "double precision floor should give the incorrect result 6.0078125");
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_dbldblarith_product);
   RUN_TEST(test_dbldblarith_sum);
   RUN_TEST(test_dbldblarith_division);
   RUN_TEST(test_dbldblarith_square);
   RUN_TEST(test_dbldblarith_sqrt);
   RUN_TEST(test_dbldblarith_floor_ceil);
   return UNITY_END();
}
