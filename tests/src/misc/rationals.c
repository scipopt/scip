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

/**@file   rationals.c
 * @brief  unit tests for rationals
 * @author Leon Eifler;
 */

#include "scip/rational.h"
#include "scip/rationalgmp.h"
#include "scip/type_clock.h"
#include "scip/scip_randnumgen.h"
#include <time.h>
#ifdef SCIP_WITH_GMP
#include <gmp.h>
#endif
#include "include/scip_test.h"


/** GLOBAL VARIABLES **/

/* corresponds to exact rational number 5223617350827361/151115727451828646838272 */
volatile double a = 0.000000034567;

/* corresponds to exact rational number 1320494817241769/562949953421312 */
volatile double b = 2.34567;

/* corresponds to exact rational number 6628035890302157/536870912 */
volatile double c = 12345678.9;

static SCIP_RATIONAL* r1;
static SCIP_RATIONAL* r2;
static SCIP_RATIONAL* rbuf;
static SCIP_RATIONALARRAY* ratar;
static SCIP_RATIONAL** ratpt;
static BMS_BLKMEM* blkmem;
static BMS_BUFMEM* buffmem;

static void setup(void)
{
   blkmem = BMScreateBlockMemory(1, 10);
   buffmem = BMScreateBufferMemory(1.2, 4, FALSE);

   (void) SCIPrationalCreateBuffer(buffmem, &rbuf);
   (void) SCIPrationalCreateBlock(blkmem, &r1);
   SCIPrationalCopyBlock(blkmem, &r2, r1);

   (void) SCIPrationalCreateBlockArray(blkmem, &ratpt, 10);
   SCIPrationalarrayCreate(&ratar, blkmem);
}

static void teardown(void)
{
   SCIPrationalFreeBlock(blkmem, &r1);
   SCIPrationalFreeBlock(blkmem, &r2);
   SCIPrationalFreeBuffer(buffmem, &rbuf);

   SCIPrationalFreeBlockArray(blkmem, &ratpt, 10);
   SCIPrationalarrayFree(&ratar, blkmem);

   BMSdestroyBlockMemory(&blkmem);
   BMSdestroyBufferMemory(&buffmem);
}

/* TEST SUITE */
void setUp(void)
{
   setup();
}

void tearDown(void)
{
   teardown();
}

#ifdef SCIP_WITH_BOOST

/* TESTS  */
void test_rationals_create_and_free(void)
{
   SCIP_RATIONAL** buffer;
   size_t nusedbuffers;
   nusedbuffers = BMSgetNUsedBufferMemory(buffmem);

   (void) SCIPrationalCreateBufferArray(buffmem, &buffer, 5);
   SCIPrationalReallocBufferArray(buffmem, &buffer, 5, 10);
   SCIPrationalFreeBufferArray(buffmem, &buffer, 10);
   TEST_ASSERT_EQUAL(nusedbuffers, BMSgetNUsedBufferMemory(buffmem));
}

/** @brief tests all the different methods to set/get rationals */
void test_rationals_setting(void)
{
   SCIP_Longint testint = 100345;
   SCIP_Real testreal = 1235.235690;
   SCIP_RATIONAL* testr;
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_BOOST)
   mpq_t gmpr;

   mpq_init(gmpr);
   mpq_set_d(gmpr, 1.2345);

   /* create some rationals with different methods*/

   (void) SCIPrationalCreateBlockGMP(blkmem, &testr, gmpr);
#endif

   /* test setter methods */
   SCIPrationalSetFraction(r1, testint, 1LL);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r1), testint, "setting from and converting back to int did not give same result");
   /* set to fp number */
   SCIPrationalSetReal(r1, testreal);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r1), testreal, "setting from and converting back to real did not give same result");
   TEST_ASSERT(SCIPrationalIsFpRepresentable(r1), "fp-rep number not detected as representable");

   /* set to string rep */
   SCIPrationalSetReal(r1, 0.1246912);
   TEST_LOG_INFO("%.17e \n", SCIPrationalGetReal(r1));
   TEST_ASSERT(SCIPrationalIsFpRepresentable(r1), "fp number 0.124691234 not fp representable");

   /* set to string rep */
   SCIPrationalSetString(r1, "1/3");
   TEST_ASSERT(!SCIPrationalIsFpRepresentable(r1), "non-fp-rep number not detected as non-representable");
   TEST_ASSERT(!SCIPrationalIsEQReal(r1, SCIPrationalGetReal(r1)), "approximation of 1/3 should not be the same");

   /* test rounding */
   TEST_LOG_INFO("printing test approx: %.17e \n", SCIPrationalGetReal(r1));
   TEST_LOG_INFO("rounding down:        %.17e \n", SCIPrationalRoundReal(r1, SCIP_R_ROUND_DOWNWARDS));
   TEST_LOG_INFO("rounding up:          %.17e \n", SCIPrationalRoundReal(r1, SCIP_R_ROUND_UPWARDS));
   TEST_LOG_INFO("rounding nearest:     %.17e \n", SCIPrationalRoundReal(r1, SCIP_R_ROUND_NEAREST));

   /* test that rounding down is lt rounding up */
   TEST_ASSERT_LESS_THAN(SCIPrationalRoundReal(r1, SCIP_R_ROUND_DOWNWARDS), SCIPrationalRoundReal(r1, SCIP_R_ROUND_UPWARDS), "rounding down should be lt rounding up");

   /* test gmp conversion */
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_BOOST)
   SCIPrationalSetGMP(r1, gmpr);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r1), mpq_get_d(gmpr), "gmp and Rational should be the same ");
   TEST_ASSERT(0 == mpq_cmp(gmpr, *SCIPrationalGetGMP(r1)));
#endif

   /* delete the rationals */
   SCIPrationalFreeBlock(blkmem, &testr);
}

/** @brief tests rational arithmetic methods */
void test_rationals_arithmetic(void)
{
   SCIP_Longint  intval;
   SCIP_Real doub;
   char buf[SCIP_MAXSTRLEN];

   doub = 12.3548933;

   /* test infinity values */
   SCIPrationalSetInfinity(r1);
   TEST_ASSERT(SCIPrationalIsInfinity(r1));
   TEST_ASSERT(!SCIPrationalIsNegInfinity(r1));
   SCIPrationalSetNegInfinity(r1);
   TEST_ASSERT(!SCIPrationalIsInfinity(r1));
   TEST_ASSERT(SCIPrationalIsNegInfinity(r1));
   TEST_ASSERT(SCIPrationalIsAbsInfinity(r1));

   /* inf printing */
   SCIPrationalSetInfinity(r1);
   SCIPrationalToString(r1, buf, SCIP_MAXSTRLEN);
   TEST_LOG_INFO("Test print inf: %s \n", buf);
   SCIPrationalSetNegInfinity(r1);
   SCIPrationalToString(r1, buf, SCIP_MAXSTRLEN);
   TEST_LOG_INFO("Test print -inf: %s \n", buf);

   /* multiplication with inf */
   SCIPrationalMultReal(r2, r1, 0);
   TEST_ASSERT(SCIPrationalIsZero(r2));
   TEST_ASSERT(!SCIPrationalIsAbsInfinity(r2));
   SCIPrationalMult(r2, r2, r1);
   TEST_ASSERT(SCIPrationalIsZero(r2));

   /* mul -inf * -1 */
   SCIPrationalSetFraction(r2, -1LL, 1LL);
   SCIPrationalMult(rbuf, r1, r2);
   TEST_ASSERT(SCIPrationalIsInfinity(rbuf));

   /* adding inf and not-inf */
   SCIPrationalSetNegInfinity(r1);
   SCIPrationalSetReal(r2, 12.3548934);
   SCIPrationalAdd(rbuf, r1, r2);

   TEST_ASSERT(SCIPrationalIsNegInfinity(rbuf));
   SCIPrationalSetInfinity(r1);
   SCIPrationalAdd(rbuf, r1, r2);
   TEST_ASSERT(SCIPrationalIsInfinity(rbuf));
   TEST_ASSERT(!SCIPrationalIsEQ(r1, r2));

   /* Difference 0.5 - - 0.5 == 1*/
   SCIPrationalSetFraction(r1, 1LL, 2LL);
   SCIPrationalSetFraction(r2, -1LL, 2LL);
   SCIPrationalDiff(rbuf, r1, r2);
   SCIPrationalSetFraction(r1, 1LL, 1LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* Diff 1 - 1.5 == -0.5 */
   SCIPrationalDiffReal(r1, r1, 1.5);
   TEST_ASSERT(SCIPrationalIsEQ(r1, r2));

   /* RelDiff -5 / -0.5 == -4.5 / 5 */
   SCIPrationalMultReal(r1, r1, 10);
   SCIPrationalRelDiff(rbuf, r1, r2);
   SCIPrationalSetFraction(r1, -9LL, 10LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* Division (-9/10) / (-0.5) == 18/10 */
   SCIPrationalDiv(rbuf, r1, r2);
   SCIPrationalSetFraction(r1, 18LL, 10LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   SCIPrationalDivReal(rbuf, r1, 2);
   SCIPrationalSetFraction(r1, 9LL, 10LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* SCIPrationalAddProd 9/10 += 9/10 * -0.5 == 9/20 */
   SCIPrationalAddProd(rbuf, r1, r2);
   SCIPrationalSetFraction(r1, 9LL, 20LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* SCIPrationalDiffProd 9/20 -= 9/20 * -0.5 == 27/40 */
   SCIPrationalDiffProd(rbuf, r1, r2);
   SCIPrationalSetFraction(r1, 27LL, 40LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* Negation */
   SCIPrationalNegate(rbuf, r1);
   SCIPrationalSetFraction(r1, -27LL, 40LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* Abs */
   SCIPrationalAbs(rbuf, r1);
   SCIPrationalSetFraction(r1, 27LL, 40LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* Invert */
   SCIPrationalInvert(rbuf, r1);
   SCIPrationalSetFraction(r1, 40LL, 27LL);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r1));

   /* min/max */
   SCIPrationalMin(r2, r1, rbuf);
   TEST_ASSERT(SCIPrationalIsEQ(rbuf, r2));
   SCIPrationalMax(r2, r1, rbuf);
   TEST_ASSERT(SCIPrationalIsEQ(r1, r2));

   /* comparisons (GT/LT/GE/LE) (only one checed, since use each other)*/
   SCIPrationalInvert(r1, r1);
   TEST_ASSERT(SCIPrationalIsGT(rbuf, r1));
   SCIPrationalSetString(r1, "inf"); /* to test an alternative method to create an infinite rational */
   TEST_ASSERT(SCIPrationalIsLT(rbuf, r1));
   SCIPrationalSetString(rbuf, "-inf"); /* to test an alternative method to create an infinite rational */
   TEST_ASSERT(SCIPrationalIsLT(rbuf, r1));
   TEST_ASSERT(SCIPrationalIsLT(rbuf, r2));

   /* is positive/negative */
   TEST_ASSERT(SCIPrationalIsPositive(r1));
   TEST_ASSERT(SCIPrationalIsNegative(rbuf));

   /* SCIPrationalIsIntegral */
   SCIPrationalSetFraction(r1, 15LL, 5LL);
   SCIPrationalSetFraction(r2, 15LL, 7LL);
   TEST_ASSERT(SCIPrationalIsIntegral(r1));
   TEST_ASSERT(!SCIPrationalIsIntegral(r2));

   /* comparing fp and rat */
   SCIPrationalSetString(r2, "123646/1215977400");
   TEST_ASSERT(!SCIPrationalIsEQReal(r2, SCIPrationalGetReal(r2)));

   SCIPrationalSetString(r1, "1/3");
   TEST_ASSERT(SCIPrationalIsLT(r2, r1));

   /* test rounding/ fp approximation */
   SCIPrationalAdd(rbuf, r1, r1);
   doub = SCIPrationalGetReal(r1);
   TEST_LOG_INFO("rounding nearest:     %.17e \n", SCIPrationalRoundReal(rbuf, SCIP_R_ROUND_NEAREST));
   TEST_LOG_INFO("rounding first:       %.17e \n", 2 * doub);
   TEST_ASSERT_LESS_OR_EQUAL(2 * doub, SCIPrationalGetReal(rbuf));
   TEST_ASSERT(!SCIPrationalIsFpRepresentable(rbuf));

   TEST_ASSERT(!SCIPrationalIsIntegral(rbuf));
   SCIPrationalMultReal(rbuf, r1, 3);
   TEST_ASSERT(SCIPrationalIsIntegral(rbuf));

   /* to string conversions (missing a few) */
   SCIPrationalToString(r1, buf, SCIP_MAXSTRLEN);
   TEST_LOG_INFO("Test print 1/3: %s \n", buf);
   SCIPrationalSetString(r1, "3/4");
   SCIPrationalRoundLong(&intval, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(intval, 0);
   SCIPrationalSetString(r1, "3/4");
   SCIPrationalRoundLong(&intval, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(intval, 1);
   SCIPrationalSetString(r1, "-5/4");
   SCIPrationalRoundLong(&intval, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(intval, -2);
   SCIPrationalSetString(r1, "-9/4");
   SCIPrationalRoundLong(&intval, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(intval, -2);
}

/** @brief test conversion methods with huge rationals */
void test_rationals_overflows(void)
{
   SCIP_Longint testlong = 1;
   SCIP_Longint reslong = 0;

   // round -1/2 up/down/nearest
   SCIPrationalSetFraction(r1, testlong, -(testlong * 2));
   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), 0);
   SCIPrationalRoundLong(&reslong, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(reslong, 0);

   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), -1);
   SCIPrationalRoundLong(&reslong, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(reslong, -1);

   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_NEAREST);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), -1);
   SCIPrationalRoundLong(&reslong, r1, SCIP_R_ROUND_NEAREST);
   TEST_ASSERT_EQUAL(reslong, -1);

   // round 1/2 up/down/nearest
   SCIPrationalSetFraction(r1, testlong, (testlong * 2));
   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), 1);
   SCIPrationalRoundLong(&reslong, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(reslong, 1);

   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), 0);
   SCIPrationalRoundLong(&reslong, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(reslong, 0);

   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_NEAREST);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), 1);
   SCIPrationalRoundLong(&reslong, r1, SCIP_R_ROUND_NEAREST);
   TEST_ASSERT_EQUAL(reslong, 1);

   // round -1/3 up/down/nearest
   SCIPrationalSetFraction(r1, testlong, -(testlong * 3));
   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_UPWARDS);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), 0);
   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), -1);
   SCIPrationalRoundInteger(r2, r1, SCIP_R_ROUND_NEAREST);
   TEST_ASSERT_EQUAL(SCIPrationalGetReal(r2), 0);
}

/** @brief tests rational array methods */
void test_rationals_arrays(void)
{
   SCIP_RATIONALARRAY* ratar2;

   TEST_LOG_INFO("testing rational array methods \n");

   // test getter and setters
   SCIPrationalSetFraction(r1, 2LL, 5LL);
   SCIPrationalarraySetVal(ratar, 5, r1);
   SCIPrationalarrayGetVal(ratar, 5, r2);
   TEST_ASSERT(SCIPrationalIsEQ(r1, r2));

   SCIPrationalarrayGetVal(ratar, 4, r2);
   TEST_ASSERT(SCIPrationalIsZero(r2));

   // test get/set
   // ensure no by-ref passing happens
   SCIPrationalarraySetVal(ratar, 10, r1);
   // array is : 0.4 0 0 0 0 0.4
   SCIPrationalarrayGetVal(ratar, 1, r2);
   SCIPrationalarrayGetVal(ratar, 10, r2);
   TEST_LOG_INFO("Two rat");
   TEST_LOG_INFO("Test arrray printing: \n ");
   // TEST_ASSERT(SCIP_OKAY == SCIPrationalarrayPrint(ratar));
   TEST_ASSERT(SCIPrationalIsEQ(r1, r2));

   // test incval
   SCIPrationalSetFraction(r2, 1LL, 5LL);
   SCIPrationalSetFraction(r1, 2LL, 5LL);
   SCIPrationalarraySetVal(ratar, 7, r2);
   SCIPrationalarrayIncVal(ratar, 7, r1);

   SCIPrationalAdd(r1, r1, r2);
   SCIPrationalarrayGetVal(ratar, 7, r2);
   TEST_ASSERT(SCIPrationalIsEQ(r1, r2));

   // test max/minidx
   TEST_ASSERT_EQUAL(10, SCIPrationalarrayGetMaxIdx(ratar));
   TEST_ASSERT_EQUAL(5, SCIPrationalarrayGetMinIdx(ratar));

   TEST_LOG_INFO("Test arrray copying: \n ");
   SCIPrationalarrayCopy(&ratar2, blkmem, ratar);
   for( int i = SCIPrationalarrayGetMinIdx(ratar); i < SCIPrationalarrayGetMaxIdx(ratar); i++ )
   {
      SCIPrationalarrayGetVal(ratar, i, r1);
      SCIPrationalarrayGetVal(ratar2, i, r2);
      TEST_ASSERT(SCIPrationalIsEQ(r1, r2));
   }

   SCIPrationalarrayFree(&ratar2, blkmem);
}
#endif

int main(void)
{
   UNITY_BEGIN();
#ifdef SCIP_WITH_BOOST
   RUN_TEST(test_rationals_create_and_free);
   RUN_TEST(test_rationals_setting);
   RUN_TEST(test_rationals_arithmetic);
   RUN_TEST(test_rationals_overflows);
   RUN_TEST(test_rationals_arrays);
#endif
   return UNITY_END();
}
