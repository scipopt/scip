/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

static SCIP_Rational* r1;
static SCIP_Rational* r2;
static SCIP_Rational* rbuf;
static SCIP_RATIONALARRAY* ratar;
static SCIP_Rational** ratpt;
static BMS_BLKMEM* blkmem;
static BMS_BUFMEM* buffmem;

static void setup(void)
{
   blkmem = BMScreateBlockMemory(1, 10);
   buffmem = BMScreateBufferMemory(1.2, 4, FALSE);

   (void) RatCreateBuffer(buffmem, &rbuf);
   (void) RatCreateBlock(blkmem, &r1);
   RatCopy(blkmem, &r2, r1);

   (void) RatCreateBlockArray(blkmem, &ratpt, 10);
   SCIPrationalarrayCreate(&ratar, blkmem);
}

static void teardown(void)
{
   RatFreeBlock(blkmem, &r1);
   RatFreeBlock(blkmem, &r2);
   RatFreeBuffer(buffmem, &rbuf);

   RatFreeBlockArray(blkmem, &ratpt, 10);
   SCIPrationalarrayFree(&ratar, blkmem);

   BMSdestroyBlockMemory(&blkmem);
   BMSdestroyBufferMemory(&buffmem);
}

/* TEST SUITE */
TestSuite(rationals, .init = setup, .fini = teardown);

#ifdef SCIP_WITH_BOOST

/* TESTS  */
Test(rationals, create_and_free)
{
   SCIP_Rational** buffer;
   size_t nusedbuffers;
   nusedbuffers = BMSgetNUsedBufferMemory(buffmem);

   (void) RatCreateBufferArray(buffmem, &buffer, 5);
   RatReallocBufferArray(buffmem, &buffer, 5, 10);
   RatFreeBufferArray(buffmem, &buffer, 10);
   cr_assert_eq(nusedbuffers, BMSgetNUsedBufferMemory(buffmem));
}

Test(rationals, setting, .description = "tests all the different methods to set/get rationals")
{
   int testint = 100345;
   SCIP_Real testreal = 1235.235690;
   SCIP_Rational* testr;
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_BOOST)
   mpq_t gmpr;

   mpq_init(gmpr);
   mpq_set_d(gmpr, 1.2345);

   /* create some rationals with different methods*/

   (void) RatCreateGMP(blkmem, &testr, gmpr);
#endif

   /* test setter methods */
   RatSetInt(r1, testint, 1);
   cr_assert_eq(RatApproxReal(r1), testint, "setting from and converting back to int did not give same result");
   /* set to fp number */
   RatSetReal(r1, testreal);
   cr_assert_eq(RatApproxReal(r1), testreal, "setting from and converting back to real did not give same result");
   cr_assert(RatIsFpRepresentable(r1), "fp-rep number not detected as representable");

   /* set to string rep */
   RatSetReal(r1, 0.1246912);
   cr_log_info("%.17e \n", RatApproxReal(r1));
   cr_assert(RatIsFpRepresentable(r1), "fp number 0.124691234 not fp representable");

   /* set to string rep */
   RatSetString(r1, "1/3");
   cr_assert(!RatIsFpRepresentable(r1), "non-fp-rep number not detected as non-representable");
   cr_assert(!RatIsEqualReal(r1, RatApproxReal(r1)), "approximation of 1/3 should not be the same");

   /* test rounding */
   cr_log_info("printing test approx: %.17e \n", RatApproxReal(r1));
   cr_log_info("rounding down:        %.17e \n", RatRoundReal(r1, SCIP_R_ROUND_DOWNWARDS));
   cr_log_info("rounding up:          %.17e \n", RatRoundReal(r1, SCIP_R_ROUND_UPWARDS));
   cr_log_info("rounding nearest:     %.17e \n", RatRoundReal(r1, SCIP_R_ROUND_NEAREST));

   /* test that rounding down is lt rounding up */
   cr_assert_lt(RatRoundReal(r1, SCIP_R_ROUND_DOWNWARDS), RatRoundReal(r1, SCIP_R_ROUND_UPWARDS), "rounding down should be lt rounding up");

   /* test gmp conversion */
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_BOOST)
   RatSetGMP(r1, gmpr);
   cr_assert_eq(RatApproxReal(r1), mpq_get_d(gmpr), "gmp and Rational should be the same ");
   cr_assert(0 == mpq_cmp(gmpr, *RatGetGMP(r1)));
#endif

   /* delete the rationals */
   RatFreeBlock(blkmem, &testr);
}

Test(rationals, arithmetic, .description = "tests rational arithmetic methods")
{
   SCIP_Longint  intval;
   SCIP_Real doub;
   char buf[SCIP_MAXSTRLEN];

   doub = 12.3548933;

   /* test infinity values */
   RatSetString(r1, "inf");
   cr_assert(RatIsInfinity(r1));
   cr_assert(!RatIsNegInfinity(r1));
   RatSetString(r1, "-inf");
   cr_assert(!RatIsInfinity(r1));
   cr_assert(RatIsNegInfinity(r1));
   cr_assert(RatIsAbsInfinity(r1));

   /* inf printing */
   RatSetString(r1, "inf");
   RatToString(r1, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print inf: %s \n", buf);
   RatSetString(r1, "-inf");
   RatToString(r1, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print -inf: %s \n", buf);

   /* multiplication with inf */
   RatMultReal(r2, r1, 0);
   cr_assert(RatIsZero(r2));
   cr_assert(!RatIsAbsInfinity(r2));
   RatMult(r2, r2, r1);
   cr_assert(RatIsZero(r2));

   /* mul -inf * -1 */
   RatSetInt(r2, -1, 1);
   RatMult(rbuf, r1, r2);
   cr_assert(RatIsInfinity(rbuf));

   /* adding inf and not-inf */
   RatSetString(r1, "-inf");
   RatSetReal(r2, 12.3548934);
   RatAdd(rbuf, r1, r2);

   cr_assert(RatIsNegInfinity(rbuf));
   RatSetString(r1, "inf");
   RatAdd(rbuf, r1, r2);
   cr_assert(RatIsInfinity(rbuf));
   cr_assert(!RatIsEqual(r1, r2));

   /* Difference 0.5 - - 0.5 == 1*/
   RatSetInt(r1, 1, 2);
   RatSetInt(r2, -1, 2);
   RatDiff(rbuf, r1, r2);
   RatSetInt(r1, 1, 1);
   cr_assert(RatIsEqual(rbuf, r1));

   /* Diff 1 - 1.5 == -0.5 */
   RatDiffReal(r1, r1, 1.5);
   cr_assert(RatIsEqual(r1, r2));

   /* RelDiff -5 / -0.5 == -4.5 / 5 */
   RatMultReal(r1, r1, 10);
   RatRelDiff(rbuf, r1, r2);
   RatSetInt(r1, -9, 10);
   cr_assert(RatIsEqual(rbuf, r1));

   /* Division (-9/10) / (-0.5) == 18/10 */
   RatDiv(rbuf, r1, r2);
   RatSetInt(r1, 18, 10);
   cr_assert(RatIsEqual(rbuf, r1));

   RatDivReal(rbuf, r1, 2);
   RatSetInt(r1, 9, 10);
   cr_assert(RatIsEqual(rbuf, r1));

   /* RatAddProd 9/10 += 9/10 * -0.5 == 9/20 */
   RatAddProd(rbuf, r1, r2);
   RatSetInt(r1, 9, 20);
   cr_assert(RatIsEqual(rbuf, r1));

   /* RatDiffProd 9/20 -= 9/20 * -0.5 == 27/40 */
   RatDiffProd(rbuf, r1, r2);
   RatSetInt(r1, 27, 40);
   cr_assert(RatIsEqual(rbuf, r1));

   /* Negation */
   RatNegate(rbuf, r1);
   RatSetInt(r1, -27, 40);
   cr_assert(RatIsEqual(rbuf, r1));

   /* Abs */
   RatAbs(rbuf, r1);
   RatSetInt(r1, 27, 40);
   cr_assert(RatIsEqual(rbuf, r1));

   /* Invert */
   RatInvert(rbuf, r1);
   RatSetInt(r1, 40, 27);
   cr_assert(RatIsEqual(rbuf, r1));

   /* min/max */
   RatMIN(r2, r1, rbuf);
   cr_assert(RatIsEqual(rbuf, r2));
   RatMAX(r2, r1, rbuf);
   cr_assert(RatIsEqual(r1, r2));

   /* comparisons (GT/LT/GE/LE) (only one checed, since use each other)*/
   RatInvert(r1, r1);
   cr_assert(RatIsGT(rbuf, r1));
   RatSetString(r1, "inf");
   cr_assert(RatIsLT(rbuf, r1));
   RatSetString(rbuf, "-inf");
   cr_assert(RatIsLT(rbuf, r1));
   cr_assert(RatIsLT(rbuf, r2));

   /* is positive/negative */
   cr_assert(RatIsPositive(r1));
   cr_assert(RatIsNegative(rbuf));

   /* RatIsIntegral */
   RatSetInt(r1, 15, 5);
   RatSetInt(r2, 15, 7);
   cr_assert(RatIsIntegral(r1));
   cr_assert(!RatIsIntegral(r2));

   /* comparing fp and rat */
   RatSetString(r2, "123646/1215977400");
   cr_assert(!RatIsEqualReal(r2, RatApproxReal(r2)));

   RatSetString(r1, "1/3");
   cr_assert(RatIsLT(r2, r1));

   /* test rounding/ fp approximation */
   RatAdd(rbuf, r1, r1);
   doub = RatApproxReal(r1);
   cr_log_info("rounding nearest:     %.17e \n", RatRoundReal(rbuf, SCIP_R_ROUND_NEAREST));
   cr_log_info("rounding first:       %.17e \n", 2 * doub);
   cr_assert_leq(2 * doub, RatApproxReal(rbuf));
   cr_assert(!RatIsFpRepresentable(rbuf));

   cr_assert(!RatIsIntegral(rbuf));
   RatMultReal(rbuf, r1, 3);
   cr_assert(RatIsIntegral(rbuf));

   /* to string conversions (missing a few) */
   RatToString(r1, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print 1/3: %s \n", buf);
   RatSetString(r1, "3/4");
   RatRoundInteger(&intval, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(intval, 0);
   RatSetString(r1, "3/4");
   RatRoundInteger(&intval, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(intval, 1);
   RatSetString(r1, "-5/4");
   RatRoundInteger(&intval, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(intval, -2);
   RatSetString(r1, "-9/4");
   RatRoundInteger(&intval, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(intval, -2);
}

Test(rationals, overflows, .description = "test conversion methods with huge rationals")
{
   SCIP_Longint testlong = 1;
   SCIP_Longint reslong = 0;

   // round -1/2 up/down/nearest
   RatSetInt(r1, testlong, -(testlong * 2));
   RatRound(r2, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(RatApproxReal(r2), 0);
   RatRoundInteger(&reslong, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(reslong, 0);

   RatRound(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(RatApproxReal(r2), -1);
   RatRoundInteger(&reslong, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(reslong, -1);

   RatRound(r2, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(RatApproxReal(r2), -1);
   RatRoundInteger(&reslong, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(reslong, -1);

   // round 1/2 up/down/nearest
   RatSetInt(r1, testlong, (testlong * 2));
   RatRound(r2, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(RatApproxReal(r2), 1);
   RatRoundInteger(&reslong, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(reslong, 1);

   RatRound(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(RatApproxReal(r2), 0);
   RatRoundInteger(&reslong, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(reslong, 0);

   RatRound(r2, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(RatApproxReal(r2), 1);
   RatRoundInteger(&reslong, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(reslong, 1);

   // round -1/3 up/down/nearest
   RatSetInt(r1, testlong, -(testlong * 3));
   RatRound(r2, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(RatApproxReal(r2), 0);
   RatRound(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(RatApproxReal(r2), -1);
   RatRound(r2, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(RatApproxReal(r2), 0);
}

Test(rationals, arrays, .description = "tests rational array methods")
{
   SCIP_RATIONALARRAY* ratar2;

   cr_log_info("testing rational array methods \n");

   // test getter and setters
   RatSetInt(r1, 2, 5);
   SCIPrationalarraySetVal(ratar, 5, r1);
   SCIPrationalarrayGetVal(ratar, 5, r2);
   cr_assert(RatIsEqual(r1, r2));

   SCIPrationalarrayGetVal(ratar, 4, r2);
   cr_assert(RatIsZero(r2));

   // test get/set
   // ensure no by-ref passing happens
   SCIPrationalarraySetVal(ratar, 10, r1);
   // array is : 0.4 0 0 0 0 0.4
   SCIPrationalarrayGetVal(ratar, 1, r2);
   SCIPrationalarrayGetVal(ratar, 10, r2);
   cr_log_info("Two rat");
   cr_log_info("Test arrray printing: \n ");
   // cr_assert(SCIP_OKAY == SCIPrationalArrayPrint(ratar));
   cr_assert(RatIsEqual(r1, r2));

   // test incval
   RatSetInt(r2, 1, 5);
   RatSetInt(r1, 2, 5);
   SCIPrationalarraySetVal(ratar, 7, r2);
   SCIPrationalarrayIncVal(ratar, 7, r1);

   RatAdd(r1, r1, r2);
   SCIPrationalarrayGetVal(ratar, 7, r2);
   cr_assert(RatIsEqual(r1, r2));

   // test max/minidx
   cr_assert_eq(10, SCIPrationalarrayGetMaxIdx(ratar));
   cr_assert_eq(5, SCIPrationalarrayGetMinIdx(ratar));

   cr_log_info("Test arrray copying: \n ");
   SCIPrationalarrayCopy(&ratar2, blkmem, ratar);
   for( int i = SCIPrationalarrayGetMinIdx(ratar); i < SCIPrationalarrayGetMaxIdx(ratar); i++ )
   {
      SCIPrationalarrayGetVal(ratar, i, r1);
      SCIPrationalarrayGetVal(ratar2, i, r2);
      cr_assert(RatIsEqual(r1, r2));
   }

   SCIPrationalarrayFree(&ratar2, blkmem);
}
#endif
