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
TestSuite(rationals, .init = setup, .fini = teardown);

#ifdef SCIP_WITH_BOOST

/* TESTS  */
Test(rationals, create_and_free)
{
   SCIP_Rational** buffer;
   size_t nusedbuffers;
   nusedbuffers = BMSgetNUsedBufferMemory(buffmem);

   (void) SCIPrationalCreateBufferArray(buffmem, &buffer, 5);
   SCIPrationalReallocBufferArray(buffmem, &buffer, 5, 10);
   SCIPrationalFreeBufferArray(buffmem, &buffer, 10);
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

   (void) SCIPrationalCreateBlockGMP(blkmem, &testr, gmpr);
#endif

   /* test setter methods */
   SCIPrationalSetInt(r1, testint, 1);
   cr_assert_eq(SCIPrationalApproxReal(r1), testint, "setting from and converting back to int did not give same result");
   /* set to fp number */
   SCIPrationalSetReal(r1, testreal);
   cr_assert_eq(SCIPrationalApproxReal(r1), testreal, "setting from and converting back to real did not give same result");
   cr_assert(SCIPrationalIsFpRepresentable(r1), "fp-rep number not detected as representable");

   /* set to string rep */
   SCIPrationalSetReal(r1, 0.1246912);
   cr_log_info("%.17e \n", SCIPrationalApproxReal(r1));
   cr_assert(SCIPrationalIsFpRepresentable(r1), "fp number 0.124691234 not fp representable");

   /* set to string rep */
   SCIPrationalSetString(r1, "1/3");
   cr_assert(!SCIPrationalIsFpRepresentable(r1), "non-fp-rep number not detected as non-representable");
   cr_assert(!SCIPrationalIsEqualReal(r1, SCIPrationalApproxReal(r1)), "approximation of 1/3 should not be the same");

   /* test rounding */
   cr_log_info("printing test approx: %.17e \n", SCIPrationalApproxReal(r1));
   cr_log_info("rounding down:        %.17e \n", SCIPrationalRoundReal(r1, SCIP_R_ROUND_DOWNWARDS));
   cr_log_info("rounding up:          %.17e \n", SCIPrationalRoundReal(r1, SCIP_R_ROUND_UPWARDS));
   cr_log_info("rounding nearest:     %.17e \n", SCIPrationalRoundReal(r1, SCIP_R_ROUND_NEAREST));

   /* test that rounding down is lt rounding up */
   cr_assert_lt(SCIPrationalRoundReal(r1, SCIP_R_ROUND_DOWNWARDS), SCIPrationalRoundReal(r1, SCIP_R_ROUND_UPWARDS), "rounding down should be lt rounding up");

   /* test gmp conversion */
#if defined(SCIP_WITH_GMP) && defined(SCIP_WITH_BOOST)
   SCIPrationalSetGMP(r1, gmpr);
   cr_assert_eq(SCIPrationalApproxReal(r1), mpq_get_d(gmpr), "gmp and Rational should be the same ");
   cr_assert(0 == mpq_cmp(gmpr, *SCIPrationalGetGMP(r1)));
#endif

   /* delete the rationals */
   SCIPrationalFreeBlock(blkmem, &testr);
}

Test(rationals, arithmetic, .description = "tests rational arithmetic methods")
{
   SCIP_Longint  intval;
   SCIP_Real doub;
   char buf[SCIP_MAXSTRLEN];

   doub = 12.3548933;

   /* test infinity values */
   SCIPrationalSetString(r1, "inf");
   cr_assert(SCIPrationalIsInfinity(r1));
   cr_assert(!SCIPrationalIsNegInfinity(r1));
   SCIPrationalSetString(r1, "-inf");
   cr_assert(!SCIPrationalIsInfinity(r1));
   cr_assert(SCIPrationalIsNegInfinity(r1));
   cr_assert(SCIPrationalIsAbsInfinity(r1));

   /* inf printing */
   SCIPrationalSetString(r1, "inf");
   SCIPrationalToString(r1, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print inf: %s \n", buf);
   SCIPrationalSetString(r1, "-inf");
   SCIPrationalToString(r1, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print -inf: %s \n", buf);

   /* multiplication with inf */
   SCIPrationalMultReal(r2, r1, 0);
   cr_assert(SCIPrationalIsZero(r2));
   cr_assert(!SCIPrationalIsAbsInfinity(r2));
   SCIPrationalMult(r2, r2, r1);
   cr_assert(SCIPrationalIsZero(r2));

   /* mul -inf * -1 */
   SCIPrationalSetInt(r2, -1, 1);
   SCIPrationalMult(rbuf, r1, r2);
   cr_assert(SCIPrationalIsInfinity(rbuf));

   /* adding inf and not-inf */
   SCIPrationalSetString(r1, "-inf");
   SCIPrationalSetReal(r2, 12.3548934);
   SCIPrationalAdd(rbuf, r1, r2);

   cr_assert(SCIPrationalIsNegInfinity(rbuf));
   SCIPrationalSetString(r1, "inf");
   SCIPrationalAdd(rbuf, r1, r2);
   cr_assert(SCIPrationalIsInfinity(rbuf));
   cr_assert(!SCIPrationalIsEqual(r1, r2));

   /* Difference 0.5 - - 0.5 == 1*/
   SCIPrationalSetInt(r1, 1, 2);
   SCIPrationalSetInt(r2, -1, 2);
   SCIPrationalDiff(rbuf, r1, r2);
   SCIPrationalSetInt(r1, 1, 1);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* Diff 1 - 1.5 == -0.5 */
   SCIPrationalDiffReal(r1, r1, 1.5);
   cr_assert(SCIPrationalIsEqual(r1, r2));

   /* RelDiff -5 / -0.5 == -4.5 / 5 */
   SCIPrationalMultReal(r1, r1, 10);
   SCIPrationalRelDiff(rbuf, r1, r2);
   SCIPrationalSetInt(r1, -9, 10);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* Division (-9/10) / (-0.5) == 18/10 */
   SCIPrationalDiv(rbuf, r1, r2);
   SCIPrationalSetInt(r1, 18, 10);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   SCIPrationalDivReal(rbuf, r1, 2);
   SCIPrationalSetInt(r1, 9, 10);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* SCIPrationalAddProd 9/10 += 9/10 * -0.5 == 9/20 */
   SCIPrationalAddProd(rbuf, r1, r2);
   SCIPrationalSetInt(r1, 9, 20);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* SCIPrationalDiffProd 9/20 -= 9/20 * -0.5 == 27/40 */
   SCIPrationalDiffProd(rbuf, r1, r2);
   SCIPrationalSetInt(r1, 27, 40);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* Negation */
   SCIPrationalNegate(rbuf, r1);
   SCIPrationalSetInt(r1, -27, 40);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* Abs */
   SCIPrationalAbs(rbuf, r1);
   SCIPrationalSetInt(r1, 27, 40);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* Invert */
   SCIPrationalInvert(rbuf, r1);
   SCIPrationalSetInt(r1, 40, 27);
   cr_assert(SCIPrationalIsEqual(rbuf, r1));

   /* min/max */
   SCIPrationalMIN(r2, r1, rbuf);
   cr_assert(SCIPrationalIsEqual(rbuf, r2));
   SCIPrationalMAX(r2, r1, rbuf);
   cr_assert(SCIPrationalIsEqual(r1, r2));

   /* comparisons (GT/LT/GE/LE) (only one checed, since use each other)*/
   SCIPrationalInvert(r1, r1);
   cr_assert(SCIPrationalIsGT(rbuf, r1));
   SCIPrationalSetString(r1, "inf");
   cr_assert(SCIPrationalIsLT(rbuf, r1));
   SCIPrationalSetString(rbuf, "-inf");
   cr_assert(SCIPrationalIsLT(rbuf, r1));
   cr_assert(SCIPrationalIsLT(rbuf, r2));

   /* is positive/negative */
   cr_assert(SCIPrationalIsPositive(r1));
   cr_assert(SCIPrationalIsNegative(rbuf));

   /* SCIPrationalIsIntegral */
   SCIPrationalSetInt(r1, 15, 5);
   SCIPrationalSetInt(r2, 15, 7);
   cr_assert(SCIPrationalIsIntegral(r1));
   cr_assert(!SCIPrationalIsIntegral(r2));

   /* comparing fp and rat */
   SCIPrationalSetString(r2, "123646/1215977400");
   cr_assert(!SCIPrationalIsEqualReal(r2, SCIPrationalApproxReal(r2)));

   SCIPrationalSetString(r1, "1/3");
   cr_assert(SCIPrationalIsLT(r2, r1));

   /* test rounding/ fp approximation */
   SCIPrationalAdd(rbuf, r1, r1);
   doub = SCIPrationalApproxReal(r1);
   cr_log_info("rounding nearest:     %.17e \n", SCIPrationalRoundReal(rbuf, SCIP_R_ROUND_NEAREST));
   cr_log_info("rounding first:       %.17e \n", 2 * doub);
   cr_assert_leq(2 * doub, SCIPrationalApproxReal(rbuf));
   cr_assert(!SCIPrationalIsFpRepresentable(rbuf));

   cr_assert(!SCIPrationalIsIntegral(rbuf));
   SCIPrationalMultReal(rbuf, r1, 3);
   cr_assert(SCIPrationalIsIntegral(rbuf));

   /* to string conversions (missing a few) */
   SCIPrationalToString(r1, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print 1/3: %s \n", buf);
   SCIPrationalSetString(r1, "3/4");
   SCIPrationalRoundInteger(&intval, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(intval, 0);
   SCIPrationalSetString(r1, "3/4");
   SCIPrationalRoundInteger(&intval, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(intval, 1);
   SCIPrationalSetString(r1, "-5/4");
   SCIPrationalRoundInteger(&intval, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(intval, -2);
   SCIPrationalSetString(r1, "-9/4");
   SCIPrationalRoundInteger(&intval, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(intval, -2);
}

Test(rationals, overflows, .description = "test conversion methods with huge rationals")
{
   SCIP_Longint testlong = 1;
   SCIP_Longint reslong = 0;

   // round -1/2 up/down/nearest
   SCIPrationalSetInt(r1, testlong, -(testlong * 2));
   SCIPrationalRound(r2, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(SCIPrationalApproxReal(r2), 0);
   SCIPrationalRoundInteger(&reslong, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(reslong, 0);

   SCIPrationalRound(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(SCIPrationalApproxReal(r2), -1);
   SCIPrationalRoundInteger(&reslong, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(reslong, -1);

   SCIPrationalRound(r2, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(SCIPrationalApproxReal(r2), -1);
   SCIPrationalRoundInteger(&reslong, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(reslong, -1);

   // round 1/2 up/down/nearest
   SCIPrationalSetInt(r1, testlong, (testlong * 2));
   SCIPrationalRound(r2, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(SCIPrationalApproxReal(r2), 1);
   SCIPrationalRoundInteger(&reslong, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(reslong, 1);

   SCIPrationalRound(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(SCIPrationalApproxReal(r2), 0);
   SCIPrationalRoundInteger(&reslong, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(reslong, 0);

   SCIPrationalRound(r2, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(SCIPrationalApproxReal(r2), 1);
   SCIPrationalRoundInteger(&reslong, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(reslong, 1);

   // round -1/3 up/down/nearest
   SCIPrationalSetInt(r1, testlong, -(testlong * 3));
   SCIPrationalRound(r2, r1, SCIP_R_ROUND_UPWARDS);
   cr_assert_eq(SCIPrationalApproxReal(r2), 0);
   SCIPrationalRound(r2, r1, SCIP_R_ROUND_DOWNWARDS);
   cr_assert_eq(SCIPrationalApproxReal(r2), -1);
   SCIPrationalRound(r2, r1, SCIP_R_ROUND_NEAREST);
   cr_assert_eq(SCIPrationalApproxReal(r2), 0);
}

Test(rationals, arrays, .description = "tests rational array methods")
{
   SCIP_RATIONALARRAY* ratar2;

   cr_log_info("testing rational array methods \n");

   // test getter and setters
   SCIPrationalSetInt(r1, 2, 5);
   SCIPrationalarraySetVal(ratar, 5, r1);
   SCIPrationalarrayGetVal(ratar, 5, r2);
   cr_assert(SCIPrationalIsEqual(r1, r2));

   SCIPrationalarrayGetVal(ratar, 4, r2);
   cr_assert(SCIPrationalIsZero(r2));

   // test get/set
   // ensure no by-ref passing happens
   SCIPrationalarraySetVal(ratar, 10, r1);
   // array is : 0.4 0 0 0 0 0.4
   SCIPrationalarrayGetVal(ratar, 1, r2);
   SCIPrationalarrayGetVal(ratar, 10, r2);
   cr_log_info("Two rat");
   cr_log_info("Test arrray printing: \n ");
   // cr_assert(SCIP_OKAY == SCIPrationalarrayPrint(ratar));
   cr_assert(SCIPrationalIsEqual(r1, r2));

   // test incval
   SCIPrationalSetInt(r2, 1, 5);
   SCIPrationalSetInt(r1, 2, 5);
   SCIPrationalarraySetVal(ratar, 7, r2);
   SCIPrationalarrayIncVal(ratar, 7, r1);

   SCIPrationalAdd(r1, r1, r2);
   SCIPrationalarrayGetVal(ratar, 7, r2);
   cr_assert(SCIPrationalIsEqual(r1, r2));

   // test max/minidx
   cr_assert_eq(10, SCIPrationalarrayGetMaxIdx(ratar));
   cr_assert_eq(5, SCIPrationalarrayGetMinIdx(ratar));

   cr_log_info("Test arrray copying: \n ");
   SCIPrationalarrayCopy(&ratar2, blkmem, ratar);
   for( int i = SCIPrationalarrayGetMinIdx(ratar); i < SCIPrationalarrayGetMaxIdx(ratar); i++ )
   {
      SCIPrationalarrayGetVal(ratar, i, r1);
      SCIPrationalarrayGetVal(ratar2, i, r2);
      cr_assert(SCIPrationalIsEqual(r1, r2));
   }

   SCIPrationalarrayFree(&ratar2, blkmem);
}
#endif
