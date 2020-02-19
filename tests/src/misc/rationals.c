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

/**@file   rationals.c
 * @brief  unit tests for double double arithmetic
 * @author Robert Lion Gottwald;
 */

#include "scip/rational.h"
#include "scip/type_clock.h"
#include "scip/scip_randnumgen.h"
#include <time.h>
#include <gmp.h>
#include "include/scip_test.h"


/** GLOBAL VARIABLES **/

/* corresponds to exact rational number 5223617350827361/151115727451828646838272 */
volatile double a = 0.000000034567;

/* corresponds to exact rational number 1320494817241769/562949953421312 */
volatile double b = 2.34567;

/* corresponds to exact rational number 6628035890302157/536870912 */
volatile double c = 12345678.9;

/* TEST SUITE */
TestSuite(rationals);

Test(rationals, creation, .description = "tests all the different methods to create rationals")
{
   SCIP_Rational* rno;
   SCIP_Rational* rstring;
   SCIP_Rational* rreal;
   SCIP_Rational* rinte;
   SCIP_Rational* rgmp;
   SCIP_Rational** rarray;
   int testint = 100345;
   SCIP_Real testreal = 1235.235690;
   mpq_t gmpr;

   BMS_BLKMEM* blkmem = BMScreateBlockMemory(1, 10);
   BMS_BUFMEM* buffmem = BMScreateBufferMemory(1.2, 4, FALSE);
   mpq_init(gmpr);
   mpq_set_d(gmpr, 1.2345);

   /* create some rationals with different methods*/
   RatCreate(&rno);
   RatCreateBlock(blkmem, &rstring);
   RatSetString(rstring, "1/3");
   RatCreateBlock(blkmem, &rreal);
   RatSetReal(rreal, 1.234236);
   RatCreateBlock(blkmem, &rinte);
   RatSetInt(rinte, 1, 3);
   RatCreateGMP(blkmem, &rgmp, gmpr);
   RatCreateBlockArray(blkmem, &rarray, 10);

   /* test setter methods */
   RatSetInt(rno, testint, 1);
   cr_assert_eq(RatApproxReal(rno), testint, "setting from and converting back to int did not give same result");
   /* set to fp number */
   RatSetReal(rno, testreal);
   cr_assert_eq(RatApproxReal(rno), testreal, "setting from and converting back to real did not give same result");
   cr_assert(RatIsFpRepresentable(rno), "fp-rep number not detected as representable");


   /* set to string rep */
   RatSetReal(rno, 0.1246912);
   cr_log_info("Test printing 0.1246912 %s", RatGetString(rno));
   cr_log_info("%.17e \n", RatApproxReal(rno));
   cr_assert(RatIsFpRepresentable(rno), "fp number 0.124691234 not fp representable");

   /* set to string rep */
   RatSetString(rno, "1/3");
   cr_assert(!RatIsFpRepresentable(rno), "non-fp-rep number not detected as non-representable");
   cr_assert(!RatIsEqualReal(rno, RatApproxReal(rno)), "approximation of 1/3 should not be the same");

   /* test rounding */
   cr_log_info("printing test approx: %.17e \n", RatApproxReal(rno));
   cr_log_info("rounding down:        %.17e \n", RatRoundReal(rno, SCIP_ROUND_DOWNWARDS));
   cr_log_info("rounding up:          %.17e \n", RatRoundReal(rno, SCIP_ROUND_UPWARDS));
   cr_log_info("rounding nearest:     %.17e \n", RatRoundReal(rno, SCIP_ROUND_NEAREST));

   /* test that rounding down is lt rounding up */
   cr_assert_lt(RatRoundReal(rno, SCIP_ROUND_DOWNWARDS), RatRoundReal(rno, SCIP_ROUND_UPWARDS), "rounding down should be lt rounding up");

   /* test gmp conversion */
   RatSetGMP(rno, gmpr);
   cr_assert_eq(RatApproxReal(rno), mpq_get_d(gmpr), "gmp and Rational should be the same ");
   cr_assert(0 == mpq_cmp(gmpr, *RatGetGMP(rno)));

   /* delete the rationals */
   RatFree(&rno);
   RatFreeBlock(blkmem, &rstring);
   RatFreeBlock(blkmem, &rreal);
   RatFreeBlock(blkmem, &rinte);
   RatFreeBlock(blkmem, &rgmp);
   RatFreeBlockArray(blkmem, &rarray, 10);
}

Test(rationals, arithmetic, .description = "tests rational arithmetic methods")
{
   SCIP_Rational* r1;
   SCIP_Rational* r2;
   SCIP_Rational* r3;
   SCIP_Rational* r4;
   SCIP_Rational* r5;
   SCIP_Rational* infpos;
   SCIP_Rational* infneg;
   long  intval;
   SCIP_Real doub;
   char buf[SCIP_MAXSTRLEN];

   BMS_BLKMEM* blkmem = BMScreateBlockMemory(1, 10);

   RatCreateString(blkmem, &infpos, "inf");
   RatCreateString(blkmem, &infneg, "-inf");

   RatCreateBlock(blkmem, &r1);
   RatCreateBlock(blkmem, &r2);
   RatCreateBlock(blkmem, &r3);
   RatCreateBlock(blkmem, &r4);
   RatCreateBlock(blkmem, &r5);

   RatSetReal(r1, 12.3548934);
   RatSetString(r2, "123646/1215977400");
   RatSetString(r3, "1/3");
   RatSetString(r4, "1/10");

   doub = 12.3548933;

   /* test infinity values */
   cr_assert(RatIsInfinity(infpos));
   cr_assert(!RatIsInfinity(infneg));
   cr_assert(RatIsNegInfinity(infneg));
   cr_assert(RatIsAbsInfinity(infneg));
   cr_assert(!RatIsNegInfinity(infpos));

   RatAdd(r5, r1, infpos);
   cr_assert(RatIsInfinity(r5));
   RatAdd(r5, r1, infneg);
   cr_assert(RatIsNegInfinity(r5));

   cr_assert(!RatIsEqual(r1, r2));
   cr_assert(!RatIsEqualReal(r2, doub));
   cr_assert(RatIsLT(r2, r3));

   RatAdd(r5, r3, r3);
   doub = RatApproxReal(r3);
   cr_log_info("rounding nearest:     %.17e \n", RatRoundReal(r5, SCIP_ROUND_NEAREST));
   cr_log_info("rounding first:       %.17e \n", 2 * doub);
   cr_assert_leq(2 * doub, RatApproxReal(r5));

   RatMultReal(infneg, infpos, 0);
   cr_assert(RatIsZero(infneg));
   cr_assert(!RatIsAbsInfinity(infneg));

   cr_assert(!RatIsFpRepresentable(r3));
   cr_assert(RatIsFpRepresentable(r1));
   cr_assert(!RatIsIntegral(r3));
   RatMultReal(r5, r3, 3);
   cr_assert(RatIsIntegral(r5));

   RatToString(r3, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print 1/3: %s \n", buf);

   RatToString(infpos, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print inf: %s \n", buf);

   RatSetString(infneg, "-inf");
   RatToString(infneg, buf, SCIP_MAXSTRLEN);
   cr_log_info("Test print -inf: %s \n", buf);

   RatSetString(r1, "3/4");
   RatRoundInteger(&intval, r1, SCIP_ROUND_DOWNWARDS);
   cr_assert_eq(intval, 0);
   RatSetString(r1, "3/4");
   RatRoundInteger(&intval, r1, SCIP_ROUND_UPWARDS);
   cr_assert_eq(intval, 1);

   RatSetString(r1, "-5/4");
   RatRoundInteger(&intval, r1, SCIP_ROUND_DOWNWARDS);
   cr_assert_eq(intval, -2);
   RatSetString(r1, "-9/4");
   RatRoundInteger(&intval, r1, SCIP_ROUND_UPWARDS);
   cr_assert_eq(intval, -2);

   RatFreeBlock(blkmem, &r1);
   RatFreeBlock(blkmem, &r2);
   RatFreeBlock(blkmem, &r3);
   RatFreeBlock(blkmem, &r4);
   RatFreeBlock(blkmem, &r5);
   RatFreeBlock(blkmem, &infpos);
   RatFreeBlock(blkmem, &infneg);
}

Test(rationals, arrays, .description = "tests rational array methods")
{
   SCIP_RATIONALARRAY* ar1;
   SCIP_RATIONALARRAY* ar2;
   SCIP_Rational* rat;
   SCIP_Rational* rat2;
   BMS_BLKMEM* blkmem = BMScreateBlockMemory(1, 10);

   cr_log_info("testing rational array methods \n");

   SCIPrationalarrayCreate(&ar1, blkmem);
   RatCreateBlock(blkmem, &rat);
   RatCreateBlock(blkmem, &rat2);

   // test getter and setters
   RatSetInt(rat, 2, 5);
   SCIPrationalarraySetVal(ar1, 5, rat);
   SCIPrationalarrayGetVal(ar1, 5, rat2);
   cr_assert(RatIsEqual(rat, rat2));

   SCIPrationalarrayGetVal(ar1, 4, rat2);
   cr_assert(RatIsZero(rat2));

   // ensure no by-ref passing happens
   SCIPrationalarraySetVal(ar1, 10, rat);
   // array is : 0.4 0 0 0 0 0.4
   SCIPrationalarrayGetVal(ar1, 1, rat2);
   SCIPrationalarrayGetVal(ar1, 10, rat2);
   cr_log_info("Two rat");
   cr_log_info("Test arrray printing: \n ");
   cr_assert(SCIP_OKAY == SCIPrationalArrayPrint(ar1));
   cr_assert(RatIsEqual(rat2, rat));

   RatSetInt(rat2, 1, 5);
   RatSetInt(rat, 2, 5);
   SCIPrationalarraySetVal(ar1, 7, rat2);
   SCIPrationalarrayIncVal(ar1, 7, rat);

   RatAdd(rat, rat, rat2);
   SCIPrationalarrayGetVal(ar1, 7, rat2);
   cr_log_info("Increased val is %s", RatGetString(rat2));
   cr_assert(RatIsEqual(rat, rat2));

   cr_assert_eq(10, SCIPrationalarrayGetMaxIdx(ar1));
   cr_assert_eq(5, SCIPrationalarrayGetMinIdx(ar1));

   cr_log_info("Test arrray copying: \n ");
   SCIPrationalarrayCopy(&ar2, blkmem, ar1);
   for( int i = SCIPrationalarrayGetMinIdx(ar1); i < SCIPrationalarrayGetMaxIdx(ar1); i++ )
   {
      SCIPrationalarrayGetVal(ar1, i, rat);
      SCIPrationalarrayGetVal(ar2, i, rat2);
      cr_assert(RatIsEqual(rat, rat2));
   }

}