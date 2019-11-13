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
   RatPrint(rno);
   printf("%.17e \n", RatApproxReal(rno));
   cr_assert(RatIsFpRepresentable(rno), "fp number 0.124691234 not fp representable");

   /* set to string rep */
   RatSetString(rno, "1/3");
   cr_assert(!RatIsFpRepresentable(rno), "non-fp-rep number not detected as non-representable");
   cr_assert(!RatIsEqualReal(rno, RatApproxReal(rno)), "approximation of 1/3 should not be the same");

   /* test rounding */
   RatPrint(rno);
   printf("printing test approx: %.17e \n", RatApproxReal(rno));
   printf("rounding down:        %.17e \n", RatRoundReal(rno, SCIP_ROUND_DOWNWARDS));
   printf("rounding up:          %.17e \n", RatRoundReal(rno, SCIP_ROUND_UPWARDS));
   printf("rounding nearest:     %.17e \n", RatRoundReal(rno, SCIP_ROUND_NEAREST));

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

Test(rationals, rounding, .description = "tests rational rounding speed")
{
   clock_t startt, endt;
   int niterations = 1000000;
   int i;
   int nrep = 0;
   double runtime = 0;
   double addval;

   srand((unsigned int)time(NULL));

   SCIP_Rational* r;  
   SCIP_Rational* r2;

   RatCreate(&r);
   RatCreate(&r2);

   printf("Testing time for performing tasks %d times\n", niterations);

   startt = clock();
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
   }
   endt = clock();
   printf(" cpu time used for setting from real: %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC);

   runtime = 0;
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      startt = clock();
      nrep += RatIsFpRepresentable(r) ? 1 : 0;
      endt = clock();
      runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   }
   endt = clock();
   printf(" cpu time used for checking fp-rep: %e \n", runtime);
   cr_assert(nrep == niterations, "error");

   runtime = 0;
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      startt = clock();
      addval += RatRoundReal(r, SCIP_ROUND_DOWNWARDS);
      addval += RatRoundReal(r, SCIP_ROUND_UPWARDS);
      endt = clock();
      runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   }
   printf(" cpu time used for rounding: %e, addval %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC, addval);


   runtime = 0;
   addval = 0;

   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      startt = clock();
      addval += RatApproxReal(r);
      endt = clock();
      runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   }
   printf(" cpu time used for apporx: %e, addval %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC, addval);

   runtime = 0;
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      RatSetReal(r2, ((float)rand())/RAND_MAX);
      startt = clock();
      RatAdd(r, r, r2);
      endt = clock();
      runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   }
   printf(" cpu time used for adding: %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC);

   runtime = 0;
   for( i = 0; i < niterations; ++i )
   {
      RatSetReal(r, ((float)rand())/RAND_MAX);
      RatSetReal(r2, ((float)rand())/RAND_MAX);
      startt = clock();
      RatMult(r, r, r2);
      endt = clock();
      runtime += ((double) (endt - startt)) / CLOCKS_PER_SEC;
   }
   printf(" cpu time used for multiplication: %e \n", ((double) (endt - startt)) / CLOCKS_PER_SEC);

   RatFree(&r);
   RatFree(&r2);
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
   RatPrint(r5);
   printf("rounding nearest:     %.17e \n", RatRoundReal(r5, SCIP_ROUND_NEAREST));
   printf("rounding first:       %.17e \n", 2 * doub);
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
   printf("Test print 1/3: %s \n", buf);

   RatToString(infpos, buf, SCIP_MAXSTRLEN);
   printf("Test print inf: %s \n", buf);

   RatSetString(infneg, "-inf");
   RatToString(infneg, buf, SCIP_MAXSTRLEN);
   printf("Test print -inf: %s \n", buf);

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