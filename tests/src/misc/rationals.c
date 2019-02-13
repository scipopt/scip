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
   rno = RcreateNoMem();
   rstring = RcreateString(blkmem, "1/3");
   rreal = RcreateReal(blkmem, 1.234236);
   rinte = RcreateInt(blkmem, 1, 3);
   rgmp = RcreateGMP(blkmem, gmpr);
   rarray = RcreateArray(blkmem, 10);

   /* test setter methods */
   RsetInt(rno, testint, 1);
   cr_assert_eq(RgetRealApprox(rno), testint, "setting from and converting back to int did not give same result");
   /* set to fp number */
   RsetReal(rno, testreal);
   cr_assert_eq(RgetRealApprox(rno), testreal, "setting from and converting back to real did not give same result");
   cr_assert(RisFpRepresentable(rno), "fp-rep number not detected as representable");

   /* set to string rep */
   RsetString(rno, "1/3");
   cr_assert(!RisFpRepresentable(rno), "non-fp-rep number not detected as non-representable");
   cr_assert(!RisEqualReal(rno, RgetRealApprox(rno)), "approximation of 1/3 should not be the same");

   /* test rounding */
   Rprint(rno);
   printf("printing test approx: %.17e \n", RgetRealApprox(rno));
   printf("rounding down:        %.17e \n", RgetRealRelax(rno, SCIP_ROUND_DOWNWARDS));
   printf("rounding up:          %.17e \n", RgetRealRelax(rno, SCIP_ROUND_UPWARDS));
   printf("rounding nearest:     %.17e \n", RgetRealRelax(rno, SCIP_ROUND_NEAREST));

   /* test that rounding down is lt rounding up */
   cr_assert_lt(RgetRealRelax(rno, SCIP_ROUND_DOWNWARDS), RgetRealRelax(rno, SCIP_ROUND_UPWARDS), "rounding down should be lt rounding up");

   /* test gmp conversion */
   RsetGMP(rno, gmpr);
   cr_assert_eq(RgetRealApprox(rno), mpq_get_d(gmpr), "gmp and Rational should be the same ");
   cr_assert(0 == mpq_cmp(gmpr, *RgetGMP(rno)));

   /* delete the rationals */
   RdeleteNoMem(&rno);
   Rdelete(blkmem, &rstring);
   Rdelete(blkmem, &rreal);
   Rdelete(blkmem, &rinte);
   Rdelete(blkmem, &rgmp);
   RdeleteArray(blkmem, &rarray, 10);

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
   SCIP_Real doub;
   char buf[SCIP_MAXSTRLEN];

   BMS_BLKMEM* blkmem = BMScreateBlockMemory(1, 10);

   infpos = RcreateString(blkmem, "inf");
   infneg = RcreateString(blkmem, "-inf");

   r1 = RcreateReal(blkmem, 12.3548934);
   r2 = RcreateString(blkmem, "123646/1215977400");
   r3 = RcreateString(blkmem, "1/3");
   r4 = RcreateString(blkmem, "1/10");
   r5 = Rcreate(blkmem);

   doub = 12.3548933;

   /* test infinity values */
   cr_assert(RisInfinity(infpos));
   cr_assert(!RisInfinity(infneg));
   cr_assert(RisNegInfinity(infneg));
   cr_assert(RisAbsInfinity(infneg));
   cr_assert(!RisNegInfinity(infpos));

   Radd(r5, r1, infpos);
   cr_assert(RisInfinity(r5));
   Radd(r5, r1, infneg);
   cr_assert(RisNegInfinity(r5));

   cr_assert(!RisEqual(r1, r2));
   cr_assert(!RisEqualReal(r2, doub));
   cr_assert(RisLT(r2, r3));

   Radd(r5, r3, r3);
   doub = RgetRealApprox(r3);
   Rprint(r5);
   printf("rounding nearest:     %.17e \n", RgetRealRelax(r5, SCIP_ROUND_NEAREST));
   printf("rounding first:       %.17e \n", 2 * doub);
   cr_assert_leq(2 * doub, RgetRealApprox(r5));

   RmultReal(infneg, infpos, 0);
   cr_assert(RisZero(infneg));
   cr_assert(!RisAbsInfinity(infneg));

   cr_assert(!RisFpRepresentable(r3));
   cr_assert(RisFpRepresentable(r1));
   cr_assert(!RisIntegral(r3));
   RmultReal(r5, r3, 3);
   cr_assert(RisIntegral(r5));

   RtoString(r3, buf);
   printf("Test print 1/3: %s \n", buf);

   RtoString(infpos, buf);
   printf("Test print inf: %s \n", buf);

   RsetString(infneg, "-inf");
   RtoString(infneg, buf);
   printf("Test print -inf: %s \n", buf);

   Rdelete(blkmem, &r1);
   Rdelete(blkmem, &r2);
   Rdelete(blkmem, &r3);
   Rdelete(blkmem, &r4);
   Rdelete(blkmem, &r5);
   Rdelete(blkmem, &infpos);
   Rdelete(blkmem, &infneg);
}