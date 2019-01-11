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

   /* create some rationals */
   rno = RcreateNoMem();
   rstring = RcreateString(blkmem, "1/3");
   rreal = RcreateReal(blkmem, 1.234236);
   rinte = RcreateInt(blkmem, 1, 3);
   rgmp = RcreateGMP(blkmem, gmpr);
   rarray = RcreateArray(blkmem, 10);

   /* test setter methods */
   RsetInt(rno, testint, 1);
   cr_assert_eq(RgetRealApprox(rno), testint, "setting from and converting back to int did not give same result");
   RsetReal(rno, testreal);
   cr_assert_eq(RgetRealApprox(rno), testreal, "setting from and converting back to real did not give same result");
   cr_assert(RisFpRepresentable(rno), "fp-rep number not detected as representable");
   RsetString(rno, "1/3");
   cr_assert(!RisFpRepresentable(rno), "non-fp-rep number not detected as non-representable");
   cr_assert(!RisEqualReal(rno, RgetRealApprox(rno)), "approximation of 1/3 should not be the same");

   Rprint(rno);
   printf("printing test approx: %.16e \n", RgetRealApprox(rno));
   printf("rounding down: %.17e \n", RgetRealRelax(rno, SCIP_ROUND_DOWNWARDS));
   printf("rounding up:   %.17e \n", RgetRealRelax(rno, SCIP_ROUND_UPWARDS));
   printf("rounding nearest: %.17e \n", RgetRealRelax(rno, SCIP_ROUND_NEAREST));

   cr_assert_lt(RgetRealRelax(rno, SCIP_ROUND_DOWNWARDS), RgetRealRelax(rno, SCIP_ROUND_UPWARDS), "rounding down should be lt rounding up");

   RsetGMP(rno, gmpr);
   cr_assert_eq(RgetRealApprox(rno), mpq_get_d(gmpr), "gmp and Rational should be the same ");
   cr_assert(0 == mpq_cmp(gmpr, *RgetGMP(rno)));

   RdeleteNoMem(&rno);
   Rdelete(blkmem, &rstring);
   Rdelete(blkmem, &rreal);
   Rdelete(blkmem, &rinte);
   Rdelete(blkmem, &rgmp);
   RdeleteArray(blkmem, &rarray, 10);

}