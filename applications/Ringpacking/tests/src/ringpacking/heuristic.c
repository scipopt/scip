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

/**@file   heuristic.c
 * @brief  unit test for testing packing heuristic
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"

#include "probdata_rpa.h"

#include "include/scip_test.h"

static SCIP* scip;

/** setup of test run */
static
void setup(void)
{
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
}

/** deinitialization method */
static
void teardown(void)
{
   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   TEST_ASSERT_NULL(scip);
   TEST_ASSERT_EQUAL(0, BMSgetMemoryUsed(), "There is a memory leak!!");
}

void setUp(void) { setup(); }
void tearDown(void) { teardown(); }

/* tests heuristic for a single circle */
void test_heuristic_single(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Bool ispacked;
   SCIP_Real x;
   SCIP_Real y;
   int elements = 0;
   int npacked;

   /* pack element into rectangle */
   SCIPpackCirclesGreedy(scip, &rext, &x, &y, -1.0, 2.0, 2.0, &ispacked, &elements, 1,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 1);
   SOFT_ASSERT(ispacked);
   SOFT_ASSERT(x == 1.0);
   SOFT_ASSERT(y == 1.0);

   /* pack element into ring */
   SCIPpackCirclesGreedy(scip, &rext, &x, &y, 2.0, -1.0, -1.0, &ispacked, &elements, 1,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 1);
   SOFT_ASSERT(ispacked);
   SOFT_ASSERT(x == -1.0);
   SOFT_ASSERT(y == 0.0);
}

/* packs two different circle types into a rectangle */
void test_heuristic_rectangle_two_types(void)
{
   SCIP_Real rexts[2] = {3.0, 0.45};
   SCIP_Real xs[5];
   SCIP_Real ys[5];
   SCIP_Bool ispacked[5];
   int elements[5] = {0, 1, 1, 1, 1};
   int npacked;
   int k;

   SCIPpackCirclesGreedy(scip, rexts, xs, ys, -1.0, 6.0, 6.0, ispacked, elements, 5,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 5);

   for( k = 0; k < 5; ++k )
      SOFT_ASSERT(ispacked[k]);
}

/* test optimal packing with two identical circles */
void test_heuristic_opt_two(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[2];
   SCIP_Real ys[2];
   SCIP_Bool ispacked[2];
   int elements[2] = {0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 2.0, -1.0, -1.0, ispacked, elements, 2,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 2);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 3.414213562373095, 3.414213562373095, ispacked, elements, 2,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 2);
}

/* test optimal packing with three identical circles */
void test_heuristic_opt_three(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[3];
   SCIP_Real ys[3];
   SCIP_Bool ispacked[3];
   int elements[3] = {0, 0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 2.1547005383792515, -1.0, -1.0, ispacked, elements, 3,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 3);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 3.9318516525781364, 3.9318516525781364, ispacked, elements, 3,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 3);
}

/* test optimal packing with four identical circles */
void test_heuristic_opt_four(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[4];
   SCIP_Real ys[4];
   SCIP_Bool ispacked[4];
   int elements[4] = {0, 0, 0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 2.414213562373095, -1.0, -1.0, ispacked, elements, 4,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 4);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 4.0, 4.0, ispacked, elements, 4,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 3); /* note that greedy fails for this example */
}

/* test optimal packing with five identical circles */
void test_heuristic_opt_five(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[5];
   SCIP_Real ys[5];
   SCIP_Bool ispacked[5];
   int elements[5] = {0, 0, 0, 0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 2.7013016167040798, -1.0, -1.0, ispacked, elements, 5,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 4); /* note that the greedy packing fails for the case of 5 circles */

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 4.82842712474619, 4.82842712474619, ispacked, elements, 5,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 4); /* note that greedy fails for this example */
}

/* test optimal packing with fix identical circles */
void test_heuristic_opt_six(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[6];
   SCIP_Real ys[6];
   SCIP_Bool ispacked[6];
   int elements[6] = {0, 0, 0, 0, 0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 3.0, -1.0, -1.0, ispacked, elements, 6,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 6);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 5.328201177351374, 5.328201177351374, ispacked, elements, 6,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 5); /* note that greedy fails for this example */
}

/* test optimal packing with seven identical circles */
void test_heuristic_opt_seven(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[7];
   SCIP_Real ys[7];
   SCIP_Bool ispacked[7];
   int elements[7] = {0, 0, 0, 0, 0, 0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 3.0, -1.0, -1.0, ispacked, elements, 7,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked == 7);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 5.732050807568877, 5.732050807568877, ispacked, elements, 7,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 6); /* note that greedy fails for this example */
}

/* test optimal packing with twenty identical circles */
void test_heuristic_opt_twenty(void)
{
   SCIP_Real rext = 1.0;
   SCIP_Real xs[20];
   SCIP_Real ys[20];
   SCIP_Bool ispacked[20];
   int elements[20] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, 5.123, -1.0, -1.0, ispacked, elements, 20,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 18);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, &rext, xs, ys, -1.0, 8.978083352821738, 8.978083352821738, ispacked, elements, 20,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 18); /* note that greedy fails for this example */
}

/* tests a more complicated example */
void test_heuristic_complex_1(void)
{
   SCIP_Real rexts[3] = {2.0, 1.0, 0.3};
   SCIP_Real xs[24];
   SCIP_Real ys[24];
   SCIP_Bool ispacked[24];
   int elements[24] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
   int npacked;

   /* pack into a ring */
   SCIPpackCirclesGreedy(scip, rexts, xs, ys, 4.0, -1.0, -1.0, ispacked, elements, 24,
      SCIP_PATTERNTYPE_CIRCULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 24);

   /* pack into a square */
   SCIPpackCirclesGreedy(scip, rexts, xs, ys, -1.0, 7.7, 7.7, ispacked, elements, 24,
      SCIP_PATTERNTYPE_RECTANGULAR, &npacked, 0);
   SOFT_ASSERT(npacked >= 24);
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_heuristic_single);
   RUN_TEST(test_heuristic_rectangle_two_types);
   RUN_TEST(test_heuristic_opt_two);
   RUN_TEST(test_heuristic_opt_three);
   RUN_TEST(test_heuristic_opt_four);
   RUN_TEST(test_heuristic_opt_five);
   RUN_TEST(test_heuristic_opt_six);
   RUN_TEST(test_heuristic_opt_seven);
   RUN_TEST(test_heuristic_opt_twenty);
   RUN_TEST(test_heuristic_complex_1);
   return UNITY_END();
}
