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

/**@file   select.c
 * @brief  unit tests for selection of unweighted and weighted median
 * @author Gregor Hendel
 */

#include "scip/pub_misc.h"
#include "scip/scip.h"

#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
#define ARRAYMEMSIZE 700
static SCIP_RANDNUMGEN* randgen;
static SCIP* scip;
static unsigned int randomseed = 42;
static int key[ARRAYMEMSIZE];
static int items[ARRAYMEMSIZE];
SCIP_Real weights[ARRAYMEMSIZE];


/* TEST SUITE */
static
void setup(void)
{
   int i;

   SCIPcreate(&scip);
   SCIPcreateRandom(scip, &randgen, randomseed, TRUE);

   for( i = 0; i < ARRAYMEMSIZE; ++i )
      items[i] = i;
}

static
void teardown(void)
{
   SCIPfreeRandom(scip, &randgen);
   SCIPfree(&scip);
}

void setUp(void)
{
   setup();
}

void tearDown(void)
{
   teardown();
}

/* TESTS  */
void test_select_create_and_free(void)
{
   /* calls setup and teardown */
}

/** @brief tests selection on a bunch of random permutations of the integers 1...n */
void test_select_random_permutation(void)
{
   int len = ARRAYMEMSIZE;
   int i;
   int j;
   /* initialize key */
   for( j = 0; j < len; ++j )
      key[j] = j;


   /* loop over all positions of the array and check whether the correct element is selected after a random permutation */
   for( i = 0; i < len; ++i )
   {
      int inputkey[ARRAYMEMSIZE];

      SCIPrandomPermuteIntArray(randgen, key, 0, len);
      /* save input permutation for debugging */
      BMScopyMemoryArray(inputkey, key, len);
      SCIPselectInt(key, i, len);

      /* the element key[i] must be the index itself */
      TEST_ASSERT_EQUAL(key[i], i, "Wrong key selected: %d ~= %d", key[i], i);

      /* check if the partial sorting correctly worked */
      for( j = 0; j < len; ++j )
      {
         int k;
         int start = j >= key[i] ? i : 0;
         int end = j >= key[i] ? len : i;

         for( k = start; k < end; ++k )
         {
            if( key[k] == j )
               break;
         }
         TEST_ASSERT_LESS_THAN(k, end, "Element %d is not in the right partition [%d,%d]\n", j, start, end);
      }
   }
}

/** @brief tests if weighted median selection always returns the minimum element when the capacity is zero */
void test_select_zero_capacity(void)
{
   int i;

   /* initialize key to be all ones */
   for( i = 0; i < ARRAYMEMSIZE; ++i )
      key[i] = 1;

   /* test if the minimum element of the key array is selected, no matter where it is at the beginning */
   for( i = 0; i < ARRAYMEMSIZE; ++i )
   {
      int medianpos;
      key[i] = 0;
      medianpos = -1;
      SCIPselectWeightedInt(key, NULL, 0.0, ARRAYMEMSIZE, &medianpos);
      TEST_ASSERT_EQUAL(medianpos, 0, "Selection process found median at pos %d, should be 0", medianpos);
      TEST_ASSERT_EQUAL(key[0], 0, "Selection process selected wrong median %d at pos 0, should be 0", key[0]);

      /* start next iteration with an all 1-s array */
      key[0] = 1;
   }
}

/** @brief tests the correct behavior if a knapsack has integral LP solution */
void test_select_exactIntegerSolution(void)
{
   int medianpos;
   SCIP_Real profits[ARRAYMEMSIZE];
   SCIP_Real capacity = 12;
   int len = 4;

   weights[0] = 6; /* 0 is the worst item and should be last */
   weights[1] = 5;
   weights[2] = 4;
   weights[3] = 3;

   /* the profits already represent the cost/weight */
   profits[0] = 1;
   profits[1] = 2;
   profits[2] = 3;
   profits[3] = 4;

   SCIPselectWeightedDownRealInt(profits, items, weights, capacity, len, &medianpos);

   TEST_ASSERT_EQUAL(medianpos, 3, "Median position %d should be the last item", medianpos);
   TEST_ASSERT_EQUAL(items[3],0, "Wrong last item %d after selection", items[3]);
}


/** @brief second test for the correct behavior for slightly larger coefficients */
void test_select_exactIntegerSolution2(void)
{
   int medianpos;
   SCIP_Real profits[ARRAYMEMSIZE];
   SCIP_Real capacity = 150;
   int len = 5;

   weights[0] = 50;
   weights[1] = 50;
   weights[2] = 50;
   weights[3] = 50;
   weights[4] = 50;

   /* the profits already represent the cost/weight */
   profits[0] = 4;
   profits[1] = 1;
   profits[2] = 3;
   profits[3] = 2;
   profits[4] = 5;

   SCIPselectWeightedDownRealInt(profits, items, weights, capacity, len, &medianpos);

   TEST_ASSERT_EQUAL(medianpos, 3, "Median position %d should be the last item", medianpos);
   TEST_ASSERT_EQUAL(items[3],3, "Wrong median item %d after selection", items[4]);
   TEST_ASSERT_EQUAL(items[4],1);
}

/** @brief realistic example where all items fit that has failed */
void test_select_everythingFits(void)
{
   SCIP_Real testweights[] = {
            2, 2, 1, 1, 3, 1, 1, 2, 2, 1,
            1, 1, 2, 2, 1, 1, 1, 6, 1, 3,
            1, 1, 3, 2, 1, 1, 6, 6, 6, 6,
            6, 6, 6, 6, 6, 6, 6, 6, 6, 6
   };
   SCIP_Real testprofits[] = {
            0.16788881410750284, 0.16788881410750284, 0.33577762821500567, 0.33577762821500567,
            0.11192587607166793, 0.33577762821500567, 0.33577762821500567, 0.16788881410750284,
            0.16788881410750284, 0.33577762821500379, 0.33577762821500567, 0.33577762821500567,
            0.16788881410750284, 0.16788881410750284, 0.33577762821500567, 0.33577762821500567,
            0.33577762821500567, 0.11465122453401815, 0.33577762821500567, 0.11192587607166855,
            0.68790734720410884, 0.33577762821500567, 0.11192587607166855, 0.16788881410750284,
            0.33577762821500567, 0.33577762821500567, 0.05596293803583427, 0.05596293803583427,
            0.05596293803583427, 0.05596293803583427, 0.05596293803583427, 0.05596293803583427,
            0.05596293803583427, 0.05596293803583427, 0.05596293803583427, 0.05596293803583427,
            0.05596293803583427, 0.05596293803583427, 0.05596293803583427, 0.05596293803583427
   };

   int len = 40;
   SCIP_Real capacity = 255;
   int medianpos;

   SCIPselectWeightedDownRealInt(testprofits, items, testweights, capacity, len, &medianpos);

   TEST_ASSERT_EQUAL(medianpos, len, "Medianposition %d != %d", medianpos, len);
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_select_create_and_free);
   RUN_TEST(test_select_random_permutation);
   RUN_TEST(test_select_zero_capacity);
   RUN_TEST(test_select_exactIntegerSolution);
   RUN_TEST(test_select_exactIntegerSolution2);
   RUN_TEST(test_select_everythingFits);
   return UNITY_END();
}
