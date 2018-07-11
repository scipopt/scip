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

/**@file   select.c
 * @brief  unit tests for selection of unweighted and weighted median
 * @author Gregor Hendel
 */

#include "scip/pub_misc.h"
#include "scip/scip.h"

#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
static SCIP_RANDNUMGEN* randgen;
static SCIP* scip;
static unsigned int randomseed = 42;

/* TEST SUITE */
static
void setup(void)
{
   SCIPcreate(&scip);
   SCIPcreateRandom(scip, &randgen, randomseed, TRUE);
}

static
void teardown(void)
{
   SCIPfreeRandom(scip, &randgen);
   SCIPfree(&scip);
}

TestSuite(select, .init = setup, .fini = teardown);

/* TESTS  */
Test(select, create_and_free)
{
   /* calls setup and teardown */
}

#define ARRAYMEMSIZE 700
Test(select, random_permutation, .description = "tests selection on a bunch of random permutations of the integers 1...n")
{
   int len = ARRAYMEMSIZE;
   int key[ARRAYMEMSIZE];
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
      cr_assert_eq(key[i], i, "Wrong key selected: %d ~= %d", key[i], i);

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
         cr_assert_lt(k, end, "Element %d is not in the right partition [%d,%d]\n", j, start, end);
      }
   }
}

/* For SCIPselectWeightedInt, we test the boundary case when the capacity is 0. In this case, the algorithm should always return
 * the first element (since the algorithm partially sorts the given array). To test this, we build an all 1 array (key) and modify
 * its i-th position to be 0, call SCIPselectWeightedInt and expect key[0] to be 0
 */
Test(select, zero_capacity, .description = "tests if weighted median selection always returns the minimum element when the capacity is zero")
{
   int i;
   int key[ARRAYMEMSIZE];

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
      cr_assert_eq(medianpos, 0, "Selection process found median at pos %d, should be 0", medianpos);
      cr_assert_eq(key[0], 0, "Selection process selected wrong median %d at pos 0, should be 0", key[0]);

      /* start next iteration with an all 1-s array */
      key[0] = 1;
   }
}
