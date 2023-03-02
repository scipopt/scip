/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   solveknapsackexactly.c
 * @brief  unit tests for exact dynamic programming algorithm for knapsack problem
 * @author Marc Pfetsch
 * @author Gregor Hendel
 */

#include "scip/cons_knapsack.h"
#include "scip/scip.h"

#include "include/scip_test.h"

#define EPS 1e-06
#define MAX_ARRAYLEN 1000

/* GLOBAL VARIABLES */
static SCIP* scip;
static SCIP_Longint weights[MAX_ARRAYLEN];
static int nitems;
static SCIP_Real profits[MAX_ARRAYLEN];
static int items[MAX_ARRAYLEN];
static int solitems[MAX_ARRAYLEN];
static int nonsolitems[MAX_ARRAYLEN];
static SCIP_Longint capacity;
static SCIP_Bool success;
static SCIP_Real solval;
static int nnonsolitems;
static int nsolitems = -1;

/* TEST SUITE */
static
void setup(void)
{
   int i;

   SCIPcreate(&scip);

   /* assign items */
   for( i = 0; i < MAX_ARRAYLEN; ++i )
      items[i] = i;
}

static
void teardown(void)
{
   SCIPfree(&scip);
}

/** run exact knapsack algorithm on the given data  */
static
SCIP_RETCODE solveKnapsack(
   void
   )
{
   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );

   return SCIP_OKAY;
}

/** is set1 contained in set2? */
static
SCIP_Bool checkSetContainment(
   int*                  set1,
   int*                  set2,
   int                   len1,
   int                   len2
   )
{
   int j;
   int sortedset2 [MAX_ARRAYLEN];

   cr_assert(len1 <= MAX_ARRAYLEN);
   cr_assert(len2 <= MAX_ARRAYLEN);

   if( len1 > len2 )
      return FALSE;

   BMScopyMemoryArray(sortedset2, set2, len2);
   SCIPsortInt(sortedset2, len2);

   /* find each item of set 1 in set 2 */
   for( j = 0; j < len1; ++j )
   {
      int pos;
      if( ! SCIPsortedvecFindInt(sortedset2, set1[j], len2, &pos) )
         return FALSE;
   }

   return TRUE;
}

TestSuite(solveknapsackexactly, .init = setup, .fini = teardown);

/* TESTS  */

/*
 * test trivial cases
 */
Test(solveknapsackexactly, test1, .description="check whether the case that all items are redundant is caught correctly")
{
   nitems = 2;
   capacity = 1LL;
   weights[0] = 2LL;
   weights[1] = 1LL;
   profits[0] = 1.0;   /* should not be taken, since capacity is too large */
   profits[1] = -1.0;  /* should not be taken, since profit is negative */

   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 0 );
   cr_assert( nnonsolitems == 2 );
   cr_assert_float_eq(solval, 0.0, EPS);
}

Test(solveknapsackexactly, test2, .description="test whether the correct items are sorted out (weight == 0 or negative profit)")
{
   /* also tests whether the case of all items fitting into the knapsack is caught correctly. */
   nitems = 5;
   capacity = 1LL;
   weights[0] = 0LL;
   weights[1] = 0LL;
   weights[2] = 2LL;
   weights[3] = 2LL;
   weights[4] = 1LL;

   profits[0] = 1.0;   /* should take item 1, since weight is 0 */
   profits[1] = -1.0;  /* should not take item 1, since profit is negative, although weight is 0 */
   profits[2] = 1.0;   /* should not take item 2, since weight is too large */
   profits[3] = -1.0;  /* should not take item 3, since weight is too large and profit is negative */
   profits[4] = 1.0;   /* should take item */

   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 3 );
   cr_assert( solitems[0] == 0 );
   cr_assert( solitems[1] == 4 );
   cr_assert( nonsolitems[0] == 1 );
   cr_assert( nonsolitems[1] == 2 );
   cr_assert( nonsolitems[2] == 3 );
   cr_assert_float_eq(solval, 2.0, EPS);
}

Test(solveknapsackexactly, test3, .description="test whether the case of all equal weights is handled correctly")
{
   nitems = 4;
   capacity = 4LL;
   weights[0] = 2LL;
   weights[1] = 2LL;
   weights[2] = 2LL;
   weights[3] = 2LL;

   profits[0] = 1.0;
   profits[1] = 2.0;
   profits[2] = 3.0;
   profits[3] = 4.0;

   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 2 );
   cr_assert( (solitems[0] == 2 && solitems[1] == 3) || (solitems[0] == 3 && solitems[1] == 2) );
   cr_assert( (nonsolitems[0] == 0 && nonsolitems[1] == 1) || (nonsolitems[0] == 1 && nonsolitems[1] == 0) );
   cr_assert_float_eq(solval, 7.0, EPS);
}

Test(solveknapsackexactly, test4, .description="test whether the case that only one item fits into the knapsack is handled correctly")
{
   nitems = 3;
   capacity = 3LL;
   weights[0] = 2LL;
   weights[1] = 3LL;
   weights[2] = 2LL;

   profits[0] = 3.0;
   profits[1] = 2.0;
   profits[2] = 1.0;

   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 1 );
   cr_assert( nnonsolitems == 2 );
   cr_assert( solitems[0] == 0 );
   cr_assert( (nonsolitems[0] == 1 && nonsolitems[1] == 2) || (nonsolitems[0] == 2 && nonsolitems[1] == 1) );
   cr_assert_float_eq(solval, 3.0, EPS);
}

Test(solveknapsackexactly, test_greedy1, .description="test greedy algorithm")
{
   nitems = 3;
   capacity = 3LL;
   weights[0] = 1LL;
   weights[1] = 2LL;
   weights[2] = 1LL;

   profits[0] = 3.0;
   profits[1] = 2.0;
   profits[2] = 0.5;

   /* Due to the weights and profits, the ordering is just 0,1,2. The greedy solution will pick {0,1}. This is also
    * optimal, since we used all the capacity and this solution is an integral optimal solution. So the greedy algorithm
    * should allow to terminate the process. Note, however, that this cannot be seen from the outside, so, e.g., a
    * debugger has to be used. */
   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 1 );
   cr_assert( (solitems[0] == 0 && solitems[1] == 1) || (solitems[0] == 1 && solitems[1] == 0) );
   cr_assert( nonsolitems[0] == 2 );
   cr_assert_float_eq(solval, 5, EPS);
}

Test(solveknapsackexactly, test_greedy2, .description="test whether greedy solution is equal to the rounded LP value")
{
   nitems = 3;
   capacity = 4LL;
   weights[0] = 1LL;
   weights[1] = 2LL;
   weights[2] = 2LL;

   profits[0] = 3.0;
   profits[1] = 1.0;
   profits[2] = 1.0;

   solveKnapsack();

   /* the solution packs two items and is optimal, since the LP optimal value is 4.5 */
   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 1 );
   cr_assert_float_eq(solval, 4.0, EPS);
}

Test(solveknapsackexactly, test_general, .description="general test")
{
   nitems = 6;
   capacity = 13LL;
   weights[0] = 7LL;
   weights[1] = 2LL;
   weights[2] = 7LL;
   weights[3] = 5LL;
   weights[4] = 1LL;
   weights[5] = 3LL;

   profits[0] = 1.0;
   profits[1] = 1.0;
   profits[2] = 1.0;
   profits[3] = 1.0;
   profits[4] = 1.0;
   profits[5] = 1.0;

   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 4 );
   cr_assert( nnonsolitems == 2 );
   cr_assert( (nonsolitems[0] == 0 && nonsolitems[1] == 2) || (nonsolitems[0] == 2 && nonsolitems[1] == 0) );
   cr_assert_float_eq(solval, 4.0, EPS);
}

/* large test */
Test(solveknapsackexactly, test_large, .description="large test")
{
   int j;

   nitems = 1000;
   capacity = nitems + 1;

   /* half of the items have a weight of 2 and a smaller profit */
   for (j = 0; j < nitems/2; ++j)
   {
      weights[j] = 2LL;
      profits[j] = j;
   }

   /* the remaining half of the items have a weight of 1 but a larger profit; all of them should be selected */
   for (j = nitems/2; j < nitems; ++j)
   {
      weights[j] = 1LL;
      profits[j] = j;
   }

   solveKnapsack();

   cr_assert( success );
   cr_assert( nsolitems == 750 ); /* the 500 latter items with higher profit, and 250 of the less profitable items */
   cr_assert( nnonsolitems == nitems - 750 );
   cr_assert( checkSetContainment(&items[nitems/2], solitems, nitems / 2, nsolitems) );
   cr_assert( checkSetContainment(&items[250], solitems, 250, nsolitems) );
   cr_assert( checkSetContainment(items, nonsolitems, 250, nnonsolitems) );
}
