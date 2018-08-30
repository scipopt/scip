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

/**@file   solveknapsackexactly.c
 * @brief  unit tests for exact dynamic programming algorithm for knapsack problem
 * @author Marc Pfetsch
 */

#include "scip/cons_knapsack.h"
#include "scip/scip.h"

#include "include/scip_test.h"

#define EPS 1e-06

/** GLOBAL VARIABLES **/
static SCIP* scip;

/* TEST SUITE */
static
void setup(void)
{
   SCIPcreate(&scip);
}

static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(solveknapsackexactly, .init = setup, .fini = teardown);

/* TESTS  */

/** test trivial cases */
Test(solveknapsackexactly, test1)
{
   int nitems;
   SCIP_Longint weights[5];
   SCIP_Real profits[5];
   int items[5];
   int solitems[5];
   int nonsolitems[5];
   SCIP_Longint capacity;
   SCIP_Bool success;
   SCIP_Real solval;
   int nnonsolitems;
   int nsolitems;

   items[0] = 0;
   items[1] = 1;
   items[2] = 2;
   items[3] = 3;
   items[4] = 4;

   /* check whether the case that all items are redundant is caught correctly */
   nitems = 2;
   capacity = 1LL;
   weights[0] = 2LL;
   weights[1] = 1LL;
   profits[0] = 1.0;   /* should not be taken, since capacity is too large */
   profits[1] = -1.0;  /* should not be taken, since profit is negative */

   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 0 );
   cr_assert( nnonsolitems == 2 );
   cr_assert_float_eq(solval, 0.0, EPS);

   /* test whether the correct items are sorted out (weight == 0 or negative profit) - also tests whether the case of
    * all items fitting into the knapsack is caught correctly. */
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
   profits[3] = -1.0;  /* should not take item 2, since weight is too large */
   profits[4] = 1.0;   /* should take item */

   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 3 );
   cr_assert( solitems[0] == 0 );
   cr_assert( solitems[1] == 4 );
   cr_assert( nonsolitems[0] == 1 );
   cr_assert( nonsolitems[1] == 2 );
   cr_assert( nonsolitems[2] == 3 );
   cr_assert_float_eq(solval, 2.0, EPS);

   /* test whether the case of all equal weights is handled correctly */
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

   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 2 );
   cr_assert( (solitems[0] == 2 && solitems[1] == 3) || (solitems[0] == 3 && solitems[1] == 2) );
   cr_assert( (nonsolitems[0] == 0 && nonsolitems[1] == 1) || (nonsolitems[0] == 1 && nonsolitems[1] == 0) );
   cr_assert_float_eq(solval, 7.0, EPS);

   /* test whether the case that only one item fits into the knapsack is handled correctly */
   nitems = 3;
   capacity = 3LL;
   weights[0] = 2LL;
   weights[1] = 3LL;
   weights[2] = 2LL;

   profits[0] = 3.0;
   profits[1] = 2.0;
   profits[2] = 1.0;

   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 1 );
   cr_assert( nnonsolitems == 2 );
   cr_assert( solitems[0] == 0 );
   cr_assert( (nonsolitems[0] == 1 && nonsolitems[1] == 2) || (nonsolitems[0] == 2 && nonsolitems[1] == 1) );
   cr_assert_float_eq(solval, 3.0, EPS);
}


/** test greedy algorithm */
Test(solveknapsackexactly, test2)
{
   int nitems;
   SCIP_Longint weights[3];
   SCIP_Real profits[3];
   int items[3];
   int solitems[3];
   int nonsolitems[3];
   SCIP_Longint capacity;
   SCIP_Bool success;
   SCIP_Real solval;
   int nnonsolitems;
   int nsolitems;

   nitems = 3;
   items[0] = 0;
   items[1] = 1;
   items[2] = 2;

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
   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 1 );
   cr_assert( (solitems[0] == 0 && solitems[1] == 1) || (solitems[0] == 1 && solitems[1] == 0) );
   cr_assert( nonsolitems[0] == 2 );
   cr_assert_float_eq(solval, 5, EPS);

   /* test whether greedy solution is equal to the rounded LP value */
   capacity = 4LL;
   weights[0] = 1LL;
   weights[1] = 2LL;
   weights[2] = 2LL;

   profits[0] = 3.0;
   profits[1] = 1.0;
   profits[2] = 1.0;

   /* the solution packs one item and is optimal, since the LP optimal value is 1.5 */
   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 2 );
   cr_assert( nnonsolitems == 1 );
   cr_assert_float_eq(solval, 4.0, EPS);
}


/** general test */
Test(solveknapsackexactly, test3)
{
   int nitems;
   SCIP_Longint weights[6];
   SCIP_Real profits[6];
   int items[6];
   int solitems[6];
   int nonsolitems[6];
   SCIP_Longint capacity;
   SCIP_Bool success;
   SCIP_Real solval;
   int nnonsolitems;
   int nsolitems;

   nitems = 6;
   items[0] = 0;
   items[1] = 1;
   items[2] = 2;
   items[3] = 3;
   items[4] = 4;
   items[5] = 5;

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

   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 4 );
   cr_assert( nnonsolitems == 2 );
   cr_assert( (nonsolitems[0] == 0 && nonsolitems[1] == 2) || (nonsolitems[0] == 2 && nonsolitems[1] == 0) );
   cr_assert_float_eq(solval, 4.0, EPS);
}


/** large test */
Test(solveknapsackexactly, test4)
{
   int nitems;
   SCIP_Longint* weights;
   SCIP_Real* profits;
   int* items;
   int* solitems;
   int* nonsolitems;
   SCIP_Longint capacity;
   SCIP_Bool success;
   SCIP_Real solval;
   int nnonsolitems;
   int nsolitems;
   int j;

   nitems = 1000;
   SCIP_CALL( SCIPallocBufferArray(scip, &items, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &profits, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, nitems) );

   for (j = 0; j < nitems/2; ++j)
   {
      items[j] = j;
      weights[j] = 2LL;
      profits[j] = j;
   }
   for (j = nitems/2; j < nitems; ++j)
   {
      items[j] = j;
      weights[j] = 1LL;
      profits[j] = j;
   }
   capacity = nitems + 1;

   SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &solval, &success) );
   cr_assert( success );
   cr_assert( nsolitems == 750 );
   cr_assert( nnonsolitems == nitems - 750 );

   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &profits);
   SCIPfreeBufferArray(scip, &weights);
   SCIPfreeBufferArray(scip, &items);
}
