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

/**@file   solveknapsackapprox.c
 * @brief  unit tests for approximate knapsack algorithms
 * @author Gregor Hendel
 */

#include "scip/cons_knapsack.h"
#include "scip/scip.h"

#include "include/scip_test.h"

/* helper methods */
/** GLOBAL VARIABLES **/

static int* items;
static SCIP_Longint* weights;
static SCIP_Real* profits;
static int* solitems;
static int* nonsolitems;
static int nsolitems;
static int nnonsolitems;
static SCIP_Real profit;
static int nitems;
static SCIP* scip;
static SCIP_Longint capacity;
static SCIP_RANDNUMGEN* randnumgen;
#define INITIALSEED 83

/** randomly initialize weights, profits, and set capacity to be the sum of the weights */
static
void randomDataInit(
   SCIP_RANDNUMGEN*     randgen,            /**< random number generator */
   SCIP_Longint         minweight           /**< minimum weight */
   )
{
   int i;
   capacity = 0;
   nsolitems = nnonsolitems = 0;
   /* initialize random weights and profits */
   for( i = 0; i < nitems; ++i )
   {
      weights[i] = (SCIP_Longint)SCIPrandomGetInt(randgen, minweight, 100);
      profits[i] = (SCIP_Real)(SCIPrandomGetInt(randgen, 1, 30));

      capacity += weights[i];
   }
}

#define MEMSIZE 1000000

/* TEST SUITE */
static
void setup(void)
{
   int i;
   items = NULL;
   weights = NULL;
   profits = NULL;
   solitems = NULL;
   nonsolitems = NULL;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &items, MEMSIZE) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &weights, MEMSIZE) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &profits, MEMSIZE) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &solitems, MEMSIZE) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nonsolitems, MEMSIZE) );
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, INITIALSEED, TRUE) );
   cr_assert_not_null(items);
   cr_assert_not_null(weights);
   cr_assert_not_null(profits);
   cr_assert_not_null(solitems);
   cr_assert_not_null(nonsolitems);

   /* initialize items to be numbered through */
   for( i = 0; i < MEMSIZE; ++i )
      items[i] = i;

   nsolitems = nnonsolitems = 0;
   capacity = 0;
   profit = 0.0;
   nitems = 0;
}

static
void testassignment(void)
{
   int i;

   SCIP_Real solitemsprofit;
   SCIP_Longint solitemsweight;

   solitemsprofit = 0;
   solitemsweight = 0;

   cr_assert_eq(nsolitems + nnonsolitems, nitems, "Mismatch in the number of solution items: %d + %d ~= %d\n", nsolitems, nnonsolitems, nitems);

   /* Dantzig's algorithm currently sorts items and other data by their cost/weight ratio, so that we have to reestablish the original sorting */
   SCIPsortIntRealLong(items, profits, weights, nitems);

   /* loop over solitems to compute weight and profit */
   for( i = 0; i < nsolitems; ++i )
   {
      solitemsweight += weights[solitems[i]];
      solitemsprofit += profits[solitems[i]];
   }

   cr_assert_leq(solitemsweight, capacity, "Capacity exceeded: %lld > %lld\n", solitemsweight, capacity);
   cr_assert_float_eq(profit, solitemsprofit, 1e-4, "Profit is different from recomputed profit: %.1f ~= %.1f\n", profit, solitemsprofit);


/*    loop over non-solution items to verify that none of them fully fits into the knapsack anymore
   for( i = 0; i < nnonsolitems; ++i )
   {
      cr_assert_gt(solitemsweight + weights[nonsolitems[i]], capacity, "item %d fits into the knapsack: %lld + %lld <= %lld\n", nonsolitems[i], solitemsweight, weights[nonsolitems[i]], capacity);
   }
*/
}

static
void callAndTestSolveKnapsackApproximately(void)
{
   SCIP_CALL( SCIPsolveKnapsackApproximately(scip, nitems, weights, profits, capacity, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, &profit));

   testassignment();

}

static
void teardown(void)
{

   SCIPfreeMemoryArray(scip, &items);
   SCIPfreeMemoryArray(scip, &weights);
   SCIPfreeMemoryArray(scip, &profits);
   SCIPfreeMemoryArray(scip, &solitems);
   SCIPfreeMemoryArray(scip, &nonsolitems);
   SCIPfreeRandom(scip, &randnumgen);

   SCIPfree(&scip);

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!!");
}

TestSuite(knapprox, .init = setup, .fini = teardown);

/* TESTS  */
Test(knapprox, create_and_free)
{
   /* calls setup and teardown */
}

Test(knapprox, single_call, .description = "tests single call of approximation algorithm")
{
   nitems = 30; /* more than 25 items to call Balas Zemel algorithm */

   randomDataInit(randnumgen, 1);
   /* half the items should fit into the knapsack */
   capacity /= 2;

   callAndTestSolveKnapsackApproximately();
}

Test(knapprox, trivialrandom, .description = "tests random data where all items fit")
{
   nitems = 30; /* more than 25 items to call Balas Zemel algorithm */

   randomDataInit(randnumgen, 1);

   callAndTestSolveKnapsackApproximately();
   cr_assert_eq(nsolitems, nitems, "Not all items were selected from trivial random data\n");

}

Test(knapprox, noitemfits, .description = "tests corner case where no single item fits")
{
   int i;
   nitems = 30; /* more than 25 items to call Balas Zemel algorithm */

   /* don't use a minimum weight of 1 because we might end up with a zero capacity */
   randomDataInit(randnumgen, 5);

   /* find the minimum weight as capacity */
   for( i = 0; i < nitems; ++i )
   {
      if( capacity > weights[i] )
         capacity = weights[i];
   }

   capacity -= 1;

   callAndTestSolveKnapsackApproximately();

   cr_assert_eq(nnonsolitems, nitems, "At least one item was selected despite too small capacity\n");
}

Test(knapprox, biginstance, .description = "tests big random data")
{
   nitems = 500; /* more than 25 items to call Balas Zemel algorithm */

   /* don't use a minimum weight of 1 because we might end up with a zero capacity */
   randomDataInit(randnumgen, 5);

   capacity /= 100;

   callAndTestSolveKnapsackApproximately();
}

Test(knapprox, equalratios, .description = "tests small instance with all ratios being equal")
{
   int i;
   nitems = 100; /* more than 25 items to call Balas Zemel algorithm */

   /* initialize decreasing profits and weights, such that all ratios profits/weight are equal to 1 */
   for( i = 0; i < nitems; ++i )
   {
      profits[i] = nitems - i;
      weights[i] = nitems - i;
   }

   capacity = 2 * nitems;

   callAndTestSolveKnapsackApproximately();
}

Test(knapprox, bigandbad, .description = "tests big instance that is already sorted (which should yield almost worst case run time)")
{
   int i;
   SCIP_Real expectedprofit;
   nitems = 999; /* more than 25 items to call Balas Zemel algorithm */

   /* initialize profits that are decreasing, but equal weights */
   for( i = 0; i < nitems; ++i )
   {
      profits[i] = (nitems - i + 1)/3.0;
      weights[i] = 2;
   }

   /* odd weight ensures that only the first 49999 items fit */
   capacity = 99;

   callAndTestSolveKnapsackApproximately();

   cr_assert_eq(nsolitems, 49, "Wrong number of solution items");

   expectedprofit = ((nitems + 1)/ 3.0) * nsolitems - (SCIP_Real)(nsolitems - 1)*nsolitems / 6.0;
   cr_assert_float_eq(profit, expectedprofit, 1e-4, "Expectedprofit %g ~= profit %g\n", expectedprofit, profit);
}

Test(knapprox, manybiginstances, .description = "tests many big instances for timing")
{
   int ntries = 20;
   int trial = 1;
   nitems = 999; /* more than 25 items to call Balas Zemel algorithm */

   do
   {
      randomDataInit(randnumgen, 10);

      capacity/= 3;

      callAndTestSolveKnapsackApproximately();

   } while( trial++ <= ntries );
}
