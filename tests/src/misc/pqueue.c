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

/**@file   pqueue.c
 * @brief  unit tests for priority queue
 * @author Gregor Hendel
 */

#include "scip/pub_misc.h"
#include "scip/scip.h"

#include "include/scip_test.h"

#define EPS 1e-04

struct TestContainer
{
   int                   val;
   int                   pos;
};

int values1[] = {
         5,
         4,
         6,
         7,
         10
};

int n1 = (int)sizeof(values1) / (int)sizeof(int);

typedef struct TestContainer TC;

static
SCIP_DECL_SORTPTRCOMP(cmpTestContainer)
{
   TC* t1 = (TC*)elem1;
   TC* t2 = (TC*)elem2;

   return t1->val - t2->val;
}

static
SCIP_DECL_PQUEUEELEMCHGPOS(elemChgPosTestContainer)
{
   TC* t = (TC*)elem;
   t->pos = newpos;
}

SCIP_PQUEUE* pqueue;
SCIP* scip;

/** initialize values of test containers */
static
void initTestContainers(
   TC*                   tcs,
   int*                  vals,
   int                   nvals
   )
{
   int i;
   /* initialize values and positions */
   for( i = 0; i < nvals; ++i )
   {
      tcs[i].val = vals[i];
      tcs[i].pos = -1;
   }
}

/** find test container with minimum value, return position */
static
int tcargmin(
   TC*                   tcs,
   int                   ntcs
   )
{
   int i;
   int minidx = 0;

   /* find minimum and store its index */
   for( i = 1; i < ntcs; ++i )
   {
      if( tcs[i].val < tcs[minidx].val )
         minidx = i;
   }

   return minidx;
}

/** insert elements into pqueue and assert that the length is correct */
static
SCIP_RETCODE insertIntoPQueue(
   SCIP_PQUEUE*          pqueue_local,
   TC*                   tcs,
   int                   ntcs
   )
{
   int i;

   /* insert each element into the queue */
   for( i = 0; i < ntcs; ++i )
   {
      int minidx;
      SCIP_CALL( SCIPpqueueInsert(pqueue_local, (void *)(&tcs[i])) );

      cr_assert_eq(SCIPpqueueNElems(pqueue_local), i + 1);
      minidx = tcargmin(tcs, i + 1);
      cr_assert_eq(tcs[minidx].val, ((TC*)SCIPpqueueFirst(pqueue_local))->val);
   }

   return SCIP_OKAY;
}

/** delete elements from pqueue by position */
static
void deleteFromPQueue(
   SCIP_PQUEUE*          pqueue_local,
   TC*                   tcs,
   int                   ntcs
   )
{
   int i;
   /* remove elements in the same order in which they were inserted */
   for( i = 0; i < ntcs; ++i )
   {
      int elempos = tcs[i].pos;
      int minidx = tcargmin(&tcs[i], ntcs - i) + i;

      /* ensure that the minimum element */
      cr_assert_eq(tcs[minidx].val, ((TC*)SCIPpqueueFirst(pqueue_local))->val);
      cr_assert(elempos < SCIPpqueueNElems(pqueue_local));
      cr_assert_eq(SCIPpqueueElems(pqueue_local)[elempos], (void*)(&tcs[i]));
      SCIPpqueueDelPos(pqueue_local, elempos);
   }
}


/* TEST SUITE */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( SCIPpqueueCreate(&pqueue, n1, 1.1, cmpTestContainer, elemChgPosTestContainer) );
}

static
void teardown(void)
{
   SCIPpqueueFree(&pqueue);

   SCIP_CALL( SCIPfree(&scip) );
}



TestSuite(pqueue, .init = setup, .fini = teardown);

/* TESTS  */
Test(pqueue, create_and_free)
{
   /* calls setup and teardown */
}


Test(pqueue, insert, .description="test insertion")
{
   int i;
   TC testc[n1];
   TC* testc2[n1];

   /* initialize values to those in values1 */
   initTestContainers(testc, values1, n1);

   /* insert test elements into pqueue */
   SCIP_CALL( insertIntoPQueue(pqueue, testc, n1) );

   /* remove elements from priority queue */
   i = 0;
   while( SCIPpqueueNElems(pqueue) )
      testc2[i++] = (TC*)SCIPpqueueRemove(pqueue);

   /* test if elements are sorted in increasing order */
   for( i = 0; i < n1 - 1; ++i )
      cr_assert(testc2[i]->val <= testc2[i+1]->val);
}

Test(pqueue, delpos, .description="test deletion using positions")
{
   TC testc[n1];

   /* initialize values to those in values1 */
   initTestContainers(testc, values1, n1);

   /* insert test elements into pqueue */
   SCIP_CALL( insertIntoPQueue(pqueue, testc, n1) );

   /* remove elements in the same order in which they were inserted */
   deleteFromPQueue(pqueue, testc, n1);
}

#define N2 50
#define NSHUF 1000
Test(pqueue, insert_and_delete_random, .description="test random value permutations to cover extreme cases")
{
   int i;
   int s;
   TC testc[N2];
   int values2[N2];

   SCIP_RANDNUMGEN* rng;

   /* initialize values array */
   for( i = 0; i < N2; ++i )
      values2[i] = i + 1;

   /* cover cases by shuffling the values array a couple of times */
   SCIP_CALL( SCIPcreateRandom(scip, &rng, 42, TRUE) );
   for( s = 0; s < NSHUF; ++s )
   {
      SCIPrandomPermuteIntArray(rng, values2, 0, N2);

      /* initialize values to those in values2 */
      initTestContainers(testc, values2, N2);

      /* insert elements into pqueue */
      SCIP_CALL( insertIntoPQueue(pqueue, testc, N2) );

      /* remove elements in the same order in which they were inserted */
      deleteFromPQueue(pqueue, testc, N2);
   }

   SCIPfreeRandom(scip, &rng);
}
