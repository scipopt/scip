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

/**@file   sort.c
 * @brief  unit tests for some sort implementations
 * @author Christoph Schubert
 */

#include<stdio.h>

#include "scip/pub_misc.h"
#include "scip/scip.h"

#include "include/scip_test.h"

/* helper methods */

/** compares two entries of an int array, returns 0 if both are equal, a negative number if the entry at position ind2
 *  is larger than the one at position ind1, and a positive number otherwise; used for sorting an array in ascending
 *  order
 */
static
SCIP_DECL_SORTINDCOMP(intComparatorAscending)
{  /*lint --e{715}*/
   int* tosort = (int*)dataptr;

   assert(tosort != NULL);

   return tosort[ind1] - tosort[ind2];
}

/** compares two entries of an int array, returns 0 if both are equal, a positive number if the entry at position ind2
 *  is larger than the one at position ind1, and a negative number otherwise; used for sorting an array in descending
 *  order
 */
static
SCIP_DECL_SORTINDCOMP(intComparatorDescending)
{  /*lint --e{715}*/
   int* tosort = (int*)dataptr;

   assert(tosort != NULL);

   return tosort[ind2] - tosort[ind1];
}

/** GLOBAL VARIABLES **/
static SCIP* scip;
static int* tosort;
static int ntosort;

/* TEST SUITE */
static
void setup(void)
{
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   ntosort = 8;
   SCIP_CALL( SCIPallocBufferArray(scip, &tosort, ntosort) );
   tosort[0] = 2;
   tosort[1] = 4;
   tosort[2] = 3;
   tosort[3] = 7;
   tosort[4] = 1;
   tosort[5] = 5;
   tosort[6] = 6;
   tosort[7] = 0;

   cr_assert_not_null(scip);
   cr_assert_not_null(tosort);
}

static
void teardown(void)
{

   SCIPfreeBufferArray(scip, &tosort);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(sort, .init = setup, .fini = teardown);

Test(sort, create_and_free)
{
   /* calls setup and teardown */
}

Test(sort, asc_sort_up, .description = "tests SCIPsort by checking that the permutation is sorted ascending")
{
   int* perm;
   int i;

   cr_assert_not_null(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsort(perm, intComparatorAscending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      cr_assert_lt(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

Test(sort, asc_sort_down, .description = "tests SCIPsortDown by checking that the permutation is sorted descending")
{
   int* perm;
   int i;

   cr_assert_not_null(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsortDown(perm, intComparatorAscending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      cr_assert_gt(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

Test(sort, desc_sort_up, .description = "tests SCIPsort by checking that the permutation is sorted descending")
{
   int* perm;
   int i;

   cr_assert_not_null(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsort(perm, intComparatorDescending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      cr_assert_gt(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}

Test(sort, desc_sort_down, .description = "tests SCIPsortDown by checking that the permutation is sorted ascending")
{
   int* perm;
   int i;

   cr_assert_not_null(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &perm, ntosort) );

   SCIPsortDown(perm, intComparatorDescending, (void*)tosort, ntosort);

   for( i = 0; i < ntosort-1; i++ )
   {
      cr_assert_lt(tosort[perm[i]], tosort[perm[i+1]]);
   }

   SCIPfreeBufferArray(scip, &perm);
}
