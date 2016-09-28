/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   select.c
 * @brief  unit tests for selection of unweighted and weighted median
 * @author Gselector Hendel
 */

#include "scip/misc.c"
#include "scip/scip.h"

#include "include/scip_test.h"

/** GLOBAL VARIABLES **/

/* TEST SUITE */
static
void setup(void)
{
}

static
void teardown(void)
{
}

TestSuite(select, .init = setup, .fini = teardown);

/* TESTS  */
Test(select, create_and_free)
{
   /* calls setup and teardown */
}

#define ARRAYMEMSIZE 30
Test(select, random_permutation, .description = "tests selection on a bunch of random permutations of the integers 1...n")
{
   int len = ARRAYMEMSIZE;
   int key[ARRAYMEMSIZE];
   unsigned int randomseed = 42;
   int i;

   /* initialize key */
   for( i = 0; i < len; ++i )
      key[i] = i;

   /* loop over all positions of the array and check whether the correct element is selected after a random permutation */
   for( i = 0; i < len; ++i )
   {
      int inputkey[ARRAYMEMSIZE];
      SCIPpermuteIntArray(key, 0, len, &randomseed);

      /* save input permutation for debugging */
      BMScopyMemoryArray(inputkey, key, len);
      SCIPselectInt(key, i, len);

      cr_assert_eq(key[i], i, "Wrong key selected: %d ~= %d", key[i], i);
   }
}
