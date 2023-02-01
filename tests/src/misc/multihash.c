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

/**@file   multihashlist.c
 * @brief  unittest for the multihash datastructure in misc.c
 * @author Merlin Viernickel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/pub_misc.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_MULTIHASH* multihash;

const int arraylen = 6;
static int myentries[] =
{
   14, 5, 92, 31, 12, 91
};
static int* myptrs[] =
{
   &myentries[0],
   &myentries[1],
   &myentries[2],
   &myentries[3],
   &myentries[4],
   &myentries[5]
};

/** get key of hash element */
static
SCIP_DECL_HASHGETKEY(getKey)
{
   return elem;
}  /*lint !e715*/


/** checks if two arrays have the same elements */
static
SCIP_DECL_HASHKEYEQ(keyEQ)
{
   int* ptr1;
   int* ptr2;

   ptr1 = (int*) key1;
   ptr2 = (int*) key2;

   return (*ptr1 == *ptr2);
}  /*lint !e715*/


/** get first entry as hash value */
static
SCIP_DECL_HASHKEYVAL(keyVal)
{
   int* ptr;

   ptr = (int*) key;

   return *ptr;
}  /*lint !e715*/

static
void setup(void)
{
   /* create scip */
   SCIP_CALL( SCIPcreate(&scip) );

   /* create multihash table */
   SCIP_CALL( SCIPmultihashCreate(&multihash, SCIPblkmem(scip), SCIPcalcMultihashSize(arraylen), getKey, keyEQ, keyVal, (void*) scip) );
}


static
void teardown(void)
{
   /* free multihash table */
   SCIPmultihashFree(&multihash);

   /* free scip */
   SCIP_CALL( SCIPfree(&scip) );
}


TestSuite(multihash, .init = setup, .fini = teardown);

Test(multihash, setup_and_teardown, .description = "test that setup and teardown work correctly")
{
}

Test(multihash, test_multihash_insertion, .description = "test that the multi hash map stores entries correctly.")
{
   int* ptr;
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPmultihashInsert(multihash, (void*) myptrs[i]) );

   cr_assert_eq(arraylen, SCIPmultihashGetNElements(multihash));
   for( i = 0; i < arraylen; i++ )
   {
      ptr = SCIPmultihashRetrieve(multihash, (void*) myptrs[i]);
      cr_assert_eq(myentries[i], *ptr);
   }
}

Test(multihash, test_multihash_remove, .description = "test that the multi hash map removes entries correctly.")
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPmultihashInsert(multihash, (void*) myptrs[i]) );

   cr_assert_eq(arraylen, SCIPmultihashGetNElements(multihash));
   for( i = 0; i < arraylen; i++ )
   {
      cr_assert(SCIPmultihashExists(multihash, (void*) myptrs[i]));
      SCIP_CALL( SCIPmultihashRemove(multihash, (void*) myptrs[i]) );
      cr_assert(! SCIPmultihashExists(multihash, (void*) myptrs[i]));
   }
}

Test(multihash, test_multihash_removeall, .description = "test that the multi hash map removes all entries at once correctly.")
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPmultihashInsert(multihash, (void*) myptrs[i]) );

   cr_assert_eq(arraylen, SCIPmultihashGetNElements(multihash));

   SCIPmultihashRemoveAll(multihash);

   for( i = 0; i < arraylen; i++ )
      cr_assert(! SCIPmultihashExists(multihash, (void*) myptrs[i]));
}

Test(multihash, test_multihash_statistics, .description = "test that the multi hash map prints statistics correctly")
{
   SCIP_MESSAGEHDLR* msghdlr;
   int i;

   msghdlr = SCIPgetMessagehdlr(scip);

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPmultihashInsert(multihash, (void*) myptrs[i]) );

   SCIPmultihashPrintStatistics(multihash, msghdlr);
}
