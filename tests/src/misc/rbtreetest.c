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
 * @brief  unit tests for red black tree datastructure
 * @author Robert Lion Gottwald
 */

#include "scip/misc.h"
#include "scip/scip.h"
#include "scip/rbtree.h"

#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
static SCIP_RANDNUMGEN* randgen;
static BMS_BLKMEM* blkmem;
static unsigned int randomseed = 42;

typedef struct SomeType
{
   SCIP_RBTREE_HOOKS;
   int key;
   SCIP_Real whateverdata;
} SOME_TYPE;

#define SOMETYPE_LT(a,b) ( a < b->key )
#define SOMETYPE_GT(a,b) ( a > b->key )

static
SCIP_DEF_RBTREE_FIND(findSomeType, int, SOME_TYPE, SOMETYPE_LT, SOMETYPE_GT)


/* TEST SUITE */
static
void setup(void)
{
   blkmem = BMScreateBlockMemory(1, 10);
   SCIP_CALL_ABORT( SCIPrandomCreate(&randgen, blkmem, randomseed) );
}

static
void teardown(void)
{
   SCIPrandomFree(&randgen, blkmem);
   BMSdestroyBlockMemory(&blkmem);
}

TestSuite(select, .init = setup, .fini = teardown);

/* TESTS  */
Test(rbtree, create_and_free)
{
   /* calls setup and teardown */
}


#define ARRAYMEMSIZE 700
Test(rbtree, rb_random_insert, .description = "tests rb tree insertion and lookup of the integers 1...n in random order",  .init = setup, .fini = teardown)
{
   int len = ARRAYMEMSIZE;
   int key[ARRAYMEMSIZE];
   int i;
   int j;
   int pos;
   SOME_TYPE* root;
   SOME_TYPE* x;
   /* initialize key */
   for( j = 0; j < len; ++j )
      key[j] = j;

   root = NULL;

   SCIPrandomPermuteIntArray(randgen, key, 0, len);
   /* allocate and insert the elements with integer keys in random order */
   for( i = 0; i < len; ++i )
   {
      SOME_TYPE* elem;
      SOME_TYPE* parent;

      BMSallocBlockMemory(blkmem, &elem);

      elem->key = key[i];
      elem->whateverdata = sqrt(key[i]);
      pos = findSomeType(root, elem->key, &parent);

      cr_assert(pos != 0);
      SCIPrbtreeInsert(&root, parent, pos, elem);
   }

   /* iterate all elements and check order */
   i = 0;
   FOR_EACH_NODE(SOME_TYPE*, node, root,
   {
      cr_assert_eq(node->key, i, "expected key %i but got %i\n", i, node->key);
      cr_assert_eq(node->whateverdata, sqrt(i));
      ++i;
   })

   /* check number of elements was correct */
   cr_assert_eq(i, len);
   /* check lookup of specific elements */
   pos = findSomeType(root, 10, &x);
   cr_assert(pos == 0);
   cr_assert(x->key == 10);
   /* delete element 10 */
   SCIPrbtreeDelete(&root, x);
   BMSfreeBlockMemory(blkmem, &x);
   /* lookup and delete 100 */
   pos = findSomeType(root, 100, &x);
   cr_assert(pos == 0);
   cr_assert(x->key == 100);
   SCIPrbtreeDelete(&root, x);
   BMSfreeBlockMemory(blkmem, &x);

   /* lookup of deleted element should not find element and return predecessor/successor */
   pos = findSomeType(root, 10, &x);
   cr_assert( (x->key == 9 && pos == -1) || (x->key == 11 && pos == 1) );

   pos = findSomeType(root, 100, &x);
   cr_assert( (x->key == 99 && pos == -1) || (x->key == 101 && pos == 1) );

   /* iterate again and check order */
   i = 0;
   FOR_EACH_NODE(SOME_TYPE*, node, root,
   {
      int k = i;
      if( k >= 10 )
         ++k;
      if( k >= 100 )
         ++k;
      cr_assert_eq(node->key, k, "expected key %i but got %i\n", k, node->key);
      cr_assert_eq(node->whateverdata, sqrt(k));
      ++i;
   })
   /* check number again (with 10 and 100 missing) */
   cr_assert_eq(i, len-2);
}
