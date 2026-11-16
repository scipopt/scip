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

/**@file   binarytree.c
 * @brief  unittest for the binary tree datastructure in misc.c
 * @author Merlin Viernickel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/pub_misc.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_BT* binarytree;
static int mydata = 4;

static
void setup(void)
{
   /* create scip */
   SCIP_CALL( SCIPcreate(&scip) );

   /* create binary tree */
   SCIP_CALL( SCIPbtCreate(&binarytree, SCIPblkmem(scip)) );
}


static
void teardown(void)
{
   /* free activity */
   SCIPbtFree(&binarytree);

   /* free scip */
   SCIP_CALL( SCIPfree(&scip) );
}

void setUp(void)
{
   setup();
}

void tearDown(void)
{
   teardown();
}

/** @brief test that setup and teardown work correctly */
void test_binarytree_setup_and_teardown(void)
{
}

/** @brief test that the binary tree checks emptiness correctly. */
void test_binarytree_empty(void)
{
   TEST_ASSERT(SCIPbtIsEmpty(binarytree));
}

/** @brief test that the binary tree adds nodes correctly. */
void test_binarytree_full(void)
{
   SCIP_BTNODE* root;
   SCIP_BTNODE* lchild;
   SCIP_BTNODE* rchild;

   /* create nodes */
   SCIP_CALL( SCIPbtnodeCreate(binarytree, &root, NULL) );
   SCIP_CALL( SCIPbtnodeCreate(binarytree, &lchild, NULL) );
   SCIP_CALL( SCIPbtnodeCreate(binarytree, &rchild, NULL) );

   /* set root */
   SCIPbtSetRoot(binarytree, root);

   /* set children */
   SCIPbtnodeSetLeftchild(root, lchild);
   SCIPbtnodeSetParent(lchild, root);
   SCIPbtnodeSetRightchild(root, rchild);
   SCIPbtnodeSetParent(rchild, root);

   /* check tree structure */
   TEST_ASSERT(SCIPbtnodeIsRoot(root));
   TEST_ASSERT(SCIPbtnodeIsLeftchild(lchild));
   TEST_ASSERT(SCIPbtnodeIsRightchild(rchild));
   TEST_ASSERT(SCIPbtnodeIsLeaf(lchild));
   TEST_ASSERT(SCIPbtnodeIsLeaf(rchild));
   TEST_ASSERT_EQUAL(root, SCIPbtGetRoot(binarytree));
   TEST_ASSERT_EQUAL(rchild, SCIPbtnodeGetSibling(lchild));
   TEST_ASSERT_EQUAL(lchild, SCIPbtnodeGetSibling(rchild));
   TEST_ASSERT_EQUAL(root, SCIPbtnodeGetParent(lchild));
   TEST_ASSERT_EQUAL(root, SCIPbtnodeGetParent(rchild));
   TEST_ASSERT_EQUAL(lchild, SCIPbtnodeGetLeftchild(root));
   TEST_ASSERT_EQUAL(rchild, SCIPbtnodeGetRightchild(root));
}

/** @brief test that the binary tree stores entry data correctly. */
void test_binarytree_data(void)
{
   SCIP_BTNODE* root;
   int* ptr;

   /* create node */
   SCIP_CALL( SCIPbtnodeCreate(binarytree, &root, NULL) );

   SCIPbtnodeSetData(root, (void*) &mydata);

   ptr = (int*) SCIPbtnodeGetData(root);
   TEST_ASSERT_EQUAL(mydata, *ptr);

   SCIPbtnodeFree(binarytree, &root);
}

int main(void)
{
   UNITY_BEGIN();
   RUN_TEST(test_binarytree_setup_and_teardown);
   RUN_TEST(test_binarytree_empty);
   RUN_TEST(test_binarytree_full);
   RUN_TEST(test_binarytree_data);
   return UNITY_END();
}
