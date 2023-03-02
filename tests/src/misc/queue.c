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

/**@file   queue.c
 * @brief  unittest for the queue datastructure in misc.c
 * @author Merlin Viernickel
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/scip.h"
#include "scip/pub_misc.h"

#include "include/scip_test.h"

static SCIP* scip;
static SCIP_QUEUE* queue;

#define arraylen 3
static SCIP_Real myrealentries[] =
{
   5.0,
   23.3,
   14.5
};
static SCIP_Real* myptrentries[] =
{
   &myrealentries[0],
   &myrealentries[1],
   &myrealentries[2]
};
static SCIP_Bool myuintentries[] =
{
   TRUE,
   FALSE,
   TRUE
};


static
void setup(void)
{
   /* create scip */
   SCIP_CALL( SCIPcreate(&scip) );

   /* create queue */
   SCIP_CALL( SCIPqueueCreate(&queue, arraylen, 2) );
}

static
void teardown(void)
{
   /* free scip */
   SCIP_CALL( SCIPfree(&scip) );

   /* free queue */
   SCIPqueueFree(&queue);
}

TestSuite(queue, .init = setup, .fini = teardown);

Test(queue, setup_and_teardown, .description = "test that setup and teardown work correctly")
{
}

Test(queue, test_queue_insertion, .description = "test that the queue stores entries correctly.")
{
   int i;

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPqueueInsert(queue, (void*) myptrentries[i]) );

   for( i = 0; i < arraylen; i++ )
   {
      cr_assert_eq(myptrentries[i], (SCIP_Real*) SCIPqueueFirst(queue));
      cr_assert_eq(myptrentries[i], (SCIP_Real*) SCIPqueueRemove(queue));
   }

   cr_assert(SCIPqueueIsEmpty(queue));
}

Test(queue, test_queue_uintinsertion, .description = "test that the queue stores unsigned integer entries correctly.")
{
   int i;

   /* create queue */
   SCIP_CALL( SCIPqueueCreate(&queue, arraylen, 2) );

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPqueueInsertUInt(queue, myuintentries[i]) );

   for( i = 0; i < arraylen; i++ )
   {
      cr_assert_eq(myuintentries[i], SCIPqueueFirstUInt(queue));
      cr_assert_eq(myuintentries[i], SCIPqueueRemoveUInt(queue));
   }

   cr_assert(SCIPqueueIsEmpty(queue));
}

Test(queue, test_queue_clear, .description = "test that the queue clears entries correctly.")
{
   int i;

   /* create queue */
   SCIP_CALL( SCIPqueueCreate(&queue, arraylen, 2) );

   for( i = 0; i < arraylen; i++ )
      SCIP_CALL( SCIPqueueInsert(queue, (void*) myptrentries[i]) );

   cr_assert_eq(arraylen, SCIPqueueNElems(queue));

   SCIPqueueClear(queue);

   cr_assert(SCIPqueueIsEmpty(queue));
}
