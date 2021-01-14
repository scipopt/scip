/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bmschunkmem.c
 * @brief  unit test for garbage collection of BMS chunk memory
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "include/scip_test.h"

/* struct for chunk elements (size must be divisible by size of pointer, i.e., 8 on 64 Bit machine) */
typedef struct _TESTSTRUCT
{
   int idx;
   int rnd;   /* unused - just to make size divisible by pointer size */
} TESTSTRUCT;

#define CHUNK_SIZE 8     /* initial chunk size (will be adapted by memory.c) */
#define CHUNK_NUM  384   /* this size creates two chunk blocks, one of 128 and one of 256 bytes */


/* global variables */
static BMS_CHKMEM* mem;

/* test suites */

/** setup of test run */
static
void setup(void)
{
   /* initialize memory allocator - no automatic garbage collection (called manually below) */
   mem = BMScreateChunkMemory(sizeof(TESTSTRUCT), CHUNK_SIZE, -1);
}

/** deinitialization method */
static
void teardown(void)
{
   /* clear all memory chunks */
   BMSclearChunkMemory(mem);

   cr_assert_eq(BMSgetChunkMemoryUsed(mem), 0, "There is a memory leak!");

   /* delete memory allocator */
   BMSdestroyChunkMemory(&mem);

   cr_assert_null(mem);
}

TestSuite(bmschunkmem, .init = setup, .fini = teardown);

/* TESTS */

/** Test that allocates memory in order to create two chunks (of size 128 and 256 bytes). Then it frees the last 257
 *  bytes, which should free the last chunk block and make the first an eager chunk.
 *
 *  This test is motivated by a false positive message from coverity:
 *    Calling "BMSfreeChunkMemory_call" frees pointer "mem->firsteager" which has already been freed.
 *  The message is a false positive, since chkmem->firsteager is corrected in unlinkEagerChunk().
 */
Test(bmschunkmem, bmsgarbagecollect)
{
   TESTSTRUCT* data[CHUNK_NUM];
   int i;

   /* allocate chunk elements */
   for (i = 0; i < CHUNK_NUM; ++i)
   {
      if ( BMSallocChunkMemory(mem, &data[i]) == NULL )
      {
         SCIPerrorMessage("No memory available.\n");
         abort();
      }
      data[i]->idx = i;
   }

   /* free last 257 of the chunks */
   for (i = CHUNK_NUM-1; i >= 127; --i)
   {
      assert( data[i] != NULL );
      BMSfreeChunkMemoryNull(mem, &data[i]);
      assert( data[i] == NULL );
   }

   /* call garbage collection */
   BMSgarbagecollectChunkMemory(mem);

   /* free remaining chunks */
   for (i = 0; i < CHUNK_NUM; ++i)
      BMSfreeChunkMemoryNull(mem, &data[i]);
}
