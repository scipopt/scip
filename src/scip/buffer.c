/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   buffer.c
 * @brief  memory buffer for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "buffer.h"
#include "memory.h"



/** creates memory buffer storage */
RETCODE SCIPbufferCreate(
   BUFFER**         buffer              /**< pointer to memory buffer storage */
   )
{
   assert(buffer != NULL);

   ALLOC_OKAY( allocMemory(buffer) );
   (*buffer)->data = NULL;
   (*buffer)->size = NULL;
   (*buffer)->used = NULL;
   (*buffer)->ndata = 0;
   (*buffer)->firstfree = 0;

   return SCIP_OKAY;
}

/** frees memory buffer */
void SCIPbufferFree(
   BUFFER**         buffer              /**< pointer to memory buffer storage */
   )
{
   int i;

   assert(buffer != NULL);

   for( i = 0; i < (*buffer)->ndata; ++i )
   {
      assert(!(*buffer)->used[i]);
      freeMemoryArrayNull(&(*buffer)->data[i]);
   }
   freeMemoryArrayNull(&(*buffer)->data);
   freeMemoryArrayNull(&(*buffer)->size);
   freeMemoryArrayNull(&(*buffer)->used);
   freeMemory(buffer);
}

/** allocates the next unused buffer */
RETCODE SCIPbufferCapture(
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   )
{
   int bufnum;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
   assert(set != NULL);
   assert(ptr != NULL);
   assert(size >= 0);
   
   /* allocate minimal 1 byte */
   if( size == 0 )
      size = 1;

   /* check, if we need additional buffers */
   if( buffer->firstfree == buffer->ndata )
   {
      int newsize;
      int i;

      /* create additional buffers */
      newsize = SCIPsetCalcMemGrowSize(set, buffer->firstfree+1);
      ALLOC_OKAY( reallocMemoryArray(&buffer->data, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&buffer->size, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&buffer->used, newsize) );
      for( i = buffer->ndata; i < newsize; ++i )
      {
         buffer->data[i] = NULL;
         buffer->size[i] = 0;
         buffer->used[i] = FALSE;
      }
      buffer->ndata = newsize;
   }
   assert(buffer->firstfree < buffer->ndata);

   /* check, if the actual buffer is large enough */
   bufnum = buffer->firstfree;
   assert(!buffer->used[bufnum]);
   if( buffer->size[bufnum] < size )
   {
      int newsize;

      /* enlarge buffer */
      newsize = SCIPsetCalcMemGrowSize(set, size);
      ALLOC_OKAY( reallocMemorySize(&buffer->data[bufnum], newsize) );
      buffer->size[bufnum] = newsize;
   }
   assert(buffer->size[bufnum] >= size);

   *ptr = buffer->data[bufnum];
   buffer->used[bufnum] = TRUE;
   buffer->firstfree++;

   debugMessage("captured buffer %d/%d at %p of size %d (required size: %d) for pointer %p\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], size, ptr);

   return SCIP_OKAY;
}

/** releases a buffer */
void SCIPbufferRelease(
   BUFFER*          buffer,             /**< memory buffer storage */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              dummysize           /**< used to get a safer define for SCIPsetReleaseBufferSize/Array */
   )
{  /*lint --e{715}*/
   int i;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
   assert(buffer->firstfree >= 1);

   /* Search the pointer in the buffer list
    * Usally, buffers are allocated and freed like a stack, such that the released pointer is
    * most likely at the end of the buffer list.
    */
   for( i = buffer->firstfree-1; i >= 0 && buffer->data[i] != *ptr; --i )
   {
   }
   assert(buffer->data[i] == *ptr);
   assert(buffer->used[i]);

   *ptr = NULL;
   buffer->used[i] = FALSE;

   while( buffer->firstfree > 0 && !buffer->used[buffer->firstfree-1] )
      buffer->firstfree--;

   debugMessage("released buffer %d/%d at %p of size %d for pointer %p, first free is %d\n", 
      i, buffer->ndata, buffer->data[i], buffer->size[i], ptr, buffer->firstfree);
}
