/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: buffer.c,v 1.12 2004/02/04 17:27:17 bzfpfend Exp $"

/**@file   buffer.c
 * @brief  methods for memory buffers for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "message.h"
#include "memory.h"
#include "set.h"
#include "buffer.h"

#include "struct_buffer.h"



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
RETCODE SCIPbufferAllocMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   )
{
   int bufnum;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
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

   debugMessage("allocated buffer %d/%d at %p of size %d (required size: %d) for pointer %p\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], size, ptr);

   return SCIP_OKAY;
}

/** allocates the next unused buffer and copies the given memory into the buffer */
RETCODE SCIPbufferDuplicateMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   void*            source,             /**< memory block to copy into the buffer */
   int              size                /**< minimal required size of the buffer */
   )
{
   assert(source != NULL);

   /* allocate a buffer of the given size */
   CHECK_OKAY( SCIPbufferAllocMem(buffer, set, ptr, size) );

   /* copy the source memory into the buffer */
   copyMemorySize(*ptr, source, size);

   return SCIP_OKAY;
}

/** reallocates the buffer to at least the given size */
RETCODE SCIPbufferReallocMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   )
{
   int bufnum;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
   assert(buffer->firstfree >= 1);
   assert(ptr != NULL);
   assert(size >= 0);

   /* Search the pointer in the buffer list
    * Usally, buffers are allocated and freed like a stack, such that the currently used pointer is
    * most likely at the end of the buffer list.
    */
   for( bufnum = buffer->firstfree-1; bufnum >= 0 && buffer->data[bufnum] != *ptr; --bufnum )
   {
   }
   assert(buffer->data[bufnum] == *ptr);
   assert(buffer->used[bufnum]);
   assert(buffer->size[bufnum] >= 1);

   /* check if the buffer has to be enlarged */
   if( size > buffer->size[bufnum] )
   {
      int newsize;

      /* enlarge buffer */
      newsize = SCIPsetCalcMemGrowSize(set, size);
      ALLOC_OKAY( reallocMemorySize(&buffer->data[bufnum], newsize) );
      buffer->size[bufnum] = newsize;
      *ptr = buffer->data[bufnum];
   }
   assert(buffer->size[bufnum] >= size);
   assert(*ptr == buffer->data[bufnum]);

   debugMessage("reallocated buffer %d/%d at %p to size %d (required size: %d) for pointer %p\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], size, ptr);

   return SCIP_OKAY;
}

/** frees a buffer */
void SCIPbufferFreeMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              dummysize           /**< used to get a safer define for SCIPsetFreeBufferSize/Array */
   )
{  /*lint --e{715}*/
   int bufnum;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
   assert(buffer->firstfree >= 1);

   /* Search the pointer in the buffer list
    * Usally, buffers are allocated and freed like a stack, such that the freed pointer is
    * most likely at the end of the buffer list.
    */
   for( bufnum = buffer->firstfree-1; bufnum >= 0 && buffer->data[bufnum] != *ptr; --bufnum )
   {
   }
   assert(buffer->data[bufnum] == *ptr);
   assert(buffer->used[bufnum]);

   *ptr = NULL;
   buffer->used[bufnum] = FALSE;

   while( buffer->firstfree > 0 && !buffer->used[buffer->firstfree-1] )
      buffer->firstfree--;

   debugMessage("freed buffer %d/%d at %p of size %d for pointer %p, first free is %d\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], ptr, buffer->firstfree);
}
