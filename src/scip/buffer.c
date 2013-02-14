/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   buffer.c
 * @brief  methods for memory buffers for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/buffer.h"

#include "scip/struct_buffer.h"



/** creates memory buffer storage */
SCIP_RETCODE SCIPbufferCreate(
   SCIP_BUFFER**         buffer              /**< pointer to memory buffer storage */
   )
{
   assert(buffer != NULL);

   SCIP_ALLOC( BMSallocMemory(buffer) );
   (*buffer)->data = NULL;
   (*buffer)->size = NULL;
   (*buffer)->used = NULL;
   (*buffer)->ndata = 0;
   (*buffer)->firstfree = 0;

   return SCIP_OKAY;
}

/** frees memory buffer */
void SCIPbufferFree(
   SCIP_BUFFER**         buffer              /**< pointer to memory buffer storage */
   )
{
   int i;

   assert(buffer != NULL);

   for( i = 0; i < (*buffer)->ndata; ++i )
   {
      assert(!(*buffer)->used[i]);
      BMSfreeMemoryArrayNull(&(*buffer)->data[i]);
   }
   BMSfreeMemoryArrayNull(&(*buffer)->data);
   BMSfreeMemoryArrayNull(&(*buffer)->size);
   BMSfreeMemoryArrayNull(&(*buffer)->used);
   BMSfreeMemory(buffer);
}

/** allocates the next unused buffer */
SCIP_RETCODE SCIPbufferAllocMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   int                   size                /**< minimal required size of the buffer */
   )
{
#ifndef SCIP_NOBUFFERMEM
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
      SCIP_ALLOC( BMSreallocMemoryArray(&buffer->data, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&buffer->size, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&buffer->used, newsize) );
      for( i = buffer->ndata; i < newsize; ++i )
      {
         buffer->data[i] = NULL;
         buffer->size[i] = 0;
         buffer->used[i] = FALSE;
      }
      buffer->ndata = newsize;
   }
   assert(buffer->firstfree < buffer->ndata);

   /* check, if the current buffer is large enough */
   bufnum = buffer->firstfree;
   assert(!buffer->used[bufnum]);
   if( buffer->size[bufnum] < size )
   {
      int newsize;

      /* enlarge buffer */
      newsize = SCIPsetCalcMemGrowSize(set, size);
      SCIP_ALLOC( BMSreallocMemorySize(&buffer->data[bufnum], newsize) );
      buffer->size[bufnum] = newsize;
   }
   assert(buffer->size[bufnum] >= size);

   *ptr = buffer->data[bufnum];
   buffer->used[bufnum] = TRUE;
   buffer->firstfree++;

   SCIPdebugMessage("allocated buffer %d/%d at %p of size %d (required size: %d) for pointer %p\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], size, (void*)ptr);

#else
   SCIP_ALLOC( BMSallocMemorySize(ptr, size) );
#endif

   return SCIP_OKAY;
}

/** allocates the next unused buffer and copies the given memory into the buffer */
SCIP_RETCODE SCIPbufferDuplicateMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   const void*           source,             /**< memory block to copy into the buffer */
   int                   size                /**< minimal required size of the buffer */
   )
{
   assert(source != NULL);

   /* allocate a buffer of the given size */
   SCIP_CALL( SCIPbufferAllocMem(buffer, set, ptr, size) );

   /* copy the source memory into the buffer */
   BMScopyMemorySize(*ptr, source, size);
   
   return SCIP_OKAY;
}

/** reallocates the buffer to at least the given size */
SCIP_RETCODE SCIPbufferReallocMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to the allocated memory buffer */
   int                   size                /**< minimal required size of the buffer */
   )
{
#ifndef SCIP_NOBUFFERMEM
   int bufnum;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
   assert(ptr != NULL);
   assert(size >= 0);

   /* if the pointer doesn't exist yet, allocate it */
   if( *ptr == NULL )
      return SCIPbufferAllocMem(buffer, set, ptr, size);

   assert(buffer->firstfree >= 1);

   /* Search the pointer in the buffer list
    * Usually, buffers are allocated and freed like a stack, such that the currently used pointer is
    * most likely at the end of the buffer list.
    */
   for( bufnum = buffer->firstfree-1; bufnum >= 0 && buffer->data[bufnum] != *ptr; --bufnum )
   {
   }
   assert(bufnum >= 0);
   assert(buffer->data[bufnum] == *ptr);
   assert(buffer->used[bufnum]);
   assert(buffer->size[bufnum] >= 1);

   /* check if the buffer has to be enlarged */
   if( size > buffer->size[bufnum] )
   {
      int newsize;

      /* enlarge buffer */
      newsize = SCIPsetCalcMemGrowSize(set, size);
      SCIP_ALLOC( BMSreallocMemorySize(&buffer->data[bufnum], newsize) );
      buffer->size[bufnum] = newsize;
      *ptr = buffer->data[bufnum];
   }
   assert(buffer->size[bufnum] >= size);
   assert(*ptr == buffer->data[bufnum]);

   SCIPdebugMessage("reallocated buffer %d/%d at %p to size %d (required size: %d) for pointer %p\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], size, (void*)ptr);

#else
   SCIP_ALLOC( BMSreallocMemorySize(ptr, size) );
#endif

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** allocates the next unused buffer; checks for integer overflow */
SCIP_RETCODE SCIPbufferAllocMemSave(
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   int                   num,                /**< number of entries to allocate */
   size_t                elemsize            /**< size of one element in the array */
   )
{
   if( ((size_t)(num)) > (UINT_MAX / elemsize) )
   {
      *ptr = NULL;
      return SCIP_NOMEMORY;
   }

   SCIP_CALL( SCIPbufferAllocMem((set)->buffer, set, (void**)(ptr), (int)(num*elemsize)) );

   return SCIP_OKAY;
}

/** allocates the next unused buffer and copies the given memory into the buffer; checks for integer overflows */
SCIP_RETCODE SCIPbufferDuplicateMemSave(
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   const void*           source,             /**< memory block to copy into the buffer */
   int                   num,                /**< number of entries to copy */
   size_t                elemsize            /**< size of one element in the array */
   )
{
   if( ((size_t)(num)) > (UINT_MAX / elemsize) )
   {
      *ptr = NULL;
      return SCIP_NOMEMORY;
   }

   SCIP_CALL( SCIPbufferDuplicateMem((set)->buffer, set, (void**)(ptr), source, (int)(num*elemsize)) );

   return SCIP_OKAY;
}

/** reallocates the buffer to at least the given size; checks for integer overflows */
SCIP_RETCODE SCIPbufferReallocMemSave(
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to the allocated memory buffer */
   int                   num,                /**< number of entries to get memory for */
   size_t                elemsize            /**< size of one element in the array */
   )
{
   if( ((size_t)(num)) > (UINT_MAX / elemsize) )
   {
      *ptr = NULL;
      return SCIP_NOMEMORY;
   }

   SCIP_CALL( SCIPbufferReallocMem((set)->buffer, set, (void**)(ptr), (int)(num*elemsize)) );

   return SCIP_OKAY;
}
#endif

/** frees a buffer */
void SCIPbufferFreeMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   void**                ptr,                /**< pointer to the allocated memory buffer */
   int                   dummysize           /**< used to get a safer define for SCIPsetFreeBufferSize/Array */
   )
{  /*lint --e{715}*/
#ifndef SCIP_NOBUFFERMEM
   int bufnum;

   assert(buffer != NULL);
   assert(buffer->firstfree <= buffer->ndata);
   assert(buffer->firstfree >= 1);
   assert(dummysize == 0);

   /* Search the pointer in the buffer list
    * Usually, buffers are allocated and freed like a stack, such that the freed pointer is
    * most likely at the end of the buffer list.
    */
   for( bufnum = buffer->firstfree-1; bufnum >= 0 && buffer->data[bufnum] != *ptr; --bufnum )
   {
   }
   assert(bufnum >= 0);
   assert(buffer->data[bufnum] == *ptr);
   assert(buffer->used[bufnum]);

   *ptr = NULL;
   buffer->used[bufnum] = FALSE;

   while( buffer->firstfree > 0 && !buffer->used[buffer->firstfree-1] )
      buffer->firstfree--;

   SCIPdebugMessage("freed buffer %d/%d at %p of size %d for pointer %p, first free is %d\n", 
      bufnum, buffer->ndata, buffer->data[bufnum], buffer->size[bufnum], (void*)ptr, buffer->firstfree);

#else
   BMSfreeMemory(ptr);
#endif
}

/** gets number of used buffers */
int SCIPbufferGetNUsed(
   SCIP_BUFFER*          buffer              /**< memory buffer storage */
   )
{
   assert(buffer != NULL);

   return buffer->firstfree;
}

/** outputs statistics about currently allocated buffers to the screen */
void SCIPbufferPrint(
   SCIP_BUFFER*          buffer              /**< memory buffer storage */
   )
{
   int totalmem;
   int i;

   assert(buffer != NULL);

   totalmem = 0;
   for( i = 0; i < buffer->ndata; ++i )
   {
      printf("[%c] %8d bytes at %p\n", buffer->used[i] ? '*' : ' ', buffer->size[i], buffer->data[i]);
      totalmem += buffer->size[i];
   }
   printf("    %8d bytes total in %d buffers\n", totalmem, buffer->ndata);
}


