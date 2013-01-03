/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the library                         */
/*          BMS --- Block Memory Shell                                       */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  BMS is distributed under the terms of the ZIB Academic License.          */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with BMS; see the file COPYING. If not email to achterberg@zib.de. */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   memory.c
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#ifdef WITH_SCIPDEF
#include "scip/def.h"
#include "scip/pub_message.h"
#endif

#include "blockmemshell/memory.h"

/*#define CHECKMEM*/


/* if we are included in SCIP, use SCIP's message output methods */
#ifdef SCIPdebugMessage
#define debugMessage SCIPdebugMessage
#define errorMessage SCIPerrorMessage
#else
#define debugMessage while( FALSE ) printf
#define errorMessage printf
#define printErrorHeader(f,l) printf("[%s:%d] ERROR: ", f, l)
#define printError printf
#endif

#define warningMessage printf
#define printInfo printf

/* define some macros (if not already defined) */
#ifndef FALSE
#define FALSE 0
#define TRUE  1
#endif
#ifndef MAX
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#endif




/*************************************************************************************
 * Standard Memory Management
 *
 * In debug mode, these methods extend malloc() and free() by logging all currently
 * allocated memory elements in an allocation list. This can be used as a simple leak
 * detection.
 *************************************************************************************/
#if !defined(NDEBUG) && defined(NPARASCIP)

typedef struct Memlist MEMLIST;         /**< memory list for debugging purposes */

/** memory list for debugging purposes */
struct Memlist
{
   const void*           ptr;                /**< pointer to allocated memory */
   size_t                size;               /**< size of memory element */
   char*                 filename;           /**< source file where the allocation was performed */
   int                   line;               /**< line number in source file where the allocation was performed */
   MEMLIST*              next;               /**< next entry in the memory list */
};

static MEMLIST*          memlist = NULL;     /**< global memory list for debugging purposes */
static long long         memused = 0;        /**< number of allocated bytes */

#ifdef CHECKMEM
/** checks, whether the number of allocated bytes match the entries in the memory list */
static
void checkMemlist(
   void
   )
{
   MEMLIST* list = memlist;
   long long used = 0;

   while( list != NULL )
   {
      used += list->size;
      list = list->next;
   }
   assert(used == memused);
}
#else
#define checkMemlist() /**/
#endif

/** adds entry to list of allocated memory */
static
void addMemlistEntry(
   const void*           ptr,                /**< pointer to allocated memory */
   size_t                size,               /**< size of memory element */
   const char*           filename,           /**< source file where the allocation was performed */
   int                   line                /**< line number in source file where the allocation was performed */
   )
{
   MEMLIST* list;

   assert(ptr != NULL && size > 0);

   list = (MEMLIST*)malloc(sizeof(MEMLIST));
   assert(list != NULL);

   list->ptr = ptr;
   list->size = size;
   list->filename = strdup(filename);
   assert(list->filename != NULL);
   list->line = line;
   list->next = memlist;
   memlist = list;
   memused += (long long)size;
   checkMemlist();
}

/** removes entry from the list of allocated memory */
static
void removeMemlistEntry(
   const void*           ptr,                /**< pointer to allocated memory */
   const char*           filename,           /**< source file where the deallocation was performed */
   int                   line                /**< line number in source file where the deallocation was performed */
   )
{
   MEMLIST* list;
   MEMLIST** listptr;

   assert(ptr != NULL);

   list = memlist;
   listptr = &memlist;
   while( list != NULL && ptr != list->ptr )
   {
      listptr = &(list->next);
      list = list->next;
   }
   if( list != NULL )
   {
      assert(ptr == list->ptr);

      *listptr = list->next;
      memused -= (long long)list->size;
      free(list->filename);
      free(list);
   }
   else
   {
      printErrorHeader(filename, line);
      printError("Tried to free unknown pointer <%p>\n", ptr);
   }
   checkMemlist();
}

/** returns the size of an allocated memory element */
size_t BMSgetPointerSize_call(
   const void*           ptr                 /**< pointer to allocated memory */
   )
{
   MEMLIST* list;

   list = memlist;
   while( list != NULL && ptr != list->ptr )
      list = list->next;
   if( list != NULL )
      return list->size;
   else
      return 0;
}

/** outputs information about currently allocated memory to the screen */
void BMSdisplayMemory_call(
   void
   )
{
   MEMLIST* list;
   long long used;

   printInfo("Allocated memory:\n");
   list = memlist;
   used = 0;
   while( list != NULL )
   {
      printInfo("%12p %8lld %s:%d\n", list->ptr, (long long)(list->size), list->filename, list->line);
      used += (long long)list->size;
      list = list->next;
   }
   printInfo("Total:    %8lld\n", memused);
   if( used != memused )
   {
      errorMessage("Used memory in list sums up to %lld instead of %lld\n", used, memused);
   }
   checkMemlist();
}

/** displays a warning message on the screen, if allocated memory exists */
void BMScheckEmptyMemory_call(
   void
   )
{
   if( memlist != NULL || memused > 0 )
   {
      warningMessage("Memory list not empty.\n");
      BMSdisplayMemory_call();
   }
}

/** returns total number of allocated bytes */
long long BMSgetMemoryUsed_call(
   void
   )
{
   return memused;
}

#else

/* these methods are implemented even in optimized mode, such that a program, that includes memory.h in debug mode
 * but links the optimized version compiles
 */

/** returns the size of an allocated memory element */
size_t BMSgetPointerSize_call(
   const void*           ptr                 /**< pointer to allocated memory */
   )
{
   return 0;
}

/** outputs information about currently allocated memory to the screen */
void BMSdisplayMemory_call(
   void
   )
{
#ifdef NPARASCIP
   printInfo("optimized version of memory shell linked - no memory diagnostics available\n");
#endif
}

/** displays a warning message on the screen, if allocated memory exists */
void BMScheckEmptyMemory_call(
   void
   )
{
#ifdef NPARASCIP
   printInfo("optimized version of memory shell linked - no memory leakage check available\n");
#endif
}

/** returns total number of allocated bytes */
long long BMSgetMemoryUsed_call(
   void
   )
{
   return 0;
}

#endif

/** allocates memory and initializes it with 0; returns NULL if memory allocation failed */
void* BMSallocClearMemory_call(
   size_t                num,                /**< number of memory element to allocate */
   size_t                size,               /**< size of one memory element to allocate */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   )
{
   void* ptr;

   debugMessage("calloc %lld elements of %lld bytes [%s:%d]\n", (long long) num, (long long)size, filename, line);

   num = MAX(num, 1);
   size = MAX(size, 1);
   ptr = calloc(num, size);

   if( ptr == NULL )
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for allocation of %lld bytes\n", ((long long) num) * ((long long) size));
   }
#if !defined(NDEBUG) && defined(NPARASCIP)
   else
      addMemlistEntry(ptr, num*size, filename, line);
#endif

   return ptr;
}

/** allocates memory; returns NULL if memory allocation failed */
void* BMSallocMemory_call(
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   )
{
   void* ptr;

   debugMessage("malloc %lld bytes [%s:%d]\n", (long long)size, filename, line);

   size = MAX(size, 1);
   ptr = malloc(size);

   if( ptr == NULL )
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for allocation of %lld bytes\n", (long long) size);
   }
#if !defined(NDEBUG) && defined(NPARASCIP)
   else
      addMemlistEntry(ptr, size, filename, line);
#endif

   return ptr;
}

/** allocates memory; returns NULL if memory allocation failed */
void* BMSreallocMemory_call(
   void*                 ptr,                /**< pointer to memory to reallocate */
   size_t                size,               /**< new size of memory element */
   const char*           filename,           /**< source file where the reallocation is performed */
   int                   line                /**< line number in source file where the reallocation is performed */
   )
{
   void* newptr;

#if !defined(NDEBUG) && defined(NPARASCIP)
   if( ptr != NULL )
      removeMemlistEntry(ptr, filename, line);
#endif

   size = MAX(size, 1);
   newptr = realloc(ptr, size);

   if( newptr == NULL )
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for reallocation of %lld bytes\n", (long long) size);
   }
#if !defined(NDEBUG) && defined(NPARASCIP)
   else
      addMemlistEntry(newptr, size, filename, line);
#endif

   /*   fprintf(stderr, "mem: %p %d\n", newptr, (int)size); */

   return newptr;
}

/** clears a memory element (i.e. fills it with zeros) */
void BMSclearMemory_call(
   void*                 ptr,                /**< pointer to memory element */
   size_t                size                /**< size of memory element */
   )
{
   if( size > 0 )
   {
      assert(ptr != NULL);
      memset(ptr, 0, size);
   }
}

/** copies the contents of one memory element into another memory element */
void BMScopyMemory_call(
   void*                 ptr,                /**< pointer to target memory element */
   const void*           source,             /**< pointer to source memory element */
   size_t                size                /**< size of memory element to copy */
   )
{
   if( size > 0 )
   {
      assert(ptr != NULL);
      assert(source != NULL);
      memcpy(ptr, source, size);
   }
}

/** moves the contents of one memory element into another memory element, should be used if both elements overlap,
 *  otherwise BMScopyMemory is faster
 */
void BMSmoveMemory_call(
   void*                 ptr,                /**< pointer to target memory element */
   const void*           source,             /**< pointer to source memory element */
   size_t                size                /**< size of memory element to copy */
   )
{
   if( size > 0 )
   {
      assert(ptr != NULL);
      assert(source != NULL);
      memmove(ptr, source, size);
   }
}

/** allocates memory and copies the contents of the given memory element into the new memory element */
void* BMSduplicateMemory_call(
   const void*           source,             /**< pointer to source memory element */
   size_t                size,               /**< size of memory element to copy */
   const char*           filename,           /**< source file where the duplication is performed */
   int                   line                /**< line number in source file where the duplication is performed */
   )
{
   void* ptr;

   assert(source != NULL || size == 0);

   ptr = BMSallocMemory_call(size, filename, line);
   if( ptr != NULL )
      BMScopyMemory_call(ptr, source, size);

   return ptr;
}

/** frees an allocated memory element */
void BMSfreeMemory_call(
   void*                 ptr,                /**< pointer to memory element */
   const char*           filename,           /**< source file where the deallocation is performed */
   int                   line                /**< line number in source file where the deallocation is performed */
   )
{
   if( ptr != NULL )
   {
#if !defined(NDEBUG) && defined(NPARASCIP)
      removeMemlistEntry(ptr, filename, line);
#endif
      free(ptr);
   }
   else
   {
      printErrorHeader(filename, line);
      printError("Tried to free null pointer\n");
   }
}




/********************************************************************
 * Chunk Memory Management
 *
 * Efficient memory management for multiple objects of the same size
 ********************************************************************/

/* 
 * block memory methods for faster memory access
 */

#define CHUNKLENGTH_MIN            1024 /**< minimal size of a chunk (in bytes) */
#define CHUNKLENGTH_MAX         1048576 /**< maximal size of a chunk (in bytes) */
#define STORESIZE_MAX              8192 /**< maximal number of elements in one chunk */
#define GARBAGE_SIZE                256 /**< size of lazy free list to start garbage collection */
#define ALIGNMENT    (sizeof(FREELIST)) /**< minimal alignment of chunks */

typedef struct Freelist FREELIST;       /**< linked list of free memory elements */
typedef struct Chunk CHUNK;             /**< chunk of memory elements */

/** linked list of free memory elements */
struct Freelist
{
   FREELIST*        next;               /**< pointer to the next free element */
};

/** chunk of memory elements */
struct Chunk
{
   void*                 store;              /**< data storage */
   void*                 storeend;           /**< points to the first byte in memory not belonging to the chunk */
   FREELIST*             eagerfree;          /**< eager free list */
   CHUNK*                nexteager;          /**< next chunk, that has a non-empty eager free list */
   CHUNK*                preveager;          /**< previous chunk, that has a non-empty eager free list */
   BMS_CHKMEM*           chkmem;             /**< chunk memory collection, this chunk belongs to */
   int                   elemsize;           /**< size of each element in the chunk */
   int                   storesize;          /**< number of elements in this chunk */
   int                   eagerfreesize;      /**< number of elements in the eager free list */
   int                   arraypos;           /**< position of chunk in the chunk header's chunkarray */
}; /* the chunk data structure must be aligned, because the storage is allocated directly behind the chunk header! */

/** collection of memory chunks of the same element size */
struct BMS_ChkMem
{
   FREELIST*             lazyfree;           /**< lazy free list of unused memory elements of all chunks of this chunk block */
   CHUNK**               chunks;             /**< array with the chunks of the chunk header */
   CHUNK*                firsteager;         /**< first chunk with a non-empty eager free list */ 
   BMS_CHKMEM*           nextchkmem;         /**< next chunk block in the block memory's hash list */
   int                   elemsize;           /**< size of each memory element in the chunk memory */
   int                   chunkssize;         /**< size of the chunks array */
   int                   nchunks;            /**< number of chunks in this chunk block (used slots of the chunk array) */
   int                   lastchunksize;      /**< number of elements in the last allocated chunk */
   int                   storesize;          /**< total number of elements in this chunk block */
   int                   lazyfreesize;       /**< number of elements in the lazy free list of the chunk block */
   int                   eagerfreesize;      /**< total number of elements of all eager free lists of the block's chunks */
   int                   initchunksize;      /**< number of elements in the first chunk */
   int                   garbagefactor;      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
#ifndef NDEBUG
   char*                 filename;           /**< source file, where this chunk block was created */
   int                   line;               /**< source line, where this chunk block was created */
   int                   ngarbagecalls;      /**< number of times, the garbage collector was called */
   int                   ngarbagefrees;      /**< number of chunks, the garbage collector freed */
#endif
};


/** aligns the given byte size corresponding to the minimal alignment */
static
void alignSize(
   size_t*               size                /**< pointer to the size to align */
   )
{
   if( *size < ALIGNMENT )
      *size = ALIGNMENT;
   else
      *size = ((*size + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT;
}

/** aligns the given byte size corresponding to the minimal alignment for chunk and block memory  */
void BMSalignMemsize(
   size_t*               size                /**< pointer to the size to align */
   )
{
   assert(ALIGNMENT == sizeof(void*));
   alignSize(size);
}

/** checks whether the given size meets the alignment conditions for chunk and block memory */
int BMSisAligned(
   size_t                size                /**< size to check for alignment */
   )
{
   assert(ALIGNMENT == sizeof(void*));
   return( size >= ALIGNMENT && size % ALIGNMENT == 0 );
}

#ifndef NDEBUG
/** checks, if the given pointer belongs to the given chunk */
static
int isPtrInChunk(
   const CHUNK*          chunk,              /**< memory chunk */
   const void*           ptr                 /**< pointer */
   )
{
   assert(chunk != NULL);
   assert(chunk->store <= chunk->storeend);

   return (ptr >= (void*)(chunk->store) && ptr < (void*)(chunk->storeend));
}
#endif

/** given a pointer, finds the chunk this pointer points to in the chunk array of the given chunk block;
 *  binary search is used;
 *  returns NULL if the pointer does not belong to the chunk block
 */
static
CHUNK* findChunk(
   const BMS_CHKMEM*     chkmem,             /**< chunk block */
   const void*           ptr                 /**< pointer */
   )
{
   CHUNK* chunk;
   int left;
   int right;
   int middle;

   assert(chkmem != NULL);
   assert(ptr != NULL);

   /* binary search for the chunk containing the ptr */
   left = 0;
   right = chkmem->nchunks-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      assert(0 <= middle && middle < chkmem->nchunks);
      chunk = chkmem->chunks[middle];
      assert(chunk != NULL);
      if( ptr < chunk->store )
         right = middle-1;
      else if( ptr >= chunk->storeend )
         left = middle+1;
      else
         return chunk;
   }

   /* ptr was not found in chunk */
   return NULL;
}

/** checks, if a pointer belongs to a chunk of the given chunk block */
static
int isPtrInChkmem(
   const BMS_CHKMEM*     chkmem,             /**< chunk block */
   const void*           ptr                 /**< pointer */
   )
{
   assert(chkmem != NULL);

   return (findChunk(chkmem, ptr) != NULL);
}



/*
 * debugging methods
 */

#ifdef CHECKMEM
/** sanity check for a memory chunk */
static
void checkChunk(
   const CHUNK*       chunk               /**< memory chunk */
   )
{
   FREELIST* eager;
   int eagerfreesize;

   assert(chunk != NULL);
   assert(chunk->store != NULL);
   assert(chunk->storeend == (void*)((char*)(chunk->store) + chunk->elemsize * chunk->storesize));
   assert(chunk->chkmem != NULL);
   assert(chunk->chkmem->elemsize == chunk->elemsize);

   if( chunk->eagerfree == NULL )
      assert(chunk->nexteager == NULL && chunk->preveager == NULL);
   else if( chunk->preveager == NULL )
      assert(chunk->chkmem->firsteager == chunk);

   if( chunk->nexteager != NULL )
      assert(chunk->nexteager->preveager == chunk);
   if( chunk->preveager != NULL )
      assert(chunk->preveager->nexteager == chunk);

   eagerfreesize = 0;
   eager = chunk->eagerfree;
   while( eager != NULL )
   {
      assert(isPtrInChunk(chunk, eager));
      eagerfreesize++;
      eager = eager->next;
   }
   assert(chunk->eagerfreesize == eagerfreesize);
}

/** sanity check for a chunk block */
static
void checkChkmem(
   const BMS_CHKMEM*     chkmem              /**< chunk block */
   )
{
   CHUNK* chunk;
   FREELIST* lazy;
   int nchunks;
   int storesize;
   int lazyfreesize;
   int eagerfreesize;
   int i;

   assert(chkmem != NULL);
   assert(chkmem->chunks != NULL || chkmem->chunkssize == 0);
   assert(chkmem->nchunks <= chkmem->chunkssize);

   nchunks = 0;    
   storesize = 0;    
   lazyfreesize = 0; 
   eagerfreesize = 0;

   for( i = 0; i < chkmem->nchunks; ++i )
   {
      chunk = chkmem->chunks[i];
      assert(chunk != NULL);

      checkChunk(chunk);
      nchunks++;
      storesize += chunk->storesize;
      eagerfreesize += chunk->eagerfreesize;
   }
   assert(chkmem->nchunks == nchunks);
   assert(chkmem->storesize == storesize);
   assert(chkmem->eagerfreesize == eagerfreesize);

   assert((chkmem->eagerfreesize == 0) ^ (chkmem->firsteager != NULL));

   if( chkmem->firsteager != NULL )
      assert(chkmem->firsteager->preveager == NULL);

   lazy = chkmem->lazyfree;
   while( lazy != NULL )
   {
      chunk = findChunk(chkmem, lazy);
      assert(chunk != NULL);
      assert(chunk->chkmem == chkmem);
      lazyfreesize++;
      lazy = lazy->next;
   }
   assert(chkmem->lazyfreesize == lazyfreesize);
}
#else
#define checkChunk(chunk) /**/
#define checkChkmem(chkmem) /**/
#endif


/** links chunk to the block's chunk array, sort it by store pointer;
 *  returns TRUE if successful, FALSE otherwise
 */
static
int linkChunk(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   CHUNK*                chunk               /**< memory chunk */
   )
{
   CHUNK* curchunk;
   int left;
   int right;
   int middle;
   int i;

   assert(chkmem != NULL);
   assert(chkmem->nchunks <= chkmem->chunkssize);
   assert(chunk != NULL);
   assert(chunk->store != NULL);

   debugMessage("linking chunk %p to chunk block %p [elemsize:%d, %d chunks]\n", 
      (void*)chunk, (void*)chkmem, chkmem->elemsize, chkmem->nchunks);

   /* binary search for the position to insert the chunk */
   left = -1;
   right = chkmem->nchunks;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(0 <= middle && middle < chkmem->nchunks);
      assert(left < middle && middle < right);
      curchunk = chkmem->chunks[middle];
      assert(curchunk != NULL);
      if( chunk->store < curchunk->store )
         right = middle;
      else
      {
         assert(chunk->store >= curchunk->storeend);
         left = middle;
      }
   }
   assert(-1 <= left && left < chkmem->nchunks);
   assert(0 <= right && right <= chkmem->nchunks);
   assert(left+1 == right);
   assert(left == -1 || chkmem->chunks[left]->storeend <= chunk->store);
   assert(right == chkmem->nchunks || chunk->storeend <= chkmem->chunks[right]->store);

   /* ensure, that chunk array can store the additional chunk */
   if( chkmem->nchunks == chkmem->chunkssize )
   {
      chkmem->chunkssize = 2*(chkmem->nchunks+1);
      BMSreallocMemoryArray(&chkmem->chunks, chkmem->chunkssize);
      if( chkmem->chunks == NULL )
         return FALSE;
   }
   assert(chkmem->nchunks < chkmem->chunkssize);
   assert(chkmem->chunks != NULL);

   /* move all chunks from 'right' to end one position to the right */
   for( i = chkmem->nchunks; i > right; --i )
   {
      chkmem->chunks[i] = chkmem->chunks[i-1];
      chkmem->chunks[i]->arraypos = i;
   }

   /* insert chunk at position 'right' */
   chunk->arraypos = right;
   chkmem->chunks[right] = chunk;
   chkmem->nchunks++;
   chkmem->storesize += chunk->storesize;

   return TRUE;
}

/** unlinks chunk from the chunk block's chunk list */
static
void unlinkChunk(
   CHUNK*                chunk               /**< memory chunk */
   )
{
   BMS_CHKMEM* chkmem;
   int i;

   assert(chunk != NULL);
   assert(chunk->eagerfree == NULL);
   assert(chunk->nexteager == NULL);
   assert(chunk->preveager == NULL);

   chkmem = chunk->chkmem;
   assert(chkmem != NULL);
   assert(chkmem->elemsize == chunk->elemsize);
   assert(0 <= chunk->arraypos && chunk->arraypos < chkmem->nchunks);
   assert(chkmem->chunks[chunk->arraypos] == chunk);
   
   debugMessage("unlinking chunk %p from chunk block %p [elemsize:%d, %d chunks]\n", 
      (void*)chunk, (void*)chkmem, chkmem->elemsize, chkmem->nchunks);

   /* remove the chunk from the chunks of the chunk block */
   for( i = chunk->arraypos; i < chkmem->nchunks-1; ++i )
   {
      chkmem->chunks[i] = chkmem->chunks[i+1];
      chkmem->chunks[i]->arraypos = i;
   }
   chkmem->nchunks--;
   chkmem->storesize -= chunk->storesize;
}

/** links chunk to the chunk block's eager chunk list */
static
void linkEagerChunk(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   CHUNK*                chunk               /**< memory chunk */
   )
{
   assert(chunk->chkmem == chkmem);
   assert(chunk->nexteager == NULL);
   assert(chunk->preveager == NULL);

   chunk->nexteager = chkmem->firsteager;
   chunk->preveager = NULL;
   if( chkmem->firsteager != NULL )
   {
      assert(chkmem->firsteager->preveager == NULL);
      chkmem->firsteager->preveager = chunk;
   }
   chkmem->firsteager = chunk;
}

/** unlinks chunk from the chunk block's eager chunk list */
static
void unlinkEagerChunk(
   CHUNK*                chunk               /**< memory chunk */
   )
{
   assert(chunk != NULL);
   assert(chunk->eagerfreesize == 0 || chunk->eagerfreesize == chunk->storesize);

   if( chunk->nexteager != NULL )
      chunk->nexteager->preveager = chunk->preveager;
   if( chunk->preveager != NULL )
      chunk->preveager->nexteager = chunk->nexteager;
   else
   {
      assert(chunk->chkmem->firsteager == chunk);
      chunk->chkmem->firsteager = chunk->nexteager;
   }
   chunk->nexteager = NULL;
   chunk->preveager = NULL;
   chunk->eagerfree = NULL;
}

/** creates a new memory chunk in the given chunk block and adds memory elements to the lazy free list;
 *  returns TRUE if successful, FALSE otherwise
 */
static
int createChunk(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   )
{
   CHUNK *newchunk;
   FREELIST *freelist;
   int i;
   int storesize;
   int retval;

   assert(chkmem != NULL);

   debugMessage("creating new chunk in chunk block %p [elemsize: %d]\n", (void*)chkmem, chkmem->elemsize);

   /* calculate store size */
   if( chkmem->nchunks == 0 )
      storesize = chkmem->initchunksize;
   else
      storesize = 2 * chkmem->lastchunksize;
   assert(storesize > 0);
   storesize = MAX(storesize, CHUNKLENGTH_MIN / chkmem->elemsize);
   storesize = MIN(storesize, CHUNKLENGTH_MAX / chkmem->elemsize);
   storesize = MIN(storesize, STORESIZE_MAX);
   storesize = MAX(storesize, 1);
   chkmem->lastchunksize = storesize;

   /* create new chunk */
   assert(BMSisAligned(sizeof(CHUNK)));
   assert( chkmem->elemsize < INT_MAX / storesize );
   assert( sizeof(CHUNK) < UINT_MAX - (size_t)(storesize * chkmem->elemsize) );
   BMSallocMemorySize(&newchunk, sizeof(CHUNK) + storesize * chkmem->elemsize);
   if( newchunk == NULL )
      return FALSE;

   /* the store is allocated directly behind the chunk header */
   newchunk->store = (void*) ((char*) newchunk + sizeof(CHUNK));
   newchunk->storeend = (void*) ((char*) newchunk->store + storesize * chkmem->elemsize);
   newchunk->eagerfree = NULL;
   newchunk->nexteager = NULL;
   newchunk->preveager = NULL;
   newchunk->chkmem = chkmem;
   newchunk->elemsize = chkmem->elemsize;
   newchunk->storesize = storesize;
   newchunk->eagerfreesize = 0;
   newchunk->arraypos = -1;

   debugMessage("allocated new chunk %p: %d elements with size %lld\n", 
      (void*)newchunk, newchunk->storesize, (long long)(newchunk->elemsize));

   /* add new memory to the lazy free list */
   for( i = 0; i < newchunk->storesize - 1; ++i )
   {
      freelist = (FREELIST*) ((char*) (newchunk->store) + i * chkmem->elemsize); /*lint !e826*/
      freelist->next = (FREELIST*) ((char*) (newchunk->store) + (i + 1) * chkmem->elemsize); /*lint !e826*/
   }

   freelist = (FREELIST*) ((char*) (newchunk->store) + (newchunk->storesize - 1) * chkmem->elemsize); /*lint !e826*/
   freelist->next = chkmem->lazyfree;
   chkmem->lazyfree = (FREELIST*) (newchunk->store);
   chkmem->lazyfreesize += newchunk->storesize;

   /* link chunk into chunk block */
   retval = linkChunk(chkmem, newchunk);

   checkChkmem(chkmem);

   return retval;
}

/** destroys a chunk without updating the chunk lists */
static
void destroyChunk(
   CHUNK*                chunk               /**< memory chunk */
   )
{
   assert(chunk != NULL);

   debugMessage("destroying chunk %p\n", (void*)chunk);

   /* free chunk header and store (allocated in one call) */
   BMSfreeMemory(&chunk);
}

/** removes a completely unused chunk, i.e. a chunk with all elements in the eager free list */
static
void freeChunk(
   CHUNK*                chunk               /**< memory chunk */
   )
{
   assert(chunk != NULL);
   assert(chunk->store != NULL);
   assert(chunk->eagerfree != NULL);
   assert(chunk->chkmem != NULL);
   assert(chunk->chkmem->chunks != NULL);
   assert(chunk->chkmem->firsteager != NULL);
   assert(chunk->eagerfreesize == chunk->storesize);

   debugMessage("freeing chunk %p of chunk block %p [elemsize: %d]\n", (void*)chunk, (void*)chunk->chkmem, chunk->chkmem->elemsize);

   /* count the deleted eager free slots */
   chunk->chkmem->eagerfreesize -= chunk->eagerfreesize;
   assert(chunk->chkmem->eagerfreesize >= 0);

   /* remove chunk from eager chunk list */
   unlinkEagerChunk(chunk);

   /* remove chunk from chunk list */
   unlinkChunk(chunk);

   /* destroy the chunk */
   destroyChunk(chunk);
}

/** returns an element of the eager free list and removes it from the list */
static
void* allocChunkElement(
   CHUNK*                chunk               /**< memory chunk */
   )
{
   FREELIST* ptr;

   assert(chunk != NULL);
   assert(chunk->eagerfree != NULL);
   assert(chunk->eagerfreesize > 0);
   assert(chunk->chkmem != NULL);

   debugMessage("allocating chunk element in chunk %p [elemsize: %d]\n", (void*)chunk, chunk->chkmem->elemsize);

   /* unlink first element in the eager free list */
   ptr = chunk->eagerfree;
   chunk->eagerfree = ptr->next;
   chunk->eagerfreesize--;
   chunk->chkmem->eagerfreesize--;

   assert((chunk->eagerfreesize == 0 && chunk->eagerfree == NULL)
      ||  (chunk->eagerfreesize != 0 && chunk->eagerfree != NULL));
   assert(chunk->chkmem->eagerfreesize >= 0);

   /* unlink chunk from eager chunk list if necessary */
   if( chunk->eagerfree == NULL )
   {
      assert(chunk->eagerfreesize == 0);
      unlinkEagerChunk(chunk);
   }

   checkChunk(chunk);

   return (void*) ptr;
}

/** puts given pointer into the eager free list and adds the chunk to the eager list of its chunk block, if necessary */
static
void freeChunkElement(
   CHUNK*                chunk,              /**< memory chunk */
   void*                 ptr                 /**< pointer */
   )
{
   assert(chunk != NULL);
   assert(chunk->chkmem != NULL);
   assert(isPtrInChunk(chunk, ptr));

   debugMessage("freeing chunk element %p of chunk %p [elemsize: %d]\n", (void*)ptr, (void*)chunk, chunk->chkmem->elemsize);

   /* link chunk to the eager chunk list if necessary */
   if( chunk->eagerfree == NULL )
   {
      assert(chunk->eagerfreesize == 0);
      linkEagerChunk(chunk->chkmem, chunk);
   }

   /* add ptr to the chunks eager free list */
   ((FREELIST*)ptr)->next = chunk->eagerfree;
   chunk->eagerfree = (FREELIST*)ptr;
   chunk->eagerfreesize++;
   chunk->chkmem->eagerfreesize++;

   checkChunk(chunk);
}

/** creates a new chunk block data structure */
static
BMS_CHKMEM* createChkmem(
   int                   size,               /**< element size of the chunk block */
   int                   initchunksize,      /**< number of elements in the first chunk of the chunk block */
   int                   garbagefactor       /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   )
{
   BMS_CHKMEM* chkmem;

   assert(BMSisAligned((size_t)size));

   BMSallocMemory(&chkmem);
   if( chkmem == NULL )
      return NULL;

   chkmem->lazyfree = NULL;
   chkmem->chunks = NULL;
   chkmem->firsteager = NULL;
   chkmem->nextchkmem = NULL;
   chkmem->elemsize = size;
   chkmem->chunkssize = 0;
   chkmem->nchunks = 0;
   chkmem->lastchunksize = 0;
   chkmem->storesize = 0;
   chkmem->lazyfreesize = 0;
   chkmem->eagerfreesize = 0;
   chkmem->initchunksize = initchunksize;
   chkmem->garbagefactor = garbagefactor;
#ifndef NDEBUG
   chkmem->filename = NULL;
   chkmem->line = 0;
   chkmem->ngarbagecalls = 0;
   chkmem->ngarbagefrees = 0;
#endif

   return chkmem;
}

/** destroys all chunks of the chunk block, but keeps the chunk block header structure */
static
void clearChkmem(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   )
{
   int i;

   assert(chkmem != NULL);

   /* destroy all chunks of the chunk block */
   for( i = 0; i < chkmem->nchunks; ++i )
      destroyChunk(chkmem->chunks[i]);

   chkmem->lazyfree = NULL;
   chkmem->firsteager = NULL;
   chkmem->nchunks = 0;
   chkmem->lastchunksize = 0;
   chkmem->storesize = 0;
   chkmem->lazyfreesize = 0;
   chkmem->eagerfreesize = 0;
}

/** deletes chunk block and frees all associated memory chunks */
static
void destroyChkmem(
   BMS_CHKMEM**          chkmem              /**< pointer to chunk block */
   )
{
   assert(chkmem != NULL);
   assert(*chkmem != NULL);

   clearChkmem(*chkmem);
   BMSfreeMemoryArrayNull(&(*chkmem)->chunks);

#ifndef NDEBUG
   BMSfreeMemoryArrayNull(&(*chkmem)->filename);
#endif

   BMSfreeMemory(chkmem);
}

/** allocates a new memory element from the chunk block */
static
void* allocChkmemElement(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   )
{
   FREELIST* ptr;

   assert(chkmem != NULL);

   /* if the lazy freelist is empty, we have to find the memory element somewhere else */
   if( chkmem->lazyfree == NULL )
   {
      assert(chkmem->lazyfreesize == 0);

      /* check for a free element in the eager freelists */
      if( chkmem->firsteager != NULL )
	 return allocChunkElement(chkmem->firsteager);

      /* allocate a new chunk */
      if( !createChunk(chkmem) )
	 return NULL;
   }

   /* now the lazy freelist should contain an element */
   assert(chkmem->lazyfree != NULL);
   assert(chkmem->lazyfreesize > 0);

   ptr = chkmem->lazyfree;
   chkmem->lazyfree = ptr->next;
   chkmem->lazyfreesize--;

   checkChkmem(chkmem);

   return (void*) ptr;
}

/** sorts the lazy free list of the chunk block into the eager free lists of the chunks, and removes completely
 *  unused chunks
 */
static
void garbagecollectChkmem(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   )
{
   CHUNK* chunk;
   CHUNK* nexteager;
   FREELIST* lazyfree;

   assert(chkmem != NULL);

   debugMessage("garbage collection for chunk block %p [elemsize: %d]\n", (void*)chkmem, chkmem->elemsize);

   /* check, if the chunk block is completely unused */
   if( chkmem->lazyfreesize + chkmem->eagerfreesize == chkmem->storesize )
   {
      clearChkmem(chkmem);
      return;
   }

#ifndef NDEBUG
   chkmem->ngarbagecalls++;
#endif

   /* put the lazy free elements into the eager free lists */
   while( chkmem->lazyfree != NULL )
   {
      /* unlink first element from the lazy free list */
      lazyfree = chkmem->lazyfree;
      chkmem->lazyfree = chkmem->lazyfree->next;
      chkmem->lazyfreesize--;

      /* identify the chunk of the element */
      chunk = findChunk(chkmem, (void*)lazyfree);
#ifndef NDEBUG
      if( chunk == NULL )
      {
         errorMessage("chunk for lazy free chunk %p not found in chunk block %p\n", (void*)lazyfree, (void*)chkmem);
      }
#endif
      assert(chunk != NULL);

      /* add the element to the chunk's eager free list */
      freeChunkElement(chunk, (void*)lazyfree);
      assert(chunk->eagerfreesize > 0);
   }
   assert(chkmem->lazyfreesize == 0);

   /* delete completely unused chunks, but keep at least one */
   chunk = chkmem->firsteager;
   while( chunk != NULL && chkmem->nchunks > 1 )
   {
      nexteager = chunk->nexteager;
      if( chunk->eagerfreesize == chunk->storesize )
      {
#ifndef NDEBUG
	 chkmem->ngarbagefrees++;
#endif
	 freeChunk(chunk);
      }
      chunk = nexteager;
   }

   checkChkmem(chkmem);
}

/** frees a memory element and returns it to the lazy freelist of the chunk block */
static
void freeChkmemElement(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   void*                 ptr,                /**< memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{  /*lint --e{715}*/
   assert(chkmem != NULL);
   assert(ptr != NULL);

#ifdef BMS_CHKMEM
   /* check, if ptr belongs to the chunk block */
   if( !isPtrInChkmem(chkmem, ptr) )
   {
      BMS_CHKMEM* correctchkmem;

      printErrorHeader(filename, line);
      printError("pointer %p does not belong to chunk block %p (size: %lld)\n", 
         ptr, chkmem, (long long)(chkmem->elemsize));
   }
#endif

   /* put ptr in lazy free list */
   ((FREELIST*)ptr)->next = chkmem->lazyfree;
   chkmem->lazyfree = (FREELIST*)ptr;
   chkmem->lazyfreesize++;

   /* check if we want to apply garbage collection */
   if( chkmem->garbagefactor >= 0 && chkmem->nchunks > 0 && chkmem->lazyfreesize >= GARBAGE_SIZE
      && chkmem->lazyfreesize + chkmem->eagerfreesize
      > chkmem->garbagefactor * (double)(chkmem->storesize) / (double)(chkmem->nchunks) )
   {
      garbagecollectChkmem(chkmem);
   }

   checkChkmem(chkmem);
}

/** creates a new chunk block data structure */
BMS_CHKMEM* BMScreateChunkMemory_call(
   size_t                size,               /**< element size of the chunk block */
   int                   initchunksize,      /**< number of elements in the first chunk of the chunk block */
   int                   garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   BMS_CHKMEM* chkmem;

   alignSize(&size);
   chkmem = createChkmem((int) size, initchunksize, garbagefactor);
   if( chkmem == NULL )
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for chunk block\n");
   }
   debugMessage("created chunk memory %p [elemsize: %d]\n", (void*)chkmem, (int)size);

   return chkmem;
}

/** clears a chunk block data structure */
void BMSclearChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   debugMessage("clearing chunk memory %p [elemsize: %d]\n", (void*)chkmem, chkmem->elemsize);

   if( chkmem != NULL )
      clearChkmem(chkmem);
   else
   {
      printErrorHeader(filename, line);
      printError("Tried to clear null chunk block\n");
   }
}

/** destroys and frees a chunk block data structure */
void BMSdestroyChunkMemory_call(
   BMS_CHKMEM**          chkmem,             /**< pointer to chunk block */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   assert(chkmem != NULL);

   debugMessage("destroying chunk memory %p [elemsize: %d]\n", (void*)*chkmem, (*chkmem)->elemsize);

   if( *chkmem != NULL )
      destroyChkmem(chkmem);
   else
   {
      printErrorHeader(filename, line);
      printError("Tried to destroy null chunk block\n");
   }
}

/** allocates a memory element of the given chunk block */
void* BMSallocChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   void* ptr;

   assert(chkmem != NULL);
   assert((int)size == chkmem->elemsize);

   /* get memory inside the chunk block */
   ptr = allocChkmemElement(chkmem);
   if( ptr == NULL )
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for new chunk\n");
   }
   debugMessage("alloced %8lld bytes in %p [%s:%d]\n", (long long)size, (void*)ptr, filename, line);

   checkChkmem(chkmem);

   return ptr;
}

/** duplicates a given memory element by allocating a new element of the same chunk block and copying the data */
void* BMSduplicateChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   const void*           source,             /**< source memory element */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   void* ptr;

   assert(chkmem != NULL);
   assert(source != NULL);
   assert((int)size == chkmem->elemsize);

   ptr = BMSallocChunkMemory_call(chkmem, size, filename, line);
   if( ptr != NULL )
      BMScopyMemorySize(ptr, source, chkmem->elemsize);

   return ptr;
}

/** frees a memory element of the given chunk block */
void BMSfreeChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   void*                 ptr,                /**< memory element to free */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   assert(chkmem != NULL);
   assert((int)size == chkmem->elemsize);
   
   debugMessage("free    %8lld bytes in %p [%s:%d]\n", (long long)chkmem->elemsize, ptr, filename, line);

   if( ptr == NULL )
   {
      printErrorHeader(filename, line);
      printError("Tried to free null block pointer\n");
      return;
   }

   /* free memory in chunk block */
   freeChkmemElement(chkmem, ptr, filename, line);
   
   checkChkmem(chkmem);
}

/** calls garbage collection of chunk block and frees chunks without allocated memory elements */
void BMSgarbagecollectChunkMemory_call(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   )
{
   debugMessage("garbage collection on chunk memory %p [elemsize: %d]\n", (void*)chkmem, chkmem->elemsize);

   garbagecollectChkmem(chkmem);
}

/** returns the number of allocated bytes in the chunk block */
long long BMSgetChunkMemoryUsed_call(
   const BMS_CHKMEM*     chkmem              /**< chunk block */
   )
{
   long long chkmemused;
   int i;

   assert(chkmem != NULL);

   chkmemused = 0;
   for( i = 0; i < chkmem->nchunks; ++i )
      chkmemused += (long long)(chkmem->chunks[i]->elemsize) * (long long)(chkmem->chunks[i]->storesize);
   
   return chkmemused;
}




/***********************************************************
 * Block Memory Management
 *
 * Efficient memory management for objects of varying sizes
 ***********************************************************/

#define CHKHASH_SIZE               1013 /**< size of chunk block hash table; should be prime */

/** collection of chunk blocks */
struct BMS_BlkMem
{
   BMS_CHKMEM*           chkmemhash[CHKHASH_SIZE]; /**< hash table with chunk blocks */
   long long             memused;            /**< total number of used bytes in the memory header */
   int                   initchunksize;      /**< number of elements in the first chunk of each chunk block */
   int                   garbagefactor;      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
};



/*
 * debugging methods
 */

#ifdef CHECKMEM
static
void checkBlkmem(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   )
{
   const BMS_CHKMEM* chkmem;
   int i;

   assert(blkmem != NULL);
   assert(blkmem->chkmemhash != NULL);

   for( i = 0; i < CHKHASH_SIZE; ++i )
   {
      chkmem = blkmem->chkmemhash[i];
      while( chkmem != NULL )
      {
         checkChkmem(chkmem);
	 chkmem = chkmem->nextchkmem;
      }
   }
}
#else
#define checkBlkmem(blkmem) /**/
#endif


/** finds the chunk block, to whick the given pointer belongs to;
 *  this could be done by selecting the chunk block of the corresponding element size, but in a case of an
 *  error (free gives an incorrect element size), we want to identify and output the correct element size
 */
static
BMS_CHKMEM* findChkmem(
   const BMS_BLKMEM*     blkmem,             /**< block memory */
   const void*           ptr                 /**< memory element to search */
   )
{
   BMS_CHKMEM* chkmem;
   int i;

   assert(blkmem != NULL);

   chkmem = NULL;
   for( i = 0; chkmem == NULL && i < CHKHASH_SIZE; ++i )
   {
      chkmem = blkmem->chkmemhash[i];
      while( chkmem != NULL && !isPtrInChkmem(chkmem, ptr) )
	 chkmem = chkmem->nextchkmem;
   }

   return chkmem;
}

/** calculates hash number of memory size */
static
int getHashNumber(
   int                   size                /**< element size */
   )
{
   assert(BMSisAligned((size_t)size));

   return (size % CHKHASH_SIZE);
}

/** creates a block memory allocation data structure */
BMS_BLKMEM* BMScreateBlockMemory_call(
   int                   initchunksize,      /**< number of elements in the first chunk of each chunk block */
   int                   garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   BMS_BLKMEM* blkmem;
   int i;

   BMSallocMemory(&blkmem);
   if( blkmem != NULL )
   {
      for( i = 0; i < CHKHASH_SIZE; ++i )
	 blkmem->chkmemhash[i] = NULL;
      blkmem->initchunksize = initchunksize;
      blkmem->garbagefactor = garbagefactor;
      blkmem->memused = 0;
   }
   else
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for block memory header\n");
   }

   return blkmem;
}

/** frees all chunk blocks in the block memory */
void BMSclearBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   BMS_CHKMEM* chkmem;
   BMS_CHKMEM* nextchkmem;
   int i;

   if( blkmem != NULL )
   {
      for( i = 0; i < CHKHASH_SIZE; ++i )
      {
	 chkmem = blkmem->chkmemhash[i];
	 while( chkmem != NULL )
	 {
	    nextchkmem = chkmem->nextchkmem;
	    destroyChkmem(&chkmem);
	    chkmem = nextchkmem;
	 }
	 blkmem->chkmemhash[i] = NULL;
      }
      blkmem->memused = 0;
   }
   else
   {
      printErrorHeader(filename, line);
      printError("Tried to clear null block memory\n");
   }
}

/** clears and deletes block memory */
void BMSdestroyBlockMemory_call(
   BMS_BLKMEM**          blkmem,             /**< pointer to block memory */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   assert(blkmem != NULL);

   if( *blkmem != NULL )
   {
      BMSclearBlockMemory_call(*blkmem, filename, line);
      BMSfreeMemory(blkmem);
      assert(*blkmem == NULL);
   }
   else
   {
      printErrorHeader(filename, line);
      printError("Tried to destroy null block memory\n");
   }
}

/** allocates memory in the block memory pool */
void* BMSallocBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   BMS_CHKMEM** chkmemptr;
   int hashnumber;
   void* ptr;

   assert( blkmem != NULL );

   /* calculate hash number of given size */
   alignSize(&size);
   hashnumber = getHashNumber((int)size);

   /* find correspoding chunk block */
   chkmemptr = &(blkmem->chkmemhash[hashnumber]);
   while( *chkmemptr != NULL && (*chkmemptr)->elemsize != (int)size )
      chkmemptr = &((*chkmemptr)->nextchkmem);

   /* create new chunk block if necessary */
   if( *chkmemptr == NULL )
   {
      *chkmemptr = createChkmem((int)size, blkmem->initchunksize, blkmem->garbagefactor);
      if( *chkmemptr == NULL )
      {
	 printErrorHeader(filename, line);
         printError("Insufficient memory for chunk block\n");
	 return NULL;
      }
#ifndef NDEBUG
      BMSduplicateMemoryArray(&(*chkmemptr)->filename, filename, strlen(filename) + 1);
      (*chkmemptr)->line = line;
#endif
   }

   /* get memory inside the chunk block */
   ptr = allocChkmemElement(*chkmemptr);
   if( ptr == NULL )
   {
      printErrorHeader(filename, line);
      printError("Insufficient memory for new chunk\n");
   }
   debugMessage("alloced %8lld bytes in %p [%s:%d]\n", (long long)size, ptr, filename, line);

   blkmem->memused += size;

   checkBlkmem(blkmem);

   return ptr;
}

/** resizes memory element in the block memory pool, and copies the data */
void* BMSreallocBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 ptr,                /**< memory element to reallocated */
   size_t                oldsize,            /**< old size of memory element */
   size_t                newsize,            /**< new size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   void* newptr;

   if( ptr == NULL )
   {
      assert(oldsize == 0);
      return BMSallocBlockMemory_call(blkmem, newsize, filename, line);
   }

   alignSize(&oldsize);
   alignSize(&newsize);
   if( oldsize == newsize )
      return ptr;

   newptr = BMSallocBlockMemory_call(blkmem, newsize, filename, line);
   if( newptr != NULL )
      BMScopyMemorySize(newptr, ptr, MIN(oldsize, newsize));
   BMSfreeBlockMemory_call(blkmem, ptr, oldsize, filename, line);

   return newptr;
}

/** duplicates memory element in the block memory pool, and copies the data */
void* BMSduplicateBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const void*           source,             /**< memory element to duplicate */
   size_t                size,               /**< size of memory elements */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   void* ptr;

   assert(source != NULL);

   ptr = BMSallocBlockMemory_call(blkmem, size, filename, line);
   if( ptr != NULL )
      BMScopyMemorySize(ptr, source, size);

   return ptr;
}

/** frees memory element in the block memory pool */
void BMSfreeBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 ptr,                /**< memory element to free */
   size_t                size,               /**< size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   )
{
   BMS_CHKMEM* chkmem;
   int hashnumber;

   assert( blkmem != NULL );

   if( ptr != NULL )
   {
      /* calculate hash number of given size */
      alignSize(&size);
      hashnumber = getHashNumber((int)size);

      debugMessage("free    %8lld bytes in %p [%s:%d]\n", (long long)size, ptr, filename, line);

      /* find correspoding chunk block */
      assert( blkmem->chkmemhash != NULL );
      chkmem = blkmem->chkmemhash[hashnumber];
      while( chkmem != NULL && chkmem->elemsize != (int)size )
	 chkmem = chkmem->nextchkmem;
      if( chkmem == NULL )
      {
	 printErrorHeader(filename, line);
         printError("Tried to free pointer <%p> in block memory <%p> of unknown size %lld\n",
            ptr, (void*)blkmem, (long long) size);
	 return;
      }
      assert(chkmem->elemsize == (int)size);

      /* free memory in chunk block */
      freeChkmemElement(chkmem, ptr, filename, line);
      
      blkmem->memused -= size;
      assert(blkmem->memused >= 0);
   }
   else if( size != 0 )
   {
      printErrorHeader(filename, line);
      printError("Tried to free null block pointer\n");
   }

   checkBlkmem(blkmem);
}

/** calls garbage collection of block memory, frees chunks without allocated memory elements, and frees
 *  chunk blocks without any chunks
 */
void BMSgarbagecollectBlockMemory_call(
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(blkmem != NULL);

   for( i = 0; i < CHKHASH_SIZE; ++i )
   {
      BMS_CHKMEM** chkmemptr;

      chkmemptr = &blkmem->chkmemhash[i];
      while( *chkmemptr != NULL )
      {
         garbagecollectChkmem(*chkmemptr);
         if( (*chkmemptr)->nchunks == 0 )
         {
            BMS_CHKMEM* nextchkmem;

            nextchkmem = (*chkmemptr)->nextchkmem;
            destroyChkmem(chkmemptr);
            *chkmemptr = nextchkmem;
         }
         else
            chkmemptr = &(*chkmemptr)->nextchkmem;
      }
   }
}

/** returns the number of allocated bytes in the block memory */
long long BMSgetBlockMemoryUsed_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);

   return blkmem->memused;
}

/** returns the size of the given memory element; returns 0, if the element is not member of the block memory */
size_t BMSgetBlockPointerSize_call(
   const BMS_BLKMEM*     blkmem,             /**< block memory */
   const void*           ptr                 /**< memory element */
   )
{
   const BMS_CHKMEM* chkmem;

   assert(blkmem != NULL);

   if( ptr == NULL )
      return 0;

   chkmem = findChkmem(blkmem, ptr);
   if( chkmem == NULL )
      return 0;

   return (size_t)(chkmem->elemsize);
}

/** outputs allocation diagnostics of block memory */
void BMSdisplayBlockMemory_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   )
{
   const BMS_CHKMEM* chkmem;
   int nblocks = 0;
   int nunusedblocks = 0;
   int totalnchunks = 0;
   int totalneagerchunks = 0;
   int totalnelems = 0;
   int totalneagerelems = 0;
   int totalnlazyelems = 0;
#ifndef NDEBUG
   int totalngarbagecalls = 0;
   int totalngarbagefrees = 0;
#endif
   long long allocedmem = 0;
   long long freemem = 0;
   int i;
   int c;

#ifndef NDEBUG
   printInfo(" ElSize #Chunk #Eag  #Elems  #EagFr  #LazFr  #GCl #GFr  Free  MBytes First Allocator\n");
#else
   printInfo(" ElSize #Chunk #Eag  #Elems  #EagFr  #LazFr  Free  MBytes\n");
#endif

   assert(blkmem != NULL);

   for( i = 0; i < CHKHASH_SIZE; ++i )
   {
      chkmem = blkmem->chkmemhash[i];
      while( chkmem != NULL )
      {
	 const CHUNK* chunk;
	 int nchunks = 0;
	 int nelems = 0;
	 int neagerchunks = 0;
	 int neagerelems = 0;

         for( c = 0; c < chkmem->nchunks; ++c )
         {
            chunk = chkmem->chunks[c];
            assert(chunk != NULL);
	    assert(chunk->elemsize == chkmem->elemsize);
	    assert(chunk->chkmem == chkmem);
	    nchunks++;
	    nelems += chunk->storesize;
	    if( chunk->eagerfree != NULL )
	    {
	       neagerchunks++;
	       neagerelems += chunk->eagerfreesize;
	    }
	 }

	 assert(nchunks == chkmem->nchunks);
	 assert(nelems == chkmem->storesize);
	 assert(neagerelems == chkmem->eagerfreesize);

	 if( nelems > 0 )
	 {
	    nblocks++;
	    allocedmem += (long long)chkmem->elemsize * (long long)nelems;
	    freemem += (long long)chkmem->elemsize * ((long long)neagerelems + (long long)chkmem->lazyfreesize);

#ifndef NDEBUG
	    printInfo("%7lld %6d %4d %7d %7d %7d %5d %4d %5.1f%% %6.1f %s:%d\n",
	       (long long)(chkmem->elemsize), nchunks, neagerchunks, nelems,
	       neagerelems, chkmem->lazyfreesize, chkmem->ngarbagecalls, chkmem->ngarbagefrees,
	       100.0 * (double) (neagerelems + chkmem->lazyfreesize) / (double) (nelems), 
               (double)chkmem->elemsize * nelems / (1024.0*1024.0),
               chkmem->filename, chkmem->line);
#else
	    printInfo("%7lld %6d %4d %7d %7d %7d %5.1f%% %6.1f\n",
	       (long long)(chkmem->elemsize), nchunks, neagerchunks, nelems,
	       neagerelems, chkmem->lazyfreesize,
	       100.0 * (double) (neagerelems + chkmem->lazyfreesize) / (double) (nelems),
               (double)chkmem->elemsize * nelems / (1024.0*1024.0));
#endif
	 }
	 else
	 {
#ifndef NDEBUG
	    printInfo("%7lld <unused>                            %5d %4d        %s:%d\n",
	       (long long)(chkmem->elemsize), chkmem->ngarbagecalls, chkmem->ngarbagefrees, 
               chkmem->filename, chkmem->line);
#else
	    printInfo("%7lld <unused>\n", (long long)(chkmem->elemsize));
#endif
	    nunusedblocks++;
	 }
         totalnchunks += nchunks;
         totalneagerchunks += neagerchunks;
         totalnelems += nelems;
         totalneagerelems += neagerelems;
         totalnlazyelems += chkmem->lazyfreesize;
#ifndef NDEBUG
         totalngarbagecalls += chkmem->ngarbagecalls;
         totalngarbagefrees += chkmem->ngarbagefrees;
#endif
	 chkmem = chkmem->nextchkmem;
      }
   }
#ifndef NDEBUG
   printInfo("  Total %6d %4d %7d %7d %7d %5d %4d %5.1f%% %6.1f\n",
      totalnchunks, totalneagerchunks, totalnelems, totalneagerelems, totalnlazyelems, 
      totalngarbagecalls, totalngarbagefrees,
      totalnelems > 0 ? 100.0 * (double) (totalneagerelems + totalnlazyelems) / (double) (totalnelems) : 0.0,
      (double)allocedmem/(1024.0*1024.0));
#else
   printInfo("  Total %6d %4d %7d %7d %7d %5.1f%% %6.1f\n",
      totalnchunks, totalneagerchunks, totalnelems, totalneagerelems, totalnlazyelems, 
      totalnelems > 0 ? 100.0 * (double) (totalneagerelems + totalnlazyelems) / (double) (totalnelems) : 0.0,
      (double)allocedmem/(1024.0*1024.0));
#endif
   printInfo("%d blocks (%d unused), %lld bytes allocated, %lld bytes free",
      nblocks + nunusedblocks, nunusedblocks, allocedmem, freemem);
   if( allocedmem > 0 )
      printInfo(" (%.1f%%)", 100.0 * (double) freemem / (double) allocedmem);
   printInfo("\n");
}

/** outputs warning messages, if there are allocated elements in the block memory */
void BMScheckEmptyBlockMemory_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   )
{
   const BMS_CHKMEM* chkmem;
   long long allocedmem = 0;
   long long freemem = 0;
   int i;
   int c;

   assert(blkmem != NULL);

   for( i = 0; i < CHKHASH_SIZE; ++i )
   {
      chkmem = blkmem->chkmemhash[i];
      while( chkmem != NULL )
      {
	 const CHUNK* chunk;
	 int nchunks = 0;
	 int nelems = 0;
	 int neagerelems = 0;

         for( c = 0; c < chkmem->nchunks; ++c )
         {
            chunk = chkmem->chunks[c];
            assert(chunk != NULL);
	    assert(chunk->elemsize == chkmem->elemsize);
	    assert(chunk->chkmem == chkmem);
	    nchunks++;
	    nelems += chunk->storesize;
	    if( chunk->eagerfree != NULL )
	       neagerelems += chunk->eagerfreesize;
	 }

	 assert(nchunks == chkmem->nchunks);
	 assert(nelems == chkmem->storesize);
	 assert(neagerelems == chkmem->eagerfreesize);

	 if( nelems > 0 )
	 {
	    allocedmem += (long long)chkmem->elemsize * (long long)nelems;
	    freemem += (long long)chkmem->elemsize * ((long long)neagerelems + (long long)chkmem->lazyfreesize);

            if( nelems != neagerelems + chkmem->lazyfreesize )
            {
#ifndef NDEBUG
               printInfo("%lld bytes (%d elements of size %lld) not freed. First Allocator: %s:%d\n",
                  (((long long)nelems - (long long)neagerelems) - (long long)chkmem->lazyfreesize)
                  * (long long)(chkmem->elemsize),
                  (nelems - neagerelems) - chkmem->lazyfreesize, (long long)(chkmem->elemsize),
                  chkmem->filename, chkmem->line);
#else
               printInfo("%lld bytes (%d elements of size %lld) not freed.\n",
                  ((nelems - neagerelems) - chkmem->lazyfreesize) * (long long)(chkmem->elemsize),
                  (nelems - neagerelems) - chkmem->lazyfreesize, (long long)(chkmem->elemsize));
#endif
            }
	 }
	 chkmem = chkmem->nextchkmem;
      }
   }

   if( allocedmem != freemem )
      printInfo("%lld bytes not freed in total.\n", allocedmem - freemem);
}
