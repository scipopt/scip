/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: memory.c,v 1.25 2003/12/01 16:14:29 bzfpfend Exp $"

/**@file   memory.c
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "memory.h"


#ifdef DEBUG
#define debugMessage                    printf("[%s:%d] debug: ", __FILE__, __LINE__); printf
#else
#define debugMessage(...)               /**/
#endif

#define TRUE     1
#define FALSE    0
#define MAX(x,y) ((x) >= (y) ? (x) : (y))     /**< returns maximum of x and y */
#define MIN(x,y) ((x) <= (y) ? (x) : (y))     /**< returns minimum of x and y */


/******************************
 * Standard Memory Management *
 ******************************/


#ifndef NOSAFEMEM

/*
 * safe memory management with leakage detection in debug mode
 */


#ifndef NDEBUG

typedef struct MemList MEMLIST;

/** memory list for debugging purposes */
struct MemList
{
   void    *ptr;
   size_t  size;
   char    *filename;
   int     line;
   MEMLIST *next;
};

static MEMLIST *memlist = NULL;
static long memused = 0;

static void
memListAdd(void *ptr, size_t size, const char *filename, int line)
{
   MEMLIST *list = (MEMLIST*)malloc(sizeof(MEMLIST));

   assert(list != NULL);
   assert(ptr != NULL && size > 0);
   list->ptr = ptr;
   list->size = size;
   list->filename = (char *) malloc((strlen(filename) + 1) * sizeof(char));
   strcpy(list->filename, filename);
   list->line = line;
   list->next = memlist;
   memlist = list;
   memused += (long)size;
}

static void
memListRemove(void *ptr, const char *filename, int line)
{
   MEMLIST *list = memlist;
   MEMLIST **listptr = &memlist;

   assert(ptr != NULL);
   while( list != NULL && ptr != list->ptr )
   {
      listptr = &(list->next);
      list = list->next;
   }
   if( list != NULL && ptr == list->ptr )
   {
      *listptr = list->next;
      memused -= (long)list->size;
      free(list->filename);
      free(list);
   }
   else
      fprintf(stderr, "[%s:%d] ERROR: Tried to free unknown pointer <%p>\n", filename, line, ptr);
}

size_t
memorySize(void *ptr)
{
   MEMLIST *list = memlist;

   while( list != NULL && ptr != list->ptr )
      list = list->next;
   if( list != NULL )
      return list->size;
   else
      return 0;
}

void
memoryDiagnostic(void)
{
   MEMLIST *list = memlist;
   long    used = 0;

   printf("Allocated memory:\n");
   while( list != NULL )
   {
      printf("%12p %8ld %s:%d\n", list->ptr, (long) (list->size),
	 list->filename, list->line);
      used += (long)list->size;
      list = list->next;
   }
   printf("Total:    %8ld\n", memused);
   if( used != memused )
      fprintf(stderr, "[%s:%d] ERROR: Used memory in list sums up to %ld instead of %ld\n",
	 __FILE__, __LINE__, used, memused);
}

void
memoryCheckEmpty(void)
{
   if( memlist != NULL || memused > 0 )
   {
      printf("Warning! Memory list not empty.\n");
      memoryDiagnostic();
   }
}
#endif

void *
allocMemory_call(size_t size, const char *filename, int line)
{
   void   *ptr = NULL;

   debugMessage("malloc %ld bytes [%s:%d]\n", (long)size, filename, line);
   size = MAX(size, 1);
   ptr = malloc(size);

   if( ptr == NULL )
      fprintf(stderr, "[%s:%d] ERROR: Insufficient memory for allocation of %ld bytes\n", filename, line, (long) size);
#ifndef NDEBUG
   else
      memListAdd(ptr, size, filename, line);
#endif
   return ptr;
}

void *
reallocMemory_call(void *ptr, size_t size, const char *filename, int line)
{
   void   *newptr = NULL;

#ifndef NDEBUG
   if( ptr != NULL )
      memListRemove(ptr, filename, line);
#endif

   size = MAX(size, 1);
   newptr = realloc(ptr, size);

   if( newptr == NULL )
      fprintf(stderr, "[%s:%d] ERROR: Insufficient memory for reallocation of %ld bytes\n", filename, line, (long) size);
#ifndef NDEBUG
   else
      memListAdd(newptr, size, filename, line);
#endif
   return newptr;
}

void
clearMemory_call(void* ptr, size_t size)
{
   assert(ptr != NULL);
   memset(ptr, 0, size);
}

void
copyMemory_call(void* ptr, const void* source, size_t size)
{
   assert(ptr != NULL);
   assert(source != NULL);
   memcpy(ptr, source, size);
}

void *
duplicateMemory_call(const void* source, size_t size, const char *filename, int line)
{
   void   *ptr = NULL;

   assert(source != NULL);
   ptr = allocMemory_call(size, filename, line);
   if( ptr != NULL )
      copyMemory_call(ptr, source, size);
   return ptr;
}

void
freeMemory_call(void **ptr, const char *filename, int line)
{
   if( *ptr != NULL )
   {
#ifndef NDEBUG
      memListRemove(*ptr, filename, line);
#endif
      free(*ptr);
      *ptr = NULL;
   }
   else
      fprintf(stderr, "[%s:%d] ERROR: Tried to free null pointer\n", filename, line);
}



#else

/*
 * standard memory management via malloc
 */

void *
duplicateMemory_call(const void* source, size_t size)
{
   void   *ptr = NULL;

   allocMemorySize(&ptr, size);
   if( ptr != NULL )
      copyMemorySize(ptr, source, size);
   return ptr;
}

#endif







/***************************
 * Block Memory Management *
 ***************************/


#ifndef NOBLOCKMEM

/* 
 * block memory methods for faster memory access
 */

#define BLOCKHASH_SIZE     1013	/* should be prime */
#define CHUNKLENGTH_MAX 1048576	/* maximal size of a chunk (in bytes) */
#define STORESIZE_MAX      8192	/* maximal number of elements in one chunk */
#define ALIGNMENT       (sizeof(FREELIST))
#define GARBAGE_SIZE        256	/* size of lazy free list to start
				   garbage collection */

typedef struct free_list FREELIST;
typedef struct chunk_header CHKHDR;
typedef struct block_header BLKHDR;

struct free_list
{
   FREELIST*        next;               /**< pointer to the next free element */
};

struct chunk_header
{
   int              elemSize;           /**< size of each element in the chunk (= size of elements in the block) */
   int              storeSize;          /**< number of elements in this chunk */
   int              eagerFreeSize;      /**< number of elements in the eager free list */
   int              arraypos;           /**< position of chunk in the blocks chunkarray */
   void*            store;              /**< data storage */
   void*            storeend;           /**< points to the first byte in memory not belonging to the chunk */
   FREELIST*        eagerFree;          /**< eager free list */
   CHKHDR*          nextEager;          /**< next chunk, that has a non-empty eager free list */
   CHKHDR*          prevEager;          /**< previous chunk, that has a non-empty eager free list */
   BLKHDR*          block;              /**< memory block, this chunk belongs to */
}; /* the chunk header must be aligned! */

struct block_header
{
   int              elemSize;           /**< size of each memory element in the block */
   int              chunkarraySize;     /**< size of the chunk array */
   int              numChunks;          /**< number of chunks in this block (used slots of the chunk array) */
   int              lastChunkSize;      /**< number of elements in the last allocated chunk */
   int              storeSize;          /**< total number of elements in this block */
   int              lazyFreeSize;       /**< number of elements in the lazy free list of the block */
   int              eagerFreeSize;      /**< total number of elements of all eager free lists of the block's chunks */
   FREELIST*        lazyFree;           /**< lazy free list */
   CHKHDR**         chunkarray;         /**< array with the chunks of the block */
   CHKHDR*          firstEager;         /**< first chunk with a non-empty eager free list */ 
   BLKHDR*          next;               /**< next memory block in the hash list */
   MEMHDR*          mem;                /**< memory header, this block belongs to */
#ifndef NDEBUG
   char*            filename;           /**< source file, where this block was created */
   int              line;               /**< source line, where this block was created */
   int              numGarbageCalls;    /**< number of times, the garbage collector was called */
   int              numGarbageFrees;    /**< number of chunks, the garbage collector freed */
#endif
};

struct memory_header
{
   BLKHDR*          blockhash[BLOCKHASH_SIZE]; /**< hash table with memory blocks */
   int              initChunkSize;      /**< number of elements in the first chunk */
   int              clearUnusedBlocks;  /**< TRUE iff unused blocks should be deleted */
   int              garbageFactor;      /**< garbage collector is called, if at least garbageFactor * avg. chunksize 
                                         *   elements are free */
   long long        memused;            /**< total number of used bytes in the memory header */
};


static void
alignSize(size_t* size)
{
   if( *size < ALIGNMENT )
      *size = ALIGNMENT;
   else
      *size = ((*size + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT;
}

void
alignMemsize(size_t* size)
{
   assert(ALIGNMENT == sizeof(void*));
   alignSize(size);
}

int
isAligned(size_t size)
{
   assert(ALIGNMENT == sizeof(void*));
   return( size >= ALIGNMENT && size % ALIGNMENT == 0 );
}

/* checks, if ptr belongs to chunk 'chk' */
static int
isPtrInChunk(const CHKHDR * chk, const void *ptr)
{
   assert(chk != NULL);
   assert(chk->store <= chk->storeend);

   return (ptr >= (void*)(chk->store) && ptr < (void*)(chk->storeend));
}

/* Given a pointer, find the chunk this pointer points to
 * in the chunk array of the block. Binary search is used.
 */
static CHKHDR *
findChunk(BLKHDR * blk, const void *ptr)
{
   CHKHDR* chk;
   int left;
   int right;
   int middle;

   assert(blk != NULL);
   assert(ptr != NULL);

   /* binary search for the chunk containing the ptr */
   left = 0;
   right = blk->numChunks-1;
   while( left <= right )
   {
      middle = (left+right)/2;
      assert(0 <= middle && middle < blk->numChunks);
      chk = blk->chunkarray[middle];
      assert(chk != NULL);
      if( ptr < chk->store )
         right = middle-1;
      else if( ptr >= chk->storeend )
         left = middle+1;
      else
         return chk;
   }

   /* ptr was not found in chunk */
   return NULL;
}

/* checks, if 'ptr' belongs to a chunk of block 'blk' */
static int
isPtrInBlock(BLKHDR * blk, const void *ptr)
{
   assert(blk != NULL);
   return (findChunk(blk, ptr) != NULL);
}

/* finds the block, to whick 'ptr' belongs to.
 * this could be done by selecting the block of the
 * corresponding element size, but in a case of an
 * error (free gives an incorrect element size), we
 * want to identify and output the correct element size
 */
static BLKHDR *
findBlock(MEMHDR *mem, const void *ptr)
{
   BLKHDR *blk = NULL;
   int     i;

   assert(mem != NULL);
   for( i = 0; blk == NULL && i < BLOCKHASH_SIZE; ++i )
   {
      blk = mem->blockhash[i];
      while( blk != NULL && !isPtrInBlock(blk, ptr) )
	 blk = blk->next;
   }
   return blk;
}

#if defined(DEBUG) && !defined(NDEBUG)
/*
 * debugging methods
 */
static void
checkChunk(CHKHDR * chk)
{
   FREELIST* eager;
   int eagerFreeSize;

   assert(chk != NULL);
   assert(chk->store != NULL);
   assert(chk->storeend == (void*)((char*)(chk->store) + chk->elemSize * chk->storeSize));
   assert(chk->block != NULL);
   assert(chk->block->elemSize == chk->elemSize);

   if( chk->eagerFree == NULL )
      assert(chk->nextEager == NULL && chk->prevEager == NULL);
   else if( chk->prevEager == NULL )
      assert(chk->block->firstEager == chk);

   if( chk->nextEager != NULL )
      assert(chk->nextEager->prevEager == chk);
   if( chk->prevEager != NULL )
      assert(chk->prevEager->nextEager == chk);

   eagerFreeSize = 0;
   eager = chk->eagerFree;
   while( eager != NULL )
   {
      assert(isPtrInChunk(chk, eager));
      eagerFreeSize++;
      eager = eager->next;
   }
   assert(chk->eagerFreeSize == eagerFreeSize);
}

static void
checkBlock(BLKHDR * blk)
{
   CHKHDR* chk;
   FREELIST* lazy;
   int numChunks;
   int storeSize;
   int lazyFreeSize;
   int eagerFreeSize;
   int i;

   assert(blk != NULL);
   assert(blk->chunkarray != NULL || blk->chunkarraySize == 0);
   assert(blk->numChunks <= blk->chunkarraySize);

   numChunks = 0;    
   storeSize = 0;    
   lazyFreeSize = 0; 
   eagerFreeSize = 0;

   for( i = 0; i < blk->numChunks; ++i )
   {
      chk = blk->chunkarray[i];
      assert(chk != NULL);

      checkChunk(chk);
      numChunks++;
      storeSize += chk->storeSize;
      eagerFreeSize += chk->eagerFreeSize;
   }
   assert(blk->numChunks == numChunks);
   assert(blk->storeSize == storeSize);
   assert(blk->eagerFreeSize == eagerFreeSize);

   assert((blk->eagerFreeSize == 0) ^ (blk->firstEager != NULL));

   if( blk->firstEager != NULL )
      assert(blk->firstEager->prevEager == NULL);

   lazy = blk->lazyFree;
   while( lazy != NULL )
   {
      chk = findChunk(blk, lazy);
      assert(chk != NULL);
      assert(chk->block == blk);
      lazyFreeSize++;
      lazy = lazy->next;
   }
   assert(blk->lazyFreeSize == lazyFreeSize);
}

static void
checkMem(MEMHDR *mem)
{
   BLKHDR* blk;
   int i;

   assert(mem != NULL);
   assert(mem->blockhash != NULL);

   for( i = 0; i < BLOCKHASH_SIZE; ++i )
   {
      blk = mem->blockhash[i];
      while( blk != NULL )
      {
	 assert(blk->mem == mem);
	 checkBlock(blk);
	 blk = blk->next;
      }
   }
}

#else

#define checkChunk(chk) /**/
#define checkBlock(blk) /**/
#define checkMem(mem) /**/

#endif


/* links chunk to the block's chunk array, sort it by store pointer
 * Returns:
 *    TRUE  if successful
 *    FALSE otherwise
 */
static int
linkChunk(BLKHDR* blk, CHKHDR* chk)
{
   CHKHDR* actchk;
   int left;
   int right;
   int middle;
   int i;

   assert(blk != NULL);
   assert(blk->numChunks <= blk->chunkarraySize);
   assert(chk != NULL);
   assert(chk->store != NULL);

   /* binary search for the position to insert the chunk */
   left = -1;
   right = blk->numChunks;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(0 <= middle && middle < blk->numChunks);
      assert(left < middle && middle < right);
      actchk = blk->chunkarray[middle];
      assert(actchk != NULL);
      if( chk->store < actchk->store )
         right = middle;
      else
      {
         assert(chk->store >= actchk->storeend);
         left = middle;
      }
   }
   assert(-1 <= left && left < blk->numChunks);
   assert(0 <= right && right <= blk->numChunks);
   assert(left+1 == right);
   assert(left == -1 || blk->chunkarray[left]->storeend <= chk->store);
   assert(right == blk->numChunks || chk->storeend <= blk->chunkarray[right]->store);

   /* ensure, that chunk array can store the additional chunk */
   if( blk->numChunks == blk->chunkarraySize )
   {
      blk->chunkarraySize = 2*(blk->numChunks+1);
      reallocMemoryArray(&blk->chunkarray, blk->chunkarraySize);
      if( blk->chunkarray == NULL )
         return FALSE;
   }
   assert(blk->numChunks < blk->chunkarraySize);
   assert(blk->chunkarray != NULL);

   /* move all chunks from 'right' to end one position to the right */
   for( i = blk->numChunks; i > right; --i )
   {
      blk->chunkarray[i] = blk->chunkarray[i-1];
      blk->chunkarray[i]->arraypos = i;
   }

   /* insert chunk at position 'right' */
   chk->arraypos = right;
   blk->chunkarray[right] = chk;
   blk->numChunks++;
   blk->storeSize += chk->storeSize;

   return TRUE;
}

/* unlinks chunk from the block's chunk list */
static void
unlinkChunk(CHKHDR* chk)
{
   BLKHDR* blk;
   int i;

   assert(chk != NULL);
   assert(chk->eagerFree == NULL);
   assert(chk->nextEager == NULL);
   assert(chk->prevEager == NULL);

   blk = chk->block;
   assert(blk != NULL);
   assert(blk->elemSize == chk->elemSize);
   assert(0 <= chk->arraypos && chk->arraypos < blk->numChunks);
   assert(blk->chunkarray[chk->arraypos] == chk);
   
   /* remove the chunk from the chunkarray of the block */
   for( i = chk->arraypos; i < blk->numChunks-1; ++i )
   {
      blk->chunkarray[i] = blk->chunkarray[i+1];
      blk->chunkarray[i]->arraypos = i;
   }
   blk->numChunks--;
   blk->storeSize -= chk->storeSize;
}

/* links chunk to the block's eager chunk list */
static void
linkEagerChunk(BLKHDR * blk, CHKHDR * chk)
{
   assert(chk->block == blk);
   assert(chk->nextEager == NULL);
   assert(chk->prevEager == NULL);

   chk->nextEager = blk->firstEager;
   chk->prevEager = NULL;
   if( blk->firstEager != NULL )
   {
      assert(blk->firstEager->prevEager == NULL);
      blk->firstEager->prevEager = chk;
   }
   blk->firstEager = chk;
}

/* unlinks chunk from the block's eager chunk list */
static void
unlinkEagerChunk(CHKHDR * chk)
{
   assert(chk != NULL);
   assert(chk->eagerFreeSize == 0 || chk->eagerFreeSize == chk->storeSize);

   if( chk->nextEager != NULL )
      chk->nextEager->prevEager = chk->prevEager;
   if( chk->prevEager != NULL )
      chk->prevEager->nextEager = chk->nextEager;
   else
   {
      assert(chk->block->firstEager == chk);
      chk->block->firstEager = chk->nextEager;
   }
   chk->nextEager = NULL;
   chk->prevEager = NULL;
   chk->eagerFree = NULL;
}

/* Allocates new memory if all chunks are full and
 * adds memory blocks to the lazy free list.
 * Parameters:
 *    blk : Pointer to headerstruct.
 * Returns:
 *    TRUE  if successful
 *    FALSE otherwise
 */
static int
createChunk(BLKHDR * blk)
{
   CHKHDR *newchunk;
   FREELIST *freeList;
   int i;
   int storeSize;
   int retval;

   assert(blk != NULL);
   assert(blk->mem != NULL);

   /* calculate store size */
   if( blk->numChunks == 0 )
      storeSize = blk->mem->initChunkSize;
   else
      storeSize = 2 * blk->lastChunkSize;
   assert(storeSize > 0);
   storeSize = MIN(storeSize, STORESIZE_MAX);
   storeSize = MIN(storeSize, CHUNKLENGTH_MAX / blk->elemSize);
   storeSize = MAX(storeSize, 1);
   blk->lastChunkSize = storeSize;

   /* create new chunk */
   assert(sizeof(CHKHDR) % ALIGNMENT == 0);
   allocMemorySize(&newchunk, sizeof(CHKHDR) + storeSize * blk->elemSize);
   if( newchunk == NULL )
      return FALSE;

   /* the store is allocated directly behind the chunk header */
   newchunk->elemSize = blk->elemSize;
   newchunk->storeSize = storeSize;
   newchunk->store = (void*) ((char*) newchunk + sizeof(CHKHDR));
   newchunk->storeend = (void*) ((char*) newchunk->store + storeSize * blk->elemSize);
   newchunk->arraypos = -1;
   newchunk->block = blk;

   debugMessage("allocated new chunk %p: %d elements with size %ld\n", 
      newchunk, newchunk->storeSize, (long)(newchunk->elemSize));

   /* add new memory to the lazy free list */
   for( i = 0; i < newchunk->storeSize - 1; ++i )
   {
      freeList = (FREELIST *) ((char *) (newchunk->store) + i * blk->elemSize); /*lint !e826*/
      freeList->next = (FREELIST *) ((char *) (newchunk->store) + (i + 1) * blk->elemSize); /*lint !e826*/
   }

   freeList = (FREELIST *) ((char *) (newchunk->store) + (newchunk->storeSize - 1) * blk->elemSize); /*lint !e826*/
   freeList->next = blk->lazyFree;
   blk->lazyFree = (FREELIST *) (newchunk->store);
   blk->lazyFreeSize += newchunk->storeSize;

   /* initialize eager free list of the chunk */
   newchunk->eagerFreeSize = 0;
   newchunk->eagerFree = NULL;

   /* initialize chunk links */
   newchunk->nextEager = NULL;
   newchunk->prevEager = NULL;

   /* link chunk into block */
   retval = linkChunk(blk, newchunk);

   checkBlock(blk);

   return retval;
}

/* destroy a chunk without updating the chunk lists
 * Parameters:
 *    chk : Pointer to chunkheader of chunk to delete.
 */
static void
destroyChunk(CHKHDR * chk)
{
   assert(chk != NULL);

   /* free chunk header and store (allocated in one call) */
   freeMemory(&chk);
}

/* remove a completely unused chunk, that is a chunk with all elements in the eager free list
 * Parameters:
 *    chk : Pointer to chunkheader of chunk to delete.
 */
static void
freeChunk(CHKHDR * chk)
{
   assert(chk != NULL);
   assert(chk->store != NULL);
   assert(chk->eagerFree != NULL);
   assert(chk->block != NULL);
   assert(chk->block->chunkarray != NULL);
   assert(chk->block->firstEager != NULL);
   assert(chk->eagerFreeSize == chk->storeSize);

   /* count the deleted eager free slots */
   chk->block->eagerFreeSize -= chk->eagerFreeSize;
   assert(chk->block->eagerFreeSize >= 0);

   /* remove chunk from eager chunk list */
   unlinkEagerChunk(chk);

   /* remove chunk from chunk list */
   unlinkChunk(chk);

   /* destroy the chunk */
   destroyChunk(chk);
}

/* returns an element of the eager free list and removes
 * it from the list
 */
static void *
allocChunkElement(CHKHDR * chk)
{
   FREELIST *ptr;

   assert(chk != NULL);
   assert(chk->eagerFree != NULL);
   assert(chk->eagerFreeSize > 0);
   assert(chk->block != NULL);

   /* unlink first element in the eager free list */
   ptr = chk->eagerFree;
   chk->eagerFree = ptr->next;
   chk->eagerFreeSize--;
   chk->block->eagerFreeSize--;

   assert((chk->eagerFreeSize == 0 && chk->eagerFree == NULL)
      ||  (chk->eagerFreeSize != 0 && chk->eagerFree != NULL));
   assert(chk->block->eagerFreeSize >= 0);

   /* unlink chunk from eager chunk list if necessary */
   if( chk->eagerFree == NULL )
   {
      assert(chk->eagerFreeSize == 0);
      unlinkEagerChunk(chk);
   }

   checkChunk(chk);
   return (void *) ptr;
}

/* puts ptr in the eager free list and adds the chunk
 * to the eager list of its block, if necessary
 */
static void
freeChunkElement(CHKHDR * chk, void *ptr)
{
   assert(chk != NULL);
   assert(chk->block != NULL);
   assert(isPtrInChunk(chk, ptr));

   /* link chunk to the eager chunk list if necessary */
   if( chk->eagerFree == NULL )
   {
      assert(chk->eagerFreeSize == 0);
      linkEagerChunk(chk->block, chk);
   }

   /* add ptr to the chunks eager free list */
   ((FREELIST *) ptr)->next = chk->eagerFree;
   chk->eagerFree = (FREELIST *) ptr;
   chk->eagerFreeSize++;
   chk->block->eagerFreeSize++;

   checkChunk(chk);
}

/*-----------------------------------------------------------------------------
 *--- Name: Create new block allocation pool in memory structure 'mem'.     ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    mem           : memory header structure.
 *    size          : Storage size per element.
 * Returns   :
 *    Pointer to block header structure.
 */
static BLKHDR *
createBlock(MEMHDR *mem, int size)
{
   BLKHDR *blk;

   assert(size >= (int)ALIGNMENT);

   allocMemory(&blk);

   if( blk != NULL )
   {
      blk->elemSize = size;
      blk->chunkarraySize = 0;
      blk->numChunks = 0;
      blk->lastChunkSize = 0;
      blk->storeSize = 0;
      blk->lazyFreeSize = 0;
      blk->eagerFreeSize = 0;
      blk->lazyFree = NULL;
      blk->chunkarray = NULL;
      blk->firstEager = NULL;
      blk->next = NULL;
      blk->mem = mem;
#ifndef NDEBUG
      blk->filename = NULL;
      blk->line = 0;
      blk->numGarbageCalls = 0;
      blk->numGarbageFrees = 0;
#endif
   }

   return blk;
}

/* destroy all chunks of the block, but keep block
 * header structure
 */
static void
clearBlock(BLKHDR * blk)
{
   int i;

   assert(blk != NULL);

   /* destroy all chunks of the block */
   for( i = 0; i < blk->numChunks; ++i )
      destroyChunk(blk->chunkarray[i]);

   blk->numChunks = 0;
   blk->lastChunkSize = 0;
   blk->storeSize = 0;
   blk->lazyFreeSize = 0;
   blk->eagerFreeSize = 0;
   blk->lazyFree = NULL;
   blk->firstEager = NULL;
}

/*-----------------------------------------------------------------------------
 *--- Name: Delete block allocation pool.                                   ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    blk : Pointer to block header of block to destroy.
 */
static void
destroyBlock(BLKHDR * blk)
{
   assert(blk != NULL);

   clearBlock(blk);

   freeMemoryArrayNull(&blk->chunkarray);

#ifndef NDEBUG
   freeMemoryArrayNull(&blk->filename);
#endif
   freeMemory(&blk);
}

/*-----------------------------------------------------------------------------
 *--- Name: Allocate a new memory element from the block.                   ---
 *-----------------------------------------------------------------------------
 * Returns:
 *    A pointer to a memory block of the block's elements size.
 */
static void *
allocBlockElement(BLKHDR * blk)
{
   FREELIST *ptr;

   assert(blk != NULL);

   /* if the lazy freelist is empty, we have to find the
      memory element somewhere else */
   if( blk->lazyFree == NULL )
   {
      assert(blk->lazyFreeSize == 0);

      /* check for a free element in the eager freelists */
      if( blk->firstEager != NULL )
	 return allocChunkElement(blk->firstEager);

      /* allocate a new chunk */
      if( createChunk(blk) == FALSE )
	 return NULL;
   }

   /* now the lazy freelist should contain an element */
   assert(blk->lazyFree != NULL);
   assert(blk->lazyFreeSize > 0);

   ptr = blk->lazyFree;
   blk->lazyFree = ptr->next;
   blk->lazyFreeSize--;

   checkBlock(blk);
   return (void *) ptr;
}

/*-----------------------------------------------------------------------------
 *--- Name: Free a block of memory and return it to the lazy freelist.      ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    blk : Pointer to the headerstruct.
 *    ptr : Pointer to the block  of memory to free.
 */
/*ARGSUSED*/
static void
freeBlockElement(BLKHDR * blk, void *ptr, const char *filename, int line)
{
   assert(blk != NULL);
   assert(ptr != NULL);

#ifdef DEBUG
   /* check, if ptr belongs to block */
   if( !isPtrInBlock(blk, ptr) )
   {
      BLKHDR *correctblk;

      fprintf(stderr, "[%s:%d] ERROR: pointer %p does not belong to block %p (size: %ld)\n", 
         filename, line, ptr, blk, (long) (blk->elemSize));

      correctblk = findBlock(blk->mem, ptr);
      if( correctblk == NULL )
	 fprintf(stderr, "[%s:%d] ERROR: -> pointer %p does not belong to memory structure %p\n", 
            filename, line, ptr, blk->mem);
      else
	 fprintf(stderr, "[%s:%d] ERROR: -> pointer %p belongs to block %p instead (size: %ld)\n", 
            filename, line, ptr, correctblk, (long) (correctblk->elemSize));
   }
#endif

   /* put ptr in lazy free list */
   ((FREELIST *) ptr)->next = blk->lazyFree;
   blk->lazyFree = (FREELIST *) ptr;
   blk->lazyFreeSize++;

   checkBlock(blk);
}

/*-----------------------------------------------------------------------------
 *--- Name: Get Hash number of memory size                                  ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    size : size of memory block.
 * Returns   :
 *    Hash number of size.
 */
static int
getHashNumber(int size)
{
   assert(size % (int)ALIGNMENT == 0);
   return ((size / (int)ALIGNMENT) % BLOCKHASH_SIZE);
}

/*-----------------------------------------------------------------------------
 *--- Name: Create a block memory allocation structure.                     ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    initChunkSize     : nr. of elements in the first chunk.
 *    clearUnusedBlocks : immedeately clear a block, if it is unused
 *    garbageFactor     : if at least garbageFactor * avg. chunksize elements
 *                        are free, call garbage collection,
 *                        a value of -1 disables garbage collection
 * Returns   :
 *    Pointer to memory header structure.
 */
MEMHDR *
createBlockMemory_call(int initChunkSize, int clearUnusedBlocks,
   int garbageFactor, const char *filename, int line)
{
   MEMHDR *mem;
   int     i;

   allocMemory(&mem);
   if( mem != NULL )
   {
      for( i = 0; i < BLOCKHASH_SIZE; ++i )
	 mem->blockhash[i] = NULL;
      mem->initChunkSize = initChunkSize;
      mem->clearUnusedBlocks = clearUnusedBlocks;
      mem->garbageFactor = garbageFactor;
      mem->memused = 0;
   }
   else
      fprintf(stderr, "[%s:%d] ERROR: Insufficient memory for memory header\n", filename, line);

   return mem;
}

/*-----------------------------------------------------------------------------
 *--- Name: frees all chunks in the memory structure.                       ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    mem : Pointer to memory header to clear.
 */
void
clearBlockMemory_call(MEMHDR *mem, const char *filename, int line)
{
   BLKHDR *blk;
   BLKHDR *nextblk;
   int     i;

   if( mem != NULL )
   {
      for( i = 0; i < BLOCKHASH_SIZE; ++i )
      {
	 blk = mem->blockhash[i];
	 while( blk != NULL )
	 {
	    nextblk = blk->next;
	    destroyBlock(blk);
	    blk = nextblk;
	 }
	 mem->blockhash[i] = NULL;
      }
      mem->memused = 0;
   }
   else
      fprintf(stderr, "[%s:%d] ERROR: Tried to clear null block\n", filename, line);
}

/*-----------------------------------------------------------------------------
 *--- Name: Delete a block allocation.                                      ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    mem : Pointer to memory header to destroy.
 */
void
destroyBlockMemory_call(MEMHDR **mem, const char *filename, int line)
{
   assert(mem != NULL);

   if( *mem != NULL )
   {
      clearBlockMemory_call(*mem, filename, line);
      freeMemory(mem);
      assert(*mem == NULL);
   }
   else
      fprintf(stderr, "[%s:%d] ERROR: Tried to destroy null block\n", filename, line);
}

/*-----------------------------------------------------------------------------
 *--- Name: Get a new block of memory.                                      ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    mem  : Pointer to memory header.
 *    size : size of requested block in bytes.
 * Returns:
 *    Pointer to a new block of memory of size "size".
 */
void *
allocBlockMemory_call(MEMHDR *mem, size_t size, const char *filename, int line)
{
   BLKHDR **blkptr;
   int     hashNumber;
   void   *ptr;

   /* calculate hash number of given size */
   alignSize(&size);
   hashNumber = getHashNumber((int)size);

   /* find correspoding block header */
   blkptr = &(mem->blockhash[hashNumber]);
   while( *blkptr != NULL && (*blkptr)->elemSize != (int)size )
      blkptr = &((*blkptr)->next);

   /* create new block header if necessary */
   if( *blkptr == NULL )
   {
      *blkptr = createBlock(mem, (int)size);
      if( *blkptr == NULL )
      {
	 fprintf(stderr, "[%s:%d] ERROR: Insufficient memory for block header\n", filename, line);
	 return NULL;
      }
#ifndef NDEBUG
      duplicateMemoryArray(&(*blkptr)->filename, filename, strlen(filename) + 1);
      (*blkptr)->line = line;
#endif
   }

   /* get memory inside the block */
   ptr = allocBlockElement(*blkptr);
   if( ptr == NULL )
      fprintf(stderr, "[%s:%d] ERROR: Insufficient memory for new chunk\n", filename, line);
   debugMessage("[%s:%4d] alloced %8ld bytes in %p\n", filename, line, (long)size, ptr);

   mem->memused += size;

   checkMem(mem);

   return ptr;
}

void *
reallocBlockMemory_call(MEMHDR *mem, void* ptr, size_t oldsize, size_t newsize, const char *filename, int line)
{
   void* newptr;

   if( ptr == NULL )
   {
      assert(oldsize == 0);
      return allocBlockMemory_call(mem, newsize, filename, line);
   }

   assert(ptr != NULL);

   alignSize(&oldsize);
   alignSize(&newsize);
   if( oldsize == newsize )
      return ptr;

   newptr = allocBlockMemory_call(mem, newsize, filename, line);
   if( newptr != NULL )
      copyMemorySize(newptr, ptr, MIN(oldsize, newsize));
   freeBlockMemory_call(mem, &ptr, oldsize, filename, line);

   return newptr;
}

void *
duplicateBlockMemory_call(MEMHDR *mem, const void* source, size_t size, const char *filename, int line)
{
   void   *ptr = NULL;

   assert(source != NULL);
   ptr = allocBlockMemory_call(mem, size, filename, line);
   if( ptr != NULL )
      copyMemorySize(ptr, source, size);

   return ptr;
}

/* sort the lazy free list into the eager free lists,
 * and remove completely unused chunks
 */
static void
garbageCollection(BLKHDR * blk)
{
   CHKHDR *chk;
   CHKHDR *nextEager;
   FREELIST *lazyFree;

   assert(blk != NULL);

#ifndef NDEBUG
   blk->numGarbageCalls++;
#endif

   /* put the lazy free elements in the eager free lists */
   while( blk->lazyFree != NULL )
   {
      /* unlink first element from the lazy free list */
      lazyFree = blk->lazyFree;
      blk->lazyFree = blk->lazyFree->next;
      blk->lazyFreeSize--;

      /* identify the chunk of the element */
      chk = findChunk(blk, (void *) lazyFree);
#ifndef NDEBUG
      if( chk == NULL )
         fprintf(stderr, "[%s:%d] ERROR: chunk for lazy free block %p not found in block %p\n",
            __FILE__, __LINE__, lazyFree, blk );
#endif
      assert(chk != NULL);

      /* add the element to the chunk's eager free list */
      freeChunkElement(chk, (void *) lazyFree);
      assert(chk->eagerFreeSize > 0);
   }
   assert(blk->lazyFreeSize == 0);

   /* delete completely unused chunks, but keep at least one */
   chk = blk->firstEager;
   while( chk != NULL && blk->numChunks > 1 )
   {
      nextEager = chk->nextEager;
      if( chk->eagerFreeSize == chk->storeSize )
      {
#ifndef NDEBUG
	 blk->numGarbageFrees++;
#endif
	 freeChunk(chk);
      }
      chk = nextEager;
   }

   checkBlock(blk);
}

/*-----------------------------------------------------------------------------
 *--- Name: Free a block of memory.                                         ---
 *-----------------------------------------------------------------------------
 * Parameters:
 *    mem  : Pointer to memory header.
 *    ptr  : Pointer to memory to free.
 *    size : size of memory block.
 */
void
freeBlockMemory_call(MEMHDR *mem, void **ptr, size_t size,
   const char *filename, int line)
{
   BLKHDR *blk;
   int     hashNumber;

   assert(ptr != NULL);

   if( *ptr != NULL )
   {
      /* calculate hash number of given size */
      alignSize(&size);
      hashNumber = getHashNumber((int)size);

      debugMessage("[%s:%4d] free    %8ld bytes in %p\n", filename, line, (long)size, *ptr );
      /* find correspoding block header */
      blk = mem->blockhash[hashNumber];
      while( blk != NULL && blk->elemSize != (int)size )
	 blk = blk->next;
      if( blk == NULL )
      {
	 fprintf(stderr, "[%s:%d] ERROR: Tried to free pointer <%p> in block memory <%p> of unknown size %ld\n",
            filename, line, *ptr, mem, (long) size);
	 return;
      }

      /* free memory in block */
      freeBlockElement(blk, *ptr, filename, line);
      *ptr = NULL;

      /* check, if we want to clear the block */
      if( mem->clearUnusedBlocks
	 && blk->lazyFreeSize + blk->eagerFreeSize == blk->storeSize )
	 clearBlock(blk);

      /* check, if we want to do garbage collection */
      if( mem->garbageFactor >= 0
         && blk->numChunks > 0
	 && blk->lazyFreeSize >= GARBAGE_SIZE
	 && blk->lazyFreeSize + blk->eagerFreeSize >
         mem->garbageFactor * (double)(blk->storeSize) / (double)(blk->numChunks) )
      {
	 garbageCollection(blk);
      }

      mem->memused -= size;
      assert(mem->memused >= 0);
   }
   else
      fprintf(stderr, "[%s:%d] ERROR: Tried to free null block pointer\n", filename, line);

   checkMem(mem);
}

long long
getBlockMemoryUsed(MEMHDR *mem)
{
   assert(mem != NULL);

   return mem->memused;
}

#ifndef NDEBUG
size_t
blockMemorySize(MEMHDR *mem, void *ptr)
{
   BLKHDR *blk;

   assert(mem != NULL);
   if( ptr == NULL )
      return 0;

   blk = findBlock(mem, ptr);
   if( blk == NULL )
      return 0;

   return (size_t)(blk->elemSize);
}

void
blockMemoryDiagnostic(MEMHDR *mem)
{
   BLKHDR *blk;
   int     numBlocks = 0;
   int     numUnusedBlocks = 0;
   int     numGarbageCalls = 0;
   int     numGarbageFrees = 0;
   long    allocedMem = 0;
   long    freeMem = 0;
   int     i;
   int     c;

   printf(" ElSize #Chk #Eag  #Elems  #EagFr  #LazFr  #GCl #GFr  Free  First Allocator\n");

   assert(mem != NULL);

   for( i = 0; i < BLOCKHASH_SIZE; ++i )
   {
      blk = mem->blockhash[i];
      while( blk != NULL )
      {
	 CHKHDR *chk;
	 int     numChunks = 0;
	 int     numElems = 0;
	 int     numEagerChunks = 0;
	 int     numEagerElems = 0;

         for( c = 0; c < blk->numChunks; ++c )
         {
            chk = blk->chunkarray[c];
            assert(chk != NULL);
	    assert(chk->elemSize == blk->elemSize);
	    assert(chk->block == blk);
	    numChunks++;
	    numElems += chk->storeSize;
	    if( chk->eagerFree != NULL )
	    {
	       numEagerChunks++;
	       numEagerElems += chk->eagerFreeSize;
	    }
	 }

	 assert(numChunks == blk->numChunks);
	 assert(numElems == blk->storeSize);
	 assert(numEagerElems == blk->eagerFreeSize);

	 if( numElems > 0 )
	 {
	    numBlocks++;
	    allocedMem += blk->elemSize * numElems;
	    freeMem += blk->elemSize * (numEagerElems + blk->lazyFreeSize);

	    printf("%7ld %4d %4d %7d %7d %7d %5d %4d %5.1f%% %s:%d\n",
	       (long) (blk->elemSize), numChunks, numEagerChunks, numElems,
	       numEagerElems, blk->lazyFreeSize, blk->numGarbageCalls,
	       blk->numGarbageFrees,
	       100.0 * (double) (numEagerElems +
		  blk->lazyFreeSize) / (double) (numElems), blk->filename,
	       blk->line);
	 }
	 else
	 {
	    printf("%7ld <unused>                          %5d %4d        %s:%d\n",
	       (long) (blk->elemSize), blk->numGarbageCalls,
	       blk->numGarbageFrees, blk->filename, blk->line);
	    numUnusedBlocks++;
	 }
	 numGarbageCalls += blk->numGarbageCalls;
	 numGarbageFrees += blk->numGarbageFrees;
	 blk = blk->next;
      }
   }

   printf("Garbage collector called %d times (freed %d chunks).\n",
      numGarbageCalls, numGarbageFrees);
   printf("%d blocks (%d unused), %ld bytes allocated, %ld bytes free",
      numBlocks + numUnusedBlocks, numUnusedBlocks, allocedMem, freeMem);
   if( allocedMem > 0 )
      printf(" (%.1f%%)", 100.0 * (double) freeMem / (double) allocedMem);
   printf("\n");
}

void
blockMemoryCheckEmpty(MEMHDR *mem)
{
   BLKHDR *blk;
   long    allocedMem = 0;
   long    freeMem = 0;
   int     i;
   int     c;

   assert(mem != NULL);

   for( i = 0; i < BLOCKHASH_SIZE; ++i )
   {
      blk = mem->blockhash[i];
      while( blk != NULL )
      {
	 CHKHDR *chk;
	 int     numChunks = 0;
	 int     numElems = 0;
	 int     numEagerElems = 0;

         for( c = 0; c < blk->numChunks; ++c )
         {
            chk = blk->chunkarray[c];
            assert(chk != NULL);
	    assert(chk->elemSize == blk->elemSize);
	    assert(chk->block == blk);
	    numChunks++;
	    numElems += chk->storeSize;
	    if( chk->eagerFree != NULL )
	       numEagerElems += chk->eagerFreeSize;
	 }

	 assert(numChunks == blk->numChunks);
	 assert(numElems == blk->storeSize);
	 assert(numEagerElems == blk->eagerFreeSize);

	 if( numElems > 0 )
	 {
	    allocedMem += blk->elemSize * numElems;
	    freeMem += blk->elemSize * (numEagerElems + blk->lazyFreeSize);

            if( numElems != numEagerElems + blk->lazyFreeSize )
               printf("%ld bytes (%d elements of size %ld) not freed. First Allocator: %s:%d\n",
                  ((numElems - numEagerElems) - blk->lazyFreeSize) * (long) (blk->elemSize),
                  (numElems - numEagerElems) - blk->lazyFreeSize, (long) (blk->elemSize),
                  blk->filename, blk->line);
	 }
	 blk = blk->next;
      }
   }

   if( allocedMem != freeMem )
      printf("%ld bytes not freed in total.\n", allocedMem - freeMem);
}

#endif




#endif
