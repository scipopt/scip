/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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

#include "memory.h"

/* check memory structures after every call */
#define DEBUG

#define BLOCKHASH_SIZE     1013	/* should be prime */
#define CHUNKLENGTH_MAX 1048576	/* maximal size of a chunk (in bytes) */
#define STORESIZE_MAX      8192	/* maximal number of elements in one chunk */
#define ALIGNMENT       (sizeof(FREELIST))
#define GARBAGE_SIZE        256	/* size of lazy free list to start
				   garbage collection */

void
errorMessage_call(const char *msg, const char *filename, int line)
{
   printf("[%s:%d] %s\n", filename, line, msg);
}

/******************************
 * Standard Memory Management *
 ******************************/

#ifndef NDEBUG

#include <string.h>

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
   MEMLIST *list = malloc(sizeof(MEMLIST));

   assert(list != NULL);
   assert(ptr != NULL && size > 0);
   list->ptr = ptr;
   list->size = size;
   list->filename = (char *) malloc((strlen(filename) + 1) * sizeof(char));
   strcpy(list->filename, filename);
   list->line = line;
   list->next = memlist;
   memlist = list;
   memused += size;
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
      memused -= list->size;
      free(list->filename);
      free(list);
   }
   else
   {
      char    s[255];

      sprintf(s, "Error! Tried to free unknown pointer <%p>.\n", ptr);
      errorMessage_call(s, filename, line);
   }
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
      used += list->size;
      list = list->next;
   }
   printf("Total:    %8ld\n", memused);
   if( used != memused )
   {
      char    s[255];

      sprintf(s, "Error! Used memory in list sums up to %ld instead of %ld\n",
	 used, memused);
      errorMessage(s);
   }
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

   /* printf( "malloc %ld bytes [%s:%d]\n", (long)size, filename, line );*/ /* ??? */
   size = MAX(size, 1);
   ptr = malloc(size);

   if( ptr == NULL )
   {
      char    s[255];

      sprintf(s, "Error! Insufficient memory for allocation of %ld bytes.\n",
	 (long) size);
      errorMessage_call(s, filename, line);
      /* abort(); */
   }
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
   {
      char    s[255];

      sprintf(s, "Error! Insufficient memory for reallocation of %ld bytes.\n", (long) size);
      errorMessage_call(s, filename, line);
      /* abort(); */
   }
#ifndef NDEBUG
   else
      memListAdd(newptr, size, filename, line);
#endif
   return newptr;
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
      errorMessage_call("Warning! Try to free null pointer.\n", filename, line);
}


/***************************
 * Block Memory Management *
 ***************************/

typedef struct free_list FREELIST;
typedef struct chunk_header CHKHDR;
typedef struct block_header BLKHDR;

struct free_list
{
   FREELIST *next;
};

struct chunk_header
{
   size_t  elem_size;
   int     storeSize;
   int     eagerFreeSize;
   void   *store;
   FREELIST *eagerFree;
   CHKHDR *next;
   CHKHDR *prev;
   CHKHDR *nextEager;
   CHKHDR *prevEager;
   BLKHDR *block;
};				/* the chunk header must be aligned! */

struct block_header
{
   size_t  elem_size;
   int     numChunks;
   int     storeSize;
   int     lazyFreeSize;
   int     eagerFreeSize;
   FREELIST *lazyFree;
   CHKHDR *chunklist;
   CHKHDR *firstEager;
   BLKHDR *next;
   MEMHDR *mem;
#ifndef NDEBUG
   char   *filename;
   int     line;
   int     numGarbageCalls;
   int     numGarbageFrees;
#endif
};

struct memory_header
{
   int     initChunkSize;
   int     clearUnusedBlocks;
   int     garbageFactor;
   BLKHDR *blockhash[BLOCKHASH_SIZE];
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
   return (ptr >= (void *) (chk->store)
      && ptr <
      (void *) ((char *) (chk->store) + chk->storeSize * chk->elem_size));
}

/* Given a pointer find the chunk this pointer points to
 * in the chunklist. The given chunk is searched first,
 * and from there, the chunk list is searched in both
 * directions to the start and to the end. */
static CHKHDR *
findChunk(CHKHDR * chk, const void *ptr)
{
   CHKHDR *leftchk;
   CHKHDR *rightchk;

   assert(ptr != NULL);

   if( chk == NULL )
      return NULL;

   leftchk = chk;
   rightchk = chk->next;

   /* search in both directions */
   while( leftchk != NULL && rightchk != NULL )
   {
      if( isPtrInChunk(leftchk, ptr) )
	 return leftchk;
      if( isPtrInChunk(rightchk, ptr) )
	 return rightchk;
      leftchk = leftchk->prev;
      rightchk = rightchk->next;
   }

   /* search in left direction */
   while( leftchk != NULL )
   {
      if( isPtrInChunk(leftchk, ptr) )
	 return leftchk;
      leftchk = leftchk->prev;
   }

   /* search in right direction */
   while( rightchk != NULL )
   {
      if( isPtrInChunk(rightchk, ptr) )
	 return rightchk;
      rightchk = rightchk->next;
   }

   /* search failed */
   return NULL;
}

/* checks, if 'ptr' belongs to a chunk of block 'blk' */
static int
isPtrInBlock(BLKHDR * blk, const void *ptr)
{
   assert(blk != NULL);
   return (findChunk(blk->chunklist, ptr) != NULL);
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
   FREELIST *eager;
   int     eagerFreeSize = 0;

   assert(chk != NULL);
   assert(chk->store != NULL);
   assert(chk->block != NULL);
   assert(chk->block->elem_size == chk->elem_size);

   if( chk->eagerFree == NULL )
      assert(chk->nextEager == NULL && chk->prevEager == NULL);
   else if( chk->prevEager == NULL )
      assert(chk->block->firstEager == chk);

   if( chk->next != NULL )
      assert(chk->next->prev == chk);
   if( chk->prev != NULL )
      assert(chk->prev->next == chk);
   if( chk->nextEager != NULL )
      assert(chk->nextEager->prevEager == chk);
   if( chk->prevEager != NULL )
      assert(chk->prevEager->nextEager == chk);

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
   CHKHDR *chk;
   FREELIST *lazy;
   int     numChunks = 0;
   int     storeSize = 0;
   int     lazyFreeSize = 0;
   int     eagerFreeSize = 0;

   assert(blk != NULL);

   chk = blk->chunklist;
   assert(chk == NULL || chk->prev == NULL);
   while( chk != NULL )
   {
      checkChunk(chk);
      numChunks++;
      storeSize += chk->storeSize;
      eagerFreeSize += chk->eagerFreeSize;
      chk = chk->next;
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
      chk = findChunk(blk->chunklist, lazy);
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
   BLKHDR *blk;
   int     i;

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
/* links chunk to the block's chunk list */
   static void
linkChunk(BLKHDR * blk, CHKHDR * chk)
{
   assert(chk->next == NULL);
   assert(chk->prev == NULL);

   chk->block = blk;
   chk->next = blk->chunklist;
   chk->prev = NULL;
   if( blk->chunklist != NULL )
   {
      assert(blk->chunklist->prev == NULL);
      blk->chunklist->prev = chk;
   }
   blk->chunklist = chk;
   blk->numChunks++;
   blk->storeSize += chk->storeSize;
}

/* unlinks chunk from the block's chunk list */
static void
unlinkChunk(CHKHDR * chk)
{
   if( chk->next != NULL )
      chk->next->prev = chk->prev;
   if( chk->prev != NULL )
      chk->prev->next = chk->next;
   else
   {
      assert(chk->block->chunklist == chk);
      chk->block->chunklist = chk->next;
   }
   chk->block->numChunks--;
   chk->block->storeSize -= chk->storeSize;
   chk->next = NULL;
   chk->prev = NULL;
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
   int     i;
   int     storeSize;

   assert(blk != NULL);
   assert(blk->mem != NULL);

   /* calculate store size */
   if( blk->chunklist == NULL )
      storeSize = blk->mem->initChunkSize;
   else
   {
      storeSize = MAX(1, CHUNKLENGTH_MAX / blk->elem_size);
      storeSize = MIN(storeSize, STORESIZE_MAX);
      storeSize = MIN(storeSize, 2 * blk->chunklist->storeSize);
   }

   /* create new chunk */
   assert(sizeof(CHKHDR) % ALIGNMENT == 0);
   allocMemorySize(newchunk, sizeof(CHKHDR) + storeSize * blk->elem_size);
   if( newchunk == NULL )
      return FALSE;

   /* the store is allocated directly behind the chunk header */
   newchunk->elem_size = blk->elem_size;
   newchunk->storeSize = storeSize;
   newchunk->store = (void *) ((char *) newchunk + sizeof(CHKHDR));

   /* printf( "allocated new chunk: %d elements with size %ld\n", 
      newchunk->storeSize, (long)(newchunk->elem_size)); ??? */

   /* add new memory to the lazy free list */
   for( i = 0; i < newchunk->storeSize - 1; ++i )
   {
      freeList =
	 (FREELIST *) ((char *) (newchunk->store) + i * blk->elem_size);
      freeList->next =
	 (FREELIST *) ((char *) (newchunk->store) + (i + 1) * blk->elem_size);
   }

   freeList =
      (FREELIST *) ((char *) (newchunk->store) + (newchunk->storeSize -
	 1) * blk->elem_size);
   freeList->next = blk->lazyFree;
   blk->lazyFree = (FREELIST *) (newchunk->store);
   blk->lazyFreeSize += newchunk->storeSize;

   /* initialize eager free list of the chunk */
   newchunk->eagerFreeSize = 0;
   newchunk->eagerFree = NULL;

   /* initialize chunk links */
   newchunk->next = NULL;
   newchunk->prev = NULL;
   newchunk->nextEager = NULL;
   newchunk->prevEager = NULL;

   /* link chunk into block */
   linkChunk(blk, newchunk);

   checkBlock(blk);
   return TRUE;
}

/* destroy a chunk without updating the chunk lists
 * Parameters:
 *    chk : Pointer to chunkheader of chunk to delete.
 */
static void
destroyChunk(CHKHDR * chk)
{
   assert(chk != NULL);

   /* printf( "destroy chunk: %d elements with size %ld\n",
      chk->storeSize, (long)(chk->block->elem_size)); ??? */

   /* free chunk header and store (allocated in one call) */
   freeMemory(chk);
}

/* remove a complely unused chunk
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
   assert(chk->block->chunklist != NULL);
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
createBlock(MEMHDR *mem, size_t size)
{
   BLKHDR *blk;

   assert(size >= ALIGNMENT);

   allocMemory(blk);

   if( blk != NULL )
   {
      blk->elem_size = size;
      blk->numChunks = 0;
      blk->storeSize = 0;
      blk->lazyFreeSize = 0;
      blk->eagerFreeSize = 0;
      blk->lazyFree = NULL;
      blk->chunklist = NULL;
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
   CHKHDR *chk;
   CHKHDR *nextchk;

   assert(blk != NULL);

   /* destroy all chunks of the block */
   chk = blk->chunklist;
   while( chk != NULL )
   {
      nextchk = chk->next;
      destroyChunk(chk);
      chk = nextchk;
   }
   blk->numChunks = 0;
   blk->storeSize = 0;
   blk->lazyFreeSize = 0;
   blk->eagerFreeSize = 0;
   blk->lazyFree = NULL;
   blk->chunklist = NULL;
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
   CHKHDR *chk;
   CHKHDR *nextchk;

   assert(blk != NULL);

   /* destroy all chunks of the block */
   chk = blk->chunklist;
   while( chk != NULL )
   {
      nextchk = chk->next;
      destroyChunk(chk);
      chk = nextchk;
   }
#ifndef NDEBUG
   freeMemoryArrayNull(blk->filename);
#endif
   freeMemory(blk);
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
   if( findChunk(blk->chunklist, ptr) == NULL )
   {
      BLKHDR *correctblk;
      char    s[255];

      sprintf(s, "pointer %p does not belong to block %p (size: %ld)", ptr,
	 blk, (long) (blk->elem_size));
      errorMessage_call(s, filename, line);

      correctblk = findBlock(blk->mem, ptr);
      if( correctblk == NULL )
      {
	 sprintf(s, "pointer %p does not belong to memory structure %p", ptr,
	    blk->mem);
	 errorMessage_call(s, filename, line);
      }
      else
      {
	 sprintf(s, "pointer %p belongs to block %p instead (size: %ld)", ptr,
	    correctblk, (long) (correctblk->elem_size));
	 errorMessage_call(s, filename, line);
      }
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
getHashNumber(size_t size)
{
   assert(size % ALIGNMENT == 0);
   return ((size / ALIGNMENT) % BLOCKHASH_SIZE);
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

   allocMemory(mem);
   if( mem != NULL )
   {
      mem->initChunkSize = initChunkSize;
      mem->clearUnusedBlocks = clearUnusedBlocks;
      mem->garbageFactor = garbageFactor;
      for( i = 0; i < BLOCKHASH_SIZE; ++i )
	 mem->blockhash[i] = NULL;
   }
   else
      errorMessage_call("Error! Insufficient memory for memory header.",
	 filename, line);

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
   }
   else
      errorMessage_call("Error! Try to clear null block.", filename, line);
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
      freeMemory(*mem);
      *mem = NULL;
   }
   else
      errorMessage_call("Error! Try to destroy null block.", filename, line);
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
   hashNumber = getHashNumber(size);

   /* find correspoding block header */
   blkptr = &(mem->blockhash[hashNumber]);
   while( *blkptr != NULL && (*blkptr)->elem_size != size )
      blkptr = &((*blkptr)->next);

   /* create new block header if necessary */
   if( *blkptr == NULL )
   {
      *blkptr = createBlock(mem, size);
      if( *blkptr == NULL )
      {
	 errorMessage_call("Error! Insufficient memory for block header.", filename, line);
	 return NULL;
      }
#ifndef NDEBUG
      allocMemoryArray((*blkptr)->filename, strlen(filename) + 1);
      assert((*blkptr)->filename != NULL);
      strcpy((*blkptr)->filename, filename);
      (*blkptr)->line = line;
#endif
   }

   /* get memory inside the block */
   ptr = allocBlockElement(*blkptr);
   if( ptr == NULL )
      errorMessage_call("Error! Insufficient memory for new chunk.", filename, line);
   /* printf( "[%s:%4d] alloced %8ld bytes in %p\n", filename, line, (long)size, ptr ); ??? */

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
      copyMemory_call(newptr, ptr, MIN(oldsize, newsize));
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
      copyMemory_call(ptr, source, size);

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
   chk = blk->chunklist;
   while( blk->lazyFree != NULL )
   {
      /* unlink first element from the lazy free list */
      lazyFree = blk->lazyFree;
      blk->lazyFree = blk->lazyFree->next;
      blk->lazyFreeSize--;

      /* identify the chunk of the element */
      chk = findChunk(chk, (void *) lazyFree);
#ifdef DEBUG
      if( chk == NULL )
      {
         char s[255];
         sprintf( s, "chunk for lazy free block %p not found in block %p",
            lazyFree, blk );
         errorMessage( s );
      }
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
      hashNumber = getHashNumber(size);

      /* printf( "[%s:%4d] free    %8ld bytes in %p\n", filename, line, (long)size, *ptr ); ??? */
      /* find correspoding block header */
      blk = mem->blockhash[hashNumber];
      while( blk != NULL && blk->elem_size != size )
	 blk = blk->next;
      if( blk == NULL )
      {
	 char    s[255];

	 sprintf(s, "Error! Tried to free pointer <%p> of unknown size %ld.\n",
	    *ptr, (long) size);
	 errorMessage_call(s, filename, line);
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
      if( mem->garbageFactor >= 0 && blk->numChunks > 0
	 && blk->lazyFreeSize >= GARBAGE_SIZE
	 && blk->lazyFreeSize + blk->eagerFreeSize >
	 mem->garbageFactor * (double)(blk->storeSize) / (double)(blk->numChunks) )
	 garbageCollection(blk);
   }
   else
      errorMessage_call("Warning! Try to free null block pointer.\n", filename,
	 line);

   checkMem(mem);
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

   return blk->elem_size;
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

   printf
      (" ElSize #Chk #Eag  #Elems  #EagFr  #LazFr  #GCl #GFr  Free  First Allocator\n");

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

	 chk = blk->chunklist;
	 while( chk != NULL )
	 {
	    assert(chk->elem_size == blk->elem_size);
	    assert(chk->block == blk);
	    numChunks++;
	    numElems += chk->storeSize;
	    if( chk->eagerFree != NULL )
	    {
	       numEagerChunks++;
	       numEagerElems += chk->eagerFreeSize;
	    }
	    chk = chk->next;
	 }

	 assert(numChunks == blk->numChunks);
	 assert(numElems == blk->storeSize);
	 assert(numEagerElems == blk->eagerFreeSize);

	 if( numElems > 0 )
	 {
	    numBlocks++;
	    allocedMem += blk->elem_size * numElems;
	    freeMem += blk->elem_size * (numEagerElems + blk->lazyFreeSize);

	    printf("%7ld %4d %4d %7d %7d %7d %5d %4d %5.1f%% %s:%d\n",
	       (long) (blk->elem_size), numChunks, numEagerChunks, numElems,
	       numEagerElems, blk->lazyFreeSize, blk->numGarbageCalls,
	       blk->numGarbageFrees,
	       100.0 * (double) (numEagerElems +
		  blk->lazyFreeSize) / (double) (numElems), blk->filename,
	       blk->line);
	 }
	 else
	 {
	    printf
	       ("%7ld <unused>                          %5d %4d        %s:%d\n",
	       (long) (blk->elem_size), blk->numGarbageCalls,
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
#endif

