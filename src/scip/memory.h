/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: memory.h,v 1.25 2005/04/12 08:48:19 bzfpfend Exp $"

/**@file   memory.h
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEMORY_H__
#define __MEMORY_H__


#include <stdlib.h>
#include <assert.h>



/*************************************************************************************
 * Standard Memory Management
 *
 * In debug mode, these methods extend malloc() and free() by logging all currently
 * allocated memory elements in an allocation list. This can be used as a simple leak
 * detection.
 *************************************************************************************/

#define allocMemory(ptr)                   (*(void**)(ptr) = allocMemory_call( sizeof(**(ptr)), __FILE__, __LINE__ ))
#define allocMemoryArray(ptr,num)          (*(void**)(ptr) = allocMemory_call( (num)*sizeof(**(ptr)), __FILE__, __LINE__ ))
#define allocMemorySize(ptr,size)          (*(void**)(ptr) = allocMemory_call( (size_t)(size), __FILE__, __LINE__ ))
#define allocMemoryCPP(size)               allocMemory_call( (size_t)(size), __FILE__, __LINE__ )
#define allocMemoryArrayCPP(num,size)      allocMemory_call( (size_t)((num)*(size)), __FILE__, __LINE__ )
#define reallocMemoryArray(ptr,num)        (*(void**)(ptr) = reallocMemory_call( *(ptr), (num)*sizeof(**(ptr)), \
                                           __FILE__, __LINE__ ))
#define reallocMemorySize(ptr,size)        (*(void**)(ptr) = reallocMemory_call( *(ptr), (size_t)(size), \
                                           __FILE__, __LINE__ ))
#define clearMemory(ptr)                   clearMemory_call( (void*)(ptr), sizeof(*(ptr)) )
#define clearMemoryArray(ptr, num)         clearMemory_call( (void*)(ptr), (num)*sizeof(*(ptr)) )
#define clearMemorySize(ptr, size)         clearMemory_call( (void*)(ptr), (size_t)(size) )
#define copyMemory(ptr, source)            copyMemory_call( (void*)(ptr), (const void*)(source), sizeof(*(ptr)) )
#define copyMemoryArray(ptr, source, num)  copyMemory_call( (void*)(ptr), (const void*)(source), (num)*sizeof(*(ptr)) )
#define copyMemorySize(ptr, source, size)  copyMemory_call( (void*)(ptr), (const void*)(source), (size_t)(size) )
#define duplicateMemory(ptr, source)       (*(void**)(ptr) = duplicateMemory_call( (const void*)(source), sizeof(**(ptr)), \
                                           __FILE__, __LINE__ ))
#define duplicateMemoryArray(ptr, source, num) \
                                           (*(void**)(ptr) = duplicateMemory_call( (const void*)(source), \
                                           (num)*sizeof(**(ptr)), __FILE__, __LINE__ ))
#define duplicateMemorySize(ptr, source, size) \
                                           (*(void**)(ptr) = duplicateMemory_call( (const void*)(source), (size_t)(size), \
                                           __FILE__, __LINE__ ))
#define freeMemory(ptr)                    { freeMemory_call( (void**)(ptr), __FILE__, __LINE__ ); \
                                             assert(*(ptr) == NULL); }
#define freeMemoryNull(ptr)                { if( *(ptr) != NULL ) freeMemory( (ptr) ); }
#define freeMemoryArray(ptr)               { freeMemory_call( (void**)(ptr), __FILE__, __LINE__ ); \
                                             assert(*(ptr) == NULL); }
#define freeMemoryArrayNull(ptr)           { if( *(ptr) != NULL ) freeMemoryArray( (ptr) ); }

#ifndef NDEBUG
#define getPointerSize(ptr)                getPointerSize_call(ptr)
#define displayMemory()                    displayMemory_call()
#define checkEmptyMemory()                 checkEmptyMemory_call()
#else
#define getPointerSize(ptr)                0
#define displayMemory()                    /**/
#define checkEmptyMemory()                 /**/
#endif


/** allocates memory; returns NULL if memory allocation failed */
extern
void* allocMemory_call(
   size_t           size,               /**< size of memory element to allocate */
   const char*      filename,           /**< source file where the allocation is performed */
   int              line                /**< line number in source file where the allocation is performed */
   );

/** allocates memory; returns NULL if memory allocation failed */
extern
void* reallocMemory_call(
   void*            ptr,                /**< pointer to memory to reallocate */
   size_t           size,               /**< new size of memory element */
   const char*      filename,           /**< source file where the reallocation is performed */
   int              line                /**< line number in source file where the reallocation is performed */
   );

/** clears a memory element (i.e. fills it with zeros) */
extern
void clearMemory_call(
   void*            ptr,                /**< pointer to memory element */
   size_t           size                /**< size of memory element */
   );

/** copies the contents of one memory element into another memory element */
extern
void copyMemory_call(
   void*            ptr,                /**< pointer to target memory element */
   const void*      source,             /**< pointer to source memory element */
   size_t           size                /**< size of memory element to copy */
   );

/** allocates memory and copies the contents of the given memory element into the new memory element */
extern
void* duplicateMemory_call(
   const void*      source,             /**< pointer to source memory element */
   size_t           size,               /**< size of memory element to copy */
   const char*      filename,           /**< source file where the duplication is performed */
   int              line                /**< line number in source file where the duplication is performed */
   );

/** frees an allocated memory element */
extern
void freeMemory_call(
   void**           ptr,                /**< address of pointer to memory element */
   const char*      filename,           /**< source file where the deallocation is performed */
   int              line                /**< line number in source file where the deallocation is performed */
   );

/** returns the size of an allocated memory element */
extern
size_t getPointerSize_call(
   const void*      ptr                 /**< pointer to allocated memory */
   );

/** outputs information about currently allocated memory to the screen */
extern
void displayMemory_call(
   void
   );

/** displays a warning message on the screen, if allocated memory exists */
extern
void checkEmptyMemory_call(
   void
   );




/********************************************************************
 * Chunk Memory Management
 *
 * Efficient memory management for multiple objects of the same size
 ********************************************************************/

typedef struct ChkMem CHKMEM;           /**< collection of memory chunks of the same element size */


#ifndef NOBLOCKMEM

#define createChunkMemory(sz,isz,gbf)      createChunkMemory_call( (sz), (isz), (gbf), __FILE__, __LINE__ )
#define clearChunkMemory(mem)              clearChunkMemory_call( (mem), __FILE__, __LINE__ )
#define clearChunkMemoryNull(mem)          { if( (mem) != NULL ) clearChunkMemory( (mem) ); }
#define destroyChunkMemory(mem)            destroyChunkMemory_call( (mem), __FILE__, __LINE__ )
#define destroyChunkMemoryNull(mem)        { if( *(mem) != NULL ) destroyChunkMemory( (mem) ); }

#define allocChunkMemory(mem,ptr)          (*(void**)(ptr) = allocChunkMemory_call((mem), sizeof(**(ptr)), \
                                           __FILE__, __LINE__))
#define duplicateChunkMemory(mem, ptr, source) \
                                           (*(void**)(ptr) = duplicateChunkMemory_call((mem), (const void*)(source), \
                                           sizeof(**(ptr)), __FILE__, __LINE__ ))
#define freeChunkMemory(mem,ptr)           { freeChunkMemory_call( (mem), (void**)(ptr), sizeof(**(ptr)), \
                                             __FILE__, __LINE__ ); \
                                             assert(*(ptr) == NULL); }
#define freeChunkMemoryNull(mem,ptr)       { if( *(ptr) != NULL ) freeChunkMemory( (mem), (ptr) ); }
#define garbagecollectChunkMemory(mem)     garbagecollectChunkMemory_call(mem)

#else

/* block memory management mapped to standard memory management */

#define createChunkMemory(sz,isz,gbf)           (void*)(0x01) /* dummy to not return a NULL pointer */
#define clearChunkMemory(mem)                   /**/
#define clearChunkMemoryNull(mem)               /**/
#define destroyChunkMemory(mem)                 /**/
#define destroyChunkMemoryNull(mem)             /**/
#define allocChunkMemory(mem,ptr)               allocMemory(ptr)
#define duplicateChunkMemory(mem, ptr, source)  duplicateMemory(ptr,source)
#define freeChunkMemory(mem,ptr)                freeMemory(ptr)
#define freeChunkMemoryNull(mem,ptr)            freeMemoryNull(ptr)
#define garbagecollectChunkMemory(mem)          /**/

#endif


/** aligns the given byte size corresponding to the minimal alignment for chunk and block memory */
extern
void alignMemsize(
   size_t*          size                /**< pointer to the size to align */
   );

/** checks whether the given size meets the alignment conditions for chunk and block memory  */
extern
int isAligned(
   size_t           size                /**< size to check for alignment */
   );

/** creates a new chunk block data structure */
extern
CHKMEM* createChunkMemory_call(
   int              size,               /**< element size of the chunk block */
   int              initchunksize,      /**< number of elements in the first chunk of the chunk block */
   int              garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                         *   elements are free (-1: disable garbage collection) */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** clears a chunk block data structure */
extern
void clearChunkMemory_call(
   CHKMEM*          chkmem,             /**< chunk block */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** destroys and frees a chunk block data structure */
extern
void destroyChunkMemory_call(
   CHKMEM**         chkmem,             /**< pointer to chunk block */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** allocates a memory element of the given chunk block */
extern
void* allocChunkMemory_call(
   CHKMEM*          chkmem,             /**< chunk block */
   size_t           size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** duplicates a given memory element by allocating a new element of the same chunk block and copying the data */
extern
void* duplicateChunkMemory_call(
   CHKMEM*          chkmem,             /**< chunk block */
   const void*      source,             /**< source memory element */
   size_t           size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** frees a memory element of the given chunk block */
extern
void freeChunkMemory_call(
   CHKMEM*          chkmem,             /**< chunk block */
   void**           ptr,                /**< pointer to memory element to free */
   size_t           size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** calls garbage collection of chunk block and frees chunks without allocated memory elements */
extern
void garbagecollectChunkMemory_call(
   CHKMEM*          chkmem              /**< chunk block */
   );




/***********************************************************
 * Block Memory Management
 *
 * Efficient memory management for objects of varying sizes
 ***********************************************************/

typedef struct BlkMem BLKMEM;           /**< block memory: collection of chunk blocks */


#ifndef NOBLOCKMEM

/* block memory methods for faster memory access */

#define createBlockMemory(csz,gbf)         createBlockMemory_call( (csz), (gbf), __FILE__, __LINE__ )
#define clearBlockMemory(mem)              clearBlockMemory_call( (mem), __FILE__, __LINE__ )
#define clearBlockMemoryNull(mem)          { if( (mem) != NULL ) clearBlockMemory( (mem) ); }
#define destroyBlockMemory(mem)            destroyBlockMemory_call( (mem), __FILE__, __LINE__ )
#define destroyBlockMemoryNull(mem)        { if( *(mem) != NULL ) destroyBlockMemory( (mem) ); }

#define allocBlockMemory(mem,ptr)          (*(void**)(ptr) = allocBlockMemory_call((mem), sizeof(**(ptr)), \
                                           __FILE__, __LINE__))
#define allocBlockMemoryArray(mem,ptr,num) (*(void**)(ptr) = allocBlockMemory_call((mem), (num)*sizeof(**(ptr)), \
                                           __FILE__, __LINE__))
#define allocBlockMemorySize(mem,ptr,size) (*(void**)(ptr) = allocBlockMemory_call((mem), (size_t)(size), \
                                           __FILE__, __LINE__))
#define reallocBlockMemoryArray(mem,ptr,oldnum,newnum) \
                                           (*(void**)(ptr) = reallocBlockMemory_call((mem), *(void**)(ptr), \
                                           (oldnum)*sizeof(**(ptr)), (newnum)*sizeof(**(ptr)), __FILE__, __LINE__))
#define reallocBlockMemorySize(mem,ptr,oldsize,newsize) \
                                           (*(void**)(ptr) = reallocBlockMemory_call((mem), *(void**)(ptr), \
                                           (size_t)(oldsize), (size_t)(newsize), __FILE__, __LINE__))
#define duplicateBlockMemory(mem, ptr, source) \
                                           (*(void**)(ptr) = duplicateBlockMemory_call((mem), (const void*)(source), \
                                           sizeof(**(ptr)), __FILE__, __LINE__ ))
#define duplicateBlockMemoryArray(mem, ptr, source, num) \
                                           (*(void**)(ptr) = duplicateBlockMemory_call( (mem), (const void*)(source), \
                                           (num)*sizeof(**(ptr)), __FILE__, __LINE__ ))
#define freeBlockMemory(mem,ptr)           { freeBlockMemory_call( (mem), (void**)(ptr), sizeof(**(ptr)), \
                                             __FILE__, __LINE__ ); \
                                             assert(*(ptr) == NULL); }
#define freeBlockMemoryNull(mem,ptr)       { if( *(ptr) != NULL ) freeBlockMemory( (mem), (ptr) ); }
#define freeBlockMemoryArray(mem,ptr,num)  { freeBlockMemory_call( (mem), (void**)(ptr), (num)*sizeof(**(ptr)), \
                                             __FILE__, __LINE__ ); \
                                             assert(*(ptr) == NULL); }
#define freeBlockMemoryArrayNull(mem,ptr,num)  { if( *(ptr) != NULL ) freeBlockMemoryArray( (mem), (ptr), (num) ); }
#define freeBlockMemorySize(mem,ptr,size)  { freeBlockMemory_call( (mem), (void**)(ptr), (size_t)(size), \
                                             __FILE__, __LINE__ ); \
                                             assert(*(ptr) == NULL); }
#define freeBlockMemorySizeNull(mem,ptr,size)  { if( *(ptr) != NULL ) freeBlockMemorySize( (mem), (ptr), (size) ); }
#define garbagecollectBlockMemory(mem)     garbagecollectBlockMemory_call(mem)
#define getBlockMemoryUsed(mem)            getBlockMemoryUsed_call(mem)
#define getBlockPointerSize(mem,ptr)       getBlockPointerSize_call((mem), (ptr))
#define displayBlockMemory(mem)            displayBlockMemory_call(mem)
#define blockMemoryCheckEmpty(mem)         checkEmptyBlockMemory_call(mem)

#else

/* block memory management mapped to standard memory management */

#define createBlockMemory(csz,gbf)                        (void*)(0x01) /* dummy to not return a NULL pointer */
#define clearBlockMemory(mem)                             /**/
#define clearBlockMemoryNull(mem)                         /**/
#define destroyBlockMemory(mem)                           /**/
#define destroyBlockMemoryNull(mem)                       /**/
#define allocBlockMemory(mem,ptr)                         allocMemory(ptr)
#define allocBlockMemoryArray(mem,ptr,num)                allocMemoryArray(ptr,num)
#define allocBlockMemorySize(mem,ptr,size)                allocMemorySize(ptr,size)
#define reallocBlockMemoryArray(mem,ptr,oldnum,newnum)    reallocMemoryArray(ptr,newnum)
#define reallocBlockMemorySize(mem,ptr,oldsize,newsize)   reallocMemorySize(ptr,newsize)
#define duplicateBlockMemory(mem, ptr, source)            duplicateMemory(ptr,source)
#define duplicateBlockMemoryArray(mem, ptr, source, num)  duplicateMemoryArray(ptr,source,num)
#define freeBlockMemory(mem,ptr)                          freeMemory(ptr)
#define freeBlockMemoryNull(mem,ptr)                      freeMemoryNull(ptr)
#define freeBlockMemoryArray(mem,ptr,num)                 freeMemoryArray(ptr)
#define freeBlockMemoryArrayNull(mem,ptr,num)             freeMemoryArrayNull(ptr)
#define freeBlockMemorySize(mem,ptr,size)                 freeMemory(ptr)
#define freeBlockMemorySizeNull(mem,ptr,size)             freeMemoryNull(ptr)
#define garbagecollectBlockMemory(mem)                    /**/
#define getBlockMemoryUsed(mem)                           0
#define getBlockPointerSize(mem,ptr)                      0
#define displayBlockMemory(mem)                           /**/
#define blockMemoryCheckEmpty(mem)                        /**/

#endif


/** creates a block memory allocation data structure */
extern
BLKMEM* createBlockMemory_call(
   int              initchunksize,      /**< number of elements in the first chunk of each chunk block */
   int              garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                         *   elements are free (-1: disable garbage collection) */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** frees all chunk blocks in the block memory */
extern
void clearBlockMemory_call(
   BLKMEM*          blkmem,             /**< block memory */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** clears and deletes block memory */
extern
void destroyBlockMemory_call(
   BLKMEM**         blkmem,             /**< pointer to block memory */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** allocates memory in the block memory pool */
void* allocBlockMemory_call(
   BLKMEM*          blkmem,             /**< block memory */
   size_t           size,               /**< size of memory element to allocate */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** resizes memory element in the block memory pool, and copies the data */
extern
void* reallocBlockMemory_call(
   BLKMEM*          blkmem,             /**< block memory */
   void*            ptr,                /**< memory element to reallocated */
   size_t           oldsize,            /**< old size of memory element */
   size_t           newsize,            /**< new size of memory element */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** duplicates memory element in the block memory pool, and copies the data */
extern
void* duplicateBlockMemory_call(
   BLKMEM*          blkmem,             /**< block memory */
   const void*      source,             /**< memory element to duplicate */
   size_t           size,               /**< size of memory elements */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** frees memory element in the block memory pool */
extern
void freeBlockMemory_call(
   BLKMEM*          blkmem,             /**< block memory */
   void**           ptr,                /**< pointer to memory element to free */
   size_t           size,               /**< size of memory element */
   const char*      filename,           /**< source file of the function call */
   int              line                /**< line number in source file of the function call */
   );

/** calls garbage collection of block memory, frees chunks without allocated memory elements, and frees
 *  chunk blocks without any chunks
 */
extern
void garbagecollectBlockMemory_call(
   BLKMEM*          blkmem              /**< block memory */
   );

/** returns the number of allocated bytes in the block memory */
long long getBlockMemoryUsed_call(
   const BLKMEM*    blkmem              /**< block memory */
   );

/** returns the size of the given memory element; returns 0, if the element is not member of the block memory */
extern
size_t getBlockPointerSize_call(
   const BLKMEM*    blkmem,             /**< block memory */
   const void*      ptr                 /**< memory element */
   );

/** outputs allocation diagnostics of block memory */
extern
void displayBlockMemory_call(
   const BLKMEM*    blkmem              /**< block memory */
   );

/** outputs warning messages, if there are allocated elements in the block memory */
extern
void checkEmptyBlockMemory_call(
   const BLKMEM*    blkmem              /**< block memory */
   );


#endif
