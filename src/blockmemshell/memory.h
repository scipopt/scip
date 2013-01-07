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

/**@file   memory.h
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BMS_MEMORY_H__
#define __BMS_MEMORY_H__

#include <limits.h>
#include <stdlib.h>

/* special thanks to Daniel Junglas for following template and macros */
#ifdef __cplusplus

template<typename T> T* docast(T*, void *v) { return reinterpret_cast<T*>(v); }

extern "C" {

#define ASSIGN(pointerstarstar, voidstarfunction) (*(pointerstarstar) = docast(*(pointerstarstar), (voidstarfunction)))

#else

#define ASSIGN(pointerstarstar, voidstarfunction) (*(pointerstarstar) = (voidstarfunction))

#endif

/*
 * Define the marco EXTERN depending if the OS is Windows or not
 */
#ifndef EXTERN

#if defined(_WIN32) || defined(_WIN64)
#define EXTERN __declspec(dllexport)
#else
#define EXTERN extern
#endif

#endif



/*************************************************************************************
 * Standard Memory Management
 *
 * In debug mode, these methods extend malloc() and free() by logging all currently
 * allocated memory elements in an allocation list. This can be used as a simple leak
 * detection.
 *************************************************************************************/

/* Check for integer overflow in allocation size */
#ifdef NDEBUG
#define BMSallocMemoryArray(ptr,num)          ASSIGN((ptr), BMSallocMemory_call( (num)*sizeof(**(ptr)), \
                                                __FILE__, __LINE__ ))
#define BMSreallocMemoryArray(ptr,num)        ASSIGN((ptr), BMSreallocMemory_call( *(ptr), (num)*sizeof(**(ptr)), \
                                                __FILE__, __LINE__ ))
#define BMSallocMemoryArrayCPP(num,size)      BMSallocMemory_call( (size_t)((num)*(size)), __FILE__, __LINE__ )
#define BMSallocClearMemoryArray(ptr,num)     ASSIGN((ptr), BMSallocClearMemory_call((size_t)(num), sizeof(**(ptr)), __FILE__, __LINE__ ))
#else
#define BMSallocMemoryArray(ptr,num)          ( ( ((size_t)(num)) > (UINT_MAX / sizeof(**(ptr))) ) \
                                                ? ( ASSIGN((ptr), NULL) )                          \
						: ( ASSIGN((ptr), BMSallocMemory_call( (num)*sizeof(**(ptr)),  \
										       __FILE__, __LINE__ )) ) )
#define BMSreallocMemoryArray(ptr,num)        ( ( ((size_t)(num)) > (UINT_MAX / sizeof(**(ptr))) ) \
                                                ? ( ASSIGN((ptr), NULL) )                          \
						: ( ASSIGN((ptr), BMSreallocMemory_call( *(ptr), (num)*sizeof(**(ptr)), \
											 __FILE__, __LINE__ )) )        )
#define BMSallocMemoryArrayCPP(num,size)      ( ( ((size_t)(num)) > (UINT_MAX / size) ) \
                                                ? ( ASSIGN((ptr), NULL) )                          \
						:( BMSallocMemory_call( (size_t)((num)*(size)), __FILE__, __LINE__ ) ) )
#define BMSallocClearMemoryArray(ptr,num)     ( ( ((size_t)(num)) > (UINT_MAX / sizeof(**(ptr))) ) \
                                                ? ( ASSIGN((ptr), NULL) )                          \
                                                : ( ASSIGN((ptr), BMSallocClearMemory_call((size_t)(num), sizeof(**(ptr)), __FILE__, __LINE__ )) ) )
#endif

#define BMSallocMemory(ptr)                   ASSIGN((ptr), BMSallocMemory_call( sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSallocMemorySize(ptr,size)          ASSIGN((ptr), BMSallocMemory_call( (size_t)(size), __FILE__, __LINE__ ))
#define BMSallocMemoryCPP(size)               BMSallocMemory_call( (size_t)(size), __FILE__, __LINE__ )
#define BMSreallocMemorySize(ptr,size)        ASSIGN((ptr), BMSreallocMemory_call( *(ptr), (size_t)(size), \
                                                __FILE__, __LINE__ ))
#define BMSclearMemory(ptr)                   BMSclearMemory_call( (void*)(ptr), sizeof(*(ptr)) )
#define BMSclearMemoryArray(ptr, num)         BMSclearMemory_call( (void*)(ptr), (num)*sizeof(*(ptr)) )
#define BMSclearMemorySize(ptr, size)         BMSclearMemory_call( (void*)(ptr), (size_t)(size) )

#define BMScopyMemory(ptr, source)            BMScopyMemory_call( (void*)(ptr), (const void*)(source), sizeof(*(ptr)) )
#define BMScopyMemoryArray(ptr, source, num)  BMScopyMemory_call( (void*)(ptr), (const void*)(source), (num)*sizeof(*(ptr)) )
#define BMScopyMemorySize(ptr, source, size)  BMScopyMemory_call( (void*)(ptr), (const void*)(source), (size_t)(size) )

#define BMSmoveMemory(ptr, source)            BMSmoveMemory_call( (void*)(ptr), (const void*)(source), sizeof(*(ptr)) )
#define BMSmoveMemoryArray(ptr, source, num)  BMSmoveMemory_call( (void*)(ptr), (const void*)(source), (num) * sizeof(*(ptr)) )
#define BMSmoveMemorySize(ptr, source, size)  BMSmoveMemory_call( (void*)(ptr), (const void*)(source), (size_t)(size) )

#define BMSduplicateMemory(ptr, source)       ASSIGN((ptr), BMSduplicateMemory_call( (const void*)(source), \
                                                sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSduplicateMemoryArray(ptr, source, num) \
                                                ASSIGN((ptr), BMSduplicateMemory_call( (const void*)(source), \
                                                (num)*sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSduplicateMemorySize(ptr, source, size) \
                                                ASSIGN((ptr), BMSduplicateMemory_call( (const void*)(source), \
                                                (size_t)(size), __FILE__, __LINE__ ))
#define BMSfreeMemory(ptr)                    { BMSfreeMemory_call( (void*)(*(ptr)), __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeMemoryNull(ptr)                { if( *(ptr) != NULL ) BMSfreeMemory( (ptr) ); }
#define BMSfreeMemoryArray(ptr)               { BMSfreeMemory_call( (void*)(*(ptr)), __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeMemoryArrayNull(ptr)           { if( *(ptr) != NULL ) BMSfreeMemoryArray( (ptr) ); }
#define BMSfreeMemorySize(ptr)                { BMSfreeMemory_call( (void*)(*(ptr)), __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeMemorySizeNull(ptr)            { if( *(ptr) != NULL ) BMSfreeMemorySize( (ptr) ); }

#ifndef NDEBUG
#define BMSgetPointerSize(ptr)                BMSgetPointerSize_call(ptr)
#define BMSdisplayMemory()                    BMSdisplayMemory_call()
#define BMScheckEmptyMemory()                 BMScheckEmptyMemory_call()
#define BMSgetMemoryUsed()                    BMSgetMemoryUsed_call()
#else
#define BMSgetPointerSize(ptr)                0
#define BMSdisplayMemory()                    /**/
#define BMScheckEmptyMemory()                 /**/
#define BMSgetMemoryUsed()                    0LL
#endif

/** allocates memory and initializes it with 0; returns NULL if memory allocation failed */
EXTERN
void* BMSallocClearMemory_call(
   size_t                num,                /**< number of memory element to allocate */
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   );

/** allocates memory; returns NULL if memory allocation failed */
EXTERN
void* BMSallocMemory_call(
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file where the allocation is performed */
   int                   line                /**< line number in source file where the allocation is performed */
   );

/** allocates memory; returns NULL if memory allocation failed */
EXTERN
void* BMSreallocMemory_call(
   void*                 ptr,                /**< pointer to memory to reallocate */
   size_t                size,               /**< new size of memory element */
   const char*           filename,           /**< source file where the reallocation is performed */
   int                   line                /**< line number in source file where the reallocation is performed */
   );

/** clears a memory element (i.e. fills it with zeros) */
EXTERN
void BMSclearMemory_call(
   void*                 ptr,                /**< pointer to memory element */
   size_t                size                /**< size of memory element */
   );

/** copies the contents of one memory element into another memory element */
EXTERN
void BMScopyMemory_call(
   void*                 ptr,                /**< pointer to target memory element */
   const void*           source,             /**< pointer to source memory element */
   size_t                size                /**< size of memory element to copy */
   );

/** moves the contents of one memory element into another memory element, should be used if both elements overlap,
 *  otherwise BMScopyMemory is faster
 */
EXTERN
void BMSmoveMemory_call(
   void*                 ptr,                /**< pointer to target memory element */
   const void*           source,             /**< pointer to source memory element */
   size_t                size                /**< size of memory element to copy */
   );

/** allocates memory and copies the contents of the given memory element into the new memory element */
EXTERN
void* BMSduplicateMemory_call(
   const void*           source,             /**< pointer to source memory element */
   size_t                size,               /**< size of memory element to copy */
   const char*           filename,           /**< source file where the duplication is performed */
   int                   line                /**< line number in source file where the duplication is performed */
   );

/** frees an allocated memory element */
EXTERN
void BMSfreeMemory_call(
   void*                 ptr,                /**< pointer to memory element */
   const char*           filename,           /**< source file where the deallocation is performed */
   int                   line                /**< line number in source file where the deallocation is performed */
   );

/** returns the size of an allocated memory element */
EXTERN
size_t BMSgetPointerSize_call(
   const void*           ptr                 /**< pointer to allocated memory */
   );

/** outputs information about currently allocated memory to the screen */
EXTERN
void BMSdisplayMemory_call(
   void
   );

/** displays a warning message on the screen, if allocated memory exists */
EXTERN
void BMScheckEmptyMemory_call(
   void
   );

/** returns total number of allocated bytes */
EXTERN
long long BMSgetMemoryUsed_call(
   void
   );




/********************************************************************
 * Chunk Memory Management
 *
 * Efficient memory management for multiple objects of the same size
 ********************************************************************/

typedef struct BMS_ChkMem BMS_CHKMEM;           /**< collection of memory chunks of the same element size */


#ifndef BMS_NOBLOCKMEM

#define BMScreateChunkMemory(sz,isz,gbf)      BMScreateChunkMemory_call( (sz), (isz), (gbf), __FILE__, __LINE__ )
#define BMSclearChunkMemory(mem)              BMSclearChunkMemory_call( (mem), __FILE__, __LINE__ )
#define BMSclearChunkMemoryNull(mem)          { if( (mem) != NULL ) BMSclearChunkMemory( (mem) ); }
#define BMSdestroyChunkMemory(mem)            BMSdestroyChunkMemory_call( (mem), __FILE__, __LINE__ )
#define BMSdestroyChunkMemoryNull(mem)        { if( *(mem) != NULL ) BMSdestroyChunkMemory( (mem) ); }

#define BMSallocChunkMemory(mem,ptr)          ASSIGN((ptr), BMSallocChunkMemory_call((mem), sizeof(**(ptr)), \
                                                __FILE__, __LINE__))
#define BMSduplicateChunkMemory(mem, ptr, source) \
                                                ASSIGN((ptr), BMSduplicateChunkMemory_call((mem), (const void*)(source), \
                                                sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSfreeChunkMemory(mem,ptr)           { BMSfreeChunkMemory_call( (mem), (void*)(*(ptr)), sizeof(**(ptr)), \
                                                  __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeChunkMemoryNull(mem,ptr)       { if( *(ptr) != NULL ) BMSfreeChunkMemory( (mem), (ptr) ); }
#define BMSgarbagecollectChunkMemory(mem)     BMSgarbagecollectChunkMemory_call(mem)
#define BMSgetChunkMemoryUsed(mem)            BMSgetChunkMemoryUsed_call(mem)

#else

/* block memory management mapped to standard memory management */

#define BMScreateChunkMemory(sz,isz,gbf)           (void*)(0x01) /* dummy to not return a NULL pointer */
#define BMSclearChunkMemory(mem)                   /**/
#define BMSclearChunkMemoryNull(mem)               /**/
#define BMSdestroyChunkMemory(mem)                 /**/
#define BMSdestroyChunkMemoryNull(mem)             /**/
#define BMSallocChunkMemory(mem,ptr)               BMSallocMemory(ptr)
#define BMSduplicateChunkMemory(mem, ptr, source)  BMSduplicateMemory(ptr,source)
#define BMSfreeChunkMemory(mem,ptr)                BMSfreeMemory(ptr)
#define BMSfreeChunkMemoryNull(mem,ptr)            BMSfreeMemoryNull(ptr)
#define BMSgarbagecollectChunkMemory(mem)          /**/
#define BMSgetChunkMemoryUsed(mem)                 0LL

#endif


/** aligns the given byte size corresponding to the minimal alignment for chunk and block memory */
EXTERN
void BMSalignMemsize(
   size_t*               size                /**< pointer to the size to align */
   );

/** checks whether the given size meets the alignment conditions for chunk and block memory  */
EXTERN
int BMSisAligned(
   size_t                size                /**< size to check for alignment */
   );

/** creates a new chunk block data structure */
EXTERN
BMS_CHKMEM* BMScreateChunkMemory_call(
   size_t                size,               /**< element size of the chunk block */
   int                   initchunksize,      /**< number of elements in the first chunk of the chunk block */
   int                   garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** clears a chunk block data structure */
EXTERN
void BMSclearChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** destroys and frees a chunk block data structure */
EXTERN
void BMSdestroyChunkMemory_call(
   BMS_CHKMEM**          chkmem,             /**< pointer to chunk block */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates a memory element of the given chunk block */
EXTERN
void* BMSallocChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** duplicates a given memory element by allocating a new element of the same chunk block and copying the data */
EXTERN
void* BMSduplicateChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   const void*           source,             /**< source memory element */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees a memory element of the given chunk block */
EXTERN
void BMSfreeChunkMemory_call(
   BMS_CHKMEM*           chkmem,             /**< chunk block */
   void*                 ptr,                /**< memory element to free */
   size_t                size,               /**< size of memory element to allocate (only needed for sanity check) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** calls garbage collection of chunk block and frees chunks without allocated memory elements */
EXTERN
void BMSgarbagecollectChunkMemory_call(
   BMS_CHKMEM*           chkmem              /**< chunk block */
   );

/** returns the number of allocated bytes in the chunk block */
EXTERN
long long BMSgetChunkMemoryUsed_call(
   const BMS_CHKMEM*     chkmem              /**< chunk block */
   );




/***********************************************************
 * Block Memory Management
 *
 * Efficient memory management for objects of varying sizes
 ***********************************************************/

typedef struct BMS_BlkMem BMS_BLKMEM;           /**< block memory: collection of chunk blocks */


#ifndef BMS_NOBLOCKMEM

/* block memory methods for faster memory access */

/* Check for integer overflow in allocation size */
#ifdef NDEBUG
#define BMSallocBlockMemoryArray(mem,ptr,num) ASSIGN((ptr), BMSallocBlockMemory_call((mem), (num)*sizeof(**(ptr)), \
                                                __FILE__, __LINE__))
#define BMSreallocBlockMemoryArray(mem,ptr,oldnum,newnum) \
                                                ASSIGN((ptr), BMSreallocBlockMemory_call((mem), (void*)(*(ptr)), \
                                                (oldnum)*sizeof(**(ptr)), (newnum)*sizeof(**(ptr)), __FILE__, __LINE__))
#else
#define BMSallocBlockMemoryArray(mem,ptr,num) ( ( ((size_t)(num)) > UINT_MAX / sizeof(**(ptr)) ) \
                                                ? ( ASSIGN((ptr), NULL) )                        \
						: ( ASSIGN((ptr), BMSallocBlockMemory_call((mem), (num)*sizeof(**(ptr)),\
											   __FILE__, __LINE__)) )       )
#define BMSreallocBlockMemoryArray(mem,ptr,oldnum,newnum) \
						( ( ((size_t)(newnum)) > UINT_MAX / sizeof(**(ptr)) )                \
                                                ? ( ASSIGN((ptr), NULL) )                        \
                                                : ( ASSIGN((ptr), BMSreallocBlockMemory_call((mem), (void*)(*(ptr)), \
                                                         (oldnum)*sizeof(**(ptr)), (newnum)*sizeof(**(ptr)), __FILE__, \
                                                         __LINE__)) )  )
#endif

#define BMScreateBlockMemory(csz,gbf)         BMScreateBlockMemory_call( (csz), (gbf), __FILE__, __LINE__ )
#define BMSclearBlockMemory(mem)              BMSclearBlockMemory_call( (mem), __FILE__, __LINE__ )
#define BMSclearBlockMemoryNull(mem)          { if( (mem) != NULL ) BMSclearBlockMemory( (mem) ); }
#define BMSdestroyBlockMemory(mem)            BMSdestroyBlockMemory_call( (mem), __FILE__, __LINE__ )
#define BMSdestroyBlockMemoryNull(mem)        { if( *(mem) != NULL ) BMSdestroyBlockMemory( (mem) ); }

#define BMSallocBlockMemory(mem,ptr)          ASSIGN((ptr), BMSallocBlockMemory_call((mem), sizeof(**(ptr)), \
                                                __FILE__, __LINE__))
#define BMSallocBlockMemorySize(mem,ptr,size) ASSIGN((ptr), BMSallocBlockMemory_call((mem), (size_t)(size), \
                                                __FILE__, __LINE__))
#define BMSreallocBlockMemorySize(mem,ptr,oldsize,newsize) \
                                                ASSIGN((ptr), BMSreallocBlockMemory_call((mem), (void*)(*(ptr)), \
                                                (size_t)(oldsize), (size_t)(newsize), __FILE__, __LINE__))
#define BMSduplicateBlockMemory(mem, ptr, source) \
                                                ASSIGN((ptr), BMSduplicateBlockMemory_call((mem), (const void*)(source), \
                                                sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSduplicateBlockMemoryArray(mem, ptr, source, num) \
                                                ASSIGN((ptr), BMSduplicateBlockMemory_call( (mem), (const void*)(source), \
                                                (num)*sizeof(**(ptr)), __FILE__, __LINE__ ))
#define BMSfreeBlockMemory(mem,ptr)           { BMSfreeBlockMemory_call( (mem), (void*)(*(ptr)), sizeof(**(ptr)), \
                                                  __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeBlockMemoryNull(mem,ptr)       { if( *(ptr) != NULL ) BMSfreeBlockMemory( (mem), (ptr) ); }
#define BMSfreeBlockMemoryArray(mem,ptr,num)  { BMSfreeBlockMemory_call( (mem), (void*)(*(ptr)), (num)*sizeof(**(ptr)), \
                                                  __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeBlockMemoryArrayNull(mem,ptr,num)  { if( *(ptr) != NULL ) BMSfreeBlockMemoryArray( (mem), (ptr), (num) ); }
#define BMSfreeBlockMemorySize(mem,ptr,size)  { BMSfreeBlockMemory_call( (mem), (void*)(*(ptr)), (size_t)(size), \
                                                  __FILE__, __LINE__ ); \
                                                  *(ptr) = NULL; }
#define BMSfreeBlockMemorySizeNull(mem,ptr,size)  { if( *(ptr) != NULL ) BMSfreeBlockMemorySize( (mem), (ptr), (size) ); }
#define BMSgarbagecollectBlockMemory(mem)     BMSgarbagecollectBlockMemory_call(mem)
#define BMSgetBlockMemoryUsed(mem)            BMSgetBlockMemoryUsed_call(mem)
#define BMSgetBlockPointerSize(mem,ptr)       BMSgetBlockPointerSize_call((mem), (ptr))
#define BMSdisplayBlockMemory(mem)            BMSdisplayBlockMemory_call(mem)
#define BMSblockMemoryCheckEmpty(mem)         BMScheckEmptyBlockMemory_call(mem)

#else

/* block memory management mapped to standard memory management */

#define BMScreateBlockMemory(csz,gbf)                        (void*)(0x01) /* dummy to not return a NULL pointer */
#define BMSclearBlockMemory(mem)                             /**/
#define BMSclearBlockMemoryNull(mem)                         /**/
#define BMSdestroyBlockMemory(mem)                           /**/
#define BMSdestroyBlockMemoryNull(mem)                       /**/
#define BMSallocBlockMemory(mem,ptr)                         BMSallocMemory(ptr)
#define BMSallocBlockMemoryArray(mem,ptr,num)                BMSallocMemoryArray(ptr,num)
#define BMSallocBlockMemorySize(mem,ptr,size)                BMSallocMemorySize(ptr,size)
#define BMSreallocBlockMemoryArray(mem,ptr,oldnum,newnum)    BMSreallocMemoryArray(ptr,newnum)
#define BMSreallocBlockMemorySize(mem,ptr,oldsize,newsize)   BMSreallocMemorySize(ptr,newsize)
#define BMSduplicateBlockMemory(mem, ptr, source)            BMSduplicateMemory(ptr,source)
#define BMSduplicateBlockMemoryArray(mem, ptr, source, num)  BMSduplicateMemoryArray(ptr,source,num)
#define BMSfreeBlockMemory(mem,ptr)                          BMSfreeMemory(ptr)
#define BMSfreeBlockMemoryNull(mem,ptr)                      BMSfreeMemoryNull(ptr)
#define BMSfreeBlockMemoryArray(mem,ptr,num)                 BMSfreeMemoryArray(ptr)
#define BMSfreeBlockMemoryArrayNull(mem,ptr,num)             BMSfreeMemoryArrayNull(ptr)
#define BMSfreeBlockMemorySize(mem,ptr,size)                 BMSfreeMemory(ptr)
#define BMSfreeBlockMemorySizeNull(mem,ptr,size)             BMSfreeMemoryNull(ptr)
#define BMSgarbagecollectBlockMemory(mem)                    /**/
#define BMSgetBlockMemoryUsed(mem)                           0LL
#define BMSgetBlockPointerSize(mem,ptr)                      0
#define BMSdisplayBlockMemory(mem)                           /**/
#define BMSblockMemoryCheckEmpty(mem)                        /**/

#endif


/** creates a block memory allocation data structure */
EXTERN
BMS_BLKMEM* BMScreateBlockMemory_call(
   int                   initchunksize,      /**< number of elements in the first chunk of each chunk block */
   int                   garbagefactor,      /**< garbage collector is called, if at least garbagefactor * avg. chunksize 
                                              *   elements are free (-1: disable garbage collection) */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees all chunk blocks in the block memory */
EXTERN
void BMSclearBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** clears and deletes block memory */
EXTERN
void BMSdestroyBlockMemory_call(
   BMS_BLKMEM**          blkmem,             /**< pointer to block memory */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** allocates memory in the block memory pool */
void* BMSallocBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   size_t                size,               /**< size of memory element to allocate */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** resizes memory element in the block memory pool, and copies the data */
EXTERN
void* BMSreallocBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 ptr,                /**< memory element to reallocated */
   size_t                oldsize,            /**< old size of memory element */
   size_t                newsize,            /**< new size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** duplicates memory element in the block memory pool, and copies the data */
EXTERN
void* BMSduplicateBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const void*           source,             /**< memory element to duplicate */
   size_t                size,               /**< size of memory elements */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** frees memory element in the block memory pool */
EXTERN
void BMSfreeBlockMemory_call(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   void*                 ptr,                /**< memory element to free */
   size_t                size,               /**< size of memory element */
   const char*           filename,           /**< source file of the function call */
   int                   line                /**< line number in source file of the function call */
   );

/** calls garbage collection of block memory, frees chunks without allocated memory elements, and frees
 *  chunk blocks without any chunks
 */
EXTERN
void BMSgarbagecollectBlockMemory_call(
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns the number of allocated bytes in the block memory */
EXTERN
long long BMSgetBlockMemoryUsed_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** returns the size of the given memory element; returns 0, if the element is not member of the block memory */
EXTERN
size_t BMSgetBlockPointerSize_call(
   const BMS_BLKMEM*     blkmem,             /**< block memory */
   const void*           ptr                 /**< memory element */
   );

/** outputs allocation diagnostics of block memory */
EXTERN
void BMSdisplayBlockMemory_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

/** outputs warning messages, if there are allocated elements in the block memory */
EXTERN
void BMScheckEmptyBlockMemory_call(
   const BMS_BLKMEM*     blkmem              /**< block memory */
   );

#ifdef __cplusplus
}
#endif

#endif
