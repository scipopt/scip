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
#pragma ident "@(#) $Id: memory.h,v 1.21 2004/12/10 14:23:01 bzfpfend Exp $"

/**@file   memory.h
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEMORY_H__
#define __MEMORY_H__


#include <stdlib.h>



/*
 * standard memory management
 */

#ifndef NOSAFEMEM

/* safe memory management with leakage detection in debug mode */

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

void*  allocMemory_call(size_t size, const char *filename, int line);
void*  reallocMemory_call(void* ptr, size_t size, const char *filename, int line);
void   clearMemory_call(void* ptr, size_t size);
void   copyMemory_call(void* ptr, const void* source, size_t size);
void*  duplicateMemory_call(const void* source, size_t size, const char *filename, int line);
void   freeMemory_call(void** ptr, const char *filename, int line);

#ifndef NDEBUG
void   memoryDiagnostic(void);
void   memoryCheckEmpty(void);
size_t memorySize(void *ptr);
#endif


#else

/* standard memory management via malloc */

#include <string.h>

#define allocMemory(ptr)                   (*(void**)(ptr) = malloc( sizeof(**(ptr)) ))
#define allocMemoryArray(ptr,num)          (*(void**)(ptr) = malloc( (num)*sizeof(**(ptr)) ))
#define allocMemorySize(ptr,size)          (*(void**)(ptr) = malloc( (size_t)(size) ))
#define reallocMemoryArray(ptr,num)        (*(void**)(ptr) = realloc( *(ptr), (num)*sizeof(**(ptr)) ))
#define reallocMemorySize(ptr,size)        (*(void**)(ptr) = realloc( *(ptr), (size_t)(size) ))
#define copyMemory(ptr, source)            memcpy( (void*)(ptr), (const void*)(source), sizeof(*(ptr)) )
#define copyMemoryArray(ptr, source, num)  memcpy( (void*)(ptr), (const void*)(source), (num)*sizeof(*(ptr)) )
#define copyMemorySize(ptr, source, size)  memcpy( (void*)(ptr), (const void*)(source), (size_t)(size) )
#define duplicateMemory(ptr, source)       (*(void**)(ptr) = duplicateMemory_call( (const void*)(source), sizeof(**(ptr))))
#define duplicateMemoryArray(ptr, source, num) (*(void**)(ptr) = duplicateMemory_call( (const void*)(source), \
                                               (num)*sizeof(**(ptr))))
#define duplicateMemorySize(ptr, source, size) (*(void**)(ptr) = duplicateMemory_call( (const void*)(source), \
                                               (size_t)(size)))
#define freeMemory(ptr)                    { free(*(ptr)); *(ptr) = NULL; }
#define freeMemoryNull(ptr)                { if( *(ptr) != NULL ) freeMemory( (ptr) ); }
#define freeMemoryArray(ptr)               freeMemory(ptr)
#define freeMemoryArrayNull(ptr)           { if( *(ptr) != NULL ) freeMemoryArray( (ptr) ); }

void *  duplicateMemory_call(const void* source, size_t size);

#ifndef NDEBUG
#define memoryDiagnostic()                 /**/
#define memoryCheckEmpty()                 /**/
#define memorySize(ptr)                    0
#endif


#endif




/*
 * block memory management
 */

typedef struct memory_header MEMHDR;


#ifndef NOBLOCKMEM

/* block memory methods for faster memory access */

#define createBlockMemory(csz,clr,gbf)     createBlockMemory_call( (csz), (clr), (gbf), __FILE__, __LINE__ )
#define freeAllBlockMemory(mem)            freeAllBlockMemory_call( (mem), __FILE__, __LINE__ )
#define freeAllBlockMemoryNull(mem)        { if( (mem) != NULL ) freeAllBlockMemory( (mem) ); }
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


void   alignMemsize(size_t* size);
int    isAligned(size_t size);

/* Create a block memory allocation structure.
 * Parameters:
 *    initChunkSize     : nr. of elements in the first chunk.
 *    clearUnusedBlocks : immedeately clear a block, if it is unused
 *    garbageFactor     : if at least garbageFactor * avg. chunksize elements
 *                        are free, call garbage collection,
 *                        a value of -1 disables garbage collection
 * Returns   :
 *    Pointer to memory header structure.
 */
MEMHDR* createBlockMemory_call(int initChunkSize, int clearUnusedBlocks, int garbageFactor,
   const char *filename, int line);

/* Free all chunks in the memory allocation structure.
 * Parameters:
 *    mem : Pointer to memory header to clear.
 */
void freeAllBlockMemory_call(MEMHDR* mem, const char *filename, int line);

/* Delete a block allocation.
 * Parameters:
 *    mem : Pointer to memory header to destroy.
 */
void destroyBlockMemory_call(MEMHDR** mem, const char *filename, int line);

/* Get a new block of memory.
 * Parameters:
 *    mem  : Pointer to memory header.
 *    size : size of requested block in bytes.
 * Returns:
 *    Pointer to a new block of memory of size "size".
 */
void* allocBlockMemory_call(MEMHDR* mem, size_t size, const char *filename, int line);

void * reallocBlockMemory_call(MEMHDR *mem, void* ptr, size_t oldsize, size_t newsize, const char *filename, int line);

void * duplicateBlockMemory_call(MEMHDR *mem, const void* source, size_t size, const char *filename, int line);

/* Free a block of memory.
 * Parameters:
 *    mem  : Pointer to memory header.
 *    ptr  : Pointer to memory to free.
 *    size : size of memory block.
 */
void freeBlockMemory_call(MEMHDR* mem, void** ptr, size_t size, const char *filename, int line);

/* get the number of used bytes in block memory */
long long getBlockMemoryUsed(MEMHDR *mem);

size_t blockMemorySize(MEMHDR *mem, void *ptr);
void   blockMemoryDiagnostic(MEMHDR *mem);
void   blockMemoryCheckEmpty(MEMHDR *mem);




#else

/* block memory management mapped to standard memory management */

#define createBlockMemory(csz,clr,gbf)       (void*)(0x01) /* dummy to not return a NULL pointer */
#define clearBlockMemory(mem)                /**/
#define clearBlockMemoryNull(mem)            /**/
#define destroyBlockMemory(mem)              /**/
#define destroyBlockMemoryNull(mem)          /**/
#define freeAllBlockMemory(mem)              /**/
#define freeAllBlockMemoryNull(mem)          /**/

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

#define getBlockMemoryUsed(mem)              0
#define blockMemorySize(mem,ptr)             0
#define blockMemoryDiagnostic(mem)           /**/
#define blockMemoryCheckEmpty(mem)           /**/


#endif




#endif
