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

/**@file   memory.h
 * @brief  memory allocation routines
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEMORY_H__
#define __MEMORY_H__

#include <stdlib.h>
#include "def.h"

typedef struct memory_header MEMHDR;

#define errorMessage(msg)                  errorMessage_call( (msg), __FILE__, __LINE__ )
#define allocMemory(ptr)                   ((ptr) = allocMemory_call( sizeof(*(ptr)), __FILE__, __LINE__ ))
#define allocMemoryArray(ptr,num)          ((ptr) = allocMemory_call( (num)*sizeof(*(ptr)), __FILE__, __LINE__ ))
#define allocMemorySize(ptr,size)          ((ptr) = allocMemory_call( (size_t)(size), __FILE__, __LINE__ ))
#define allocMemoryCPP(size)               allocMemory_call( (size_t)(size), __FILE__, __LINE__ )
#define allocMemoryArrayCPP(num,size)      allocMemory_call( (size_t)((num)*(size)), __FILE__, __LINE__ )
#define reallocMemoryArray(ptr,num)        reallocMemory_call( (ptr), (num)*sizeof(*(ptr)), __FILE__, __LINE__ )
#define reallocMemorySize(ptr,size)        reallocMemory_call( (ptr), (size_t)(size), __FILE__, __LINE__ )
#define copyMemory(ptr, source)            copyMemory_call( (void*)(ptr), (const void*)source, sizeof(*(ptr)) )
#define copyMemoryArray(ptr, source, num)  copyMemory_call( (void*)(ptr), (const void*)source, (num)*sizeof(*(ptr)) )
#define copyMemorySize(ptr, source, size)  copyMemory_call( (void*)(ptr), (const void*)source, (size_t)(size) )
#define duplicateMemory(ptr, source)       ((ptr) = duplicateMemory_call( (const void*)source, sizeof(*(ptr)), \
                                           __FILE__, __LINE__ ))
#define duplicateMemoryArray(ptr, source, num) \
                                           ((ptr) = duplicateMemory_call( (const void*)source, (num)*sizeof(*(ptr)), \
                                           __FILE__, __LINE__ ))
#define duplicateMemorySize(ptr, source, size) \
                                           ((ptr) = duplicateMemory_call( (const void*)source, (size_t)(size), \
                                           __FILE__, __LINE__ ))
#define freeMemory(ptr)                    freeMemory_call( (void**)(&(ptr)), __FILE__, __LINE__ )
#define freeMemoryNull(ptr)                if( (ptr) != NULL ) freeMemory( (ptr) )
#define freeMemoryArray(ptr)               freeMemory_call( (void**)(&(ptr)), __FILE__, __LINE__ )
#define freeMemoryArrayNull(ptr)           if( (ptr) != NULL ) freeMemoryArray( (ptr) )

#define createBlockMemory(csz,clr,gbf)     createBlockMemory_call( (csz), (clr), (gbf), __FILE__, __LINE__ )
#define clearBlockMemory(mem)              clearBlockMemory_call( (mem), __FILE__, __LINE__ )
#define destroyBlockMemory(mem)            destroyBlockMemory_call( (MEMHDR**)(&(mem)), __FILE__, __LINE__ )
#define destroyBlockMemoryNull(mem)        if( (mem) != NULL ) destroyBlockMemory( (mem) )

#define allocBlockMemory(mem,ptr)          ((ptr) = allocBlockMemory_call((mem),sizeof(*(ptr)),__FILE__,__LINE__))
#define allocBlockMemoryArray(mem,ptr,num) ((ptr) = allocBlockMemory_call((mem),(num)*sizeof(*(ptr)),__FILE__,__LINE__))
#define duplicateBlockMemory(mem, ptr, source) \
                                           ((ptr) = duplicateBlockMemory_call((mem), (const void*)source, \
                                           sizeof(*(ptr)), __FILE__, __LINE__ ))
#define duplicateBlockMemoryArray(mem, ptr, source, num) \
                                           ((ptr) = duplicateBlockMemory_call( (mem), (const void*)source, \
                                           (num)*sizeof(*(ptr)), __FILE__, __LINE__ ))
#define freeBlockMemory(mem,ptr)           freeBlockMemory_call( (mem), (void**)(&(ptr)), sizeof(*(ptr)), \
                                           __FILE__, __LINE__ )
#define freeBlockMemoryNull(mem,ptr)       if( (ptr) != NULL ) freeBlockMemory( (mem), (ptr) )
#define freeBlockMemoryArray(mem,ptr,num)  freeBlockMemory_call( (mem), (void**)(&(ptr)), (num)*sizeof(*(ptr)), \
                                           __FILE__, __LINE__ )
#define freeBlockMemoryArrayNull(mem,ptr,num)  if( (ptr) != NULL ) freeBlockMemoryArray( (mem), (ptr), (num) )


#ifndef NDEBUG
void   memoryDiagnostic( void );
void   memoryCheckEmpty( void );
size_t memorySize( void *ptr );
size_t blockMemorySize( MEMHDR *mem, void *ptr );
void   blockMemoryDiagnostic( MEMHDR *mem );
#endif

void   errorMessage_call( const char *msg, const char *filename, int line );
void*  allocMemory_call( size_t size, const char *filename, int line );
void*  reallocMemory_call( void* ptr, size_t size, const char *filename, int line );
void   copyMemory_call(void* ptr, const void* source, size_t size);
void*  duplicateMemory_call(const void* source, size_t size, const char *filename, int line);
void   freeMemory_call( void** ptr, const char *filename, int line );

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
MEMHDR* createBlockMemory_call( int initChunkSize, int clearUnusedBlocks,
				int garbageFactor,
				const char *filename, int line );

/* Free all chunks in the memory allocation structure.
 * Parameters:
 *    mem : Pointer to memory header to clear.
 */
void clearBlockMemory_call( MEMHDR* mem, 
			    const char *filename, int line );

/* Delete a block allocation.
 * Parameters:
 *    mem : Pointer to memory header to destroy.
 */
void destroyBlockMemory_call( MEMHDR** mem,
			      const char *filename, int line );

/* Get a new block of memory.
 * Parameters:
 *    mem  : Pointer to memory header.
 *    size : size of requested block in bytes.
 * Returns:
 *    Pointer to a new block of memory of size "size".
 */
void* allocBlockMemory_call( MEMHDR* mem, size_t size,
			     const char *filename, int line );

void * duplicateBlockMemory_call(MEMHDR *mem, const void* source, size_t size, const char *filename, int line);

/* Free a block of memory.
 * Parameters:
 *    mem  : Pointer to memory header.
 *    ptr  : Pointer to memory to free.
 *    size : size of memory block.
 */
void freeBlockMemory_call( MEMHDR* mem, void** ptr, size_t size,
			   const char *filename, int line );



#endif /* __BLKALLOC_H__ */
