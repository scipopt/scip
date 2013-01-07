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

/**@file   buffer.h
 * @brief  internal methods for memory buffers for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BUFFER_H__
#define __SCIP_BUFFER_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_buffer.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates memory buffer storage */
extern
SCIP_RETCODE SCIPbufferCreate(
   SCIP_BUFFER**         buffer              /**< pointer to memory buffer */
   );

/** frees memory buffer storage */
extern
void SCIPbufferFree(
   SCIP_BUFFER**         buffer              /**< pointer to memory buffer */
   );

/** allocates the next unused buffer */
extern
SCIP_RETCODE SCIPbufferAllocMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   int                   size                /**< minimal required size of the buffer */
   );

/** allocates the next unused buffer and copies the given memory into the buffer */
extern
SCIP_RETCODE SCIPbufferDuplicateMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   const void*           source,             /**< memory block to copy into the buffer */
   int                   size                /**< minimal required size of the buffer */
   );

/** reallocates the buffer to at least the given size */
extern
SCIP_RETCODE SCIPbufferReallocMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to the allocated memory buffer */
   int                   size                /**< minimal required size of the buffer */
   );

#ifndef NDEBUG
/** allocates the next unused buffer; checks for integer overflow */
extern
SCIP_RETCODE SCIPbufferAllocMemSave(
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   int                   num,                /**< number of entries to allocate */
   size_t                elemsize            /**< size of one element in the array */
   );

/** allocates the next unused buffer and copies the given memory into the buffer; checks for integer overflows */
extern
SCIP_RETCODE SCIPbufferDuplicateMemSave(
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to store the allocated memory buffer */
   const void*           source,             /**< memory block to copy into the buffer */
   int                   num,                /**< number of entries to copy */
   size_t                elemsize            /**< size of one element in the array */
   );

/** reallocates the buffer to at least the given size; checks for integer overflows */
extern
SCIP_RETCODE SCIPbufferReallocMemSave(
   SCIP_SET*             set,                /**< global SCIP settings */
   void**                ptr,                /**< pointer to the allocated memory buffer */
   int                   num,                /**< number of entries to get memory for */
   size_t                elemsize            /**< size of one element in the array */
   );
#endif

/** frees a buffer */
extern
void SCIPbufferFreeMem(
   SCIP_BUFFER*          buffer,             /**< memory buffer storage */
   void**                ptr,                /**< pointer to the allocated memory buffer */
   int                   dummysize           /**< used to get a safer define for SCIPsetFreeBufferSize/Array */
   );

/** gets number of used buffers */
extern
int SCIPbufferGetNUsed(
   SCIP_BUFFER*          buffer              /**< memory buffer storage */
   );

/** outputs statistics about currently allocated buffers to the screen */
extern
void SCIPbufferPrint(
   SCIP_BUFFER*          buffer              /**< memory buffer storage */
   );

#ifdef __cplusplus
}
#endif

#endif
