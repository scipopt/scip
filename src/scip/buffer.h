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
#pragma ident "@(#) $Id: buffer.h,v 1.13 2005/02/14 13:35:39 bzfpfend Exp $"

/**@file   buffer.h
 * @brief  internal methods for memory buffers for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BUFFER_H__
#define __BUFFER_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_buffer.h"


/** creates memory buffer storage */
extern
RETCODE SCIPbufferCreate(
   BUFFER**         buffer              /**< pointer to memory buffer */
   );

/** frees memory buffer storage */
extern
void SCIPbufferFree(
   BUFFER**         buffer              /**< pointer to memory buffer */
   );

/** allocates the next unused buffer */
extern
RETCODE SCIPbufferAllocMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   SET*             set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   );

/** allocates the next unused buffer and copies the given memory into the buffer */
extern
RETCODE SCIPbufferDuplicateMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   SET*             set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   void*            source,             /**< memory block to copy into the buffer */
   int              size                /**< minimal required size of the buffer */
   );

/** reallocates the buffer to at least the given size */
extern
RETCODE SCIPbufferReallocMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   SET*             set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   );

/** frees a buffer */
extern
void SCIPbufferFreeMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              dummysize           /**< used to get a safer define for SCIPsetFreeBufferSize/Array */
   );

/** gets number of used buffers */
extern
int SCIPbufferGetNUsed(
   BUFFER*          buffer              /**< memory buffer storage */
   );

#endif
