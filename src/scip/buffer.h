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
#pragma ident "@(#) $Id: buffer.h,v 1.5 2003/11/21 10:35:32 bzfpfend Exp $"

/**@file   buffer.h
 * @brief  memory buffer for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __BUFFER_H__
#define __BUFFER_H__


typedef struct Buffer BUFFER;


#include "def.h"
#include "retcode.h"
#include "set.h"


/** memory buffer storage for temporary objects */
struct Buffer
{
   void**           data;               /**< allocated memory chunks for arbitrary data */
   int*             size;               /**< sizes of buffers in bytes */
   Bool*            used;               /**< TRUE iff corresponding buffer is in use */
   int              ndata;              /**< number of memory chunks */
   int              firstfree;          /**< first unused memory chunk */
};



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
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   );

/** allocates the next unused buffer and copies the given memory into the buffer */
extern
RETCODE SCIPbufferDuplicateMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   void*            source,             /**< memory block to copy into the buffer */
   int              size                /**< minimal required size of the buffer */
   );

/** reallocates the buffer to at least the given size */
extern
RETCODE SCIPbufferReallocMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   );

/** frees a buffer */
extern
void SCIPbufferFreeMem(
   BUFFER*          buffer,             /**< memory buffer storage */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              dummysize           /**< used to get a safer define for SCIPsetReleaseBufferSize/Array */
   );

#endif
