/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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



extern
RETCODE SCIPbufferCreate(               /**< creates memory buffer storage */
   BUFFER**         buffer              /**< pointer to memory buffer */
   );

extern
void SCIPbufferFree(                    /**< frees memory buffer storage */
   BUFFER**         buffer              /**< pointer to memory buffer */
   );

extern
RETCODE SCIPbufferCapture(              /**< allocates the next unused buffer */
   BUFFER*          buffer,             /**< memory buffer storage */
   const SET*       set,                /**< global SCIP settings */
   void**           ptr,                /**< pointer to store the allocated memory buffer */
   int              size                /**< minimal required size of the buffer */
   );

extern
void SCIPbufferRelease(                 /**< releases a buffer */
   BUFFER*          buffer,             /**< memory buffer storage */
   void**           ptr,                /**< pointer to the allocated memory buffer */
   int              dummysize           /**< used to get a safer define for SCIPsetReleaseBufferSize/Array */
   );

#endif
