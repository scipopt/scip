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

/**@file   mem.h
 * @brief  block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEM_H__
#define __MEM_H__


typedef struct Mem MEM;


#include "def.h"
#include "memory.h"
#include "retcode.h"
#include "set.h"


/** various block memory buffers */
struct Mem
{
   MEMHDR*          treemem;            /**< ptr to memory blocks for the tree */
   MEMHDR*          statemem;           /**< ptr to memory blocks for LP states */
   MEMHDR*          lpmem;              /**< ptr to memory blocks for LP data */
   MEMHDR*          dommem;             /**< ptr to memory blocks for domains of variables */
   MEMHDR*          consmem;            /**< ptr to memory blocks for constraint data */
   MEMHDR*          primalmem;          /**< ptr to memory blocks for primal solutions */
   MEMHDR*          tempmem;            /**< ptr to memory blocks for short living objects */
   void**           ptrbuf;             /**< buffer for storing temporary pointer arrays */
   char*            charbuf;            /**< buffer for storing temporary char arrays */
   int*             intbuf;             /**< buffer for storing temporary int arrays */
   Real*            realbuf;            /**< buffer for storing temporary real arrays */
   int              ptrbufsize;         /**< size of pointer array buffer */
   int              charbufsize;        /**< size of char array buffer */
   int              intbufsize;         /**< size of int array buffer */
   int              realbufsize;        /**< size of real array buffer */
};



extern
RETCODE SCIPmemCreate(                  /**< creates block memory structures */
   MEM**            mem                 /**< pointer to block memory structure */
   );

extern
RETCODE SCIPmemGetPtrbuf(               /**< returns buffer for storing pointer array */
   void***          ptrbuf,             /**< pointer to a pointer array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of pointer buffer */
   );

extern
RETCODE SCIPmemGetCharbuf(              /**< returns buffer for storing char array */
   char**           charbuf,            /**< pointer to char array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of char buffer */
   );

extern
RETCODE SCIPmemGetIntbuf(               /**< returns buffer for storing int array */
   int**            intbuf,             /**< pointer to int array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of int buffer */
   );

extern
RETCODE SCIPmemGetRealbuf(              /**< returns buffer for storing Real array */
   Real**           realbuf,            /**< pointer to Real array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of real buffer */
   );

#endif
