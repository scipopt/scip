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

/**@file   mem.c
 * @brief  block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "mem.h"


static
RETCODE memEnsurePtrbufSize(            /**< ensures, that pointer array buffer can store at least num entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(mem != NULL);
   
   if( num > mem->ptrbufsize )
   {
      int newsize;

      newsize = SCIPcalcBufGrowSize(set, num);
      freeBlockMemoryArrayNull(mem->tempmem, mem->ptrbuf, mem->ptrbufsize);
      ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, mem->ptrbuf, newsize) );
      mem->ptrbufsize = newsize;
   }
   assert(num <= mem->ptrbufsize);

   return SCIP_OKAY;
}

static
RETCODE memEnsureCharbufSize(           /**< ensures, that char array buffer can store at least num entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(mem != NULL);
   
   if( num > mem->charbufsize )
   {
      int newsize;

      newsize = SCIPcalcBufGrowSize(set, num);
      freeBlockMemoryArrayNull(mem->tempmem, mem->charbuf, mem->charbufsize);
      ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, mem->charbuf, newsize) );
      mem->charbufsize = newsize;
   }
   assert(num <= mem->charbufsize);

   return SCIP_OKAY;
}

static
RETCODE memEnsureIntbufSize(            /**< ensures, that int array buffer can store at least num entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(mem != NULL);
   
   if( num > mem->intbufsize )
   {
      int newsize;

      newsize = SCIPcalcBufGrowSize(set, num);
      freeBlockMemoryArrayNull(mem->tempmem, mem->intbuf, mem->intbufsize);
      ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, mem->intbuf, newsize) );
      mem->intbufsize = newsize;
   }
   assert(num <= mem->intbufsize);

   return SCIP_OKAY;
}

static
RETCODE memEnsureRealbufSize(           /**< ensures, that real array buffer can store at least num entries */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(mem != NULL);
   
   if( num > mem->realbufsize )
   {
      int newsize;

      newsize = SCIPcalcBufGrowSize(set, num);
      freeBlockMemoryArrayNull(mem->tempmem, mem->realbuf, mem->realbufsize);
      ALLOC_OKAY( allocBlockMemoryArray(mem->tempmem, mem->realbuf, newsize) );
      mem->realbufsize = newsize;
   }
   assert(num <= mem->realbufsize);

   return SCIP_OKAY;
}

void** SCIPmemGetPtrbuf(                /**< returns buffer for storing pointer array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of pointer buffer */
   )
{
   CHECK_NULL( memEnsurePtrbufSize(mem, set, size) );
   return mem->ptrbuf;
}

char* SCIPmemGetCharbuf(                /**< returns buffer for storing char array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of char buffer */
   )
{
   CHECK_NULL( memEnsureCharbufSize(mem, set, size) );
   return mem->charbuf;
}

int* SCIPmemGetIntbuf(                  /**< returns buffer for storing int array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of int buffer */
   )
{
   CHECK_NULL( memEnsureIntbufSize(mem, set, size) );
   return mem->intbuf;
}

Real* SCIPmemGetRealbuf(                /**< returns buffer for storing Real array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of real buffer */
   )
{
   CHECK_NULL( memEnsureRealbufSize(mem, set, size) );
   return mem->realbuf;
}
