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



/*
 * dynamically sized arrays
 */

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



/*
 * external methods
 */

RETCODE SCIPmemCreate(                  /**< creates block memory structures */
   MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   ALLOC_OKAY( allocMemory(*mem) );

   ALLOC_OKAY( (*mem)->probmem = createBlockMemory(1, TRUE, 10) );
   ALLOC_OKAY( (*mem)->solvemem = createBlockMemory(1, FALSE, 10) );
#if 0 /* ??? */
   ALLOC_OKAY( (*mem)->treemem = createBlockMemory(1, FALSE, 10) );
   ALLOC_OKAY( (*mem)->statemem = createBlockMemory(1, FALSE, 10) );
   ALLOC_OKAY( (*mem)->lpmem = createBlockMemory(1, FALSE, 10) );
   ALLOC_OKAY( (*mem)->dommem = createBlockMemory(1, FALSE, 10) );
   ALLOC_OKAY( (*mem)->consmem = createBlockMemory(1, FALSE, 10) );
   ALLOC_OKAY( (*mem)->primalmem = createBlockMemory(1, FALSE, 10) );
#endif
   ALLOC_OKAY( (*mem)->tempmem = createBlockMemory(1, FALSE, 10) );
   (*mem)->ptrbuf = NULL;
   (*mem)->charbuf = NULL;
   (*mem)->intbuf = NULL;
   (*mem)->realbuf = NULL;
   (*mem)->ptrbufsize = 0;
   (*mem)->charbufsize = 0;
   (*mem)->intbufsize = 0;
   (*mem)->realbufsize = 0;

   return SCIP_OKAY;
}

RETCODE SCIPmemFree(                    /**< frees block memory structures */
   MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   destroyBlockMemory((*mem)->probmem);
   destroyBlockMemory((*mem)->solvemem);
#if 0 /* ??? */
   destroyBlockMemory((*mem)->treemem);
   destroyBlockMemory((*mem)->statemem);
   destroyBlockMemory((*mem)->lpmem);
   destroyBlockMemory((*mem)->dommem);
   destroyBlockMemory((*mem)->consmem);
   destroyBlockMemory((*mem)->primalmem);
#endif
   destroyBlockMemory((*mem)->tempmem);

   freeMemory(*mem);

   return SCIP_OKAY;
}
   
RETCODE SCIPmemGetPtrbuf(               /**< returns buffer for storing pointer array */
   void***          ptrbuf,             /**< pointer to a pointer array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of pointer buffer */
   )
{
   assert(ptrbuf != NULL);

   CHECK_OKAY( memEnsurePtrbufSize(mem, set, size) );
   *ptrbuf = mem->ptrbuf;

   return SCIP_OKAY;
}

RETCODE SCIPmemGetCharbuf(              /**< returns buffer for storing char array */
   char**           charbuf,            /**< pointer to char array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of char buffer */
   )
{
   assert(charbuf != NULL);

   CHECK_OKAY( memEnsureCharbufSize(mem, set, size) );
   *charbuf = mem->charbuf;

   return SCIP_OKAY;
}

RETCODE SCIPmemGetIntbuf(               /**< returns buffer for storing int array */
   int**            intbuf,             /**< pointer to int array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of int buffer */
   )
{
   assert(intbuf != NULL);

   CHECK_OKAY( memEnsureIntbufSize(mem, set, size) );
   *intbuf = mem->intbuf;

   return SCIP_OKAY;
}

RETCODE SCIPmemGetRealbuf(              /**< returns buffer for storing Real array */
   Real**           realbuf,            /**< pointer to Real array */
   MEM*             mem,                /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   int              size                /**< minimal size of real buffer */
   )
{
   assert(realbuf != NULL);

   CHECK_OKAY( memEnsureRealbufSize(mem, set, size) );
   *realbuf = mem->realbuf;

   return SCIP_OKAY;
}
