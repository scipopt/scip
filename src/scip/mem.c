/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   mem.c
 * @brief  block memory pools and memory buffers
 * @author Tobias Achterberg
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/mem.h"
#include "scip/pub_message.h"



/** creates block and buffer memory structures */
SCIP_RETCODE SCIPmemCreate(
   SCIP_MEM**            mem                 /**< pointer to block and buffer memory structure */
   )
{
   assert(mem != NULL);

   SCIP_ALLOC( BMSallocMemory(mem) );

   /* alloc block memory */
   SCIP_ALLOC( (*mem)->setmem = BMScreateBlockMemory(1, 10) );
   SCIP_ALLOC( (*mem)->probmem = BMScreateBlockMemory(1, 10) );

   /* alloc memory buffers */
   SCIP_ALLOC( (*mem)->buffer = BMScreateBufferMemory(SCIP_DEFAULT_MEM_ARRAYGROWFAC, SCIP_DEFAULT_MEM_ARRAYGROWINIT, FALSE) );
   SCIP_ALLOC( (*mem)->cleanbuffer = BMScreateBufferMemory(SCIP_DEFAULT_MEM_ARRAYGROWFAC, SCIP_DEFAULT_MEM_ARRAYGROWINIT, TRUE) );

   SCIPdebugMessage("created setmem   block memory at <%p>\n", (void*)(*mem)->setmem);
   SCIPdebugMessage("created probmem  block memory at <%p>\n", (void*)(*mem)->probmem);

   SCIPdebugMessage("created       buffer memory at <%p>\n", (void*)(*mem)->buffer);
   SCIPdebugMessage("created clean buffer memory at <%p>\n", (void*)(*mem)->cleanbuffer);

   return SCIP_OKAY;
}

/** frees block and buffer memory structures */
SCIP_RETCODE SCIPmemFree(
   SCIP_MEM**            mem                 /**< pointer to block and buffer memory structure */
   )
{
   assert(mem != NULL);

   /* free memory buffers */
   BMSdestroyBufferMemory(&(*mem)->cleanbuffer);
   BMSdestroyBufferMemory(&(*mem)->buffer);

   /* free block memory */
   BMSdestroyBlockMemory(&(*mem)->probmem);
   BMSdestroyBlockMemory(&(*mem)->setmem);

   BMSfreeMemory(mem);

   return SCIP_OKAY;
}

/** returns the total number of bytes used in block and buffer memory */
SCIP_Longint SCIPmemGetUsed(
   SCIP_MEM*             mem                 /**< pointer to block and buffer memory structure */
   )
{
   assert(mem != NULL);

   return BMSgetBlockMemoryUsed(mem->setmem) + BMSgetBlockMemoryUsed(mem->probmem)
      + BMSgetBufferMemoryUsed(mem->buffer) + BMSgetBufferMemoryUsed(mem->cleanbuffer);
}
