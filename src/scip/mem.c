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

/**@file   mem.c
 * @brief  block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/mem.h"
#include "scip/pub_message.h"



/** creates block memory structures */
SCIP_RETCODE SCIPmemCreate(
   SCIP_MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   SCIP_ALLOC( BMSallocMemory(mem) );

   SCIP_ALLOC( (*mem)->setmem = BMScreateBlockMemory(1, 10) );
   SCIP_ALLOC( (*mem)->probmem = BMScreateBlockMemory(1, 10) );

   SCIPdebugMessage("created setmem   block memory at <%p>\n", (void*)(*mem)->setmem);
   SCIPdebugMessage("created probmem  block memory at <%p>\n", (void*)(*mem)->probmem);

   return SCIP_OKAY;
}

/** frees block memory structures */
SCIP_RETCODE SCIPmemFree(
   SCIP_MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   BMSdestroyBlockMemory(&(*mem)->probmem);
   BMSdestroyBlockMemory(&(*mem)->setmem);

   BMSfreeMemory(mem);

   return SCIP_OKAY;
}

/** returns the total number of bytes used in block memory */
SCIP_Longint SCIPmemGetUsed(
   SCIP_MEM*             mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   return BMSgetBlockMemoryUsed(mem->setmem) + BMSgetBlockMemoryUsed(mem->probmem);
}
