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



RETCODE SCIPmemCreate(                  /**< creates block memory structures */
   MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   ALLOC_OKAY( allocMemory(*mem) );

   ALLOC_OKAY( (*mem)->probmem = createBlockMemory(1, TRUE, 10 ) );
   ALLOC_OKAY( (*mem)->solvemem = createBlockMemory(1, FALSE, 10) );

   debugMessage("created probmem  block memory at <%p>\n", (*mem)->probmem);
   debugMessage("created solvemem block memory at <%p>\n", (*mem)->solvemem);

   return SCIP_OKAY;
}

RETCODE SCIPmemFree(                    /**< frees block memory structures */
   MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   destroyBlockMemory((*mem)->probmem);
   destroyBlockMemory((*mem)->solvemem);

   freeMemory(*mem);

   return SCIP_OKAY;
}

