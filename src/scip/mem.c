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
#pragma ident "@(#) $Id: mem.c,v 1.12 2003/11/21 10:35:37 bzfpfend Exp $"

/**@file   mem.c
 * @brief  block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "mem.h"



/** creates block memory structures */
RETCODE SCIPmemCreate(
   MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   ALLOC_OKAY( allocMemory(mem) );

   ALLOC_OKAY( (*mem)->setmem = createBlockMemory(1, TRUE, 10 ) );
   ALLOC_OKAY( (*mem)->probmem = createBlockMemory(1, TRUE, 10 ) );
   ALLOC_OKAY( (*mem)->solvemem = createBlockMemory(1, FALSE, 10) );

   debugMessage("created setmem   block memory at <%p>\n", (*mem)->setmem);
   debugMessage("created probmem  block memory at <%p>\n", (*mem)->probmem);
   debugMessage("created solvemem block memory at <%p>\n", (*mem)->solvemem);

   return SCIP_OKAY;
}

/** frees block memory structures */
RETCODE SCIPmemFree(
   MEM**            mem                 /**< pointer to block memory structure */
   )
{
   assert(mem != NULL);

   destroyBlockMemory(&(*mem)->solvemem);
   destroyBlockMemory(&(*mem)->probmem);
   destroyBlockMemory(&(*mem)->setmem);

   freeMemory(mem);

   return SCIP_OKAY;
}

