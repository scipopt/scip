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

/**@file   stat.c
 * @brief  problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "memory.h"
#include "stat.h"


RETCODE SCIPstatCreate(                 /**< creates problem statistics data */
   STAT**           stat                /**< pointer to problem statistics data */
   )
{
   assert(stat != NULL);

   ALLOC_OKAY( allocMemory(*stat) );
   (*stat)->numvaridx = 0;
   (*stat)->numcolidx = 0;
   (*stat)->numrowidx = 0;

   return SCIP_OKAY;
}

RETCODE SCIPstatFree(                   /**< frees problem statistics data */
   STAT**           stat                /**< pointer to problem statistics data */
   )
{
   assert(stat != NULL);
   assert(*stat != NULL);

   freeMemory(*stat);

   return SCIP_OKAY;
}

