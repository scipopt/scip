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
   (*stat)->marked_nvaridx = 0;
   (*stat)->marked_ncolidx = 0;
   (*stat)->marked_nrowidx = 0;
   SCIPstatReset(*stat);

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

void SCIPstatMark(                      /**< marks statistics to be able to reset them when solving process is freed */
   STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);
   assert(stat->marked_nvaridx == -1);
   assert(stat->marked_ncolidx == -1);
   assert(stat->marked_nrowidx == -1);
   assert(stat->nlp == 0);
   assert(stat->nprimallp == 0);
   assert(stat->nduallp == 0);

   stat->marked_nvaridx = stat->nvaridx;
   stat->marked_ncolidx = stat->ncolidx;
   stat->marked_nrowidx = stat->nrowidx;
}

void SCIPstatReset(                     /**< reset statistics to the data before solving started */
   STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);
   assert(stat->marked_nvaridx >= 0);
   assert(stat->marked_ncolidx >= 0);
   assert(stat->marked_nrowidx >= 0);
   
   stat->nvaridx = stat->marked_nvaridx;
   stat->ncolidx = stat->marked_ncolidx;
   stat->nrowidx = stat->marked_nrowidx;
   stat->nlp = 0;
   stat->nprimallp = 0;
   stat->nduallp = 0;
   stat->nnodes = 0;
   stat->lastdispnode = 0;

   stat->marked_nvaridx = -1;
   stat->marked_ncolidx = -1;
   stat->marked_nrowidx = -1;
}
