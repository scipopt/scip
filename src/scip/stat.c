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

/**@file   stat.c
 * @brief  problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "memory.h"
#include "stat.h"


/** creates problem statistics data */
RETCODE SCIPstatCreate(
   STAT**           stat                /**< pointer to problem statistics data */
   )
{
   assert(stat != NULL);

   ALLOC_OKAY( allocMemory(stat) );
   (*stat)->marked_nvaridx = 0;
   (*stat)->marked_ncolidx = 0;
   (*stat)->marked_nrowidx = 0;
   SCIPstatReset(*stat);

   return SCIP_OKAY;
}

/** frees problem statistics data */
RETCODE SCIPstatFree(
   STAT**           stat                /**< pointer to problem statistics data */
   )
{
   assert(stat != NULL);
   assert(*stat != NULL);

   freeMemory(stat);

   return SCIP_OKAY;
}

/** marks statistics to be able to reset them when solving process is freed */
void SCIPstatMark(
   STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);
   assert(stat->marked_nvaridx == -1);
   assert(stat->marked_ncolidx == -1);
   assert(stat->marked_nrowidx == -1);
   assert(stat->lpcount == 0);
   assert(stat->nlps == 0);
   assert(stat->nprimallps == 0);
   assert(stat->nduallps == 0);

   stat->marked_nvaridx = stat->nvaridx;
   stat->marked_ncolidx = stat->ncolidx;
   stat->marked_nrowidx = stat->nrowidx;
}

/** reset statistics to the data before solving started */
void SCIPstatReset(
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
   stat->lpcount = 0;
   stat->nlps = 0;
   stat->nprimallps = 0;
   stat->nduallps = 0;
   stat->nlpiterations = 0;
   stat->nprimallpiterations = 0;
   stat->nduallpiterations = 0;
   stat->nstrongbranch = 0;
   stat->nseparounds = 0;
   stat->nnodes = 0;
   stat->nboundchanges = 0;
   stat->lastdispnode = 0;
   stat->ndisplines = 0;
   stat->maxdepth = -1;
   stat->plungedepth = 0;

   stat->marked_nvaridx = -1;
   stat->marked_ncolidx = -1;
   stat->marked_nrowidx = -1;
}
