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
   STAT**           stat,               /**< pointer to problem statistics data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(stat != NULL);
   assert(set != NULL);

   ALLOC_OKAY( allocMemory(stat) );

   CHECK_OKAY( SCIPclockCreate(&(*stat)->solvingtime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->presolvingtime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->primallptime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->duallptime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->strongbranchtime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->lppricingtime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->lpsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->pseudosoltime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->redcoststrtime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*stat)->nodeactivationtime, SCIP_CLOCKTYPE_DEFAULT) );

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

   SCIPclockFree(&(*stat)->solvingtime);
   SCIPclockFree(&(*stat)->presolvingtime);
   SCIPclockFree(&(*stat)->primallptime);
   SCIPclockFree(&(*stat)->duallptime);
   SCIPclockFree(&(*stat)->strongbranchtime);
   SCIPclockFree(&(*stat)->lppricingtime);
   SCIPclockFree(&(*stat)->lpsoltime);
   SCIPclockFree(&(*stat)->pseudosoltime);
   SCIPclockFree(&(*stat)->redcoststrtime);
   SCIPclockFree(&(*stat)->nodeactivationtime);

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
   
   SCIPclockReset(stat->solvingtime);
   SCIPclockReset(stat->presolvingtime);
   SCIPclockReset(stat->primallptime);
   SCIPclockReset(stat->duallptime);
   SCIPclockReset(stat->strongbranchtime);
   SCIPclockReset(stat->lppricingtime);
   SCIPclockReset(stat->lpsoltime);
   SCIPclockReset(stat->pseudosoltime);
   SCIPclockReset(stat->redcoststrtime);
   SCIPclockReset(stat->nodeactivationtime);

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
   stat->ndivinglpiterations = 0;
   stat->nsblpiterations = 0;
   stat->nstrongbranch = 0;
   stat->nseparounds = 0;
   stat->nlppricings = 0;
   stat->nlppricingvars = 0;
   stat->nredcoststrcalls = 0;
   stat->nredcoststrfound = 0;
   stat->nnodes = 0;
   stat->nboundchanges = 0;
   stat->nlpsolsfound = 0;
   stat->npssolsfound = 0;
   stat->lastdispnode = 0;
   stat->ndisplines = 0;
   stat->maxdepth = -1;
   stat->plungedepth = 0;
   stat->memsavemode = FALSE;

   stat->marked_nvaridx = -1;
   stat->marked_ncolidx = -1;
   stat->marked_nrowidx = -1;
}

/** resets display statistics, such that a new header line is displayed before the next display line */
void SCIPstatResetDisplay(
   STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->lastdispnode = 0;
   stat->ndisplines = 0;
}

/** depending on the current memory usage, switches mode flag to standard or memory saving mode */
void SCIPstatUpdateMemsaveMode(
   STAT*            stat,               /**< problem statistics data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(stat != NULL);
   assert(set != NULL);

   if( set->memlimit >= 0 && SCIPsetIsLT(set, set->memsavefac, 1.0) )
   {
      Longint memused;
      char s[MAXSTRLEN];

      memused = SCIPgetMemUsed(set->scip);
      if( !stat->memsavemode && memused >= set->memsavefac * set->memlimit )
      {
         /* switch to memory saving mode */
         sprintf(s, "(node %lld) switching to memory saving mode (mem: %lld/%lld)", 
            stat->nnodes, memused, set->memlimit);
         infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL, s);
         stat->memsavemode = TRUE;
      }
      else if( stat->memsavemode && memused < 0.5 * set->memsavefac * set->memlimit )
      {
         /* switch to standard mode */
         sprintf(s, "(node %lld) switching to standard mode (mem: %lld/%lld)", 
            stat->nnodes, memused, set->memlimit);
         infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL, s);
         stat->memsavemode = FALSE;
      }
   }
   else
      stat->memsavemode = FALSE;
}
