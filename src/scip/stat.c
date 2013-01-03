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

/**@file   stat.c
 * @brief  methods for problem statistics
 * @author Tobias Achterberg
 * @author Stefan Heinz
 * @author Gregor Hendel
 * @author Gerald Gamrath
 * @author Marc Pfetsch
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/vbc.h"
#include "scip/mem.h"
#include "scip/history.h"



/** creates problem statistics data */
SCIP_RETCODE SCIPstatCreate(
   SCIP_STAT**           stat,               /**< pointer to problem statistics data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   assert(stat != NULL);
   assert(set != NULL);

   SCIP_ALLOC( BMSallocMemory(stat) );

   SCIP_CALL( SCIPclockCreate(&(*stat)->solvingtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->presolvingtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->primallptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->duallptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->lexduallptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->barrierlptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->divinglptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->strongbranchtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->conflictlptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->lpsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->pseudosoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->nodeactivationtime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->nlpsoltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*stat)->copyclock, SCIP_CLOCKTYPE_DEFAULT) );

   SCIP_CALL( SCIPhistoryCreate(&(*stat)->glbhistory, blkmem) );
   SCIP_CALL( SCIPhistoryCreate(&(*stat)->glbhistorycrun, blkmem) );
   SCIP_CALL( SCIPvbcCreate(&(*stat)->vbc, messagehdlr) );

   (*stat)->status = SCIP_STATUS_UNKNOWN;
   (*stat)->marked_nvaridx = 0;
   (*stat)->marked_ncolidx = 0;
   (*stat)->marked_nrowidx = 0;
   (*stat)->userinterrupt = FALSE;
   (*stat)->userrestart = FALSE;
   (*stat)->inrestart = FALSE;
   (*stat)->collectvarhistory = TRUE;
   (*stat)->subscipdepth = 0;

   SCIPstatReset(*stat);

   return SCIP_OKAY;
}

/** frees problem statistics data */
SCIP_RETCODE SCIPstatFree(
   SCIP_STAT**           stat,               /**< pointer to problem statistics data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(stat != NULL);
   assert(*stat != NULL);

   SCIPclockFree(&(*stat)->solvingtime);
   SCIPclockFree(&(*stat)->presolvingtime);
   SCIPclockFree(&(*stat)->primallptime);
   SCIPclockFree(&(*stat)->duallptime);
   SCIPclockFree(&(*stat)->lexduallptime);
   SCIPclockFree(&(*stat)->barrierlptime);
   SCIPclockFree(&(*stat)->divinglptime);
   SCIPclockFree(&(*stat)->strongbranchtime);
   SCIPclockFree(&(*stat)->conflictlptime);
   SCIPclockFree(&(*stat)->lpsoltime);
   SCIPclockFree(&(*stat)->pseudosoltime);
   SCIPclockFree(&(*stat)->nodeactivationtime);
   SCIPclockFree(&(*stat)->nlpsoltime);
   SCIPclockFree(&(*stat)->copyclock);

   SCIPhistoryFree(&(*stat)->glbhistory, blkmem);
   SCIPhistoryFree(&(*stat)->glbhistorycrun, blkmem);
   SCIPvbcFree(&(*stat)->vbc);

   BMSfreeMemory(stat);

   return SCIP_OKAY;
}

/** diables the collection of any statistic for a variable */
void SCIPstatDisableVarHistory(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->collectvarhistory = FALSE;
}

/** enables the collection of statistics for a variable */
void SCIPstatEnableVarHistory(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->collectvarhistory = TRUE;
}

/** marks statistics to be able to reset them when solving process is freed */
void SCIPstatMark(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->marked_nvaridx = stat->nvaridx;
   stat->marked_ncolidx = stat->ncolidx;
   stat->marked_nrowidx = stat->nrowidx;
}

/** reset statistics to the data before solving started */
void SCIPstatReset(
   SCIP_STAT*            stat                /**< problem statistics data */
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
   SCIPclockReset(stat->lexduallptime);
   SCIPclockReset(stat->barrierlptime);
   SCIPclockReset(stat->divinglptime);
   SCIPclockReset(stat->strongbranchtime);
   SCIPclockReset(stat->conflictlptime);
   SCIPclockReset(stat->lpsoltime);
   SCIPclockReset(stat->pseudosoltime);
   SCIPclockReset(stat->nodeactivationtime);
   SCIPclockReset(stat->nlpsoltime);
   SCIPclockReset(stat->copyclock);

   SCIPhistoryReset(stat->glbhistory);

   stat->vsidsweight = 1.0;
   stat->nlpiterations = 0;
   stat->nrootlpiterations = 0;
   stat->nprimallpiterations = 0;
   stat->nduallpiterations = 0;
   stat->nlexduallpiterations = 0;
   stat->nbarrierlpiterations = 0;
   stat->nprimalresolvelpiterations = 0;
   stat->ndualresolvelpiterations = 0;
   stat->nlexdualresolvelpiterations = 0;
   stat->nnodelpiterations = 0;
   stat->ninitlpiterations = 0;
   stat->ndivinglpiterations = 0;
   stat->nsblpiterations = 0;
   stat->nrootsblpiterations = 0;
   stat->nconflictlpiterations = 0;
   stat->ntotalnodes = 0;
   stat->ncreatednodes = 0;
   stat->nlpsolsfound = 0;
   stat->npssolsfound = 0;
   stat->domchgcount = 0;
   stat->nboundchgs = 0;
   stat->nholechgs = 0;
   stat->nprobboundchgs = 0;
   stat->nprobholechgs = 0;
   stat->nruns = 0;
   stat->nconfrestarts = 0;
   stat->nrootboundchgs = 0;
   stat->nrootintfixings = 0;
   stat->prevrunnvars = 0;
   stat->nvaridx = stat->marked_nvaridx;
   stat->ncolidx = stat->marked_ncolidx;
   stat->nrowidx = stat->marked_nrowidx;
   stat->lpcount = 0;
   stat->nlps = 0;
   stat->nrootlps = 0;
   stat->nprimallps = 0;
   stat->nprimalzeroitlps = 0;
   stat->nduallps = 0;
   stat->ndualzeroitlps = 0;
   stat->nlexduallps = 0;
   stat->nbarrierlps = 0;
   stat->nbarrierzeroitlps = 0;
   stat->nprimalresolvelps = 0;
   stat->ndualresolvelps = 0;
   stat->nlexdualresolvelps = 0;
   stat->nnodelps = 0;
   stat->ninitlps = 0;
   stat->ndivinglps = 0;
   stat->nstrongbranchs = 0;
   stat->nrootstrongbranchs = 0;
   stat->nconflictlps = 0;
   stat->nnlps = 0;
   stat->maxtotaldepth = -1;
   stat->nactiveconss = 0;
   stat->nenabledconss = 0;
   stat->solindex = 0;
   stat->memsavemode = FALSE;
   stat->nnodesbeforefirst = -1;
   stat->nrunsbeforefirst = -1;
   stat->firstprimalheur = NULL; 
   stat->firstprimaltime = SCIP_DEFAULT_INFINITY;
   stat->firstprimalbound = SCIP_DEFAULT_INFINITY;
   stat->primalzeroittime = 0.0;
   stat->dualzeroittime = 0.0;
   stat->barrierzeroittime = 0.0;
   stat->maxcopytime = SCIP_REAL_MIN;
   stat->mincopytime = SCIP_REAL_MAX;
   stat->ncopies = 0;
   stat->marked_nvaridx = -1;
   stat->marked_ncolidx = -1;
   stat->marked_nrowidx = -1;

   SCIPstatResetImplications(stat);
   SCIPstatResetPresolving(stat);
}

/** reset implication counter */
void SCIPstatResetImplications(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->nimplications = 0;
}

/** reset presolving and current run specific statistics */
void SCIPstatResetPresolving(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->npresolrounds = 0;
   stat->npresolfixedvars = 0;
   stat->npresolaggrvars = 0;
   stat->npresolchgvartypes = 0;
   stat->npresolchgbds = 0;
   stat->npresoladdholes = 0;
   stat->npresoldelconss = 0;
   stat->npresoladdconss = 0;
   stat->npresolupgdconss = 0;
   stat->npresolchgcoefs = 0;
   stat->npresolchgsides = 0;

   SCIPstatResetCurrentRun(stat);
}

/** reset current branch and bound run specific statistics */
void SCIPstatResetCurrentRun(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->nnodes = 0;
   stat->ncreatednodesrun = 0;
   stat->nactivatednodes = 0;
   stat->ndeactivatednodes = 0;
   stat->nbacktracks = 0;
   stat->ndelayedcutoffs = 0;
   stat->nreprops = 0;
   stat->nrepropboundchgs = 0;
   stat->nrepropcutoffs = 0;
   stat->lastdivenode = 0;
   stat->lastconflictnode = 0;
   stat->bestsolnode = 0;
   stat->rootlowerbound = SCIP_REAL_MIN;
   stat->lastbranchvar = NULL;
   stat->status = SCIP_STATUS_UNKNOWN;
   stat->lastbranchdir = SCIP_BRANCHDIR_DOWNWARDS;
   stat->nrootboundchgsrun = 0;
   stat->nrootintfixingsrun = 0;
   stat->npricerounds = 0;
   stat->nseparounds = 0;
   stat->maxdepth = -1;
   stat->plungedepth = 0;

   SCIPhistoryReset(stat->glbhistorycrun);

   SCIPstatResetDisplay(stat);
}

/** resets display statistics, such that a new header line is displayed before the next display line */
void SCIPstatResetDisplay(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->lastdispnode = 0;
   stat->ndisplines = 0;
}

/** increases LP count, such that all lazy updates depending on the LP are enforced again */
void SCIPstatEnforceLPUpdates(
   SCIP_STAT*            stat                /**< problem statistics data */
   )
{
   assert(stat != NULL);

   stat->lpcount++;
}

/** depending on the current memory usage, switches mode flag to standard or memory saving mode */
void SCIPstatUpdateMemsaveMode(
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_MEM*             mem                 /**< block memory pools */
   )
{
   assert(stat != NULL);
   assert(set != NULL);

   if( SCIPsetIsLT(set, set->mem_savefac, 1.0) )
   {
      SCIP_Longint memused;

      memused = SCIPmemGetUsed(mem);
      if( !stat->memsavemode && memused >= set->mem_savefac * set->limit_memory * 1024.0 * 1024.0 )
      {
         /* switch to memory saving mode */
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(node %"SCIP_LONGINT_FORMAT") switching to memory saving mode (mem: %.1fM/%.1fM)\n",
            stat->nnodes, (SCIP_Real)memused/(1024.0*1024.0), set->limit_memory);
         stat->memsavemode = TRUE;
         set->nodesel = NULL;
      }
      else if( stat->memsavemode && memused < 0.5 * set->mem_savefac * set->limit_memory * 1024.0 * 1024.0 )
      {
         /* switch to standard mode */
         SCIPmessagePrintVerbInfo(messagehdlr, set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(node %"SCIP_LONGINT_FORMAT") switching to standard mode (mem: %.1fM/%.1fM)\n",
            stat->nnodes, (SCIP_Real)memused/(1024.0*1024.0), set->limit_memory);
         stat->memsavemode = FALSE;
         set->nodesel = NULL;
      }
   }
   else
      stat->memsavemode = FALSE;
}
