/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol.c,v 1.25 2005/02/07 14:08:25 bzfpfend Exp $"

/**@file   presol.c
 * @brief  methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "memory.h"
#include "set.h"
#include "clock.h"
#include "paramset.h"
#include "scip.h"
#include "presol.h"

#include "struct_presol.h"



/*
 * presolver methods
 */

/** compares two presolvers w. r. to their priority */
DECL_SORTPTRCOMP(SCIPpresolComp)
{  /*lint --e{715}*/
   return ((PRESOL*)elem2)->priority - ((PRESOL*)elem1)->priority;
}

/** method to call, when the priority of a presolver was changed */
static
DECL_PARAMCHGD(paramChgdPresolPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPresolPriority() to mark the presols unsorted */
   CHECK_OKAY( SCIPsetPresolPriority(scip, (PRESOL*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a presolver */
RETCODE SCIPpresolCreate(
   PRESOL**         presol,             /**< pointer to store presolver */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of presolver */
   const char*      desc,               /**< description of presolver */
   int              priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int              maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(presol != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   ALLOC_OKAY( allocMemory(presol) );
   ALLOC_OKAY( duplicateMemoryArray(&(*presol)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*presol)->desc, desc, strlen(desc)+1) );
   (*presol)->priority = priority;
   (*presol)->presolfree = presolfree;
   (*presol)->presolinit = presolinit;
   (*presol)->presolexit = presolexit;
   (*presol)->presolinitpre = presolinitpre;
   (*presol)->presolexitpre = presolexitpre;
   (*presol)->presolexec = presolexec;
   (*presol)->presoldata = presoldata;
   CHECK_OKAY( SCIPclockCreate(&(*presol)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*presol)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "presolving/%s/priority", name);
   sprintf(paramdesc, "priority of presolver <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*presol)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdPresolPriority, (PARAMDATA*)(*presol)) ); /*lint !e740*/
   sprintf(paramname, "presolving/%s/maxrounds", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname,
         "maximal number of presolving rounds the presolver participates in (-1: no limit)",
         &(*presol)->maxrounds, maxrounds, -1, INT_MAX, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** frees memory of presolver */   
RETCODE SCIPpresolFree(
   PRESOL**         presol,             /**< pointer to presolver data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(*presol != NULL);
   assert(!(*presol)->initialized);
   assert(set != NULL);

   /* call destructor of presolver */
   if( (*presol)->presolfree != NULL )
   {
      CHECK_OKAY( (*presol)->presolfree(set->scip, *presol) );
   }

   SCIPclockFree(&(*presol)->clock);
   freeMemoryArray(&(*presol)->name);
   freeMemoryArray(&(*presol)->desc);
   freeMemory(presol);

   return SCIP_OKAY;
}

/** initializes presolver */
RETCODE SCIPpresolInit(
   PRESOL*          presol,             /**< presolver */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   if( presol->initialized )
   {
      errorMessage("presolver <%s> already initialized\n", presol->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(presol->clock);

   presol->lastnfixedvars = 0;
   presol->lastnaggrvars = 0;
   presol->lastnchgvartypes = 0;
   presol->lastnchgbds = 0;
   presol->lastnaddholes = 0;
   presol->lastndelconss = 0;
   presol->lastnupgdconss = 0;
   presol->lastnchgcoefs = 0;
   presol->lastnchgsides = 0;
   presol->nfixedvars = 0;
   presol->naggrvars = 0;
   presol->nchgvartypes = 0;
   presol->nchgbds = 0;
   presol->naddholes = 0;
   presol->ndelconss = 0;
   presol->nupgdconss = 0;
   presol->nchgcoefs = 0;
   presol->nchgsides = 0;

   /* call initialization method of presolver */
   if( presol->presolinit != NULL )
   {
      CHECK_OKAY( presol->presolinit(set->scip, presol) );
   }
   presol->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes presolver */
RETCODE SCIPpresolExit(
   PRESOL*          presol,             /**< presolver */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(presol != NULL);
   assert(set != NULL);

   if( !presol->initialized )
   {
      errorMessage("presolver <%s> not initialized\n", presol->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of presolver */
   if( presol->presolexit != NULL )
   {
      CHECK_OKAY( presol->presolexit(set->scip, presol) );
   }
   presol->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs presolver that the presolving process is being started */
RETCODE SCIPpresolInitpre(
   PRESOL*          presol,             /**< presolver */
   SET*             set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(presol != NULL);
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   presol->lastnfixedvars = 0;
   presol->lastnaggrvars = 0;
   presol->lastnchgvartypes = 0;
   presol->lastnchgbds = 0;
   presol->lastnaddholes = 0;
   presol->lastndelconss = 0;
   presol->lastnupgdconss = 0;
   presol->lastnchgcoefs = 0;
   presol->lastnchgsides = 0;

   /* call presolving initialization method of presolver */
   if( presol->presolinitpre != NULL )
   {
      CHECK_OKAY( presol->presolinitpre(set->scip, presol, result) );

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_FEASIBLE )
      {
         errorMessage("presolving initialization method of presolver <%s> returned invalid result <%d>\n", 
            presol->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** informs presolver that the presolving process is finished */
RETCODE SCIPpresolExitpre(
   PRESOL*          presol,             /**< presolver */
   SET*             set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(presol != NULL);
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   /* call presolving deinitialization method of presolver */
   if( presol->presolexitpre != NULL )
   {
      CHECK_OKAY( presol->presolexitpre(set->scip, presol, result) );

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_UNBOUNDED
         && *result != SCIP_FEASIBLE )
      {
         errorMessage("presolving deinitialization method of presolver <%s> returned invalid result <%d>\n", 
            presol->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** executes presolver */
RETCODE SCIPpresolExec(
   PRESOL*          presol,             /**< presolver */
   SET*             set,                /**< global SCIP settings */
   int              nrounds,            /**< number of presolving rounds already done */
   int*             nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*             naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*             nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*             nchgbds,            /**< pointer to total number of variable bounds tightend of all presolvers */
   int*             naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*             ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*             nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*             nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*             nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   int nnewfixedvars;
   int nnewaggrvars;
   int nnewchgvartypes;
   int nnewchgbds;
   int nnewholes;
   int nnewdelconss;
   int nnewupgdconss;
   int nnewchgcoefs;
   int nnewchgsides;

   assert(presol != NULL);
   assert(presol->presolexec != NULL);
   assert(set != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nchgvartypes != NULL);
   assert(nchgbds != NULL);
   assert(naddholes != NULL);
   assert(ndelconss != NULL);
   assert(nupgdconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check number of presolving rounds */
   if( presol->maxrounds >= 0 && nrounds >= presol->maxrounds )
      return SCIP_OKAY;

   debugMessage("calling presolver <%s>\n", presol->name);

   /* calculate the number of changes since last call */
   nnewfixedvars = *nfixedvars - presol->lastnfixedvars;
   nnewaggrvars = *naggrvars - presol->lastnaggrvars;
   nnewchgvartypes = *nchgvartypes - presol->lastnchgvartypes;
   nnewchgbds = *nchgbds - presol->lastnchgbds;
   nnewholes = *naddholes - presol->lastnaddholes;
   nnewdelconss = *ndelconss - presol->lastndelconss;
   nnewupgdconss = *nupgdconss - presol->lastnupgdconss;
   nnewchgcoefs = *nchgcoefs - presol->lastnchgcoefs;
   nnewchgsides = *nchgsides - presol->lastnchgsides;
   assert(nnewfixedvars >= 0);
   assert(nnewaggrvars >= 0);
   assert(nnewchgvartypes >= 0);
   assert(nnewchgbds >= 0);
   assert(nnewholes >= 0);
   assert(nnewdelconss >= 0);
   assert(nnewupgdconss >= 0);
   assert(nnewchgcoefs >= 0);
   assert(nnewchgsides >= 0);

   /* remember the old number of changes */
   presol->lastnfixedvars = *nfixedvars;
   presol->lastnaggrvars = *naggrvars;
   presol->lastnchgvartypes = *nchgvartypes;
   presol->lastnchgbds = *nchgbds;
   presol->lastnaddholes = *naddholes;
   presol->lastndelconss = *ndelconss;
   presol->lastnupgdconss = *nupgdconss;
   presol->lastnchgcoefs = *nchgcoefs;
   presol->lastnchgsides = *nchgsides;

   /* start timing */
   SCIPclockStart(presol->clock, set);

   /* call external method */
   CHECK_OKAY( presol->presolexec(set->scip, presol, nrounds,
                  nnewfixedvars, nnewaggrvars, nnewchgvartypes, nnewchgbds, nnewholes,
                  nnewdelconss, nnewupgdconss, nnewchgcoefs, nnewchgsides,
                  nfixedvars, naggrvars, nchgvartypes, nchgbds, naddholes,
                  ndelconss, nupgdconss, nchgcoefs, nchgsides, result) );

   /* stop timing */
   SCIPclockStop(presol->clock, set);

   /* count the new changes */
   presol->nfixedvars += *nfixedvars - presol->lastnfixedvars;
   presol->naggrvars += *naggrvars - presol->lastnaggrvars;
   presol->nchgvartypes += *nchgvartypes - presol->lastnchgvartypes;
   presol->nchgbds += *nchgbds - presol->lastnchgbds;
   presol->naddholes += *naddholes - presol->lastnaddholes;
   presol->ndelconss += *ndelconss - presol->lastndelconss;
   presol->nupgdconss += *nupgdconss - presol->lastnupgdconss;
   presol->nchgcoefs += *nchgcoefs - presol->lastnchgcoefs;
   presol->nchgsides += *nchgsides - presol->lastnchgsides;

   /* check result code of callback method */
   if( *result != SCIP_CUTOFF
      && *result != SCIP_UNBOUNDED
      && *result != SCIP_SUCCESS
      && *result != SCIP_DIDNOTFIND
      && *result != SCIP_DIDNOTRUN )
   {
      errorMessage("presolver <%s> returned invalid result <%d>\n", presol->name, *result);
      return SCIP_INVALIDRESULT;
   }

   return SCIP_OKAY;
}

/** gets user data of presolver */
PRESOLDATA* SCIPpresolGetData(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->presoldata;
}

/** sets user data of presolver; user has to free old data in advance! */
void SCIPpresolSetData(
   PRESOL*          presol,             /**< presolver */
   PRESOLDATA*      presoldata          /**< new presolver user data */
   )
{
   assert(presol != NULL);

   presol->presoldata = presoldata;
}

/** gets name of presolver */
const char* SCIPpresolGetName(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->name;
}

/** gets description of presolver */
const char* SCIPpresolGetDesc(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->desc;
}

/** gets priority of presolver */
int SCIPpresolGetPriority(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->priority;
}

/** sets priority of presolver */
void SCIPpresolSetPriority(
   PRESOL*          presol,             /**< presolver */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the presolver */
   )
{
   assert(presol != NULL);
   assert(set != NULL);
   
   presol->priority = priority;
   set->presolssorted = FALSE;
}

/** is presolver initialized? */
Bool SCIPpresolIsInitialized(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->initialized;
}

/** gets time in seconds used in this presolver */
Real SCIPpresolGetTime(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return SCIPclockGetTime(presol->clock);
}

/** gets number of variables fixed in presolver */
int SCIPpresolGetNFixedVars(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nfixedvars;
}

/** gets number of variables aggregated in presolver */
int SCIPpresolGetNAggrVars(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->naggrvars;
}

/** gets number of variable types changed in presolver */
int SCIPpresolGetNVarTypes(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgvartypes;
}

/** gets number of bounds changed in presolver */
int SCIPpresolGetNChgBds(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgbds;
}

/** gets number of holes added to domains of variables in presolver */
int SCIPpresolGetNAddHoles(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->naddholes;
}

/** gets number of constraints deleted in presolver */
int SCIPpresolGetNDelConss(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->ndelconss;
}

/** gets number of constraints upgraded in presolver */
int SCIPpresolGetNUpgdConss(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nupgdconss;
}

/** gets number of coefficients changed in presolver */
int SCIPpresolGetNChgCoefs(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgcoefs;
}

/** gets number of constraint sides changed in presolver */
int SCIPpresolGetNChgSides(
   PRESOL*          presol              /**< presolver */
   )
{
   assert(presol != NULL);

   return presol->nchgsides;
}

