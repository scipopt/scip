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
#pragma ident "@(#) $Id: presol.c,v 1.10 2003/11/21 10:35:38 bzfpfend Exp $"

/**@file   presol.c
 * @brief  methods and datastructures for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "presol.h"
#include "clock.h"


/** presolver */
struct Presol
{
   char*            name;               /**< name of presolver */
   char*            desc;               /**< description of presolver */
   int              priority;           /**< priority of the presolver */
   DECL_PRESOLFREE  ((*presolfree));    /**< destructor of presolver */
   DECL_PRESOLINIT  ((*presolinit));    /**< initialize presolver */
   DECL_PRESOLEXIT  ((*presolexit));    /**< deinitialize presolver */
   DECL_PRESOLEXEC  ((*presolexec));    /**< presolving execution method */
   PRESOLDATA*      presoldata;         /**< presolver data */
   CLOCK*           clock;              /**< presolving time */
   unsigned int     initialized:1;      /**< is presolver initialized? */
   int              lastnfixedvars;     /**< number of variables fixed before the last call to the presolver */
   int              lastnaggrvars;      /**< number of variables aggregated before the last call to the presolver */
   int              lastnchgvartypes;   /**< number of variable type changes before the last call to the presolver */
   int              lastnchgbds;        /**< number of variable bounds tightend before the last call to the presolver */
   int              lastnaddholes;      /**< number of domain holes added before the last call to the presolver */
   int              lastndelconss;      /**< number of deleted constraints before the last call to the presolver */
   int              lastnupgdconss;     /**< number of upgraded constraints before the last call to the presolver */
   int              lastnchgcoefs;      /**< number of changed coefficients before the last call to the presolver */
   int              lastnchgsides;      /**< number of changed left or right hand sides before the last call */
   int              nfixedvars;         /**< total number of variables fixed by this presolver */
   int              naggrvars;          /**< total number of variables aggregated by this presolver */
   int              nchgvartypes;       /**< total number of variable type changes by this presolver */
   int              nchgbds;            /**< total number of variable bounds tightend by this presolver */
   int              naddholes;          /**< total number of domain holes added by this presolver */
   int              ndelconss;          /**< total number of deleted constraints by this presolver */
   int              nupgdconss;         /**< total number of upgraded constraints by this presolver */
   int              nchgcoefs;          /**< total number of changed coefficients by this presolver */
   int              nchgsides;          /**< total number of changed left or right hand sides by this presolver */
};




/*
 * presolver methods
 */

/** creates a presolver */
RETCODE SCIPpresolCreate(
   PRESOL**         presol,             /**< pointer to store presolver */
   const char*      name,               /**< name of presolver */
   const char*      desc,               /**< description of presolver */
   int              priority,           /**< priority of the presolver */
   DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver */
   DECL_PRESOLINIT  ((*presolinit)),    /**< initialize presolver */
   DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialize presolver */
   DECL_PRESOLEXEC  ((*presolexec)),    /**< presolving execution method */
   PRESOLDATA*      presoldata          /**< presolver data */
   )
{
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
   (*presol)->presolexec = presolexec;
   (*presol)->presoldata = presoldata;
   CHECK_OKAY( SCIPclockCreate(&(*presol)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*presol)->initialized = FALSE;

   return SCIP_OKAY;
}

/** frees memory of presolver */   
RETCODE SCIPpresolFree(
   PRESOL**         presol,             /**< pointer to presolver data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(presol != NULL);
   assert(*presol != NULL);
   assert(!(*presol)->initialized);

   /* call destructor of presolver */
   if( (*presol)->presolfree != NULL )
   {
      CHECK_OKAY( (*presol)->presolfree(scip, *presol) );
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
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(presol != NULL);
   assert(scip != NULL);

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

   if( presol->presolinit != NULL )
   {
      CHECK_OKAY( presol->presolinit(scip, presol) );
   }
   presol->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes presolver */
RETCODE SCIPpresolExit(
   PRESOL*          presol,             /**< presolver */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(presol != NULL);
   assert(scip != NULL);

   if( !presol->initialized )
   {
      errorMessage("presolver <%s> not initialized\n", presol->name);
      return SCIP_INVALIDCALL;
   }

   if( presol->presolexit != NULL )
   {
      CHECK_OKAY( presol->presolexit(scip, presol) );
   }
   presol->initialized = FALSE;

   return SCIP_OKAY;
}

/** executes presolver */
RETCODE SCIPpresolExec(
   PRESOL*          presol,             /**< presolver */
   const SET*       set,                /**< global SCIP settings */
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

   debugMessage("calling presolver <%s>\n", presol->name);

   *result = SCIP_DIDNOTRUN;

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

