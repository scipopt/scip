/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: sepa.c,v 1.31 2004/02/04 17:27:41 bzfpfend Exp $"

/**@file   sepa.c
 * @brief  methods and datastructures for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "paramset.h"
#include "sepastore.h"
#include "scip.h"
#include "sepa.h"

#include "struct_sepa.h"



/** compares two separators w. r. to their priority */
DECL_SORTPTRCOMP(SCIPsepaComp)
{  /*lint --e{715}*/
   return ((SEPA*)elem2)->priority - ((SEPA*)elem1)->priority;
}

/** method to call, when the priority of a separator was changed */
static
DECL_PARAMCHGD(paramChgdSepaPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSepaPriority() to mark the sepas unsorted */
   CHECK_OKAY( SCIPsetSepaPriority(scip, (SEPA*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a separator */
RETCODE SCIPsepaCreate(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of the separator */
   int              freq,               /**< frequency for calling separator */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SEPADATA*        sepadata            /**< separator data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(sepa != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(sepaexec != NULL);

   ALLOC_OKAY( allocMemory(sepa) );
   ALLOC_OKAY( duplicateMemoryArray(&(*sepa)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*sepa)->desc, desc, strlen(desc)+1) );
   (*sepa)->priority = priority;
   (*sepa)->freq = freq;
   (*sepa)->sepafree = sepafree;
   (*sepa)->sepainit = sepainit;
   (*sepa)->sepaexit = sepaexit;
   (*sepa)->sepaexec = sepaexec;
   (*sepa)->sepadata = sepadata;
   CHECK_OKAY( SCIPclockCreate(&(*sepa)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*sepa)->lastsepanode = -1;
   (*sepa)->ncalls = 0;
   (*sepa)->ncutsfound = 0;
   (*sepa)->ncallsatnode = 0;
   (*sepa)->ncutsfoundatnode = 0;
   (*sepa)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "separating/%s/priority", name);
   sprintf(paramdesc, "priority of separator <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*sepa)->priority, priority, INT_MIN, INT_MAX, 
                  paramChgdSepaPriority, (PARAMDATA*)(*sepa)) ); /*lint !e740*/
   sprintf(paramname, "separating/%s/freq", name);
   sprintf(paramdesc, "frequency for calling separator <%s> (-1: never, 0: only in root node)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*sepa)->freq, freq, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of separator */
RETCODE SCIPsepaFree(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(sepa != NULL);
   assert(*sepa != NULL);
   assert(!(*sepa)->initialized);
   assert(scip != NULL);

   /* call destructor of separator */
   if( (*sepa)->sepafree != NULL )
   {
      CHECK_OKAY( (*sepa)->sepafree(scip, *sepa) );
   }

   SCIPclockFree(&(*sepa)->clock);
   freeMemoryArray(&(*sepa)->name);
   freeMemoryArray(&(*sepa)->desc);
   freeMemory(sepa);

   return SCIP_OKAY;
}

/** initializes separator */
RETCODE SCIPsepaInit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(sepa != NULL);
   assert(scip != NULL);

   if( sepa->initialized )
   {
      errorMessage("separator <%s> already initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(sepa->clock);

   sepa->lastsepanode = -1;
   sepa->ncalls = 0;
   sepa->ncutsfound = 0;
   sepa->ncallsatnode = 0;
   sepa->ncutsfoundatnode = 0;

   if( sepa->sepainit != NULL )
   {
      CHECK_OKAY( sepa->sepainit(scip, sepa) );
   }
   sepa->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of separator */
RETCODE SCIPsepaExit(
   SEPA*            sepa,               /**< separator */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(sepa != NULL);
   assert(scip != NULL);

   if( !sepa->initialized )
   {
      errorMessage("separator <%s> not initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   if( sepa->sepaexit != NULL )
   {
      CHECK_OKAY( sepa->sepaexit(scip, sepa) );
   }
   sepa->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls execution method of separator */
RETCODE SCIPsepaExec(
   SEPA*            sepa,               /**< separator */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              actdepth,           /**< depth of active node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(sepa != NULL);
   assert(sepa->sepaexec != NULL);
   assert(sepa->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(actdepth >= 0);
   assert(result != NULL);

   if( (actdepth == 0 && sepa->freq == 0) || (sepa->freq > 0 && actdepth % sepa->freq == 0) )
   {
      int oldncutsfound;
      int ncutsfound;

      debugMessage("executing separator <%s>\n", sepa->name);

      oldncutsfound = SCIPsepastoreGetNCutsFound(sepastore);

      /* reset the statistics for current node */
      if( sepa->lastsepanode != stat->nnodes )
      {
         sepa->lastsepanode = stat->nnodes;
         sepa->ncallsatnode = 0;
         sepa->ncutsfoundatnode = 0;
      }

      /* start timing */
      SCIPclockStart(sepa->clock, set);

      /* call external separation method */
      CHECK_OKAY( sepa->sepaexec(set->scip, sepa, result) );

      /* stop timing */
      SCIPclockStop(sepa->clock, set);

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_SEPARATED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_CONSADDED
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         errorMessage("execution method of separator <%s> returned invalid result <%d>\n", 
            sepa->name, *result);
         return SCIP_INVALIDRESULT;
      }
      if( *result != SCIP_DIDNOTRUN )
      {
         sepa->ncalls++;
         sepa->ncallsatnode++;
      }

      ncutsfound = SCIPsepastoreGetNCutsFound(sepastore) - oldncutsfound;
      sepa->ncutsfound += ncutsfound;
      sepa->ncutsfoundatnode += ncutsfound;
   }
   else
      *result = SCIP_DIDNOTRUN;
   
   return SCIP_OKAY;
}

/** gets user data of separator */
SEPADATA* SCIPsepaGetData(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->sepadata;
}

/** sets user data of separator; user has to free old data in advance! */
void SCIPsepaSetData(
   SEPA*            sepa,               /**< separator */
   SEPADATA*        sepadata            /**< new separator user data */
   )
{
   assert(sepa != NULL);

   sepa->sepadata = sepadata;
}

/** gets name of separator */
const char* SCIPsepaGetName(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->name;
}

/** gets description of separator */
const char* SCIPsepaGetDesc(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->desc;
}

/** gets priority of separator */
int SCIPsepaGetPriority(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->priority;
}

/** sets priority of separator */
void SCIPsepaSetPriority(
   SEPA*            sepa,               /**< primal sepaistic */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the separator */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);
   
   sepa->priority = priority;
   set->sepassorted = FALSE;
}

/** gets frequency of separator */
int SCIPsepaGetFreq(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->freq;
}

/** gets time in seconds used in this separator */
Real SCIPsepaGetTime(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return SCIPclockGetTime(sepa->clock);
}

/** gets the total number of times, the separator was called */
Longint SCIPsepaGetNCalls(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncalls;
}

/** gets the number of times, the separator was called at the current node */
int SCIPsepaGetNCallsAtNode(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncallsatnode;
}

/** gets the total number of cutting planes found by this separator */
Longint SCIPsepaGetNCutsFound(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfound;
}

/** gets the number of cutting planes found by this separator at the current node */
Longint SCIPsepaGetNCutsFoundAtNode(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfoundatnode;
}

/** is separator initialized? */
Bool SCIPsepaIsInitialized(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->initialized;
}
