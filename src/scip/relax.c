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
#pragma ident "@(#) $Id: relax.c,v 1.7 2005/02/07 14:08:26 bzfpfend Exp $"

/**@file   relax.c
 * @brief  methods and datastructures for relaxators
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
#include "scip.h"
#include "relax.h"

#include "struct_relax.h"



/** compares two relaxators w. r. to their priority */
DECL_SORTPTRCOMP(SCIPrelaxComp)
{  /*lint --e{715}*/
   return ((RELAX*)elem2)->priority - ((RELAX*)elem1)->priority;
}

/** method to call, when the priority of a relaxator was changed */
static
DECL_PARAMCHGD(paramChgdRelaxPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetRelaxPriority() to mark the relaxs unsorted */
   CHECK_OKAY( SCIPsetRelaxPriority(scip, (RELAX*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a relaxator */
RETCODE SCIPrelaxCreate(
   RELAX**          relax,              /**< pointer to relaxator data structure */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of relaxator */
   const char*      desc,               /**< description of relaxator */
   int              priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
   int              freq,               /**< frequency for calling relaxator */
   DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxator */
   DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxator */
   DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxator */
   DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxator */
   DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxator */
   DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxator */
   RELAXDATA*       relaxdata           /**< relaxator data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(relax != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(relaxexec != NULL);

   ALLOC_OKAY( allocMemory(relax) );
   ALLOC_OKAY( duplicateMemoryArray(&(*relax)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*relax)->desc, desc, strlen(desc)+1) );
   (*relax)->priority = priority;
   (*relax)->freq = freq;
   (*relax)->relaxfree = relaxfree;
   (*relax)->relaxinit = relaxinit;
   (*relax)->relaxexit = relaxexit;
   (*relax)->relaxinitsol = relaxinitsol;
   (*relax)->relaxexitsol = relaxexitsol;
   (*relax)->relaxexec = relaxexec;
   (*relax)->relaxdata = relaxdata;
   CHECK_OKAY( SCIPclockCreate(&(*relax)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*relax)->ncalls = 0;
   (*relax)->lastsolvednode = -1;
   (*relax)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "relaxing/%s/priority", name);
   sprintf(paramdesc, "priority of relaxator <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*relax)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdRelaxPriority, (PARAMDATA*)(*relax)) ); /*lint !e740*/
   sprintf(paramname, "relaxing/%s/freq", name);
   sprintf(paramdesc, "frequency for calling relaxator <%s> (-1: never, 0: only in root node)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*relax)->freq, freq, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of relaxator */
RETCODE SCIPrelaxFree(
   RELAX**          relax,              /**< pointer to relaxator data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(*relax != NULL);
   assert(!(*relax)->initialized);
   assert(set != NULL);

   /* call destructor of relaxator */
   if( (*relax)->relaxfree != NULL )
   {
      CHECK_OKAY( (*relax)->relaxfree(set->scip, *relax) );
   }

   SCIPclockFree(&(*relax)->clock);
   freeMemoryArray(&(*relax)->name);
   freeMemoryArray(&(*relax)->desc);
   freeMemory(relax);

   return SCIP_OKAY;
}

/** initializes relaxator */
RETCODE SCIPrelaxInit(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   if( relax->initialized )
   {
      errorMessage("relaxator <%s> already initialized\n", relax->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(relax->clock);
   relax->ncalls = 0;

   if( relax->relaxinit != NULL )
   {
      CHECK_OKAY( relax->relaxinit(set->scip, relax) );
   }
   relax->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of relaxator */
RETCODE SCIPrelaxExit(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   if( !relax->initialized )
   {
      errorMessage("relaxator <%s> not initialized\n", relax->name);
      return SCIP_INVALIDCALL;
   }

   if( relax->relaxexit != NULL )
   {
      CHECK_OKAY( relax->relaxexit(set->scip, relax) );
   }
   relax->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs relaxator that the branch and bound process is being started */
RETCODE SCIPrelaxInitsol(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   /* call solving process initialization method of relaxator */
   if( relax->relaxinitsol != NULL )
   {
      CHECK_OKAY( relax->relaxinitsol(set->scip, relax) );
   }

   return SCIP_OKAY;
}

/** informs relaxator that the branch and bound process data is being freed */
RETCODE SCIPrelaxExitsol(
   RELAX*           relax,              /**< relaxator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of relaxator */
   if( relax->relaxexitsol != NULL )
   {
      CHECK_OKAY( relax->relaxexitsol(set->scip, relax) );
   }

   return SCIP_OKAY;
}

/** calls execution method of relaxator */
RETCODE SCIPrelaxExec(
   RELAX*           relax,              /**< relaxator */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth,              /**< depth of current node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(relax != NULL);
   assert(relax->relaxexec != NULL);
   assert(relax->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* check, if the relaxation is already solved */
   if( relax->lastsolvednode == stat->ntotalnodes )
      return SCIP_OKAY;

   relax->lastsolvednode = stat->ntotalnodes;

   if( (depth == 0 && relax->freq == 0) || (relax->freq > 0 && depth % relax->freq == 0) )
   {
      debugMessage("executing relaxator <%s>\n", relax->name);

      /* start timing */
      SCIPclockStart(relax->clock, set);

      /* call external relaxration method */
      CHECK_OKAY( relax->relaxexec(set->scip, relax, result) );

      /* stop timing */
      SCIPclockStop(relax->clock, set);

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_SUCCESS
         && *result != SCIP_SUSPENDED
         && *result != SCIP_DIDNOTRUN )
      {
         errorMessage("execution method of relaxator <%s> returned invalid result <%d>\n", 
            relax->name, *result);
         return SCIP_INVALIDRESULT;
      }
      if( *result != SCIP_DIDNOTRUN )
      {
         relax->ncalls++;
         if( *result == SCIP_SUSPENDED )
            SCIPrelaxMarkUnsolved(relax);
      }
   }

   return SCIP_OKAY;
}

/** gets user data of relaxator */
RELAXDATA* SCIPrelaxGetData(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->relaxdata;
}

/** sets user data of relaxator; user has to free old data in advance! */
void SCIPrelaxSetData(
   RELAX*           relax,              /**< relaxator */
   RELAXDATA*       relaxdata           /**< new relaxator user data */
   )
{
   assert(relax != NULL);

   relax->relaxdata = relaxdata;
}

/** gets name of relaxator */
const char* SCIPrelaxGetName(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->name;
}

/** gets description of relaxator */
const char* SCIPrelaxGetDesc(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->desc;
}

/** gets priority of relaxator */
int SCIPrelaxGetPriority(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->priority;
}

/** sets priority of relaxator */
void SCIPrelaxSetPriority(
   RELAX*           relax,              /**< relaxator */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the relaxator */
   )
{
   assert(relax != NULL);
   assert(set != NULL);
   
   relax->priority = priority;
   set->relaxssorted = FALSE;
}

/** gets frequency of relaxator */
int SCIPrelaxGetFreq(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->freq;
}

/** gets time in seconds used in this relaxator */
Real SCIPrelaxGetTime(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return SCIPclockGetTime(relax->clock);
}

/** gets the total number of times, the relaxator was called */
Longint SCIPrelaxGetNCalls(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->ncalls;
}

/** is relaxator initialized? */
Bool SCIPrelaxIsInitialized(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->initialized;
}

/** returns whether the relaxation was completely solved at the current node */
Bool SCIPrelaxIsSolved(
   RELAX*           relax,              /**< relaxator */
   STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(relax != NULL);
   assert(stat != NULL);

   return (relax->lastsolvednode == stat->ntotalnodes);
}

/** marks the current relaxation unsolved, s.t. the relaxator is called again in the next solving round */
void SCIPrelaxMarkUnsolved(
   RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   relax->lastsolvednode = -1;
}

