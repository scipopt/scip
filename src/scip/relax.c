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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: relax.c,v 1.11 2005/08/24 17:26:55 bzfpfend Exp $"

/**@file   relax.c
 * @brief  methods and datastructures for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/paramset.h"
#include "scip/scip.h"
#include "scip/relax.h"

#include "scip/struct_relax.h"



/** compares two relaxators w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPrelaxComp)
{  /*lint --e{715}*/
   return ((SCIP_RELAX*)elem2)->priority - ((SCIP_RELAX*)elem1)->priority;
}

/** method to call, when the priority of a relaxator was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdRelaxPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetRelaxPriority() to mark the relaxs unsorted */
   SCIP_CALL( SCIPsetRelaxPriority(scip, (SCIP_RELAX*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a relaxator */
SCIP_RETCODE SCIPrelaxCreate(
   SCIP_RELAX**          relax,              /**< pointer to relaxator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of relaxator */
   const char*           desc,               /**< description of relaxator */
   int                   priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxator */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxator */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxator */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxator */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxator */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxator */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxator */
   SCIP_RELAXDATA*       relaxdata           /**< relaxator data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(relax != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(relaxexec != NULL);

   SCIP_ALLOC( BMSallocMemory(relax) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*relax)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*relax)->desc, desc, strlen(desc)+1) );
   (*relax)->priority = priority;
   (*relax)->freq = freq;
   (*relax)->relaxfree = relaxfree;
   (*relax)->relaxinit = relaxinit;
   (*relax)->relaxexit = relaxexit;
   (*relax)->relaxinitsol = relaxinitsol;
   (*relax)->relaxexitsol = relaxexitsol;
   (*relax)->relaxexec = relaxexec;
   (*relax)->relaxdata = relaxdata;
   SCIP_CALL( SCIPclockCreate(&(*relax)->relaxclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*relax)->ncalls = 0;
   (*relax)->lastsolvednode = -1;
   (*relax)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "relaxing/%s/priority", name);
   sprintf(paramdesc, "priority of relaxator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*relax)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdRelaxPriority, (SCIP_PARAMDATA*)(*relax)) ); /*lint !e740*/
   sprintf(paramname, "relaxing/%s/freq", name);
   sprintf(paramdesc, "frequency for calling relaxator <%s> (-1: never, 0: only in root node)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*relax)->freq, freq, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** calls destructor and frees memory of relaxator */
SCIP_RETCODE SCIPrelaxFree(
   SCIP_RELAX**          relax,              /**< pointer to relaxator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(*relax != NULL);
   assert(!(*relax)->initialized);
   assert(set != NULL);

   /* call destructor of relaxator */
   if( (*relax)->relaxfree != NULL )
   {
      SCIP_CALL( (*relax)->relaxfree(set->scip, *relax) );
   }

   SCIPclockFree(&(*relax)->relaxclock);
   BMSfreeMemoryArray(&(*relax)->name);
   BMSfreeMemoryArray(&(*relax)->desc);
   BMSfreeMemory(relax);

   return SCIP_OKAY;
}

/** initializes relaxator */
SCIP_RETCODE SCIPrelaxInit(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   if( relax->initialized )
   {
      SCIPerrorMessage("relaxator <%s> already initialized\n", relax->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(relax->relaxclock);
   relax->ncalls = 0;

   if( relax->relaxinit != NULL )
   {
      SCIP_CALL( relax->relaxinit(set->scip, relax) );
   }
   relax->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of relaxator */
SCIP_RETCODE SCIPrelaxExit(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   if( !relax->initialized )
   {
      SCIPerrorMessage("relaxator <%s> not initialized\n", relax->name);
      return SCIP_INVALIDCALL;
   }

   if( relax->relaxexit != NULL )
   {
      SCIP_CALL( relax->relaxexit(set->scip, relax) );
   }
   relax->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs relaxator that the branch and bound process is being started */
SCIP_RETCODE SCIPrelaxInitsol(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   /* call solving process initialization method of relaxator */
   if( relax->relaxinitsol != NULL )
   {
      SCIP_CALL( relax->relaxinitsol(set->scip, relax) );
   }

   return SCIP_OKAY;
}

/** informs relaxator that the branch and bound process data is being freed */
SCIP_RETCODE SCIPrelaxExitsol(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(relax != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of relaxator */
   if( relax->relaxexitsol != NULL )
   {
      SCIP_CALL( relax->relaxexitsol(set->scip, relax) );
   }

   return SCIP_OKAY;
}

/** calls execution method of relaxator */
SCIP_RETCODE SCIPrelaxExec(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
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
      SCIPdebugMessage("executing relaxator <%s>\n", relax->name);

      /* start timing */
      SCIPclockStart(relax->relaxclock, set);

      /* call external relaxration method */
      SCIP_CALL( relax->relaxexec(set->scip, relax, result) );

      /* stop timing */
      SCIPclockStop(relax->relaxclock, set);

      /* evaluate result */
      if( *result != SCIP_CUTOFF
         && *result != SCIP_CONSADDED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_SUCCESS
         && *result != SCIP_SUSPENDED
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("execution method of relaxator <%s> returned invalid result <%d>\n", 
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
SCIP_RELAXDATA* SCIPrelaxGetData(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->relaxdata;
}

/** sets user data of relaxator; user has to free old data in advance! */
void SCIPrelaxSetData(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_RELAXDATA*       relaxdata           /**< new relaxator user data */
   )
{
   assert(relax != NULL);

   relax->relaxdata = relaxdata;
}

/** gets name of relaxator */
const char* SCIPrelaxGetName(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->name;
}

/** gets description of relaxator */
const char* SCIPrelaxGetDesc(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->desc;
}

/** gets priority of relaxator */
int SCIPrelaxGetPriority(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->priority;
}

/** sets priority of relaxator */
void SCIPrelaxSetPriority(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the relaxator */
   )
{
   assert(relax != NULL);
   assert(set != NULL);
   
   relax->priority = priority;
   set->relaxssorted = FALSE;
}

/** gets frequency of relaxator */
int SCIPrelaxGetFreq(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->freq;
}

/** gets time in seconds used in this relaxator */
SCIP_Real SCIPrelaxGetTime(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return SCIPclockGetTime(relax->relaxclock);
}

/** gets the total number of times, the relaxator was called */
SCIP_Longint SCIPrelaxGetNCalls(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->ncalls;
}

/** is relaxator initialized? */
SCIP_Bool SCIPrelaxIsInitialized(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   return relax->initialized;
}

/** returns whether the relaxation was completely solved at the current node */
SCIP_Bool SCIPrelaxIsSolved(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   )
{
   assert(relax != NULL);
   assert(stat != NULL);

   return (relax->lastsolvednode == stat->ntotalnodes);
}

/** marks the current relaxation unsolved, s.t. the relaxator is called again in the next solving round */
void SCIPrelaxMarkUnsolved(
   SCIP_RELAX*           relax               /**< relaxator */
   )
{
   assert(relax != NULL);

   relax->lastsolvednode = -1;
}

