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
#pragma ident "@(#) $Id: sepa.c,v 1.48 2005/08/22 18:35:47 bzfpfend Exp $"

/**@file   sepa.c
 * @brief  methods and datastructures for separators
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
#include "scip/sepastore.h"
#include "scip/scip.h"
#include "scip/sepa.h"

#include "scip/struct_sepa.h"



/** compares two separators w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPsepaComp)
{  /*lint --e{715}*/
   return ((SCIP_SEPA*)elem2)->priority - ((SCIP_SEPA*)elem1)->priority;
}

/** method to call, when the priority of a separator was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdSepaPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSepaPriority() to mark the sepas unsorted */
   SCIP_CALL( SCIPsetSepaPriority(scip, (SCIP_SEPA*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a separator */
SCIP_RETCODE SCIPsepaCreate(
   SCIP_SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXEC    ((*sepaexec)),      /**< execution method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(sepa != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(sepaexec != NULL);

   SCIP_ALLOC( BMSallocMemory(sepa) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*sepa)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*sepa)->desc, desc, strlen(desc)+1) );
   (*sepa)->priority = priority;
   (*sepa)->freq = freq;
   (*sepa)->sepafree = sepafree;
   (*sepa)->sepainit = sepainit;
   (*sepa)->sepaexit = sepaexit;
   (*sepa)->sepainitsol = sepainitsol;
   (*sepa)->sepaexitsol = sepaexitsol;
   (*sepa)->sepaexec = sepaexec;
   (*sepa)->sepadata = sepadata;
   SCIP_CALL( SCIPclockCreate(&(*sepa)->sepaclock, SCIP_CLOCKTYPE_DEFAULT) );
   (*sepa)->lastsepanode = -1;
   (*sepa)->ncalls = 0;
   (*sepa)->ncutsfound = 0;
   (*sepa)->ncallsatnode = 0;
   (*sepa)->ncutsfoundatnode = 0;
   (*sepa)->wasdelayed = FALSE;
   (*sepa)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "separating/%s/priority", name);
   sprintf(paramdesc, "priority of separator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*sepa)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdSepaPriority, (SCIP_PARAMDATA*)(*sepa)) ); /*lint !e740*/

   sprintf(paramname, "separating/%s/freq", name);
   sprintf(paramdesc, "frequency for calling separator <%s> (-1: never, 0: only in root node)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*sepa)->freq, freq, -1, INT_MAX, NULL, NULL) );

   sprintf(paramname, "separating/%s/delay", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, blkmem, paramname,
         "should separator be delayed, if other separators found cuts?",
         &(*sepa)->delay, delay, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of separator */
SCIP_RETCODE SCIPsepaFree(
   SCIP_SEPA**           sepa,               /**< pointer to separator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(*sepa != NULL);
   assert(!(*sepa)->initialized);
   assert(set != NULL);

   /* call destructor of separator */
   if( (*sepa)->sepafree != NULL )
   {
      SCIP_CALL( (*sepa)->sepafree(set->scip, *sepa) );
   }

   SCIPclockFree(&(*sepa)->sepaclock);
   BMSfreeMemoryArray(&(*sepa)->name);
   BMSfreeMemoryArray(&(*sepa)->desc);
   BMSfreeMemory(sepa);

   return SCIP_OKAY;
}

/** initializes separator */
SCIP_RETCODE SCIPsepaInit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   if( sepa->initialized )
   {
      SCIPerrorMessage("separator <%s> already initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(sepa->sepaclock);

   sepa->lastsepanode = -1;
   sepa->ncalls = 0;
   sepa->ncutsfound = 0;
   sepa->ncallsatnode = 0;
   sepa->ncutsfoundatnode = 0;

   if( sepa->sepainit != NULL )
   {
      SCIP_CALL( sepa->sepainit(set->scip, sepa) );
   }
   sepa->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of separator */
SCIP_RETCODE SCIPsepaExit(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   if( !sepa->initialized )
   {
      SCIPerrorMessage("separator <%s> not initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   if( sepa->sepaexit != NULL )
   {
      SCIP_CALL( sepa->sepaexit(set->scip, sepa) );
   }
   sepa->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs separator that the branch and bound process is being started */
SCIP_RETCODE SCIPsepaInitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   /* call solving process initialization method of separator */
   if( sepa->sepainitsol != NULL )
   {
      SCIP_CALL( sepa->sepainitsol(set->scip, sepa) );
   }

   return SCIP_OKAY;
}

/** informs separator that the branch and bound process data is being freed */
SCIP_RETCODE SCIPsepaExitsol(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of separator */
   if( sepa->sepaexitsol != NULL )
   {
      SCIP_CALL( sepa->sepaexitsol(set->scip, sepa) );
   }

   return SCIP_OKAY;
}

/** calls execution method of separator */
SCIP_RETCODE SCIPsepaExec(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_SEPASTORE*       sepastore,          /**< separation storage */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(sepa != NULL);
   assert(sepa->sepaexec != NULL);
   assert(sepa->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( (depth == 0 && sepa->freq == 0) || (sepa->freq > 0 && depth % sepa->freq == 0) || sepa->wasdelayed )
   {
      if( !sepa->delay || execdelayed )
      {
         int oldncutsstored;
         int ncutsfound;

         SCIPdebugMessage("executing separator <%s>\n", sepa->name);

         oldncutsstored = SCIPsepastoreGetNCutsStored(sepastore);

         /* reset the statistics for current node */
         if( sepa->lastsepanode != stat->ntotalnodes )
         {
            sepa->ncallsatnode = 0;
            sepa->ncutsfoundatnode = 0;
         }

         /* start timing */
         SCIPclockStart(sepa->sepaclock, set);

         /* call external separation method */
         SCIP_CALL( sepa->sepaexec(set->scip, sepa, result) );

         /* stop timing */
         SCIPclockStop(sepa->sepaclock, set);

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
         {
            sepa->ncalls++;
            sepa->ncallsatnode++;
            sepa->lastsepanode = stat->ntotalnodes;
         }
         ncutsfound = SCIPsepastoreGetNCutsStored(sepastore) - oldncutsstored;
         sepa->ncutsfound += ncutsfound;
         sepa->ncutsfoundatnode += ncutsfound;

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("execution method of separator <%s> returned invalid result <%d>\n", 
               sepa->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPdebugMessage("separator <%s> was delayed\n", sepa->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether separator was delayed */
      sepa->wasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;
   
   return SCIP_OKAY;
}

/** gets user data of separator */
SCIP_SEPADATA* SCIPsepaGetData(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->sepadata;
}

/** sets user data of separator; user has to free old data in advance! */
void SCIPsepaSetData(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata            /**< new separator user data */
   )
{
   assert(sepa != NULL);

   sepa->sepadata = sepadata;
}

/** gets name of separator */
const char* SCIPsepaGetName(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->name;
}

/** gets description of separator */
const char* SCIPsepaGetDesc(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->desc;
}

/** gets priority of separator */
int SCIPsepaGetPriority(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->priority;
}

/** sets priority of separator */
void SCIPsepaSetPriority(
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the separator */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);
   
   sepa->priority = priority;
   set->sepassorted = FALSE;
}

/** gets frequency of separator */
int SCIPsepaGetFreq(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->freq;
}

/** gets time in seconds used in this separator */
SCIP_Real SCIPsepaGetTime(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return SCIPclockGetTime(sepa->sepaclock);
}

/** gets the total number of times, the separator was called */
SCIP_Longint SCIPsepaGetNCalls(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncalls;
}

/** gets the number of times, the separator was called at the current node */
int SCIPsepaGetNCallsAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncallsatnode;
}

/** gets the total number of cutting planes found by this separator */
SCIP_Longint SCIPsepaGetNCutsFound(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfound;
}

/** gets the number of cutting planes found by this separator at the current node */
SCIP_Longint SCIPsepaGetNCutsFoundAtNode(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->ncutsfoundatnode;
}

/** should separator be delayed, if other separators found cuts? */
SCIP_Bool SCIPsepaIsDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->delay;
}

/** was separator delayed at the last call? */
SCIP_Bool SCIPsepaWasDelayed(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->wasdelayed;
}

/** is separator initialized? */
SCIP_Bool SCIPsepaIsInitialized(
   SCIP_SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->initialized;
}
