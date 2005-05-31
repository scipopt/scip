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
#pragma ident "@(#) $Id: sepa.c,v 1.47 2005/05/31 17:20:20 bzfpfend Exp $"

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
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of separator */
   const char*      desc,               /**< description of separator */
   int              priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int              freq,               /**< frequency for calling separator */
   Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
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
   (*sepa)->sepainitsol = sepainitsol;
   (*sepa)->sepaexitsol = sepaexitsol;
   (*sepa)->sepaexec = sepaexec;
   (*sepa)->sepadata = sepadata;
   CHECK_OKAY( SCIPclockCreate(&(*sepa)->clock, SCIP_CLOCKTYPE_DEFAULT) );
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
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*sepa)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdSepaPriority, (PARAMDATA*)(*sepa)) ); /*lint !e740*/

   sprintf(paramname, "separating/%s/freq", name);
   sprintf(paramdesc, "frequency for calling separator <%s> (-1: never, 0: only in root node)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*sepa)->freq, freq, -1, INT_MAX, NULL, NULL) );

   sprintf(paramname, "separating/%s/delay", name);
   CHECK_OKAY( SCIPsetAddBoolParam(set, blkmem, paramname,
         "should separator be delayed, if other separators found cuts?",
         &(*sepa)->delay, delay, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of separator */
RETCODE SCIPsepaFree(
   SEPA**           sepa,               /**< pointer to separator data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(*sepa != NULL);
   assert(!(*sepa)->initialized);
   assert(set != NULL);

   /* call destructor of separator */
   if( (*sepa)->sepafree != NULL )
   {
      CHECK_OKAY( (*sepa)->sepafree(set->scip, *sepa) );
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
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

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
      CHECK_OKAY( sepa->sepainit(set->scip, sepa) );
   }
   sepa->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of separator */
RETCODE SCIPsepaExit(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   if( !sepa->initialized )
   {
      errorMessage("separator <%s> not initialized\n", sepa->name);
      return SCIP_INVALIDCALL;
   }

   if( sepa->sepaexit != NULL )
   {
      CHECK_OKAY( sepa->sepaexit(set->scip, sepa) );
   }
   sepa->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs separator that the branch and bound process is being started */
RETCODE SCIPsepaInitsol(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   /* call solving process initialization method of separator */
   if( sepa->sepainitsol != NULL )
   {
      CHECK_OKAY( sepa->sepainitsol(set->scip, sepa) );
   }

   return SCIP_OKAY;
}

/** informs separator that the branch and bound process data is being freed */
RETCODE SCIPsepaExitsol(
   SEPA*            sepa,               /**< separator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(sepa != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of separator */
   if( sepa->sepaexitsol != NULL )
   {
      CHECK_OKAY( sepa->sepaexitsol(set->scip, sepa) );
   }

   return SCIP_OKAY;
}

/** calls execution method of separator */
RETCODE SCIPsepaExec(
   SEPA*            sepa,               /**< separator */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   SEPASTORE*       sepastore,          /**< separation storage */
   int              depth,              /**< depth of current node */
   Bool             execdelayed,        /**< execute separator even if it is marked to be delayed */
   RESULT*          result              /**< pointer to store the result of the callback method */
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

         debugMessage("executing separator <%s>\n", sepa->name);

         oldncutsstored = SCIPsepastoreGetNCutsStored(sepastore);

         /* reset the statistics for current node */
         if( sepa->lastsepanode != stat->ntotalnodes )
         {
            sepa->ncallsatnode = 0;
            sepa->ncutsfoundatnode = 0;
         }

         /* start timing */
         SCIPclockStart(sepa->clock, set);

         /* call external separation method */
         CHECK_OKAY( sepa->sepaexec(set->scip, sepa, result) );

         /* stop timing */
         SCIPclockStop(sepa->clock, set);

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
            errorMessage("execution method of separator <%s> returned invalid result <%d>\n", 
               sepa->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         debugMessage("separator <%s> was delayed\n", sepa->name);
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
   SEPA*            sepa,               /**< separator */
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

/** should separator be delayed, if other separators found cuts? */
Bool SCIPsepaIsDelayed(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->delay;
}

/** was separator delayed at the last call? */
Bool SCIPsepaWasDelayed(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->wasdelayed;
}

/** is separator initialized? */
Bool SCIPsepaIsInitialized(
   SEPA*            sepa                /**< separator */
   )
{
   assert(sepa != NULL);

   return sepa->initialized;
}
