/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop.c
 * @brief  methods and datastructures for propagators
 * @author Tobias Achterberg
 * @author Timo Berthold
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
#include "scip/var.h"
#include "scip/scip.h"
#include "scip/prop.h"
#include "scip/pub_misc.h"

#include "scip/struct_prop.h"



/** compares two propagators w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPpropComp)
{  /*lint --e{715}*/
   return ((SCIP_PROP*)elem2)->priority - ((SCIP_PROP*)elem1)->priority;
}

/** method to call, when the priority of a propagator was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdPropPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPropPriority() to mark the props unsorted */
   SCIP_CALL( SCIPsetPropPriority(scip, (SCIP_PROP*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given propagator to a new scip */
SCIP_RETCODE SCIPpropCopyInclude(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(prop != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( prop->propcopy != NULL )
   {
      SCIPdebugMessage("including propagator %s in subscip %p\n", SCIPpropGetName(prop), (void*)set->scip);
      SCIP_CALL( prop->propcopy(set->scip, prop) );
   }
   return SCIP_OKAY;
}

/** creates a propagator */
SCIP_RETCODE SCIPpropCreate(
   SCIP_PROP**           prop,               /**< pointer to propagator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into subscips */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(prop != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(propexec != NULL);

   SCIP_ALLOC( BMSallocMemory(prop) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*prop)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*prop)->desc, desc, strlen(desc)+1) );
   (*prop)->priority = priority;
   (*prop)->freq = freq;
   (*prop)->propcopy = propcopy;
   (*prop)->propfree = propfree;
   (*prop)->propinit = propinit;
   (*prop)->propexit = propexit;
   (*prop)->propinitsol = propinitsol;
   (*prop)->propexitsol = propexitsol;
   (*prop)->propexec = propexec;
   (*prop)->propresprop = propresprop;
   (*prop)->propdata = propdata;
   SCIP_CALL( SCIPclockCreate(&(*prop)->proptime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*prop)->resproptime, SCIP_CLOCKTYPE_DEFAULT) );
   (*prop)->ncalls = 0;
   (*prop)->nrespropcalls = 0;
   (*prop)->ncutoffs = 0;
   (*prop)->ndomredsfound = 0;
   (*prop)->wasdelayed = FALSE;
   (*prop)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of propagator <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*prop)->priority, TRUE, priority, INT_MIN/4, INT_MAX/4, 
         paramChgdPropPriority, (SCIP_PARAMDATA*)(*prop)) ); /*lint !e740*/

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/freq", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "frequency for calling propagator <%s> (-1: never, 0: only in root node)", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*prop)->freq, FALSE, freq, -1, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "propagating/%s/delay", name);
   SCIP_CALL( SCIPsetAddBoolParam(set, blkmem, paramname,
         "should propagator be delayed, if other propagators found reductions?",
         &(*prop)->delay, TRUE, delay, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of propagator */
SCIP_RETCODE SCIPpropFree(
   SCIP_PROP**           prop,               /**< pointer to propagator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(*prop != NULL);
   assert(!(*prop)->initialized);
   assert(set != NULL);

   /* call destructor of propagator */
   if( (*prop)->propfree != NULL )
   {
      SCIP_CALL( (*prop)->propfree(set->scip, *prop) );
   }

   SCIPclockFree(&(*prop)->proptime);
   SCIPclockFree(&(*prop)->resproptime);
   BMSfreeMemoryArray(&(*prop)->name);
   BMSfreeMemoryArray(&(*prop)->desc);
   BMSfreeMemory(prop);

   return SCIP_OKAY;
}

/** initializes propagator */
SCIP_RETCODE SCIPpropInit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   if( prop->initialized )
   {
      SCIPerrorMessage("propagator <%s> already initialized\n", prop->name);
      return SCIP_INVALIDCALL;
   }

   if( set->misc_resetstat )
   {
      SCIPclockReset(prop->proptime);
      SCIPclockReset(prop->resproptime);

      prop->ncalls = 0;
      prop->nrespropcalls = 0;
      prop->ncutoffs = 0;
      prop->ndomredsfound = 0;
   }

   if( prop->propinit != NULL )
   {
      SCIP_CALL( prop->propinit(set->scip, prop) );
   }
   prop->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of propagator */
SCIP_RETCODE SCIPpropExit(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   if( !prop->initialized )
   {
      SCIPerrorMessage("propagator <%s> not initialized\n", prop->name);
      return SCIP_INVALIDCALL;
   }

   if( prop->propexit != NULL )
   {
      SCIP_CALL( prop->propexit(set->scip, prop) );
   }
   prop->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs propagator that the branch and bound process is being started */
SCIP_RETCODE SCIPpropInitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call solving process initialization method of propagator */
   if( prop->propinitsol != NULL )
   {
      SCIP_CALL( prop->propinitsol(set->scip, prop) );
   }

   return SCIP_OKAY;
}

/** informs propagator that the branch and bound process data is being freed */
SCIP_RETCODE SCIPpropExitsol(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of propagator */
   if( prop->propexitsol != NULL )
   {
      SCIP_CALL( prop->propexitsol(set->scip, prop) );
   }

   return SCIP_OKAY;
}

/** calls execution method of propagator */
SCIP_RETCODE SCIPpropExec(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Bool             execdelayed,        /**< execute propagator even if it is marked to be delayed */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(prop != NULL);
   assert(prop->propexec != NULL);
   assert(prop->freq >= -1);
   assert(set != NULL);
   assert(set->scip != NULL);
   assert(stat != NULL);
   assert(depth >= 0);
   assert(result != NULL);

   if( (depth == 0 && prop->freq == 0) || (prop->freq > 0 && depth % prop->freq == 0) || prop->wasdelayed )
   {
      if( !prop->delay || execdelayed )
      {
         SCIP_Longint oldndomchgs;
         SCIP_Longint oldnprobdomchgs;
         
         SCIPdebugMessage("executing propagator <%s>\n", prop->name);
         
         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnprobdomchgs = stat->nprobboundchgs + stat->nprobholechgs;
         
         /* start timing */
         SCIPclockStart(prop->proptime, set);
         
         /* call external propagation method */
         SCIP_CALL( prop->propexec(set->scip, prop, result) );
         
         /* stop timing */
         SCIPclockStop(prop->proptime, set);
         
         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            prop->ncalls++;
         if( *result == SCIP_CUTOFF )
            prop->ncutoffs++;

         /* update domain reductions; therefore remove the domain
          * reduction counts which were generated in probing mode */
         prop->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
         prop->ndomredsfound -= (stat->nprobboundchgs + stat->nprobholechgs - oldnprobdomchgs);

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            SCIPerrorMessage("execution method of propagator <%s> returned invalid result <%d>\n", 
               prop->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         SCIPdebugMessage("propagator <%s> was delayed\n", prop->name);
         *result = SCIP_DELAYED;
      }

      /* remember whether propagator was delayed */
      prop->wasdelayed = (*result == SCIP_DELAYED);
   }
   else
      *result = SCIP_DIDNOTRUN;
   
   return SCIP_OKAY;
}

/** resolves the given conflicting bound, that was deduced by the given propagator, by putting all "reason" bounds
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictLb() and SCIPaddConflictUb()
 */
SCIP_RETCODE SCIPpropResolvePropagation(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int                   inferinfo,          /**< user inference information attached to the bound change */
   SCIP_BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(prop != NULL);
   assert((inferboundtype == SCIP_BOUNDTYPE_LOWER
         && SCIPvarGetLbAtIndex(infervar, bdchgidx, TRUE) > SCIPvarGetLbGlobal(infervar))
      || (inferboundtype == SCIP_BOUNDTYPE_UPPER
         && SCIPvarGetUbAtIndex(infervar, bdchgidx, TRUE) < SCIPvarGetUbGlobal(infervar)));
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   if( prop->propresprop != NULL )
   {

      /* start timing */
      SCIPclockStart(prop->resproptime, set);

      SCIP_CALL( prop->propresprop(set->scip, prop, infervar, inferinfo, inferboundtype, bdchgidx,
            result) );
      
      /* stop timing */
      SCIPclockStop(prop->resproptime, set);

      /* update statistic */
      prop->nrespropcalls++;
      
      /* check result code */
      if( *result != SCIP_SUCCESS && *result != SCIP_DIDNOTFIND )
      {
         SCIPerrorMessage("propagation conflict resolving method of propagator <%s> returned invalid result <%d>\n", 
            prop->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   else
   {
      SCIPerrorMessage("propagation conflict resolving method of propagator <%s> is not implemented\n", prop->name);
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** gets user data of propagator */
SCIP_PROPDATA* SCIPpropGetData(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->propdata;
}

/** sets user data of propagator; user has to free old data in advance! */
void SCIPpropSetData(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_PROPDATA*        propdata            /**< new propagator user data */
   )
{
   assert(prop != NULL);

   prop->propdata = propdata;
}

/** gets name of propagator */
const char* SCIPpropGetName(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->name;
}

/** gets description of propagator */
const char* SCIPpropGetDesc(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->desc;
}

/** gets priority of propagator */
int SCIPpropGetPriority(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->priority;
}

/** sets priority of propagator */
void SCIPpropSetPriority(
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the propagator */
   )
{
   assert(prop != NULL);
   assert(set != NULL);
   
   prop->priority = priority;
   set->propssorted = FALSE;
}

/** gets frequency of propagator */
int SCIPpropGetFreq(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->freq;
}

/** gets time in seconds used in this propagator for propagation */
SCIP_Real SCIPpropGetTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->proptime);
}

/** gets time in seconds used in this propagator for resolve propagation */
SCIP_Real SCIPpropGetRespropTime(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->resproptime);
}

/** gets the total number of times, the propagator was called */
SCIP_Longint SCIPpropGetNCalls(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ncalls;
}

/** gets the total number of times, the propagator was called for resolving a propagation */
SCIP_Longint SCIPpropGetNRespropCalls(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->nrespropcalls;
}

/** gets total number of times, this propagator detected a cutoff */
SCIP_Longint SCIPpropGetNCutoffs(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ncutoffs;
}

/** gets total number of domain reductions found by this propagator */
SCIP_Longint SCIPpropGetNDomredsFound(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ndomredsfound;
}

/** should propagator be delayed, if other propagators found reductions? */
SCIP_Bool SCIPpropIsDelayed(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->delay;
}

/** was propagator delayed at the last call? */
SCIP_Bool SCIPpropWasDelayed(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->wasdelayed;
}

/** is propagator initialized? */
SCIP_Bool SCIPpropIsInitialized(
   SCIP_PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->initialized;
}
