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
#pragma ident "@(#) $Id: prop.c,v 1.10 2005/02/14 13:35:47 bzfpfend Exp $"

/**@file   prop.c
 * @brief  methods and datastructures for propagators
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
#include "scip/var.h"
#include "scip/scip.h"
#include "scip/prop.h"

#include "scip/struct_prop.h"



/** compares two propagators w. r. to their priority */
DECL_SORTPTRCOMP(SCIPpropComp)
{  /*lint --e{715}*/
   return ((PROP*)elem2)->priority - ((PROP*)elem1)->priority;
}

/** method to call, when the priority of a propagator was changed */
static
DECL_PARAMCHGD(paramChgdPropPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPropPriority() to mark the props unsorted */
   CHECK_OKAY( SCIPsetPropPriority(scip, (PROP*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a propagator */
RETCODE SCIPpropCreate(
   PROP**           prop,               /**< pointer to propagator data structure */
   SET*             set,                /**< global SCIP settings */
   BLKMEM*          blkmem,             /**< block memory for parameter settings */
   const char*      name,               /**< name of propagator */
   const char*      desc,               /**< description of propagator */
   int              priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int              freq,               /**< frequency for calling propagator */
   Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   PROPDATA*        propdata            /**< propagator data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(prop != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(freq >= -1);
   assert(propexec != NULL);

   ALLOC_OKAY( allocMemory(prop) );
   ALLOC_OKAY( duplicateMemoryArray(&(*prop)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*prop)->desc, desc, strlen(desc)+1) );
   (*prop)->priority = priority;
   (*prop)->freq = freq;
   (*prop)->propfree = propfree;
   (*prop)->propinit = propinit;
   (*prop)->propexit = propexit;
   (*prop)->propinitsol = propinitsol;
   (*prop)->propexitsol = propexitsol;
   (*prop)->propexec = propexec;
   (*prop)->propresprop = propresprop;
   (*prop)->propdata = propdata;
   CHECK_OKAY( SCIPclockCreate(&(*prop)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*prop)->ncalls = 0;
   (*prop)->ncutoffs = 0;
   (*prop)->ndomredsfound = 0;
   (*prop)->wasdelayed = FALSE;
   (*prop)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "propagating/%s/priority", name);
   sprintf(paramdesc, "priority of propagator <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*prop)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdPropPriority, (PARAMDATA*)(*prop)) ); /*lint !e740*/

   sprintf(paramname, "propagating/%s/freq", name);
   sprintf(paramdesc, "frequency for calling propagator <%s> (-1: never, 0: only in root node)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*prop)->freq, freq, -1, INT_MAX, NULL, NULL) );

   sprintf(paramname, "propagating/%s/delay", name);
   CHECK_OKAY( SCIPsetAddBoolParam(set, blkmem, paramname,
         "should propagator be delayed, if other propagators found reductions?",
         &(*prop)->delay, delay, NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of propagator */
RETCODE SCIPpropFree(
   PROP**           prop,               /**< pointer to propagator data structure */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(*prop != NULL);
   assert(!(*prop)->initialized);
   assert(set != NULL);

   /* call destructor of propagator */
   if( (*prop)->propfree != NULL )
   {
      CHECK_OKAY( (*prop)->propfree(set->scip, *prop) );
   }

   SCIPclockFree(&(*prop)->clock);
   freeMemoryArray(&(*prop)->name);
   freeMemoryArray(&(*prop)->desc);
   freeMemory(prop);

   return SCIP_OKAY;
}

/** initializes propagator */
RETCODE SCIPpropInit(
   PROP*            prop,               /**< propagator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   if( prop->initialized )
   {
      errorMessage("propagator <%s> already initialized\n", prop->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(prop->clock);

   prop->ncalls = 0;
   prop->ncutoffs = 0;
   prop->ndomredsfound = 0;

   if( prop->propinit != NULL )
   {
      CHECK_OKAY( prop->propinit(set->scip, prop) );
   }
   prop->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of propagator */
RETCODE SCIPpropExit(
   PROP*            prop,               /**< propagator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   if( !prop->initialized )
   {
      errorMessage("propagator <%s> not initialized\n", prop->name);
      return SCIP_INVALIDCALL;
   }

   if( prop->propexit != NULL )
   {
      CHECK_OKAY( prop->propexit(set->scip, prop) );
   }
   prop->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs propagator that the branch and bound process is being started */
RETCODE SCIPpropInitsol(
   PROP*            prop,               /**< propagator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call solving process initialization method of propagator */
   if( prop->propinitsol != NULL )
   {
      CHECK_OKAY( prop->propinitsol(set->scip, prop) );
   }

   return SCIP_OKAY;
}

/** informs propagator that the branch and bound process data is being freed */
RETCODE SCIPpropExitsol(
   PROP*            prop,               /**< propagator */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(prop != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of propagator */
   if( prop->propexitsol != NULL )
   {
      CHECK_OKAY( prop->propexitsol(set->scip, prop) );
   }

   return SCIP_OKAY;
}

/** calls execution method of propagator */
RETCODE SCIPpropExec(
   PROP*            prop,               /**< propagator */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   int              depth,              /**< depth of current node */
   Bool             execdelayed,        /**< execute propagator even if it is marked to be delayed */
   RESULT*          result              /**< pointer to store the result of the callback method */
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
         Longint oldndomchgs;
         
         debugMessage("executing propagator <%s>\n", prop->name);
         
         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         
         /* start timing */
         SCIPclockStart(prop->clock, set);
         
         /* call external propagation method */
         CHECK_OKAY( prop->propexec(set->scip, prop, result) );
         
         /* stop timing */
         SCIPclockStop(prop->clock, set);
         
         /* update statistics */
         if( *result != SCIP_DIDNOTRUN && *result != SCIP_DELAYED )
            prop->ncalls++;
         if( *result == SCIP_CUTOFF )
            prop->ncutoffs++;
         prop->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;

         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_DIDNOTFIND
            && *result != SCIP_DIDNOTRUN
            && *result != SCIP_DELAYED )
         {
            errorMessage("execution method of propagator <%s> returned invalid result <%d>\n", 
               prop->name, *result);
            return SCIP_INVALIDRESULT;
         }
      }
      else
      {
         debugMessage("propagator <%s> was delayed\n", prop->name);
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
RETCODE SCIPpropResolvePropagation(
   PROP*            prop,               /**< propagator */
   SET*             set,                /**< global SCIP settings */
   VAR*             infervar,           /**< variable whose bound was deduced by the constraint */
   int              inferinfo,          /**< user inference information attached to the bound change */
   BOUNDTYPE        inferboundtype,     /**< bound that was deduced (lower or upper bound) */
   BDCHGIDX*        bdchgidx,           /**< bound change index, representing the point of time where change took place */
   RESULT*          result              /**< pointer to store the result of the callback method */
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
      CHECK_OKAY( prop->propresprop(set->scip, prop, infervar, inferinfo, inferboundtype, bdchgidx,
            result) );
      
      /* check result code */
      if( *result != SCIP_SUCCESS && *result != SCIP_DIDNOTFIND )
      {
         errorMessage("propagation conflict resolving method of propagator <%s> returned invalid result <%d>\n", 
            prop->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   else
   {
      errorMessage("propagation conflict resolving method of propagator <%s> is not implemented\n", prop->name);
      return SCIP_PLUGINNOTFOUND;
   }

   return SCIP_OKAY;
}

/** gets user data of propagator */
PROPDATA* SCIPpropGetData(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->propdata;
}

/** sets user data of propagator; user has to free old data in advance! */
void SCIPpropSetData(
   PROP*            prop,               /**< propagator */
   PROPDATA*        propdata            /**< new propagator user data */
   )
{
   assert(prop != NULL);

   prop->propdata = propdata;
}

/** gets name of propagator */
const char* SCIPpropGetName(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->name;
}

/** gets description of propagator */
const char* SCIPpropGetDesc(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->desc;
}

/** gets priority of propagator */
int SCIPpropGetPriority(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->priority;
}

/** sets priority of propagator */
void SCIPpropSetPriority(
   PROP*            prop,               /**< propagator */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the propagator */
   )
{
   assert(prop != NULL);
   assert(set != NULL);
   
   prop->priority = priority;
   set->propssorted = FALSE;
}

/** gets frequency of propagator */
int SCIPpropGetFreq(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->freq;
}

/** gets time in seconds used in this propagator */
Real SCIPpropGetTime(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return SCIPclockGetTime(prop->clock);
}

/** gets the total number of times, the propagator was called */
Longint SCIPpropGetNCalls(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ncalls;
}

/** gets total number of times, this propagator detected a cutoff */
Longint SCIPpropGetNCutoffs(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ncutoffs;
}

/** gets total number of domain reductions found by this propagator */
Longint SCIPpropGetNDomredsFound(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->ndomredsfound;
}

/** should propagator be delayed, if other propagators found reductions? */
Bool SCIPpropIsDelayed(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->delay;
}

/** was propagator delayed at the last call? */
Bool SCIPpropWasDelayed(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->wasdelayed;
}

/** is propagator initialized? */
Bool SCIPpropIsInitialized(
   PROP*            prop                /**< propagator */
   )
{
   assert(prop != NULL);

   return prop->initialized;
}
