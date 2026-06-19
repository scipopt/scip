/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_sync.c
 * @ingroup DEFPLUGINS_PROP
 * @brief  propagator for applying global bound changes that were communicated by other
 *         concurrent solvers
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/concsolver.h"
#include "scip/concurrent.h"
#include "scip/prop_sync.h"
#include "scip/pub_message.h"
#include "scip/pub_prop.h"
#include "scip/pub_var.h"
#include "scip/scip_concurrent.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_probing.h"
#include "scip/scip_prop.h"
#include "scip/scip_var.h"
#include "scip/scip_message.h"
#include "scip/syncstore.h"
#include "tpi/tpi.h"

/* fundamental propagator properties */
#define PROP_NAME              "sync"
#define PROP_DESC              "propagator for synchronization of bound changes"
#define PROP_PRIORITY                     (INT_MAX/4) /**< propagator priority */
#define PROP_FREQ                                  -1 /**< propagator frequency */
#define PROP_DELAY                              FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING            SCIP_PROPTIMING_ALWAYS /**< propagation timing mask */

#define PROP_PRESOL_PRIORITY          (INT_MAX/4) /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_ALWAYS /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_VAR**            bndvar;        /**< array of variables with a bound change */
   SCIP_Real*            bndval;        /**< array of new bound values */
   SCIP_BOUNDTYPE*       bndtype;       /**< array of bound types */
   int                   nbnds;         /**< number of boundchanges */
   int                   bndsize;       /**< current size of bound change array */
   SCIP_Longint          ntightened;    /**< number of tightened bounds */
   SCIP_Longint          ntightenedint; /**< number of tightened bounds of integer variables */
   SCIP_Bool             useboundpool;  /**< should tightened global bounds be drained from the shared board every node? */
   SCIP_VAR**            boardvars;     /**< this solver's variables in the communication variable order (not owned) */
   int                   nboardvars;    /**< number of variables in boardvars */
   int                   lastboundversion; /**< board change counter at the last drain, to skip the scan when unchanged */
   SCIP_CONCSOLVER*      concsolver;    /**< the concurrent solver this propagator belongs to, for the shared statistics */
};


/*
 * Local methods
 */


/** apply the stored bound changes */
static
SCIP_RETCODE applyBoundChanges(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        data,               /**< propagator data */
   SCIP_RESULT*          result,             /**< result of propagations */
   int*                  ntightened,         /**< pointer to store the number of tightened bounds */
   int*                  ntightenedint       /**< pointer to store the number of tightened integer bounds */
   )
{
   int i;

   assert(data != NULL);
   assert(result != NULL);
   assert(ntightened != NULL);
   assert(ntightenedint  != NULL);

   *ntightened = 0;
   *ntightenedint = 0;

   SCIPdisableConcurrentBoundStorage(scip);
   *result = SCIP_DIDNOTFIND;

   for( i = 0; i < data->nbnds; ++i )
   {
      SCIP_Bool infeas;
      SCIP_Bool tightened;

      SCIP_CALL( SCIPvarGetProbvarBound(&data->bndvar[i], &data->bndval[i], &data->bndtype[i]) );

      /* cannot change bounds of multi-aggregated variables so skip this bound-change */
      if( SCIPvarGetStatus(data->bndvar[i]) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      if( data->bndtype[i] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, data->bndvar[i], data->bndval[i], FALSE, &infeas, &tightened) );
      }
      else
      {
         assert(data->bndtype[i] == SCIP_BOUNDTYPE_UPPER);
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, data->bndvar[i], data->bndval[i], FALSE, &infeas, &tightened) );
      }

      if( tightened )
      {
         ++(*ntightened);
         if( SCIPvarIsNonimpliedIntegral(data->bndvar[i]) )
            ++(*ntightenedint);
      }

      if( infeas )
      {
#ifndef NDEBUG
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "sync propagator found cutoff in thread %i.\n", SCIPtpiGetThreadNum());
#endif
         *result = SCIP_CUTOFF;
         break;
      }
   }

   data->nbnds = 0;
   SCIPenableConcurrentBoundStorage(scip);

   return SCIP_OKAY;
}

/** applies one shared global bound from the bound board to the local problem, resolving the variable
 *  through aggregations and skipping multi-aggregated variables; updates the tightening counters and
 *  the propagation result. The tightening is idempotent, so re-applying an already known bound is a no-op.
 */
static
SCIP_RETCODE applyBoardBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        data,               /**< propagator data */
   SCIP_VAR*             var,                /**< this solver's variable in communication order */
   SCIP_Real             val,                /**< shared bound value in communication space */
   SCIP_BOUNDTYPE        boundtype,          /**< whether val is a lower or an upper bound */
   SCIP_RESULT*          result              /**< pointer to update the propagation result */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(var != NULL);

   /* resolve to the active problem variable, adjusting the value and bound type through aggregations */
   SCIP_CALL( SCIPvarGetProbvarBound(&var, &val, &boundtype) );

   /* cannot change bounds of multi-aggregated variables so skip this bound */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, val, FALSE, &infeas, &tightened) );
   }
   else
   {
      SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, val, FALSE, &infeas, &tightened) );
   }

   if( tightened )
   {
      ++data->ntightened;
      if( SCIPvarIsNonimpliedIntegral(var) )
         ++data->ntightenedint;
   }

   if( infeas )
      *result = SCIP_CUTOFF;
   else if( tightened && *result != SCIP_CUTOFF )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/*
 * Callback methods of propagator
 */

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSync)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeBlockMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitSync)
{  /*lint --e{715}*/
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   data->bndsize = 0;
   data->nbnds = 0;
   data->bndvar = NULL;
   data->bndval = NULL;
   data->bndtype = NULL;
   data->ntightened = 0;
   data->ntightenedint = 0;

   return SCIP_OKAY;
}

/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitSync)
{  /*lint --e{715}*/
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &data->bndvar, data->bndsize);
   SCIPfreeBlockMemoryArrayNull(scip, &data->bndval, data->bndsize);
   SCIPfreeBlockMemoryArrayNull(scip, &data->bndtype, data->bndsize);

   return SCIP_OKAY;
}

/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolSync)
{  /*lint --e{715}*/
   SCIP_PROPDATA*  data;
   int ntightened;
   int ntightenedint;

   assert(prop != NULL);
   assert(result != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   *result = SCIP_DIDNOTRUN;

   if( data->nbnds == 0 || SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* remember number of tightened bounds before applying new bound tightenings */
   SCIP_CALL( applyBoundChanges(scip, data, result, &ntightened, &ntightenedint) );

   /* add number of tightened bounds to the total number of presolving boundchanges */
   if( ntightened > 0 )
   {
      *nchgbds += ntightened;
      data->ntightened += ntightened;
      data->ntightenedint += ntightened;
      if( *result != SCIP_CUTOFF )
         *result = SCIP_SUCCESS;
   }

   SCIPpropSetFreq(prop, -1);

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSync)
{  /*lint --e{715}*/
   SCIP_PROPDATA*  data;
   int ntightened;
   int ntightenedint;

   assert(prop != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   *result = SCIP_DIDNOTRUN;

   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   SCIP_CALL( applyBoundChanges(scip, data, result, &ntightened, &ntightenedint) );

   if( ntightened > 0 )
   {
      data->ntightened += ntightened;
      data->ntightenedint += ntightenedint;
      if( *result != SCIP_CUTOFF )
         *result = SCIP_REDUCEDDOM;
   }

   /* drain the shared bound board: apply the tightest global bounds contributed by the other concurrent
    * solvers. The board is monotone and the tightening is idempotent, so no ordering or barrier is needed
    * and a stale read only delays a tightening to the next node. The change counter gates the scan, so the
    * common no-change case costs a single read; only when a bound was shared since the last drain is the
    * board scanned. The version is read before the scan, so a change landing during the scan is caught next time.
    */
   if( data->useboundpool && *result != SCIP_CUTOFF )
   {
      SCIP_SYNCSTORE* syncstore;
      int version;

      syncstore = SCIPgetSyncstore(scip);
      version = SCIPsyncstoreGetBoundVersion(syncstore);

      if( version != data->lastboundversion )
      {
         SCIP_Longint nbnds0;
         SCIP_Longint nintbnds0;
         int v;

         nbnds0 = data->ntightened;
         nintbnds0 = data->ntightenedint;

         for( v = 0; v < data->nboardvars && *result != SCIP_CUTOFF; ++v )
         {
            SCIP_VAR* var;
            SCIP_Real lb;
            SCIP_Real ub;

            var = data->boardvars[v];
            if( var == NULL )
               continue;

            SCIPsyncstoreGetBound(syncstore, v, &lb, &ub);

            if( !SCIPisInfinity(scip, -lb) )
               SCIP_CALL( applyBoardBound(scip, data, var, lb, SCIP_BOUNDTYPE_LOWER, result) );

            if( *result != SCIP_CUTOFF && !SCIPisInfinity(scip, ub) )
               SCIP_CALL( applyBoardBound(scip, data, var, ub, SCIP_BOUNDTYPE_UPPER, result) );
         }

         /* report the bounds applied from the board as received tighter bounds, so the concurrent
          * solver statistics reflect the immediate sharing the same way the synchronization-point
          * path does (the propagator's own ntightened counters are updated in applyBoardBound)
          */
         if( data->concsolver != NULL && data->ntightened > nbnds0 )
            SCIPconcsolverAddNTighterBnds(data->concsolver, data->ntightened - nbnds0, data->ntightenedint - nintbnds0);

         data->lastboundversion = version;
      }
   }

   /* with the bound board active the propagator stays awake at every node to drain new bounds; otherwise
    * it sleeps and is woken again when bounds are pushed via SCIPpropSyncAddBndchg
    */
   if( data->useboundpool )
      SCIPpropSetFreq(prop, 1);
   else
      SCIPpropSetFreq(prop, -1);

   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

/** creates the sync propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropSync(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata = NULL;
   SCIP_PROP* prop = NULL;

   SCIP_CALL( SCIPallocBlockMemory(scip, &propdata) );
   propdata->useboundpool = FALSE;
   propdata->boardvars = NULL;
   propdata->nboardvars = 0;
   propdata->lastboundversion = 0;
   propdata->concsolver = NULL;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecSync, propdata) );
   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSync) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitSync) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSync) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolSync, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );

   return SCIP_OKAY;
}

/** enables the bound-board drain in the sync propagator; it then runs at every node and applies the
 *  tightened global variable bounds published immediately by the other concurrent solvers
 */
void SCIPpropSyncEnableBoundPool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< the sync propagator */
   SCIP_CONCSOLVER*      concsolver,         /**< the concurrent solver this SCIP instance belongs to */
   SCIP_VAR**            vars,               /**< variables of this SCIP in the communication variable order */
   int                   nvars               /**< number of variables */
   )
{
   SCIP_PROPDATA* data;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(concsolver != NULL);
   assert(vars != NULL);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   data->useboundpool = TRUE;
   data->boardvars = vars;
   data->nboardvars = nvars;
   data->lastboundversion = 0;
   data->concsolver = concsolver;
   SCIPpropSetFreq(prop, 1);
}


/** adds a boundchange to the sync propagator */
SCIP_RETCODE SCIPpropSyncAddBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< sync propagator */
   SCIP_VAR*             var,                /**< variable for bound */
   SCIP_Real             val,                /**< value of bound */
   SCIP_BOUNDTYPE        bndtype             /**< type of bound */
   )
{
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   SCIP_STRINGEQ( SCIPpropGetName(prop), PROP_NAME, SCIP_INVALIDCALL );

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   if( data->nbnds + 1 > data->bndsize )
   {
      int newsize;
      newsize = SCIPcalcMemGrowSize(scip, data->nbnds+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->bndvar, data->bndsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->bndval, data->bndsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->bndtype, data->bndsize, newsize) );
      data->bndsize = newsize;
   }

   data->bndvar[data->nbnds] = var;
   data->bndval[data->nbnds] = val;
   data->bndtype[data->nbnds] = bndtype;

   if( data->nbnds == 0 )
   {
      SCIPpropSetFreq(prop, 1);
   }
   ++data->nbnds;

   return SCIP_OKAY;
}

/** returns the total number of tightened bounds found by the sync propagator */
SCIP_Longint SCIPpropSyncGetNTightenedBnds(
   SCIP_PROP*            prop                /**< sync propagator */
   )
{
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   return data->ntightened;
}

/** returns the total number of tightened bounds for integer variables found by the sync propagator */
SCIP_Longint SCIPpropSyncGetNTightenedIntBnds(
   SCIP_PROP*            prop                /**< sync propagator */
   )
{
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   return data->ntightenedint;
}
