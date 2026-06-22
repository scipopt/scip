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

/**@file   syncstore.c
 * @ingroup PARALLEL
 * @brief  the function definitions of the synchronization store
 * @author Leona Gottwald
 * @author Stephen J. Maher
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/pub_message.h"
#include "scip/concsolver.h"
#include "scip/struct_concsolver.h"
#include "scip/prob.h"
#include "scip/scip.h"
#include "blockmemshell/memory.h"
#include "tpi/tpi.h"
#include "scip/struct_syncstore.h"
#include "scip/concurrent.h"
#include "scip/syncstore.h"
#include "scip/boundstore.h"


/** computes the size of the array of synchronization data, such that
 *  it cannot ever happen that a synchronization data is reused while still
 *  not read by any thread */
static
int getNSyncdata(
   SCIP*                 scip                /**< SCIP main data structure */
   )
{
   int maxnsyncdelay;

   SCIP_CALL_ABORT( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &maxnsyncdelay) );

   return 2 * (maxnsyncdelay + 1);
}

/** creates and captures a new synchronization store */
SCIP_RETCODE SCIPsyncstoreCreate(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to return the created synchronization store */
   )
{
   assert(syncstore != NULL);

   SCIPdebugMessage("SCIPsyncstoreCreate()\n");

   SCIP_ALLOC( BMSallocMemory(syncstore) );

   (*syncstore)->mode = SCIP_PARA_DETERMINISTIC;                      /* initializing the mode */
   (*syncstore)->initialized = FALSE;
   (*syncstore)->syncdata = NULL;
   (*syncstore)->stopped = FALSE;
   (*syncstore)->nuses = 1;
   (*syncstore)->bestminobj = SCIP_REAL_MAX;
   (*syncstore)->bestdualbound = -SCIP_REAL_MAX;
   (*syncstore)->winnerid = -1;
   (*syncstore)->winnerstatus = SCIP_STATUS_UNKNOWN;
   (*syncstore)->printincumbents = FALSE;
   (*syncstore)->usesolpool = FALSE;
   (*syncstore)->poolsols = NULL;
   (*syncstore)->poolobjs = NULL;
   (*syncstore)->poolsource = NULL;
   (*syncstore)->poolnvals = NULL;
   (*syncstore)->npoolsols = 0;
   (*syncstore)->poolsolssize = 0;

   SCIP_CALL( SCIPtpiInitLock(&(*syncstore)->lock) );

   return SCIP_OKAY;
}

/** releases a synchronization store */
SCIP_RETCODE SCIPsyncstoreRelease(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to the synchronization store */
   )
{
   int references;

   assert(syncstore != NULL);
   if( *syncstore == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPtpiAcquireLock((*syncstore)->lock) );
   (*syncstore)->nuses -= 1;
   references = (*syncstore)->nuses;
   SCIP_CALL( SCIPtpiReleaseLock((*syncstore)->lock) );

   if( references == 0 )
   {
      if( (*syncstore)->initialized )
      {
         SCIP_CALL( SCIPsyncstoreExit(*syncstore) );
      }

      assert(!(*syncstore)->initialized);
      SCIPtpiDestroyLock(&(*syncstore)->lock);
      BMSfreeMemory(syncstore);
   }
   else
   {
      *syncstore = NULL;
   }

   return SCIP_OKAY;
}

/** captures a synchronization store */
SCIP_RETCODE SCIPsyncstoreCapture(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   SCIP_CALL( SCIPtpiAcquireLock(syncstore->lock) );

   ++(syncstore->nuses);

   SCIP_CALL( SCIPtpiReleaseLock(syncstore->lock) );

   return SCIP_OKAY;
}

/** initialize the syncstore for the given SCIP instance */
SCIP_RETCODE SCIPsyncstoreInit(
   SCIP*                 scip                /**< SCIP main data structure */
   )
{
   SCIP_SYNCSTORE* syncstore;
   int i;
   int j;
   int paramode;

   assert(scip != NULL);
   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);
   syncstore->mainscip = scip;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/gap", &syncstore->limit_gap) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/absgap", &syncstore->limit_absgap) );
   syncstore->lastsync = NULL;
   syncstore->nsolvers = SCIPgetNConcurrentSolvers(scip);

   syncstore->ninitvars = SCIPgetNVars(scip) + SCIPgetNFixedVars(scip);
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsols", &syncstore->maxnsols) );
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &syncstore->maxnsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/minsyncdelay", &syncstore->minsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqinit", &syncstore->syncfreqinit) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqmax", &syncstore->syncfreqmax) );
   syncstore->nsyncdata = getNSyncdata(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &(syncstore->syncdata), syncstore->nsyncdata) );

   for( i = 0; i < syncstore->nsyncdata; ++i )
   {
      syncstore->syncdata[i].syncnum = -1;
      SCIP_CALL( SCIPboundstoreCreate(syncstore->mainscip, &syncstore->syncdata[i].boundstore, syncstore->ninitvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solobj, syncstore->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solsource, syncstore->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols, syncstore->maxnsols) );

      for( j = 0; j < syncstore->maxnsols; ++j )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols[j], syncstore->ninitvars) );
      }

      SCIP_CALL( SCIPtpiInitLock(&(syncstore->syncdata[i].lock)) );
      SCIP_CALL( SCIPtpiInitCondition(&(syncstore->syncdata[i].allsynced)) );
   }

   syncstore->initialized = TRUE;
   syncstore->stopped = FALSE;
   syncstore->bestminobj = SCIP_REAL_MAX;
   syncstore->bestdualbound = -SCIP_REAL_MAX;
   syncstore->winnerid = -1;
   syncstore->winnerstatus = SCIP_STATUS_UNKNOWN;
   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/printincumbents", &syncstore->printincumbents) );

   SCIP_CALL( SCIPgetIntParam(scip, "parallel/mode", &paramode) );
   syncstore->mode = (SCIP_PARALLELMODE) paramode;

   /* the solution pool bypasses the synchronization-point barrier and is therefore only
    * available in opportunistic mode; deterministic mode keeps the regular protocol
    */
   SCIP_CALL( SCIPgetBoolParam(scip, "concurrent/solpool", &syncstore->usesolpool) );
   syncstore->usesolpool = syncstore->usesolpool && syncstore->mode == SCIP_PARA_OPPORTUNISTIC;
   syncstore->npoolsols = 0;

   SCIP_CALL( SCIPtpiInit(syncstore->nsolvers, INT_MAX, FALSE) );
   SCIP_CALL( SCIPautoselectDisps(scip) );

   if( syncstore->mode == SCIP_PARA_DETERMINISTIC )
   {
      /* in deterministic mode use the number of non-zeros and the number of variables to get a good
       * syncdelay and maximum syncfreq
       */
      syncstore->minsyncdelay *= 0.01 * (SCIPgetNNZs(scip) * SCIPgetNVars(scip)); /*lint !e790*/
      syncstore->syncfreqmax *= 0.01 * (SCIPgetNNZs(scip) * SCIPgetNVars(scip));  /*lint !e790*/
   }

   return SCIP_OKAY;
}

/** deinitializes the synchronization store */
SCIP_RETCODE SCIPsyncstoreExit(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   int i;
   int j;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   SCIP_CALL( SCIPtpiExit() );

   for( i = 0; i < syncstore->nsyncdata; ++i )
   {
      SCIPtpiDestroyLock(&(syncstore->syncdata[i].lock));
      SCIPtpiDestroyCondition(&(syncstore->syncdata[i].allsynced));
      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solobj, syncstore->maxnsols);
      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solsource, syncstore->maxnsols);
      SCIPboundstoreFree(syncstore->mainscip,  &syncstore->syncdata[i].boundstore);

      for( j = 0; j < syncstore->maxnsols; ++j )
      {
         SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols[j], syncstore->ninitvars);
      }

      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols, syncstore->maxnsols);
   }

   SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata, syncstore->nsyncdata);

   /* the solution pool is allocated with standard memory because it is filled from the
    * solver threads, where the main SCIP's block memory must not be used
    */
   for( i = 0; i < syncstore->npoolsols; ++i )
   {
      BMSfreeMemoryArray(&syncstore->poolsols[i]);
   }
   BMSfreeMemoryArrayNull(&syncstore->poolsols);
   BMSfreeMemoryArrayNull(&syncstore->poolobjs);
   BMSfreeMemoryArrayNull(&syncstore->poolsource);
   BMSfreeMemoryArrayNull(&syncstore->poolnvals);
   syncstore->npoolsols = 0;
   syncstore->poolsolssize = 0;
   syncstore->usesolpool = FALSE;

   syncstore->initialized = FALSE;
   syncstore->stopped = FALSE;

   return SCIP_OKAY;
}

/** checks whether the solve-is-stopped flag in the syncstore has been set by any thread */
SCIP_Bool SCIPsyncstoreSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   SCIP_Bool stopped;

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(syncstore->lock) );

   stopped = syncstore->stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(syncstore->lock) );

   return stopped;
}

/** sets the solve-is-stopped flag in the syncstore so that subsequent calls to
 *  SCIPsyncstoreSolveIsStopped will return the given value in any thread
 */
void SCIPsyncstoreSetSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Bool             stopped             /**< flag if the solve is stopped */
   )
{
   SCIP_CALL_ABORT( SCIPtpiAcquireLock(syncstore->lock) );

   syncstore->stopped = stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(syncstore->lock) );

   /* wake up any solvers blocked in SCIPsyncstoreEnsureAllSynced so they can
    * check the stopped flag and exit instead of waiting for all syncs to complete */
   if( stopped && syncstore->initialized )
   {
      int i;

      for( i = 0; i < syncstore->nsyncdata; ++i )
      {
         SCIP_CALL_ABORT( SCIPtpiBroadcastCondition(syncstore->syncdata[i].allsynced) );
      }
   }
}

/** appends an improving solution to the solution pool; the caller must hold the syncstore lock.
 *  Entries are append-only and immutable once published, and the size counter is incremented
 *  only after the entry is fully written, so the lock-free size hint read by the sync
 *  heuristics never exposes a partial entry.
 */
static
void syncstoreAddPoolSol(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Real             minobj,             /**< objective value in minimization-normalized original objective space */
   int                   ownerid,            /**< index of the concurrent solver that found the solution */
   SCIP_Real*            solvals,            /**< solution values in the communication variable order */
   int                   nsolvals            /**< number of solution values */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->usesolpool);
   assert(solvals != NULL);

   if( syncstore->npoolsols == syncstore->poolsolssize )
   {
      int newsize;

      newsize = MAX(64, 2 * syncstore->poolsolssize);
      SCIP_ALLOC_ABORT( BMSreallocMemoryArray(&syncstore->poolsols, newsize) );
      SCIP_ALLOC_ABORT( BMSreallocMemoryArray(&syncstore->poolobjs, newsize) );
      SCIP_ALLOC_ABORT( BMSreallocMemoryArray(&syncstore->poolsource, newsize) );
      SCIP_ALLOC_ABORT( BMSreallocMemoryArray(&syncstore->poolnvals, newsize) );
      syncstore->poolsolssize = newsize;
   }

   SCIP_ALLOC_ABORT( BMSduplicateMemoryArray(&syncstore->poolsols[syncstore->npoolsols], solvals, nsolvals) );
   syncstore->poolobjs[syncstore->npoolsols] = minobj;
   syncstore->poolsource[syncstore->npoolsols] = ownerid;
   syncstore->poolnvals[syncstore->npoolsols] = nsolvals;

   ++syncstore->npoolsols;
}

/** prints a line for a new globally best solution found by a concurrent solver; the caller must
 *  hold the syncstore lock so that the printed objective values are monotone and the lines of
 *  different solvers do not interleave
 */
static
void syncstorePrintIncumbent(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Real             minobj,             /**< objective value in minimization-normalized original objective space */
   const char*           solvername          /**< name of the concurrent solver that found the solution */
   )
{
   SCIP_Real origobj;

   assert(syncstore != NULL);

   /* the communicated objective lives in the space of the main SCIP's transformed
    * problem, of which the concurrent solvers are copies; that space may be scaled
    * and shifted relative to the user's original problem (e.g. objective scaling by
    * a gcd during the transformation), so map the value back before printing
    */
   origobj = SCIPretransformObj(syncstore->mainscip, minobj);
   SCIPverbMessage(syncstore->mainscip, SCIP_VERBLEVEL_NORMAL, NULL, "%.2f I.SOL %.9g  (%s)\n",
      SCIPgetSolvingTime(syncstore->mainscip), origobj, solvername);
}

/** updates the minimization-normalized objective value of the best solution found by any
 *  concurrent solver; used for printing new incumbents and as the improvement filter of the
 *  solution pool. If the solution pool is enabled and solution values are passed, the
 *  improving solution is published in the pool so that the other concurrent solvers can
 *  install it as an incumbent at their next drain. If published is not NULL, it returns whether
 *  the solution was actually added to the pool, so the caller can count it as shared.
 */
void SCIPsyncstoreUpdateBestMinObj(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Real             minobj,             /**< objective value in minimization-normalized original objective space */
   const char*           solvername,         /**< name of the concurrent solver that found the solution */
   int                   ownerid,            /**< index of the concurrent solver that found the solution */
   SCIP_Real*            solvals,            /**< solution values in the communication variable order, or NULL to
                                              *   only communicate the objective value */
   int                   nsolvals,           /**< number of solution values */
   SCIP_Bool*            published           /**< pointer to return whether the solution was added to the pool, or NULL */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   if( published != NULL )
      *published = FALSE;

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(syncstore->lock) );

   if( minobj < syncstore->bestminobj )
   {
      syncstore->bestminobj = minobj;

      /* publish the improving solution in the pool so the other solvers can install it (pool mode only) */
      if( syncstore->usesolpool && solvals != NULL )
      {
         syncstoreAddPoolSol(syncstore, minobj, ownerid, solvals, nsolvals);

         if( published != NULL )
            *published = TRUE;
      }

      /* optionally log the new global incumbent; independent of the solution pool and active in both modes */
      if( syncstore->printincumbents )
         syncstorePrintIncumbent(syncstore, minobj, solvername);
   }

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(syncstore->lock) );
}

/** updates the minimization-normalized best global dual bound; the concurrent solvers call this
 *  immediately while solving so that SCIPgetConcurrentDualbound stays up to date in solution-pool mode
 */
void SCIPsyncstoreUpdateBestDualbound(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Real             dualbound           /**< dual bound in minimization-normalized original objective space */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(syncstore->lock) );

   if( dualbound > syncstore->bestdualbound )
      syncstore->bestdualbound = dualbound;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(syncstore->lock) );
}

/** returns whether the solution pool for immediate solution sharing is enabled */
SCIP_Bool SCIPsyncstoreSolPoolEnabled(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);

   return syncstore->usesolpool;
}

/** returns the main SCIP that initialized the synchronization store; the bounds reported in
 *  solution-pool mode live in this SCIP's objective space, so it is used to convert them
 */
SCIP* SCIPsyncstoreGetMainScip(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);

   return syncstore->mainscip;
}

/** gets the current number of solutions in the solution pool; reads the counter without
 *  locking, so it may lag behind concurrent pushes, which is fine for its use as a
 *  cheap size hint that is polled at every node
 */
int SCIPsyncstoreGetNPoolSols(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->npoolsols;
}

/** gets the solution with the given index from the solution pool; entries are immutable
 *  once published, so the returned values can be read without holding any lock
 */
void SCIPsyncstoreGetPoolSol(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   int                   idx,                /**< index of the pooled solution, must be smaller than the size hint */
   SCIP_Real**           solvals,            /**< pointer to return the solution values in communication variable order */
   int*                  nsolvals,           /**< pointer to return the number of solution values */
   int*                  ownerid             /**< pointer to return the index of the contributing concurrent solver */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(solvals != NULL);
   assert(nsolvals != NULL);
   assert(ownerid != NULL);

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(syncstore->lock) );

   assert(idx >= 0 && idx < syncstore->npoolsols);
   *solvals = syncstore->poolsols[idx];
   *nsolvals = syncstore->poolnvals[idx];
   *ownerid = syncstore->poolsource[idx];

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(syncstore->lock) );
}

/** gets the minimization-normalized objective value of the best solution found by any
 *  concurrent solver, or SCIP_REAL_MAX if no solution was registered yet
 */
SCIP_Real SCIPsyncstoreGetBestMinObj(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   SCIP_Real minobj;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(syncstore->lock) );

   minobj = syncstore->bestminobj;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(syncstore->lock) );

   return minobj;
}

/** gets the upperbound from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastUpperbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   /* in solution-pool mode lastsync is not advanced during solving; report the immediately-tracked global
    * primal bound (bestminobj), which the getter then converts to the caller's objective space and sense
    */
   if( syncstore->usesolpool )
      return syncstore->bestminobj < SCIP_REAL_MAX ? syncstore->bestminobj : SCIPinfinity(syncstore->mainscip);

   return syncstore->lastsync == NULL ? SCIPinfinity(syncstore->mainscip) : syncstore->lastsync->bestupperbound;
}

/** gets the lowerbound from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastLowerbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   /* in solution-pool mode lastsync is not advanced during solving; report the immediately-tracked global
    * dual bound, which the getter then converts to the caller's objective space and sense
    */
   if( syncstore->usesolpool )
      return syncstore->bestdualbound > -SCIP_REAL_MAX ? syncstore->bestdualbound : -SCIPinfinity(syncstore->mainscip);

   return syncstore->lastsync == NULL ? -SCIPinfinity(syncstore->mainscip) : syncstore->lastsync->bestlowerbound;
}

/** gets the number of solutions from the last synchronization */
int SCIPsyncstoreGetLastNSols(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0 : syncstore->lastsync->nsols;
}

/** gets the number of boundchanges from the last synchronization */
int SCIPsyncstoreGetLastNBounds(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0 : SCIPboundstoreGetNChgs(syncstore->lastsync->boundstore);
}

/** gets total memory used by all solvers from the last synchronization */
SCIP_Longint SCIPsyncstoreGetLastMemTotal(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0 : syncstore->lastsync->memtotal;
}

/** gets the synchronization frequency from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastSyncfreq(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->lastsync == NULL ? 0.0 : syncstore->lastsync->syncfreq;
}

/** get synchronization data with given number. It is the responsibility of the caller
 *  to only ask for a synchronization number that still exists, which is checked
 *  with an assert in debug mode. */
SCIP_SYNCDATA* SCIPsyncstoreGetSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum             /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   )
{
   int j;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   j = (int) syncnum % syncstore->nsyncdata;

   /* check if requested syncnumber still exists if in debug mode */
   assert(syncstore->syncdata[j].syncnum == syncnum);

   return &syncstore->syncdata[j];
}

/** get the next synchronization data that should be read and
 *  adjust the delay. Returns NULL if no more data should be read due to minimum delay */
SCIP_SYNCDATA* SCIPsyncstoreGetNextSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq,           /**< the current synchronization frequency */
   SCIP_Longint          writenum,           /**< number of synchronizations the solver has written to */
   SCIP_Real*            delay               /**< pointer holding the current synchronization delay */
   )
{
   SCIP_Real newdelay;
   SCIP_Longint nextsyncnum;

   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(delay != NULL);

   if( syncdata == NULL )
   {
      nextsyncnum = 0;
   }
   else
   {
      if( syncdata->status != SCIP_STATUS_UNKNOWN )
         return NULL;

      nextsyncnum = syncdata->syncnum + 1;
   }

   if( nextsyncnum == writenum )
      return NULL;

   newdelay = *delay - syncfreq;

   /* if the delay would get too small we do not want to read the next syncdata.
    * But due to the limited length of the syncdata array we might need to
    * read this synchronization data anyways which is checked by the second part
    * of the if condition
    */
   if( newdelay < syncstore->minsyncdelay && nextsyncnum >= writenum - syncstore->maxnsyncdelay )
      return NULL;

   *delay = newdelay;
   assert(syncstore->syncdata[nextsyncnum % syncstore->nsyncdata].syncnum == nextsyncnum);

   return &syncstore->syncdata[nextsyncnum % syncstore->nsyncdata];
}

/** ensures that the given synchronization data has been written by
 *  all solvers upon return of this function and blocks the caller if necessary. */
SCIP_RETCODE SCIPsyncstoreEnsureAllSynced(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   /* check if waiting is required, make sure to hold the lock */
   SCIP_CALL( SCIPtpiAcquireLock(syncdata->lock) );

   while( syncdata->syncedcount < syncstore->nsolvers && !syncstore->stopped )
   {
      /* yes, so wait on the condition variable
       * (automatically releases the lock and reacquires it after the waiting)
       */
      SCIP_CALL( SCIPtpiWaitCondition(syncdata->allsynced, syncdata->lock) );
   }

   SCIP_CALL( SCIPtpiReleaseLock(syncdata->lock) );

   return SCIP_OKAY;
}

/** Tries to acquire the synchronization data with the given number for non-blocking reading.
 *  On success the synchronization data is returned in locked state and must be released with
 *  SCIPsyncstoreUnlockSyncdata after reading. If the data is not ready, NULL is returned and
 *  lost indicates why: if the slot was overwritten by a newer synchronization the data is
 *  lost and the reader should skip this number, otherwise the synchronization has not been
 *  written by all solvers yet and the reader should retry later. Never blocks the caller.
 */
SCIP_RETCODE SCIPsyncstoreTryLockCompleteSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum,            /**< the number of the synchronization to read */
   SCIP_SYNCDATA**       syncdata,           /**< pointer to return the locked synchronization data, or NULL */
   SCIP_Bool*            lost                /**< pointer to return whether the data was overwritten and is lost */
   )
{
   SCIP_SYNCDATA* data;

   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);
   assert(lost != NULL);

   data = &syncstore->syncdata[syncnum % syncstore->nsyncdata];

   SCIP_CALL( SCIPtpiAcquireLock(data->lock) );

   if( data->syncnum == syncnum && data->syncedcount == syncstore->nsolvers )
   {
      *syncdata = data;
      *lost = FALSE;
      return SCIP_OKAY;
   }

   /* a different synchronization number means the slot was recycled and the data of this
    * synchronization can never be read anymore; the same number with an incomplete synced
    * count means some solvers have not written it yet and the reader should retry later
    */
   *lost = data->syncnum != syncnum;
   *syncdata = NULL;

   SCIP_CALL( SCIPtpiReleaseLock(data->lock) );

   return SCIP_OKAY;
}

/** releases the lock of a synchronization data acquired with SCIPsyncstoreTryLockCompleteSyncdata */
SCIP_RETCODE SCIPsyncstoreUnlockSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data to unlock */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);

   SCIP_CALL( SCIPtpiReleaseLock(syncdata->lock) );

   return SCIP_OKAY;
}

/** Start synchronization for the given concurrent solver.
 *  Needs to be followed by a call to SCIPsyncstoreFinishSync if
 *  the syncdata that is returned is not NULL
 */
SCIP_RETCODE SCIPsyncstoreStartSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum,            /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   SCIP_SYNCDATA**       syncdata            /**< pointer to return the synchronization data */
   )
{
   int i;

   assert(syncdata != NULL);
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   if( SCIPsyncstoreSolveIsStopped(syncstore) )
   {
      *syncdata = NULL;
      return SCIP_OKAY;
   }

   i = syncnum % syncstore->nsyncdata; /*lint !e712*/
   *syncdata = &syncstore->syncdata[i];
   assert(*syncdata != NULL);

   SCIP_CALL( SCIPtpiAcquireLock((*syncdata)->lock) );

   if( (*syncdata)->syncnum != syncnum )
   {
      /* recycle freely: in solution-pool mode the winner is on the scalar channel, not this ring, so
       * never skipping the write lets a finishing solver record its terminal status and stop the portfolio */
      SCIPboundstoreClear((*syncdata)->boundstore);
      (*syncdata)->nsols = 0;
      (*syncdata)->memtotal = SCIPgetMemTotal(syncstore->mainscip);
      (*syncdata)->syncedcount = 0;
      (*syncdata)->bestupperbound = SCIPinfinity(syncstore->mainscip);
      (*syncdata)->bestlowerbound = -(*syncdata)->bestupperbound;
      (*syncdata)->status = SCIP_STATUS_UNKNOWN;
      (*syncdata)->winner = 0;
      (*syncdata)->syncnum = syncnum;
      (*syncdata)->syncfreq = 0.0;
   }

   return SCIP_OKAY;
}

/** whether terminal status (newstatus,newid) beats (curstatus,curid): closer to OPTIMAL wins, UNKNOWN
 *  worst, ties by smaller id; matches the ranking in SCIPsyncdataSetStatus */
static
SCIP_Bool syncstoreStatusIsBetter(
   SCIP_STATUS           newstatus,          /**< candidate terminal status */
   int                   newid,              /**< candidate solver id */
   SCIP_STATUS           curstatus,          /**< current best terminal status */
   int                   curid               /**< current best solver id */
   )
{
   if( curstatus == SCIP_STATUS_UNKNOWN )
      return newstatus != SCIP_STATUS_UNKNOWN;

   if( newstatus == SCIP_STATUS_UNKNOWN )
      return FALSE;

   /* both are real statuses; the smaller enum value is closer to SCIP_STATUS_OPTIMAL and thus better */
   if( newstatus != curstatus )
      return newstatus < curstatus;

   return newid < curid;
}

/** finishes synchronization for the synchronization data */
SCIP_RETCODE SCIPsyncstoreFinishSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA**       syncdata            /**< the synchronization data */
   )
{
   SCIP_Bool printline = FALSE;

   assert(syncdata != NULL);
   assert((*syncdata) != NULL);
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   ++(*syncdata)->syncedcount;

   /* record the globally best terminal status on the scalar winner channel; in solution-pool mode the
    * ring is recycled without the barrier, so lastsync may be overwritten or hold a worse, later status */
   if( (*syncdata)->status != SCIP_STATUS_UNKNOWN )
   {
      /* shared scalar but only the syncdata lock is held here; guard it (syncdata->syncstore ordering
       * matches the all-synced branch's SetSolveIsStopped) */
      SCIP_CALL( SCIPtpiAcquireLock(syncstore->lock) );

      if( syncstore->winnerstatus == SCIP_STATUS_UNKNOWN ||
         syncstoreStatusIsBetter((*syncdata)->status, (*syncdata)->winner, syncstore->winnerstatus, syncstore->winnerid) )
      {
         syncstore->winnerstatus = (*syncdata)->status;
         syncstore->winnerid = (*syncdata)->winner;
      }

      SCIP_CALL( SCIPtpiReleaseLock(syncstore->lock) );

      if( syncstore->lastsync != *syncdata )
      {
         syncstore->lastsync = *syncdata;
         printline = TRUE;
      }
   }

   if( (*syncdata)->syncedcount == syncstore->nsolvers )
   {
      if( (*syncdata)->status != SCIP_STATUS_UNKNOWN ||
         (SCIPgetConcurrentGap(syncstore->mainscip) <= syncstore->limit_gap) ||
         (SCIPgetNLimSolsFound(syncstore->mainscip) > 0 && REALABS(SCIPgetConcurrentPrimalbound(syncstore->mainscip) - SCIPgetConcurrentDualbound(syncstore->mainscip)) <= syncstore->limit_absgap) )
         SCIPsyncstoreSetSolveIsStopped(syncstore, TRUE);

      syncstore->lastsync = *syncdata;
      printline = TRUE;

      SCIP_CALL( SCIPtpiBroadcastCondition((*syncdata)->allsynced) );
   }

   SCIP_CALL( SCIPtpiReleaseLock((*syncdata)->lock) );

   if( printline )
   {
      SCIP_CALL( SCIPprintDisplayLine(syncstore->mainscip, NULL, SCIP_VERBLEVEL_HIGH, TRUE) );
   }

   *syncdata = NULL;

   return SCIP_OKAY;
}

/** gets status in synchronization data */
SCIP_STATUS SCIPsyncdataGetStatus(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->status;
}

/** gets the solver that had the best status, or -1 if solve is not stopped yet */
int SCIPsyncstoreGetWinner(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   /* solution-pool mode: the ring/lastsync is recycled, read the scalar winner */
   if( syncstore->usesolpool )
      return syncstore->winnerstatus == SCIP_STATUS_UNKNOWN ? -1 : syncstore->winnerid;

   if( syncstore->lastsync == NULL || syncstore->lastsync->status == SCIP_STATUS_UNKNOWN )
      return -1;

   return syncstore->lastsync->winner;
}

/** how many solvers have already finished synchronizing on this synchronization data */
int SCIPsyncdataGetNSynced(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->syncedcount;
}

/** how many solvers have are running concurrently */
int SCIPsyncstoreGetNSolvers(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   return syncstore->nsolvers;
}

/** read amount total memory used from synchronization data */
SCIP_Longint SCIPsyncdataGetMemTotal(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->memtotal;
}

/** read the synchronization frequency from a synchronization data */
SCIP_Real SCIPsyncdataGetSyncFreq(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->syncfreq;
}

/** read the upperbound stored in a synchronization data */
SCIP_Real SCIPsyncdataGetUpperbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->bestupperbound;
}

/** read the lowerbound stored in a synchronization data */
SCIP_Real SCIPsyncdataGetLowerbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->bestlowerbound;
}

/** read the solutions stored in a synchronization data */
void SCIPsyncdataGetSolutions(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real***          solvalues,          /**< array of buffers containing the solution values */
   int**                 solowner,           /**< array of ownerids of solutions */
   int*                  nsols               /**< pointer to return number of solutions */
   )
{
   assert(syncdata != NULL);
   assert(solvalues != NULL);
   assert(solowner != NULL);
   assert(nsols != NULL);

   *solvalues = syncdata->sols;
   *solowner = syncdata->solsource;
   *nsols = syncdata->nsols;
}

/** read bound changes stored in the synchronization data */
SCIP_BOUNDSTORE* SCIPsyncdataGetBoundChgs(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   )
{
   assert(syncdata != NULL);

   return syncdata->boundstore;
}

/** write the synchronization frequency to a synchronization data */
void SCIPsyncdataSetSyncFreq(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq            /**< the synchronization frequency */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);

   syncdata->syncfreq = MIN(syncfreq, syncstore->syncfreqmax);
}

/** set status in the synchronization data */
void SCIPsyncdataSetStatus(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_STATUS           status,             /**< the status */
   int                   solverid            /**< identifier of the solver that has this status */
   )
{
   assert(syncdata != NULL);

   /* check if status is better than current one (closer to SCIP_STATUS_OPTIMAL assumed to be followed by
    * SCIP_STATUS_INFEASIBLE and SCIP_STATUS_UNBOUNDED) and break ties by the solverid; remember the solver with the
    * best status so that the winner will be selected deterministically
    */
   if( syncdata->winner < 0 )
   {
      syncdata->status = status;
      syncdata->winner = solverid;
   }
   else if( syncdata->status < SCIP_STATUS_OPTIMAL )
   {
      if( status > syncdata->status || (status == syncdata->status && solverid < syncdata->winner) )
      {
         syncdata->status = status;
         syncdata->winner = solverid;
      }
   }
   else if( syncdata->status > SCIP_STATUS_OPTIMAL && status >= SCIP_STATUS_OPTIMAL )
   {
      if( status < syncdata->status || (status == syncdata->status && solverid < syncdata->winner) )
      {
         syncdata->status = status;
         syncdata->winner = solverid;
      }
   }
}

/** adds memory used to the synchronization data */
void SCIPsyncdataAddMemTotal(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Longint          memtotal            /**< the number of bytes used */
   )
{
   assert(syncdata != NULL);

   syncdata->memtotal += memtotal;
}

/** set upperbound to the synchronization data */
void SCIPsyncdataSetUpperbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_Real             upperbound          /**< the upperbound */
   )
{
   assert(syncdata != NULL);

   syncdata->bestupperbound = MIN(syncdata->bestupperbound, upperbound);
}

/** set lowerbound to the synchronization data */
void SCIPsyncdataSetLowerbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the lowerbound should be added to */
   SCIP_Real             lowerbound          /**< the lowerbound */
   )
{
   assert(syncdata != NULL);

   syncdata->bestlowerbound = MAX(syncdata->bestlowerbound, lowerbound);
}

/** gives a buffer to store the solution values, or NULL if solution should not be stored
 *  because there are already better solutions stored.
 */
void SCIPsyncdataGetSolutionBuffer(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Real             solobj,             /**< the objective value of the solution */
   int                   ownerid,            /**< an identifier for the owner of the solution, e.g. the thread number */
   SCIP_Real**           buffer              /**< pointer to return a buffer for the solution values, which must be set
                                              *   if the buffer is not NULL */
   )
{
   int pos;
   int i;

   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);
   assert(buffer != NULL);

   for( pos = 0; pos < syncdata->nsols; ++pos )
   {
      if( syncdata->solobj[pos] < solobj || (syncdata->solobj[pos] == solobj && ownerid < syncdata->solsource[pos]) ) /*lint !e777*/
         break;
   }

   if( syncdata->nsols < syncstore->maxnsols )
   {
      for( i = syncdata->nsols; i > pos; --i )
      {
         syncdata->solobj[i] = syncdata->solobj[i - 1];
         syncdata->solsource[i] = syncdata->solsource[i - 1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i - 1]);
      }

      ++syncdata->nsols;
   }
   else
   {
      --pos;

      for( i = 0; i < pos; ++i )
      {
         syncdata->solobj[i] = syncdata->solobj[i + 1];
         syncdata->solsource[i] = syncdata->solsource[i + 1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i + 1]);
      }
   }

   if( pos >= 0 )
   {
      syncdata->solobj[pos] = solobj;
      syncdata->solsource[pos] = ownerid;
      *buffer = syncdata->sols[pos];
   }
   else
   {
      *buffer = NULL;
   }
}

/** adds bound changes to the synchronization data */
SCIP_RETCODE SCIPsyncdataAddBoundChanges(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_BOUNDSTORE*      boundstore          /**< bound store containing the bounds to add */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);
   assert(syncdata != NULL);
   assert(boundstore != NULL);

   SCIP_CALL( SCIPboundstoreMerge(syncstore->mainscip, syncdata->boundstore, boundstore) );

   return SCIP_OKAY;
}

/** is synchronization store initialized */
SCIP_Bool SCIPsyncstoreIsInitialized(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);

   return syncstore->initialized;
}

/** returns the mode of the synchronization store */
SCIP_PARALLELMODE SCIPsyncstoreGetMode(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   )
{
   assert(syncstore != NULL);

   return syncstore->mode;
}
