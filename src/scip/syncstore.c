/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   parallel.c
 * @ingroup PARAINTERFACE
 * @brief  the interface functions for tinycthread
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
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

static
int getNSyncdata(
   SCIP*       scip
   )
{
   int maxnsyncdelay;
   SCIP_CALL_ABORT( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &maxnsyncdelay) );

   return 2 * (maxnsyncdelay + 1);
}

/** creates and captures a new synchronization store */
SCIP_RETCODE SCIPsyncstoreCreate(
   SCIP_SYNCSTORE**         syncstore            /**< pointer to return the created synchronization store */
   )
{
   assert(syncstore != NULL);

   SCIPdebugMessage("SCIPsyncstoreCreate()\n");

   SCIP_ALLOC( BMSallocMemory(syncstore) );

   (*syncstore)->mode = SCIP_PARA_DETERMINISTIC;                      /* initialising the mode */
   (*syncstore)->initialized = FALSE;
   (*syncstore)->syncdata = NULL;
   (*syncstore)->stopped = FALSE;
   (*syncstore)->references = 1;
   SCIP_CALL( SCIPtpiInitLock(&(*syncstore)->lock) );

   return SCIP_OKAY;
}

/** releases a synchronization store */
SCIP_RETCODE SCIPsyncstoreRelease(
   SCIP_SYNCSTORE**         syncstore            /**< pointer to the synchronization store */
   )
{
   int references;

   assert(syncstore != NULL);
   assert(*syncstore != NULL);

   SCIP_CALL( SCIPtpiAcquireLock(&(*syncstore)->lock) );
   (*syncstore)->references -= 1;
   references = (*syncstore)->references;
   SCIP_CALL( SCIPtpiReleaseLock(&(*syncstore)->lock) );

   if( references == 0 )
   {
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
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   SCIP_CALL( SCIPtpiAcquireLock(&syncstore->lock) );

   ++(syncstore->references);

   SCIP_CALL( SCIPtpiReleaseLock(&syncstore->lock) );

   return SCIP_OKAY;
}

/** initialize the syncstore for the given SCIP instance */
SCIP_RETCODE SCIPsyncstoreInit(
   SCIP*                    scip                 /**< SCIP main datastructure */
   )
{
   SCIP_SYNCSTORE* syncstore;
   int nvars;
   int i;
   int j;
   int paramode;

   assert(scip != NULL);
   syncstore = SCIPgetSyncstore(scip);
   assert(syncstore != NULL);
   syncstore->mainscip = scip;
   syncstore->lastsync = NULL;
   syncstore->nsolvers = SCIPgetNConcurrentSolvers(scip);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsols", &syncstore->maxnsols) );
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &syncstore->maxnsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/minsyncdelay", &syncstore->minsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqinit", &syncstore->syncfreqinit) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqmax", &syncstore->syncfreqmax) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqfactor", &syncstore->syncfreqfactor) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/targetprogress", &syncstore->targetprogress) );
   syncstore->nsyncdata = getNSyncdata(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &(syncstore->syncdata), syncstore->nsyncdata) );
   for( i = 0; i < syncstore->nsyncdata; ++i )
   {
      syncstore->syncdata[i].syncnum = -1;
      SCIP_CALL( SCIPboundstoreCreate(syncstore->mainscip, &syncstore->syncdata[i].boundstore, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solobj, syncstore->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].solsource, syncstore->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols, syncstore->maxnsols) );
      for( j = 0; j < syncstore->maxnsols; ++j )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols[j], nvars) );
      }
      SCIP_CALL( SCIPtpiInitLock(&(syncstore->syncdata[i].lock)) );
      SCIP_CALL( SCIPtpiInitCondition(&(syncstore->syncdata[i].allsynced)) );
   }

   syncstore->initialized = TRUE;
   syncstore->stopped = FALSE;

   SCIP_CALL( SCIPgetIntParam(scip, "parallel/mode", &paramode) );
   syncstore->mode = (SCIP_PARALLELMODE) paramode;

   SCIP_CALL( SCIPtpiInit(syncstore->nsolvers, INT_MAX, FALSE) );
   SCIP_CALL( SCIPautoselectDisps(scip) );

   return SCIP_OKAY;
}

/** deinitializes the synchronization store */
SCIP_RETCODE SCIPsyncstoreExit(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
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

      for(j = 0; j < syncstore->maxnsols; ++j)
      {
         SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols[j], SCIPgetNVars(syncstore->mainscip));
      }
      SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata[i].sols, syncstore->maxnsols);
   }

   SCIPfreeBlockMemoryArray(syncstore->mainscip, &syncstore->syncdata, syncstore->nsyncdata);

   syncstore->initialized = FALSE;
   syncstore->stopped = FALSE;

   return SCIP_OKAY;
}

/** initialize the synchronization timing parameters for the first synchronization */
void SCIPsyncstoreInitSyncTiming(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_Real                time                 /**< the time the solver spent before the first synchronization */
   )
{
   assert(syncstore->syncdata[0].syncnum == 0);

   syncstore->syncdata[0].syncfreq = MAX(syncstore->syncdata[0].syncfreq, time);

   if( syncstore->syncdata[0].syncedcount == syncstore->nsolvers - 1 )
   {
      syncstore->minsyncdelay *= syncstore->syncdata[0].syncfreq;
      syncstore->syncfreqmax *= syncstore->syncdata[0].syncfreq;
      syncstore->minsyncdelay *= syncstore->syncdata[0].syncfreq;
      syncstore->syncdata[0].syncfreq *= syncstore->syncfreqinit;
   }
}

/** checks whether the solve-is-stopped flag in the syncstore has been set by any thread */
SCIP_Bool SCIPsyncstoreSolveIsStopped(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   SCIP_Bool stopped;

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&syncstore->lock) );

   stopped = syncstore->stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&syncstore->lock) );

   return stopped;
}

/** sets the solve-is-stopped flag in the SPI so that subsequent calls to
 *  SCIPsyncstoreSolveIsStopped will return the given value in any thread
 */
void SCIPsyncstoreSetSolveIsStopped(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_Bool                stopped              /**< flag if the solve is stopped */
   )
{
   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&syncstore->lock) );

   syncstore->stopped = stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&syncstore->lock) );
}

/** gets the upperbound from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastUpperbound(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->lastsync == NULL ? SCIPinfinity(syncstore->mainscip) : syncstore->lastsync->bestupperbound;
}

/** gets the lowerbound from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastLowerbound(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->lastsync == NULL ? -SCIPinfinity(syncstore->mainscip) : syncstore->lastsync->bestlowerbound;
}

/** gets the number of solutions from the last synchronization */
int SCIPsyncstoreGetLastNSols(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->lastsync == NULL ? 0 : syncstore->lastsync->nsols;
}

/** gets the number of boundchanges from the last synchronization */
int SCIPsyncstoreGetLastNBounds(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->lastsync == NULL ? 0 : SCIPboundstoreGetNChgs(syncstore->lastsync->boundstore);
}

/** gets total memory used by all solvers from the last synchronization */
SCIP_Longint SCIPsyncstoreGetLastMemTotal(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->lastsync == NULL ? 0 : syncstore->lastsync->memtotal;
}

/** gets the synchronization frequency from the last synchronization */
SCIP_Real SCIPsyncstoreGetLastSyncfreq(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->lastsync == NULL ? 0.0 : syncstore->lastsync->syncfreq;
}

/** get synchronization data with given number. It is the responsibility of the caller
 *  to only ask for a synchronization number that still exists, which is checked
 *  with an assert in debug mode. */
SCIP_SYNCDATA* SCIPsyncstoreGetSyncdata(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_Longint             syncnum              /**< the number of the synchronization to start, which
                                                  *   must be increasing between calls of the same thread */
   )
{
   int j = syncnum % syncstore->nsyncdata;

   /* check if requested syncnumber still exists if in debug mode */
   assert( syncstore->syncdata[j].syncnum == syncnum );

   return &syncstore->syncdata[j];
}

/** get the next synchronization data that should be read and
 *  adjust the delay. Returns NULL if no more data should be read due to minimum delay */
SCIP_SYNCDATA* SCIPsyncstoreGetNextSyncdata(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data */
   SCIP_Longint             writenum,            /**< number of synchronizations the solver has written to */
   SCIP_Real*               delay                /**< pointer holding the current synchronization delay */
   )
{
   SCIP_Real newdelay;
   SCIP_Longint nextsyncnum;

   nextsyncnum = syncdata->syncnum + 1;
   newdelay = *delay - syncdata->syncfreq;

   /* if the delay would get too small we dont want to read the next syncdata.
    * But due to the limited length of the syncdata array we might need to
    * read this synchronization data anyways which is checked by the second part
    * of the if condition
    */
   if( newdelay < syncstore->minsyncdelay && nextsyncnum >= writenum - syncstore->maxnsyncdelay )
      return NULL;

   *delay = newdelay;

   return &syncstore->syncdata[nextsyncnum % syncstore->nsyncdata];
}

/** ensures that the given synchronization data has been written by
 *  all solvers upon return of this function and blocks the caller if necessary. */
SCIP_RETCODE SCIPsyncstoreEnsureAllSynced(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   /* check if waiting is required, make sure to hold the lock */
   SCIP_CALL( SCIPtpiAcquireLock(&syncdata->lock) );
   while( syncdata->syncedcount < syncstore->nsolvers )
   {
      /* yes, so wait on the condition variable
       * (automatically releases the lock and reacquires it after the waiting)
       */
      SCIP_CALL( SCIPtpiWaitCondition(&syncdata->allsynced, &syncdata->lock) );
   }
   SCIP_CALL( SCIPtpiReleaseLock(&syncdata->lock) );

   return SCIP_OKAY;
}

/** Start synchronization for the given concurrent solver.
 *  Needs to be followed by a call to SCIPsyncstoreFinishSync if
 *  the syncdata that is returned is not NULL
 */
SCIP_RETCODE SCIPsyncstoreStartSync(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_Longint             syncnum,             /**< the number of the synchronization to start, which
                                                  *   must be increasing between calls of the same thread */
   SCIP_SYNCDATA**          syncdata             /**< pointer to return the synchronization data */
   )
{
   int i;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   if( SCIPsyncstoreSolveIsStopped(syncstore) )
   {
      *syncdata = NULL;
      return SCIP_OKAY;
   }

   i = syncnum % syncstore->nsyncdata;
   *syncdata = &syncstore->syncdata[i];
   assert( *syncdata != NULL );

   SCIP_CALL( SCIPtpiAcquireLock(&(*syncdata)->lock) );

   if((*syncdata)->syncnum != syncnum )
   {
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

/** finishes synchronization for the synchronization data */
SCIP_RETCODE SCIPsyncstoreFinishSync(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_SYNCDATA**          syncdata             /**< the synchronization data */
   )
{
   assert(syncstore != NULL);
   assert(syncstore->initialized);

   ++(*syncdata)->syncedcount;

   if( (*syncdata)->syncedcount == syncstore->nsolvers )
   {
      if( (*syncdata)->status != SCIP_STATUS_UNKNOWN )
         SCIPsyncstoreSetSolveIsStopped(syncstore, TRUE);
      syncstore->lastsync = *syncdata;
      SCIP_CALL( SCIPtpiBroadcastCondition(&(*syncdata)->allsynced) );
   }
   SCIP_CALL( SCIPtpiReleaseLock(&(*syncdata)->lock) );

   if( *syncdata == syncstore->lastsync )
   {
      SCIP_CALL( SCIPprintDisplayLine(syncstore->mainscip, NULL, SCIP_VERBLEVEL_HIGH, TRUE) );
   }
   *syncdata = NULL;

   return SCIP_OKAY;
}

/** gets status in synchronization data */
SCIP_STATUS SCIPsyncdataGetStatus(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->status;
}

/** gets the solver that had the best status, or -1 if solve is not stopped yet */
int SCIPsyncstoreGetWinner(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   if( syncstore->lastsync == NULL || syncstore->lastsync->status == SCIP_STATUS_UNKNOWN )
      return -1;

   return syncstore->lastsync->winner;
}

/** how many solvers have already finished synchronizing on this sychronization data */
int SCIPsyncdataGetNSynced(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->syncedcount;
}

/** how many solvers have are running concurrently */
int SCIPsyncstoreGetNSolvers(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->nsolvers;
}


/** read amount total memory used from synchronization data */
SCIP_Longint SCIPsyncdataGetMemTotal(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->memtotal;
}

/** read the synchronization frequency from a synchronization data */
SCIP_Real SCIPsyncdataGetSyncFreq(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->syncfreq;
}

/** read the upperbound stored in a synchronization data */
SCIP_Real SCIPsyncdataGetUpperbound(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->bestupperbound;
}

/** read the lowerbound stored in a synchronization data */
SCIP_Real SCIPsyncdataGetLowerbound(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->bestlowerbound;
}

/** read the solutions stored in a synchronization data */
void SCIPsyncdataGetSolutions(
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data */
   SCIP_Real***             solvalues,           /**< array of buffers containing the solution values */
   int**                    solowner,            /**< array of ownerids of solutions */
   int*                     nsols                /**< pointer to return number of solutions */
   )
{
   *solvalues = syncdata->sols;
   *solowner = syncdata->solsource;
   *nsols = syncdata->nsols;
}

/** read bound changes stored in the synchronization data */
SCIP_BOUNDSTORE* SCIPsyncdataGetBoundChgs(
   SCIP_SYNCDATA*           syncdata             /**< the synchronization data */
   )
{
   return syncdata->boundstore;
}

/** write the synchronization frequency to a synchronization data */
void SCIPsyncdataSetSyncFreq(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data */
   SCIP_Real                syncfreq             /**< the synchronization frequency */
   )
{
   syncdata->syncfreq = MIN(syncfreq, syncstore->syncfreqmax);
}

/** set status in the synchronization data */
void SCIPsyncdataSetStatus(
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data the upperbound should be added to */
   SCIP_STATUS              status,              /**< the status */
   int                      solverid             /**< identifier of te solver that has this status */
   )
{
   /* check if status is better than current one (closer to SCIP_STATUS_OPTIMAL),
    * break ties by the solverid, and remember the solver wit the best status
    * so that the winner will be selected deterministically */
   if( syncdata->status < SCIP_STATUS_OPTIMAL )
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
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data the solution should be added to */
   SCIP_Longint             memtotal             /**< the number of bytes used */
   )
{
   syncdata->memtotal += memtotal;
}

/** set upperbound to the synchronization data */
void SCIPsyncdataSetUpperbound(
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data the upperbound should be added to */
   SCIP_Real                upperbound           /**< the upperbound */
   )
{
   syncdata->bestupperbound = MIN(syncdata->bestupperbound, upperbound);
}

/** set lowerbound to the synchronization data */
void SCIPsyncdataSetLowerbound(
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data the lowerbound should be added to */
   SCIP_Real                lowerbound           /**< the lowerbound */
   )
{
   syncdata->bestlowerbound = MAX(syncdata->bestlowerbound, lowerbound);
}

/** gives a buffer to store the solution values, or NULL if solution should not be stored
 *  because there are already better solutions stored.
 */
void SCIPsyncdataGetSolutionBuffer(
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data the solution should be added to */
   SCIP_Real                solobj,              /**< the objective value of the solution */
   int                      ownerid,             /**< an identifier for the owner of the solution, e.g. the thread number */
   SCIP_Real**              buffer               /**< pointer to return a buffer for the solution values, which must be set
                                                  *   if the buffer is not NULL */
   )
{
   int                  pos;
   int                  i;

   assert(syncstore != NULL);
   assert(syncstore->initialized);

   for( pos = 0; pos < syncdata->nsols; ++pos )
   {
      if(syncdata->solobj[pos] < solobj || (syncdata->solobj[pos] == solobj && ownerid < syncdata->solsource[pos]  ) )
         break;
   }

   if(syncdata->nsols < syncstore->maxnsols)
   {
      for( i = syncdata->nsols; i > pos; --i )
      {
         syncdata->solobj[i] = syncdata->solobj[i-1];
         syncdata->solsource[i] = syncdata->solsource[i-1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i-1]);
      }
      ++syncdata->nsols;
   }
   else
   {
      --pos;
      for( i = 0; i < pos; ++i )
      {
         syncdata->solobj[i] = syncdata->solobj[i+1];
         syncdata->solsource[i] = syncdata->solsource[i+1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i+1]);
      }
   }

   if(pos >= 0)
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
   SCIP_SYNCSTORE*          syncstore,           /**< the synchronization store */
   SCIP_SYNCDATA*           syncdata,            /**< the synchronization data */
   SCIP_BOUNDSTORE*         boundstore           /**< bound store containing the bounds to add */
   )
{
   SCIP_CALL( SCIPboundstoreMerge(syncstore->mainscip, syncdata->boundstore, boundstore) );

   return SCIP_OKAY;
}

/** is synchronization store initialized */
SCIP_Bool SCIPsyncstoreIsInitialized(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->initialized;
}

/** returns the mode of the synchronization store */
SCIP_PARALLELMODE SCIPsyncstoreGetMode(
   SCIP_SYNCSTORE*          syncstore            /**< the synchronization store */
   )
{
   return syncstore->mode;
}




#if 0




/** creates and captures a syncstore object */
SCIP_RETCODE SCIPsyncstoreCreate(
   SCIP_SYNCSTORE**                 spi            /**< pointer to an parallel interface structure */
   )
{
   assert(spi != NULL);

   SCIPdebugMessage("SCIPsyncstoreCreate()\n");

   SCIP_ALLOC( BMSallocMemory(spi) );

   (*spi)->mode = SCIP_PARA_DETERMINISTIC;                      /* initialising the mode */
   (*spi)->initialized = FALSE;
   (*spi)->syncdata = NULL;
   (*spi)->stopped = FALSE;
   (*spi)->references = 1;
   SCIP_CALL( SCIPtpiInitLock(&(*spi)->lock) );

   return SCIP_OKAY;
}



/** deletes the parallel interface object */
SCIP_RETCODE SCIPsyncstoreRelease(
   SCIP_SYNCSTORE**                 spi            /**< pointer to an parallel interface structure */
   )
{
   int references;
   assert(spi != NULL);
   assert(*spi != NULL);

   SCIP_CALL( SCIPtpiAcquireLock(&(*spi)->lock) );
   (*spi)->references -= 1;
   references = (*spi)->references;
   SCIP_CALL( SCIPtpiReleaseLock(&(*spi)->lock) );

   if( references == 0 )
   {
      assert(!(*spi)->initialized);
      SCIPtpiDestroyLock(&(*spi)->lock);
      BMSfreeMemory(spi);
   }
   else
   {
      *spi = NULL;
   }

   return SCIP_OKAY;
}

/** copies pointer to spi */
SCIP_RETCODE SCIPsyncstoreCapture(
   SCIP_SYNCSTORE*                 syncstore       /**< syncstore to capture */
   )
{
   SCIP_CALL( SCIPtpiAcquireLock(&syncstore->lock) );

   ++(syncstore->references);

   SCIP_CALL( SCIPtpiReleaseLock(&syncstore->lock) );

   return SCIP_OKAY;
}

/** checks whether the solve-is-stopped flag in the SPI has been set by any thread */
SCIP_Bool SCIPsyncstoreSolveIsStopped(
   SCIP_SYNCSTORE*                 spi             /**< pointer to an parallel interface structure */
   )
{
   SCIP_Bool stopped;

   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&spi->lock) );

   stopped = spi->stopped;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&spi->lock) );

   return stopped;
}

/** sets the solve-is-stopped flag in the SPI so that subsequent calls to
 *  SCIPsyncstoreSolveIsStopped will return TRUE in any thread */
void SCIPsyncstoreSignalSolveIsStopped(
   SCIP_SYNCSTORE*                 spi             /**< pointer to a parallel interface structure */
   )
{
   SCIP_CALL_ABORT( SCIPtpiAcquireLock(&spi->lock) );

   spi->stopped = TRUE;

   SCIP_CALL_ABORT( SCIPtpiReleaseLock(&spi->lock) );
}

/** gets the upperbound in the last synchronization */
SCIP_Real SCIPsyncstoreGetLastUpperbound(
   SCIP_SYNCSTORE*                 spi             /**< pointer to a parallel interface structure */
   )
{
   return spi->lastsync == NULL ? SCIPinfinity(spi->mainscip) : spi->lastsync->bestprimalbound;
}

/** gets the lowerbound in the last synchronization */
SCIP_Real SCIPsyncstoreGetLastLowerbound(
   SCIP_SYNCSTORE*                 spi             /**< pointer to a parallel interface structure */
   )
{
   return spi->lastsync == NULL ? -SCIPinfinity(spi->mainscip) : spi->lastsync->bestdualbound;
}

/** gets the number of solutions in the last synchronization */
int SCIPsyncstoreGetLastNSols(
   SCIP_SYNCSTORE*    spi                          /**< pointer to a parallel interface structure */
   )
{
   return spi->lastsync == NULL ? 0 : spi->lastsync->nsols;
}

/** gets the number of boundchanges in the last synchronization */
int SCIPsyncstoreGetLastNBounds(
   SCIP_SYNCSTORE*    spi                          /**< pointer to a parallel interface structure */
   )
{
   return spi->lastsync == NULL ? 0 : SCIPboundstoreGetNChgs(spi->lastsync->boundstore);
}

/** gets total memory used by all solvers at the last synchronization */
SCIP_Longint SCIPsyncstoreGetLastMemUsed(
   SCIP_SYNCSTORE*    spi                          /**< pointer to a parallel interface structure */
   )
{
   return spi->lastsync == NULL ? 0 : spi->lastsync->memused;
}

/** gets the synchronization frequency in the last synchronization */
SCIP_Real SCIPsyncstoreGetLastSyncfreq(
   SCIP_SYNCSTORE*    spi                          /**< pointer to a parallel interface structure */
   )
{
   return spi->lastsync == NULL ? 0.0 : spi->lastsync->syncfreq;
}

/** returns the current synchronization frequency */
SCIP_Longint SCIPsyncstoreGetMemUsed(
   SCIP_CONCSOLVER*          concsolver      /**< pointer to a concurrent solver */
   )
{
   return concsolver->syncdata == NULL ? 0 : concsolver->syncdata->memused;
}

/** returns the current synchronization frequency */
SCIP_Real SCIPsyncstoreGetSyncFreq(
   SCIP_CONCSOLVER*          concsolver      /**< pointer to a concurrent solver */
   )
{
   return concsolver->syncdata == NULL ? 0.0 : concsolver->syncdata->syncfreq;
}

static SCIP_RETCODE execConcsolver(void *args)
{
   SCIP* scip = (SCIP*) args;
   SCIP_CALL( SCIPexecConcurrentSolver(scip, SCIPtpiGetThreadNum()) );
   return SCIP_OKAY;
}

/** start solving in parallel using the given set of concurrent solvers */
SCIP_RETCODE SCIPsyncstoreSolvePortfolio(
   SCIP*                scip                 /**< pointer to scip datastructure */
   )
{
   SCIP_SYNCSTORE* spi;
   int idx;
   SCIP_STATUS beststatus;
   int jobid;
   int i;
   SCIP_RETCODE retcode;
   SCIP_CONCSOLVER*       * concsolvers;

   spi = SCIPgetSyncstore(scip);

   concsolvers = SCIPgetConcurrentSolvers(scip);
   spi->stopped = FALSE;
   jobid = SCIPtpiGetNewJobID();
   SPI_PARA
   {
      SPI_MASTER
      {
         for( i = 0; i < spi->nsolvers; ++i )
         {
            SCIP_JOB*         job;
            SCIP_SUBMITSTATUS status;

            SCIP_CALL_ABORT( SCIPtpiCreateJob(&job, jobid, execConcsolver, scip ) );

            SCIP_CALL_ABORT( SCIPtpiSumbitJob(job, &status) );
            if( status != SCIP_SUBMIT_SUCCESS )
               retcode = SCIP_ERROR;
         }
         assert(spi != NULL);
      }
   }

   retcode = SCIPtpiCollectJobs(jobid);
   idx = 0;
   beststatus = SCIP_STATUS_UNKNOWN;
   for( i = 0; i < spi->nsolvers; ++i )
   {
      SCIP_STATUS solverstatus;
      solverstatus = SCIPconcsolverGetStatus(concsolvers[i]);
      if( (beststatus < SCIP_STATUS_OPTIMAL && solverstatus > beststatus) ||
          (beststatus > SCIP_STATUS_OPTIMAL && solverstatus >= SCIP_STATUS_OPTIMAL && solverstatus < beststatus) )
      {
         idx = i;
         beststatus = solverstatus;
      }
   }
   assert(beststatus != SCIP_STATUS_UNKNOWN);

   SCIP_CALL( SCIPconcsolverGetSolvingData(concsolvers[idx], scip) );

   return retcode;
}

/** Start synchronization for the given concurrent solver.
 * Needs to be followed by a call to SCIPsyncstoreFinishSync if
 * the syncdata that is returned is not NULL */
SCIP_RETCODE SCIPsyncstoreStartSync(
   SCIP_SYNCSTORE*            spi,                 /**< pointer to an parallel interface structure */
   SCIP_SYNCDATA** syncdata,           /**< pointer to store the synchronization data to use */
   SCIP_CONCSOLVER*             concsolver            /**< pointer to concurrent solver */
   )
{
   SCIP_Longint syncnum;
   int i;

   assert(spi != NULL);
   assert(spi->initialized);

   if( SCIPsyncstoreSolveIsStopped(spi) )
   {
      *syncdata = NULL;
      return SCIP_OKAY;
   }

   syncnum = (concsolver->nsyncs)++;
   i = syncnum % spi->nsyncdata;
   *syncdata = &spi->syncdata[i];
   assert( *syncdata != NULL );
   concsolver->syncdelay += concsolver->timesincelastsync;

   SCIP_CALL( SCIPtpiAcquireLock(&(*syncdata)->lock) );

   if((*syncdata)->syncnum != syncnum )
   {
      SCIPboundstoreClear((*syncdata)->boundstore);
      (*syncdata)->nsols = 0;
      (*syncdata)->memused = SCIPgetMemUsed(spi->mainscip);
      (*syncdata)->syncedcount = 0;
      (*syncdata)->bestprimalbound = SCIPconcsolverGetPrimalbound(concsolver);
      (*syncdata)->bestdualbound = SCIPconcsolverGetDualbound(concsolver);
      (*syncdata)->status = SCIP_STATUS_UNKNOWN;
      (*syncdata)->syncnum = syncnum;
      if( syncnum == 0 )
      {
         (*syncdata)->syncfreq = concsolver->timesincelastsync;
      }
   }
   else
   {
      (*syncdata)->bestdualbound = MAX(SCIPconcsolverGetDualbound(concsolver), (*syncdata)->bestdualbound);
      (*syncdata)->bestprimalbound = MIN(SCIPconcsolverGetPrimalbound(concsolver), (*syncdata)->bestprimalbound);

      if( syncnum == 0 )
      {
         (*syncdata)->syncfreq = MAX(concsolver->timesincelastsync, (*syncdata)->syncfreq);
      }
   }

   return SCIP_OKAY;
}

/** Finishes synchronization for the given concurrent solver */
SCIP_RETCODE SCIPsyncstoreFinishSync(
   SCIP_SYNCSTORE*            spi,                 /**< pointer to an parallel interface structure */
   SCIP_SYNCDATA** syncdata,           /**< pointer to the synchronization data used */
   SCIP_CONCSOLVER*     concsolver           /**< pointer to concurrent solver */
   )
{
   SCIP_Longint          syncnum;
   SCIP_SYNCDATA*   currentsyncdata;
   int                   i;
   SCIP_Longint          s;
   int                   j;
   SCIP_Real             newdelay;
   SCIP_STATUS           concsolverstatus;

   assert(spi != NULL);
   assert(spi->initialized);

   concsolverstatus = SCIPconcsolverGetStatus(concsolver);

   (*syncdata)->status = MAX(concsolverstatus, (*syncdata)->status);

   if( concsolver->syncdata )
      syncnum = concsolver->syncdata->syncnum+1;
   else
      syncnum = 0;

   ++(*syncdata)->syncedcount;

   if( (*syncdata)->syncedcount == spi->nsolvers )
   {
      if( (*syncdata)->syncnum == 0 )
      {
         spi->minsyncdelay *= (*syncdata)->syncfreq;
         spi->syncfreqmax *= (*syncdata)->syncfreq;
         (*syncdata)->syncfreq *= spi->syncfreqinit;
         SCIPdebugMessage("chose initial syncfreq %g and delay %g\n", (*syncdata)->syncfreq, spi->minsyncdelay);
      }
      else
      {
         SCIP_Bool dbok;
         SCIP_Bool pbok;
         SCIP_Real progress;
         SCIP_Real prevpb;
         SCIP_Real prevdb;
         SCIP_Real newpb;
         SCIP_Real newdb;
         SCIP_Real freqfactor;

         j = ((*syncdata)->syncnum - 1) % spi->nsyncdata;
         prevpb = spi->syncdata[j].bestprimalbound;
         prevdb = spi->syncdata[j].bestdualbound;
         newpb = (*syncdata)->bestprimalbound;
         newdb = (*syncdata)->bestdualbound;
         dbok = !(SCIPisInfinity(spi->mainscip, prevdb) || SCIPisInfinity(spi->mainscip, -prevdb));
         pbok = !(SCIPisInfinity(spi->mainscip, prevpb) || SCIPisInfinity(spi->mainscip, -prevpb));
         if( dbok && pbok )
            progress = REALABS(SCIPrelDiff(prevpb - prevdb, newpb - newdb));
         else if( dbok )
            progress = REALABS(SCIPrelDiff(prevdb, newdb));
         else if( pbok )
            progress = REALABS(SCIPrelDiff(prevpb, newpb));
         else if( (SCIPisInfinity(spi->mainscip, newdb) || SCIPisInfinity(spi->mainscip, -newdb)) &&
                  (SCIPisInfinity(spi->mainscip, newpb) || SCIPisInfinity(spi->mainscip, -newpb)) )
            progress = 0.0;
         else
            progress = 1.0;

         /*TODO: targetprogress to spi and freqfactor too */
         if( progress < 0.5 * spi->targetprogress )
            freqfactor = spi->syncfreqfactor;
         else if( progress > 2 * spi->targetprogress )
            freqfactor = 0.5 + 0.5/spi->syncfreqfactor;
         else
            freqfactor = 1.0;

         (*syncdata)->syncfreq = MIN(spi->syncdata[j].syncfreq * freqfactor, spi->syncfreqmax);
         SCIPdebugMessage("new syncfreq is %g\n", (*syncdata)->syncfreq);
      }

      if( (*syncdata)->status != SCIP_STATUS_UNKNOWN )
         SCIPsyncstoreSignalSolveIsStopped(spi);
      spi->lastsync = *syncdata;
      SCIP_CALL( SCIPtpiBroadcastCondition(&(*syncdata)->allsynced) );
   }
   else if(  (*syncdata)->status != SCIP_STATUS_UNKNOWN && concsolverstatus == SCIP_STATUS_UNKNOWN )
   {
      SCIP_CALL( SCIPconcsolverStop(concsolver) );
   }
   SCIP_CALL( SCIPtpiReleaseLock(&(*syncdata)->lock) );

   if( *syncdata == spi->lastsync )
   {
      SCIP_CALL( SCIPprintDisplayLine(spi->mainscip, NULL, SCIP_VERBLEVEL_HIGH, TRUE) );
   }
   *syncdata = NULL;

   if( concsolverstatus != SCIP_STATUS_UNKNOWN || SCIPsyncstoreSolveIsStopped(spi) )
      return SCIP_OKAY;

   for( s = syncnum; s < concsolver->nsyncs; ++s )
   {
      newdelay = concsolver->syncdelay;
      if( s == 0 )
      {
         newdelay -= concsolver->timesincelastsync;
      }
      else
      {
         newdelay -= spi->syncdata[((s-1) % spi->nsyncdata)].syncfreq;
      }
      if( concsolver->syncdata != NULL && newdelay < spi->minsyncdelay && s >= concsolver->nsyncs - spi->maxnsyncdelay )
         return SCIP_OKAY;

      j = s % spi->nsyncdata;
      currentsyncdata = &(spi->syncdata[j]);
      concsolver->syncdata = currentsyncdata;
      assert(currentsyncdata->syncnum == s);
      concsolver->syncdelay = newdelay;

      /* check if waiting in required */
      SCIP_CALL( SCIPtpiAcquireLock(&currentsyncdata->lock) );
      while( currentsyncdata->syncedcount < spi->nsolvers )
      {
         SCIP_CALL( SCIPtpiWaitCondition(&currentsyncdata->allsynced, &currentsyncdata->lock) );
      }
      SCIP_CALL( SCIPtpiReleaseLock(&currentsyncdata->lock) );

      if( SCIPsyncstoreSolveIsStopped(spi) )
      {
         /*SCIP_CALL( SCIPconcsolverStop(concsolver) );*/
         return SCIP_OKAY;
      }

      for( i = 0; i < currentsyncdata->nsols; ++i )
      {
         if( SCIPconcsolverGetIdx(concsolver) == currentsyncdata->solsource[i] )
            continue;

         SCIP_CALL( SCIPconcsolverAddSol(concsolver, currentsyncdata->sols[i]) );
      }

      if( SCIPboundstoreGetNChgs(currentsyncdata->boundstore) > 0 )
         SCIP_CALL( SCIPconcsolverAddVarBounds(concsolver, currentsyncdata->boundstore) );

      if( currentsyncdata->status != SCIP_STATUS_UNKNOWN )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** adds memory used to the synchronization data */
void SCIPsyncstoreSyncAddMemUsed(
   SCIP_SYNCDATA*      syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Real                memused              /**< the number of bytes used */
   )
{
   syncdata->memused += memused;
}

/** Add a solution to the synchronization data */
SCIP_RETCODE SCIPsyncstoreSyncAddSolution(
   SCIP_SYNCSTORE*                spi,                 /**< pointer to an parallel interface structure */
   SCIP_SYNCDATA*      syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_CONCSOLVER*         concsolver,          /**< the concurrent solver the solution is coming from */
   SCIP_CONCSOLVERSOL*      concsol              /**< solution of the concurrent solver that should be added to synchronization data */
   )
{
   SCIP_Real            solobj;
   int                  pos;
   int                  i;
   int                  idx;

   assert(spi != NULL);
   assert(spi->initialized);

   solobj = SCIPconcsolGetObj(concsolver, concsol);
   idx = SCIPconcsolverGetIdx(concsolver);
   for( pos = 0; pos < syncdata->nsols; ++pos )
   {
      if(syncdata->solobj[pos] < solobj || (syncdata->solobj[pos] == solobj && idx < syncdata->solsource[pos]  ) )
         break;
   }

   if(syncdata->nsols < spi->maxnsols)
   {
      for( i = syncdata->nsols; i > pos; --i )
      {
         syncdata->solobj[i] = syncdata->solobj[i-1];
         syncdata->solsource[i] = syncdata->solsource[i-1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i-1]);
      }
      ++syncdata->nsols;
   }
   else
   {
      --pos;
      for( i = 0; i < pos; ++i )
      {
         syncdata->solobj[i] = syncdata->solobj[i+1];
         syncdata->solsource[i] = syncdata->solsource[i+1];
         SCIPswapPointers((void**) &syncdata->sols[i], (void**) &syncdata->sols[i+1]);
      }
   }

   if(pos >= 0)
   {
      syncdata->solobj[pos] = solobj;
      syncdata->solsource[pos] = idx;
      SCIP_CALL( SCIPconcsolGetVals(concsolver, concsol, syncdata->sols[pos]) );
   }

   return SCIP_OKAY;
}

/** Adds bound changes to the synchronization data */
SCIP_RETCODE SCIPsyncstoreSyncAddBoundChanges(
   SCIP_SYNCSTORE*              spi,                 /**< pointer to an parallel interface structure */
   SCIP_SYNCDATA*    syncdata,           /**< the synchronization data the bounds should be added to */
   SCIP_SYNCSTORE_BOUNDSTORAGE* boundstore         /**< bound storage containing the bounds to add */
   )
{
   SCIP_CALL( SCIPboundStorageMerge(spi->mainscip, syncdata->boundstore, boundstore) );

   return SCIP_OKAY;
}

/** Is synchronization enabled, i.e. has the SPI been initialized */
SCIP_Bool SCIPsyncstoreIsInitialized(
   SCIP_SYNCSTORE*              spi                  /**< pointer to an parallel interface structure */
   )
{
   return spi->initialized;
}

/** initialize the spi for the given number of threads */
SCIP_RETCODE SCIPsyncstoreInit(
   SCIP*                  scip                 /**< pointer to scip datastructure */
   )
{
   SCIP_SYNCSTORE* spi;
   int nvars;
   int i;
   int j;
   int paramode;
   SCIP_Bool blockwhenfull;

   assert(scip != NULL);
   spi = SCIPgetSPI(scip);
   assert(spi != NULL);
   spi->mainscip = scip;
   spi->lastsync = NULL;
   spi->nsolvers = SCIPgetNConcurrentSolvers(scip);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsols", &spi->maxnsols) );
   SCIP_CALL( SCIPgetIntParam(scip, "concurrent/sync/maxnsyncdelay", &spi->maxnsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/minsyncdelay", &spi->minsyncdelay) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqinit", &spi->syncfreqinit) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqmax", &spi->syncfreqmax) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/freqfactor", &spi->syncfreqfactor) );
   SCIP_CALL( SCIPgetRealParam(scip, "concurrent/sync/targetprogress", &spi->targetprogress) );
   spi->nsyncdata = getNSyncdata(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(spi->mainscip, &(spi->syncdata), spi->nsyncdata) );
   for( i = 0; i < spi->nsyncdata; ++i )
   {
      spi->syncdata[i].syncnum = -1;
      SCIP_CALL( SCIPboundStorageCreate(spi->mainscip, &spi->syncdata[i].boundstore, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(spi->mainscip, &spi->syncdata[i].solobj, spi->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(spi->mainscip, &spi->syncdata[i].solsource, spi->maxnsols) );
      SCIP_CALL( SCIPallocBlockMemoryArray(spi->mainscip, &spi->syncdata[i].sols, spi->maxnsols) );
      for( j = 0; j < spi->maxnsols; ++j )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(spi->mainscip, &spi->syncdata[i].sols[j], nvars) );
      }
      SCIP_CALL( SCIPtpiInitLock(&(spi->syncdata[i].lock)) );
      SCIP_CALL( SCIPtpiInitCondition(&(spi->syncdata[i].allsynced)) );
   }

   spi->initialized = TRUE;
   spi->stopped = FALSE;

   SCIP_CALL( SCIPgetIntParam(scip, "parallel/mode", &paramode) );
   SCIP_CALL( SCIPgetIntParam(scip, "parallel/queuesize", &spi->queuesize) );
   SCIP_CALL( SCIPgetBoolParam(scip, "parallel/blockqueuewhenfull", &blockwhenfull) );
   spi->mode = (SCIP_PARALLELMODE) paramode;

   SCIP_CALL( SCIPtpiInit(spi->nsolvers, spi->queuesize, blockwhenfull) );
   SCIP_CALL( SCIPautoselectDisps(scip) );

   return SCIP_OKAY;
}

/** deinitializes the spi */
SCIP_RETCODE SCIPsyncstoreExit(
   SCIP_SYNCSTORE*    spi                          /**< pointer to an parallel interface structure */
   )
{
   int i;
   int j;

   assert(spi != NULL);
   assert(spi->initialized);

   SCIP_CALL( SCIPtpiExit() );

   for( i = 0; i < spi->nsyncdata; ++i )
   {
      SCIPtpiDestroyLock(&(spi->syncdata[i].lock));
      SCIPtpiDestroyCondition(&(spi->syncdata[i].allsynced));
      SCIPfreeBlockMemoryArray(spi->mainscip, &spi->syncdata[i].solobj, spi->maxnsols);
      SCIPfreeBlockMemoryArray(spi->mainscip, &spi->syncdata[i].solsource, spi->maxnsols);
      SCIPboundStorageFree(spi->mainscip,  &spi->syncdata[i].boundstore);

      for(j = 0; j < spi->maxnsols; ++j)
      {
         SCIPfreeBlockMemoryArray(spi->mainscip, &spi->syncdata[i].sols[j], SCIPgetNVars(spi->mainscip));
      }
      SCIPfreeBlockMemoryArray(spi->mainscip, &spi->syncdata[i].sols, spi->maxnsols);
   }

   SCIPfreeBlockMemoryArray(spi->mainscip, &spi->syncdata, spi->nsyncdata);

   spi->initialized = FALSE;
   spi->stopped = FALSE;

   return SCIP_OKAY;
}


SCIP_SYNCSTOREMODE SCIPsyncstoreGetMode(
   SCIP_SYNCSTORE*    spi
   )
{
   return spi->mode;
}

#endif
