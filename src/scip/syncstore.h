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

/**@file   syncstore.h
 * @ingroup PARALLEL
 * @brief  the function declarations for the synchronization store
 * @author Leona Gottwald
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SYNCSTORE_H__
#define __SYNCSTORE_H__

#include "scip/def.h"
#include "scip/type_syncstore.h"
#include "scip/type_scip.h"
#include "scip/type_retcode.h"

/** creates and captures a new synchronization store */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreCreate(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to return the created synchronization store */
   );

/** releases a synchronization store */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreRelease(
   SCIP_SYNCSTORE**      syncstore           /**< pointer to the synchronization store */
   );

/** captures a synchronization store */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreCapture(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** initialize the syncstore for the given SCIP instance */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreInit(
   SCIP*                 scip                /**< SCIP main data structure */
   );

/** deinitializes the synchronization store */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreExit(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** checks whether the solve-is-stopped flag in the syncstore has been set by any thread */
SCIP_EXPORT
SCIP_Bool SCIPsyncstoreSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** sets the solve-is-stopped flag in the syncstore so that subsequent calls to
 *  SCIPsyncstoreSolveIsStopped will return the given value in any thread
 */
SCIP_EXPORT
void SCIPsyncstoreSetSolveIsStopped(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Bool             stopped             /**< flag if the solve is stopped */
   );

/** updates the minimization-normalized objective value of the best solution found by any
 *  concurrent solver; used for printing new incumbents and as the improvement filter of the
 *  solution pool. If the solution pool is enabled and solution values are passed, the
 *  improving solution is published in the pool so that the other concurrent solvers can
 *  install it as an incumbent at their next drain. If published is not NULL, it returns whether
 *  the solution was actually added to the pool, so the caller can count it as shared.
 */
SCIP_EXPORT
void SCIPsyncstoreUpdateBestMinObj(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Real             minobj,             /**< objective value in minimization-normalized original objective space */
   const char*           solvername,         /**< name of the concurrent solver that found the solution */
   int                   ownerid,            /**< index of the concurrent solver that found the solution */
   SCIP_Real*            solvals,            /**< solution values in the communication variable order, or NULL to
                                              *   only communicate the objective value */
   int                   nsolvals,           /**< number of solution values */
   SCIP_Bool*            published           /**< pointer to return whether the solution was added to the pool, or NULL */
   );

/** updates the minimization-normalized best global dual bound; the concurrent solvers call this
 *  immediately while solving so that SCIPgetConcurrentDualbound stays up to date in solution-pool mode
 */
SCIP_EXPORT
void SCIPsyncstoreUpdateBestDualbound(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Real             dualbound           /**< dual bound in minimization-normalized original objective space */
   );

/** returns whether the solution pool for immediate solution sharing is enabled */
SCIP_EXPORT
SCIP_Bool SCIPsyncstoreSolPoolEnabled(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** returns the main SCIP that initialized the synchronization store; the bounds reported in
 *  solution-pool mode live in this SCIP's objective space, so it is used to convert them
 */
SCIP_EXPORT
SCIP* SCIPsyncstoreGetMainScip(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the current number of solutions in the solution pool; reads the counter without
 *  locking, so it may lag behind concurrent pushes, which is fine for its use as a
 *  cheap size hint that is polled at every node
 */
SCIP_EXPORT
int SCIPsyncstoreGetNPoolSols(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the solution with the given index from the solution pool; entries are immutable
 *  once published, so the returned values can be read without holding any lock
 */
SCIP_EXPORT
void SCIPsyncstoreGetPoolSol(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   int                   idx,                /**< index of the pooled solution, must be smaller than the size hint */
   SCIP_Real**           solvals,            /**< pointer to return the solution values in communication variable order */
   int*                  nsolvals,           /**< pointer to return the number of solution values */
   int*                  ownerid             /**< pointer to return the index of the contributing concurrent solver */
   );

/** Tries to acquire the synchronization data with the given number for non-blocking reading.
 *  On success the synchronization data is returned in locked state and must be released with
 *  SCIPsyncstoreUnlockSyncdata after reading. If the data is not ready, NULL is returned and
 *  lost indicates why: if the slot was overwritten by a newer synchronization the data is
 *  lost and the reader should skip this number, otherwise the synchronization has not been
 *  written by all solvers yet and the reader should retry later. Never blocks the caller.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreTryLockCompleteSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum,            /**< the number of the synchronization to read */
   SCIP_SYNCDATA**       syncdata,           /**< pointer to return the locked synchronization data, or NULL */
   SCIP_Bool*            lost                /**< pointer to return whether the data was overwritten and is lost */
   );

/** releases the lock of a synchronization data acquired with SCIPsyncstoreTryLockCompleteSyncdata */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreUnlockSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data to unlock */
   );

/** gets the minimization-normalized objective value of the best solution found by any
 *  concurrent solver, or SCIP_REAL_MAX if no solution was registered yet
 */
SCIP_EXPORT
SCIP_Real SCIPsyncstoreGetBestMinObj(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the upperbound from the last synchronization */
SCIP_EXPORT
SCIP_Real SCIPsyncstoreGetLastUpperbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the lowerbound from the last synchronization */
SCIP_EXPORT
SCIP_Real SCIPsyncstoreGetLastLowerbound(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the number of solutions from the last synchronization */
SCIP_EXPORT
int SCIPsyncstoreGetLastNSols(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the number of boundchanges from the last synchronization */
SCIP_EXPORT
int SCIPsyncstoreGetLastNBounds(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets total memory used by all solvers from the last synchronization */
SCIP_EXPORT
SCIP_Longint SCIPsyncstoreGetLastMemTotal(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** gets the synchronization frequency from the last synchronization */
SCIP_EXPORT
SCIP_Real SCIPsyncstoreGetLastSyncfreq(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** get synchronization data with given number. It is the responsibility of the caller
 *  to only ask for a synchronization number that still exists. */
SCIP_EXPORT
SCIP_SYNCDATA* SCIPsyncstoreGetSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum             /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   );

/** get the next synchronization data that should be read and
 *  adjust the delay. Returns NULL if no more data should be read due to minimum delay */
SCIP_EXPORT
SCIP_SYNCDATA* SCIPsyncstoreGetNextSyncdata(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq,           /**< the current synchronization frequency */
   SCIP_Longint          writenum,           /**< number of synchronizations the solver has written to */
   SCIP_Real*            delay               /**< pointer holding the current synchronization delay */
   );

/** ensures that the given synchronization data has been written by
 *  all solvers upon return of this function and blocks the caller if necessary. */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreEnsureAllSynced(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** Start synchronization for the given concurrent solver.
 *  Needs to be followed by a call to SCIPsyncstoreFinishSync if
 *  the syncdata that is returned is not NULL
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreStartSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_Longint          syncnum,            /**< the number of the synchronization to start, which
                                              *   must be increasing between calls of the same thread */
   SCIP_SYNCDATA**       syncdata            /**< pointer to return the synchronization data */
   );

/** finishes synchronization for the synchronization data */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncstoreFinishSync(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA**       syncdata            /**< the synchronization data */
   );

/** gets status in synchronization data */
SCIP_EXPORT
SCIP_STATUS SCIPsyncdataGetStatus(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** gets the solver that had the best status, or -1 if solve is not stopped yet */
SCIP_EXPORT
int SCIPsyncstoreGetWinner(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** how many solvers have already finished synchronizing on this synchronization data */
SCIP_EXPORT
int SCIPsyncdataGetNSynced(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** how many solvers have are running concurrently */
SCIP_EXPORT
int SCIPsyncstoreGetNSolvers(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** read amount of memory used from synchronization data */
SCIP_EXPORT
SCIP_Longint SCIPsyncdataGetMemTotal(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the synchronization frequency from a synchronization data */
SCIP_EXPORT
SCIP_Real SCIPsyncdataGetSyncFreq(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the upperbound stored in a synchronization data */
SCIP_EXPORT
SCIP_Real SCIPsyncdataGetUpperbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the lowerbound stored in a synchronization data */
SCIP_EXPORT
SCIP_Real SCIPsyncdataGetLowerbound(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** read the solutions stored in a synchronization data */
SCIP_EXPORT
void SCIPsyncdataGetSolutions(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real***          solvalues,          /**< pointer to return array of buffers containing the solution values */
   int**                 solowner,           /**< pointer to return array of ownerids of solutions */
   int*                  nsols               /**< pointer to return number of solutions */
   );

/** read bound changes stored in the synchronization data */
SCIP_EXPORT
SCIP_BOUNDSTORE* SCIPsyncdataGetBoundChgs(
   SCIP_SYNCDATA*        syncdata            /**< the synchronization data */
   );

/** write the synchronization frequency to a synchronization data */
SCIP_EXPORT
void SCIPsyncdataSetSyncFreq(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_Real             syncfreq            /**< the synchronization frequency */
   );

/** set status in the synchronization data */
SCIP_EXPORT
void SCIPsyncdataSetStatus(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_STATUS           status,             /**< the status */
   int                   solverid            /**< identifier of the solver that has this status */
   );

/** adds memory used to the synchronization data */
SCIP_EXPORT
void SCIPsyncdataAddMemTotal(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Longint          memtotal            /**< the number of bytes used */
   );

/** set upperbound to the synchronization data */
SCIP_EXPORT
void SCIPsyncdataSetUpperbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the upperbound should be added to */
   SCIP_Real             upperbound          /**< the upperbound */
   );

/** set lowerbound to the synchronization data */
SCIP_EXPORT
void SCIPsyncdataSetLowerbound(
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the lowerbound should be added to */
   SCIP_Real             lowerbound          /**< the lowerbound */
   );

/** gives a buffer to store the solution values, or NULL if solution should not be stored
 *  because there are already better solutions stored.
 */
SCIP_EXPORT
void SCIPsyncdataGetSolutionBuffer(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data the solution should be added to */
   SCIP_Real             solobj,             /**< the objective value of the solution */
   int                   ownerid,            /**< an identifier for the owner of the solution, e.g. the thread number */
   SCIP_Real**           buffer              /**< pointer to return a buffer for the solution values, which must be set
                                              *   if the buffer is not NULL */
   );

/** adds bound changes to the synchronization data */
SCIP_EXPORT
SCIP_RETCODE SCIPsyncdataAddBoundChanges(
   SCIP_SYNCSTORE*       syncstore,          /**< the synchronization store */
   SCIP_SYNCDATA*        syncdata,           /**< the synchronization data */
   SCIP_BOUNDSTORE*      boundstore          /**< bound store containing the bounds to add */
   );

/** is synchronization store initialized */
SCIP_EXPORT
SCIP_Bool SCIPsyncstoreIsInitialized(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

/** returns the mode of the synchronization store */
SCIP_EXPORT
SCIP_PARALLELMODE SCIPsyncstoreGetMode(
   SCIP_SYNCSTORE*       syncstore           /**< the synchronization store */
   );

#endif
