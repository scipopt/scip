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

/**@file   struct_syncstore.h
 * @ingroup PARALLEL
 * @brief  the struct definitions for the synchronization store
 * @author Stephen J. Maher
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SYNCSTORE_H__
#define __STRUCT_SYNCSTORE_H__

#include "scip/type_syncstore.h"
#include "tpi/type_tpi.h"
#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct SCIP_SyncStore
{
   int                   nuses;              /**< number of uses of the synchronization store */
   SCIP_PARALLELMODE     mode;               /**< the mode for the parallel solving */
   SCIP_Bool             initialized;        /**< flag to indicate whether the syncstore has been initialized */
   int                   ninitvars;          /**< number of variables it has been initialized for */
   SCIP_SYNCDATA*        syncdata;           /**< array of size nsyncdata, containing the synchronization data
                                              *   for each active synchronization */
   SCIP_SYNCDATA*        lastsync;           /**< pointer to the last synchronization data that has been synchronized
                                              *   by all threads */

   SCIP*                 mainscip;           /**< the SCIP instance that was used for initializing the syncstore */
   SCIP_Real             limit_gap;          /**< relative gap limit in main SCIP */
   SCIP_Real             limit_absgap;       /**< absolute gap limit in main SCIP */
   SCIP_Bool             stopped;            /**< flag to indicate if the solving is stopped */
   SCIP_LOCK*            lock;               /**< lock to protect the syncstore data structure from data races */

   int                   nsyncdata;          /**< the size of the synchronization data array */
   SCIP_Real             minsyncdelay;       /**< the minimum delay before a synchronization data may be read */
   int                   maxnsyncdelay;      /**< maximum number of synchronizations before the reading of the next
                                              *   synchronization data is enforced regardless of the minimal synchronization
                                              *   delay */
   SCIP_Real             syncfreqinit;       /**< the initial synchronization frequency which is read from the settings
                                              *   of the main SCIP when the syncstore is initialized */
   SCIP_Real             syncfreqmax;        /**< the maximum synchronization frequency */
   int                   maxnsols;           /**< maximum number of solutions that can be shared in one synchronization */
   int                   nsolvers;           /**< number of solvers synchronizing with this syncstore */

   SCIP_Bool             printincumbents;    /**< should a line be printed for every new globally best solution found
                                              *   by a concurrent solver? */
   SCIP_Real             bestminobj;         /**< minimization-normalized original objective value of the best solution found
                                              *   by any concurrent solver; updated immediately whenever a solver finds a
                                              *   new incumbent, in both parallel modes */

   /* The following members implement the immediate sharing of solutions, dual bounds, and global variable
    * bounds outside of the regular synchronization points. They are only used in opportunistic mode, when
    * the solution pool respectively the bound pool is enabled.
    */
   SCIP_Bool             usesolpool;         /**< should new globally best solutions be shared immediately through
                                              *   the solution pool instead of only at synchronization points? */
   SCIP_Real**           poolsols;           /**< solution value arrays of the pooled solutions; entries are
                                              *   append-only and immutable once published */
   SCIP_Real*            poolobjs;           /**< minimization-normalized original objective values of the pooled solutions */
   int*                  poolsource;         /**< index of the concurrent solver that contributed each pooled solution */
   int*                  poolnvals;          /**< number of values stored for each pooled solution */
   int                   npoolsols;          /**< current number of solutions in the pool; incremented only after the
                                              *   entry is fully written, so it doubles as a lock-free size hint */
   int                   poolsolssize;       /**< allocated capacity of the pool arrays */
   SCIP_Real             bestdualbound;      /**< minimization-normalized best global dual bound, updated immediately by
                                              *   the concurrent solvers when the solution pool is enabled and read by
                                              *   SCIPgetConcurrentDualbound; the primal counterpart is bestminobj */
   int                   winnerid;           /**< solver id of the best terminal status reported so far, or -1; kept
                                              *   directly in the syncstore because when the solution pool is enabled
                                              *   the syncdata slots may be reused before the main SCIP reads them */
   SCIP_STATUS           winnerstatus;       /**< the best terminal status reported so far, or SCIP_STATUS_UNKNOWN */
   SCIP_Bool             useboundpool;       /**< should tightened global variable bounds be shared immediately through
                                              *   the bound pool instead of only at the synchronization points? */
   SCIP_Real*            boundpoollb;        /**< per communication variable, the tightest global lower bound contributed
                                              *   by any concurrent solver, or -infinity; monotonically non-decreasing */
   SCIP_Real*            boundpoolub;        /**< per communication variable, the tightest global upper bound contributed
                                              *   by any concurrent solver, or +infinity; monotonically non-increasing */
   int                   boundpoolsize;      /**< number of communication variables the bound pool is sized for */
   int                   boundpoolversion;   /**< counter bumped whenever a pooled bound is tightened; incremented only
                                              *   after the new bound is written, so it is a safe lock-free change hint */
};


struct SCIP_SyncData
{
   SCIP_Real*            solobj;             /**< array with the objective value of all stored solutions */
   SCIP_Real**           sols;               /**< array with the solution values of each variable for all stored solutions */
   int*                  solsource;          /**< the solverid of the solution came from */
   int                   nsols;              /**< number of solutions currently stored in the synchronization data */
   SCIP_Real             bestlowerbound;     /**< largest lower bound on the objective value that was stored in this
                                              *   synchronization data */
   SCIP_Real             bestupperbound;     /**< smallest upper bound on the objective value that was stored in this
                                              *   synchronization data */
   SCIP_Longint          syncnum;            /**< the synchronization number of this synchronization data */
   int                   winner;             /**< the solverid of the solver with the best status */
   SCIP_STATUS           status;             /**< the best status that was stored in this synchronization data */
   SCIP_LOCK*            lock;               /**< a lock to protect this synchronization data */
   int                   syncedcount;        /**< a counter of how many solvers have finished writing to this synchronization data */
   SCIP_CONDITION*       allsynced;          /**< a condition variable to signal when the last solver has finished writing to this
                                              *   synchronization data */
   SCIP_BOUNDSTORE*      boundstore;         /**< a boundstore for storing all the bound changes that were added to this
                                              *   synchronization data */
   SCIP_Real             syncfreq;           /**< the synchronization frequency that was set in this synchronization data */
   SCIP_Longint          memtotal;           /**< the total amount of memory used by all solvers including the main SCIP */
};

/** struct for storing the position of variables lower and upper bound in the boundstore */
typedef struct
{
   int                   pos[2];             /**< stores at pos[SCIP_BOUNDTYPE_LOWER] the position of the lowerbound and
                                              *   at pos[SCIP_BOUNDTYPE_UPPER] the position of the upperbound */
} BoundPos;

/** struct for storing a single boundchange in the boundstore */
typedef struct
{
   int                   varidx;             /**< the variables position in the variable array of the main SCIP */
   SCIP_Real             newbound;           /**< the variables new bound */
   SCIP_BOUNDTYPE        boundtype;          /**< the type of the variables new bound */
} BoundChg;

struct SCIP_BoundStore
{
   int                   nvars;              /**< the number of variables to store bounds for */
   BoundPos*             bndpos;             /**< array of size nvars to store the positions for all the bound changes
                                              *   stored in this boundstore */
   BoundChg*             bndchg;             /**< array of boundchanges */
   int                   nbndchg;            /**< the number of boundchanges stored in this bound store */
   int                   bndchgsize;         /**< the size of the bound change array */
};

#ifdef __cplusplus
}
#endif

#endif
