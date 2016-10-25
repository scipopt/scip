/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   conflictstore.c
 * @brief  methods for storing conflicts
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/event.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/prob.h"
#include "scip/scip.h"
#include "scip/def.h"
#include "scip/cons_linear.h"
#include "scip/struct_conflictstore.h"


#define DEFAULT_CONFLICTSTORE_DUALSIZE     100 /* default size of conflict storage */
#define DEFAULT_CONFLICTSTORE_MINSIZE     2000 /* default minimal size of a dynamic conflict storage */
#define DEFAULT_CONFLICTSTORE_MAXSIZE    30000 /* maximal size of a dynamic conflict storage (multiplied by 3) */
#define DEFAULT_CONFLICTSTORE_SIZE       10000 /* default size of conflict storage */

/* event handler properties */
#define EVENTHDLR_NAME         "ConflictStore"
#define EVENTHDLR_DESC         "Solution event handler for conflict store."


/* exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecConflictstore)
{
   assert(eventhdlr != NULL);
   assert(eventdata != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_BESTSOLFOUND);

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_CALL( SCIPcleanConflictStoreNewIncumbant(scip, event) );
   }

   return SCIP_OKAY;
}

/* initialize the conflict storage */
static
SCIP_RETCODE initConflictstore(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(conflictstore != NULL);

   /* calculate the maximal size of the conflict storage */
   if( conflictstore->maxstoresize == -1 )
   {
      SCIP_CALL( SCIPsetGetIntParam(set, "conflict/graph/maxstoresize", &conflictstore->maxstoresize) );

      /* the size should be dynamic wrt number of variables after presolving */
      if( conflictstore->maxstoresize == -1 )
      {
         int nconss;
         int nvars;

         nconss = SCIPprobGetNConss(transprob);
         nvars = SCIPprobGetNVars(transprob);

         conflictstore->initstoresize = DEFAULT_CONFLICTSTORE_MINSIZE;
         conflictstore->initstoresize += 2*nconss;

         if( nvars/2 <= 500 )
            conflictstore->initstoresize += (int) DEFAULT_CONFLICTSTORE_MAXSIZE/100;
         else if( nvars/2 <= 5000 )
            conflictstore->initstoresize += (int) DEFAULT_CONFLICTSTORE_MAXSIZE/10;
         else
            conflictstore->initstoresize += DEFAULT_CONFLICTSTORE_MAXSIZE/2;

         conflictstore->initstoresize = MIN(conflictstore->initstoresize, DEFAULT_CONFLICTSTORE_MAXSIZE);
         conflictstore->storesize = conflictstore->initstoresize;
         conflictstore->maxstoresize = (int)(3.0 * MIN(conflictstore->initstoresize, DEFAULT_CONFLICTSTORE_MAXSIZE));
      }
      else
      {
         conflictstore->initstoresize = conflictstore->maxstoresize;
         conflictstore->storesize = conflictstore->maxstoresize;
      }

#ifdef NDEBUG
      if( conflictstore->maxstoresize == 0 )
         SCIPdebugMessage("usage of conflict pool is disabled.\n");
      else
         SCIPdebugMessage("[init,max] size of conflict pool is [%d,%d].\n",
               conflictstore->initstoresize, conflictstore->maxstoresize);
#endif

      if( set->conf_cleanbnddepend )
      {
         /* add solution event to the eventfilter */
         SCIP_CALL( SCIPeventfilterAdd(eventfilter, blkmem, set, SCIP_EVENTTYPE_BESTSOLFOUND, conflictstore->eventhdlr,
               (SCIP_EVENTDATA*)conflictstore, NULL) );
      }
   }

   return SCIP_OKAY;
}

/** resizes conflict and primal bound arrays to be able to store at least num entries */
static
SCIP_RETCODE conflictstoreEnsureMem(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflictstore != NULL);
   assert(set != NULL);

   /* we do not allocate more memory as allowed */
   if( conflictstore->conflictsize == conflictstore->maxstoresize )
      return SCIP_OKAY;

   if( num > conflictstore->conflictsize )
   {
      int newsize;
      int i;

      /* initialize the complete data structure */
      if( conflictstore->conflictsize == 0 )
      {
         newsize = MIN(conflictstore->storesize, DEFAULT_CONFLICTSTORE_SIZE);
         SCIP_CALL( SCIPqueueCreate(&conflictstore->slotqueue, newsize, 2) );
         SCIP_CALL( SCIPqueueCreate(&conflictstore->orderqueue, newsize, 2) );
         SCIP_ALLOC( BMSallocMemoryArray(&conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSallocMemoryArray(&conflictstore->primalbounds, newsize) );
      }
      else
      {
         newsize = MIN(conflictstore->storesize, SCIPsetCalcMemGrowSize(set, num));
         SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->primalbounds, newsize) );
      }

      /* add all new slots (oldsize,...,newsize-1) with a shift of +1 to the slotqueue */
      for( i = conflictstore->conflictsize; i < newsize; i++ )
      {
         conflictstore->conflicts[i] = NULL;
         conflictstore->primalbounds[i] = -SCIPsetInfinity(set);
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (i+1)) );
      }
      conflictstore->conflictsize = newsize;
   }
   assert(num <= conflictstore->conflictsize);

   assert(SCIPqueueNElems(conflictstore->slotqueue)+conflictstore->nconflicts == conflictstore->conflictsize );
   assert(SCIPqueueNElems(conflictstore->slotqueue)+SCIPqueueNElems(conflictstore->orderqueue) == conflictstore->conflictsize);

   return SCIP_OKAY;
}

/* increase the dynamic storage if we could not delete enough conflicts */
static
void adjustStorageSize(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   ndelconfs           /**< number of deleted conflicts */
   )
{
   assert(conflictstore != NULL);

   /* increase storage */
   if( conflictstore->storesize - (conflictstore->nconflicts - ndelconfs) <= set->conf_maxconss
      && conflictstore->storesize < conflictstore->maxstoresize )
   {
      conflictstore->storesize += MIN(set->conf_maxconss, (int)(ceil(0.001 * conflictstore->storesize)));
      conflictstore->storesize = MIN(conflictstore->storesize, conflictstore->maxstoresize);
   }

   return;
}

/** remove all deleted conflicts from the storage */
static
SCIP_RETCODE cleanDeletedConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int*                  ndelconfs           /**< pointer to store the number of deleted conflicts */
   )
{
   int firstidx;

   assert(conflictstore != NULL);
   assert(SCIPqueueNElems(conflictstore->slotqueue)+conflictstore->nconflicts == conflictstore->conflictsize );
   assert(SCIPqueueNElems(conflictstore->slotqueue)+SCIPqueueNElems(conflictstore->orderqueue) == conflictstore->conflictsize);

   (*ndelconfs) = 0;
   firstidx = -1;

   while( firstidx != ((int) (size_t) SCIPqueueFirst(conflictstore->orderqueue)-1) )
   {
      SCIP_CONS* conflict;
      int idx;

      assert(!SCIPqueueIsEmpty(conflictstore->orderqueue));
      idx = ((int) (size_t) SCIPqueueRemove(conflictstore->orderqueue)) - 1;
      assert(idx >= 0 && idx < conflictstore->conflictsize);

      if( conflictstore->conflicts[idx] == NULL )
         continue;

      /* get the oldest conflict */
      conflict = conflictstore->conflicts[idx];

      /* check whether the constraint is already marked as deleted */
      if( SCIPconsIsDeleted(conflict) )
      {
         /* release the constraint */
         SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );

         conflictstore->ncbconflicts -= (SCIPsetIsInfinity(set, REALABS(conflictstore->primalbounds[idx])) ? 0 : 1);

         /* clean the conflict and primal bound array */
         conflictstore->conflicts[idx] = NULL;
         conflictstore->primalbounds[idx] = -SCIPsetInfinity(set);

         /* add the id shifted by +1 to the queue of empty slots */
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );

         ++(*ndelconfs);
      }
      else
      {
         /* remember the first conflict that is not deleted */
         if( firstidx == -1 )
            firstidx = idx;

         SCIP_CALL( SCIPqueueInsert(conflictstore->orderqueue, (void*) (size_t) (idx+1)) );
      }
   }

   SCIPdebugMessage("removed %d/%d as deleted marked conflicts.\n", *ndelconfs, conflictstore->nconflicts);

   assert(SCIPqueueNElems(conflictstore->slotqueue)+conflictstore->nconflicts-(*ndelconfs) == conflictstore->conflictsize );
   assert(SCIPqueueNElems(conflictstore->slotqueue)+SCIPqueueNElems(conflictstore->orderqueue) == conflictstore->conflictsize);

   return SCIP_OKAY;
}

/** clean up the storage */
static
SCIP_RETCODE conflictstoreCleanUpStorage(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_CONS* conflict;
   SCIP_Real maxage;
   int idx;
   int ndelconfs;
   int ndelconfstmp;
   int nseenconfs;
   int tmpidx;
   int nimpr;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);

   assert(SCIPqueueNElems(conflictstore->slotqueue)+conflictstore->nconflicts == conflictstore->conflictsize );
   assert(SCIPqueueNElems(conflictstore->slotqueue)+SCIPqueueNElems(conflictstore->orderqueue) == conflictstore->conflictsize);

   /* the storage is empty  */
   if( conflictstore->nconflicts == 0 )
   {
      assert(SCIPqueueNElems(conflictstore->slotqueue) == conflictstore->conflictsize);
      return SCIP_OKAY;
   }
   assert(conflictstore->nconflicts >= 1);

   /* increase the number of clean up */
   ++conflictstore->ncleanups;

   ndelconfs = 0;
   ndelconfstmp = 0;
   nseenconfs = 0;

   /* remove all as deleted marked conflicts */
   SCIP_CALL( cleanDeletedConflicts(conflictstore, set, blkmem, &ndelconfstmp) );
   ndelconfs += ndelconfstmp;
   ndelconfstmp = 0;

   assert(SCIPqueueNElems(conflictstore->slotqueue)+conflictstore->nconflicts-ndelconfs == conflictstore->conflictsize );
   assert(SCIPqueueNElems(conflictstore->slotqueue)+SCIPqueueNElems(conflictstore->orderqueue) == conflictstore->conflictsize);

   /* return if at least one conflict could be deleted */
   if( ndelconfs > 0 )
      goto TERMINATE;

   /* only clean up the storage if it is filled enough */
   if( conflictstore->nconflicts < conflictstore->conflictsize )
      goto TERMINATE;

   assert(!SCIPqueueIsEmpty(conflictstore->orderqueue));

   nimpr = 0;
   tmpidx = -1;
   maxage = -SCIPsetInfinity(set);

   /* find a conflict with a locally maximum age */
   while( nseenconfs < conflictstore->nconflicts-ndelconfs )
   {
      /* get the oldest conflict */
      assert(!SCIPqueueIsEmpty(conflictstore->orderqueue));
      idx = ((int) (size_t) SCIPqueueRemove(conflictstore->orderqueue)) - 1;
      assert(idx >= 0 && idx < conflictstore->conflictsize);

      if( conflictstore->conflicts[idx] == NULL )
      {
         /* add the id shifted by +1 to the queue of empty slots */
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );
         continue;
      }

      /* get the oldest conflict */
      conflict = conflictstore->conflicts[idx];
      assert(!SCIPconsIsDeleted(conflict));

      ++nseenconfs;

      /* check if the age of the conflict is positive and larger than maxage; do nothing we have seen enough improvements
       *
       * note: this can be expensive if the store is to large. but experiments have been shown, that removing all
       *       conflicts older than a local maximum can be worse.
       */
      if( SCIPsetIsGT(set, SCIPconsGetAge(conflict), 0.0) && SCIPsetIsLT(set, maxage, SCIPconsGetAge(conflict)) )
      {
         maxage = SCIPconsGetAge(conflict);
         tmpidx = idx;
         ++nimpr;
      }

      /* reinsert the id */
      SCIP_CALL( SCIPqueueInsert(conflictstore->orderqueue, (void*) (size_t) (idx+1)) );
   }

   /* no conflict was chosen because all conflicts have age 0 */
   assert(tmpidx >= 0 || SCIPsetIsInfinity(set, -maxage));
   assert(!SCIPqueueIsEmpty(conflictstore->orderqueue));
   maxage = tmpidx == -1 ? 0 : maxage;

   /* iterate over all conflicts and remove those with an age larger or equal the local maximum maxage */
   nseenconfs = 0;
   ndelconfstmp = 0;
   while( nseenconfs < conflictstore->nconflicts-ndelconfs )
   {
      /* get the oldest conflict */
      assert(!SCIPqueueIsEmpty(conflictstore->orderqueue));
      idx = ((int) (size_t) SCIPqueueRemove(conflictstore->orderqueue)) - 1;
      assert(idx >= 0 && idx < conflictstore->conflictsize);

      if( conflictstore->conflicts[idx] == NULL )
      {
         /* add the id shifted by +1 to the queue of empty slots */
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );
         continue;
      }

      conflict = conflictstore->conflicts[idx];
      ++nseenconfs;
      assert(conflict != NULL);
      assert(!SCIPconsIsDeleted(conflict));

      if( SCIPsetIsLT(set, SCIPconsGetAge(conflict), maxage) )
      {
         SCIP_CALL( SCIPqueueInsert(conflictstore->orderqueue, (void*) (size_t) (idx+1)) );
         continue;
      }

      /* mark the constraint as deleted */
      SCIP_CALL( SCIPconsDelete(conflict, blkmem, set, stat, transprob) );
      SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );

      /* clean the conflict and primal bound array */
      conflictstore->conflicts[idx] = NULL;
      conflictstore->primalbounds[idx] = -SCIPsetInfinity(set);

      /* add the id shifted by +1 to the queue of empty slots */
      SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );

      ++ndelconfstmp;
      SCIPdebugMessage("-> removed conflict at pos=%d with age=%g\n", idx, maxage);

      /* all conflicts have age 0, we delete the oldest conflicts */
      if( SCIPsetIsEQ(set, maxage, 0.0) )
      {
         assert(tmpidx == -1);
         break;
      }
   }

   assert(SCIPqueueNElems(conflictstore->orderqueue) <= conflictstore->maxstoresize);
   ndelconfs += ndelconfstmp;

   /* adjust size of the storage if we use a dynamic store */
   if( set->conf_maxstoresize == -1 )
      adjustStorageSize(conflictstore, set, ndelconfs);
   assert(conflictstore->initstoresize <= conflictstore->storesize);
   assert(conflictstore->storesize <= conflictstore->maxstoresize);

  TERMINATE:
   SCIPdebugMessage("clean-up #%lld: removed %d/%d conflicts, %d depending on cutoff bound\n",
         conflictstore->ncleanups, ndelconfs, conflictstore->nconflicts, conflictstore->ncbconflicts);
   conflictstore->nconflicts -= ndelconfs;

   assert(SCIPqueueNElems(conflictstore->slotqueue)+conflictstore->nconflicts == conflictstore->conflictsize );
   assert(SCIPqueueNElems(conflictstore->slotqueue)+SCIPqueueNElems(conflictstore->orderqueue) == conflictstore->conflictsize);

   return SCIP_OKAY;
}


/** creates conflict storage */
SCIP_RETCODE SCIPconflictstoreCreate(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict storage */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflictstore != NULL);

   SCIP_ALLOC( BMSallocMemory(conflictstore) );

   (*conflictstore)->conflicts = NULL;
   (*conflictstore)->dualrays = NULL;
   (*conflictstore)->primalbounds = NULL;
   (*conflictstore)->slotqueue = NULL;
   (*conflictstore)->orderqueue = NULL;
   (*conflictstore)->conflictsize = 0;
   (*conflictstore)->nconflicts = 0;
   (*conflictstore)->ndualrays = 0;
   (*conflictstore)->ncbconflicts = 0;
   (*conflictstore)->nconflictsfound = 0;
   (*conflictstore)->initstoresize = -1;
   (*conflictstore)->storesize = -1;
   (*conflictstore)->maxstoresize = -1;
   (*conflictstore)->ncleanups = 0;
   (*conflictstore)->lastnodenum = -1;

   /* create event handler for LP events */
   SCIP_CALL( SCIPeventhdlrCreate(&(*conflictstore)->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, eventExecConflictstore, NULL) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(set, (*conflictstore)->eventhdlr) );

   if( (*conflictstore)->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for conflictstore not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* initialize the conflict handler */
   SCIP_CALL( SCIPeventhdlrInit((*conflictstore)->eventhdlr, set) );

   return SCIP_OKAY;
}

/** frees conflict storage */
SCIP_RETCODE SCIPconflictstoreFree(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_CONS* conflict;

   assert(conflictstore != NULL);
   assert(*conflictstore != NULL);

   if( (*conflictstore)->nconflictsfound > 0 && set->conf_cleanbnddepend )
   {
      /* remove solution event from eventfilter */
      SCIP_CALL( SCIPeventfilterDel(eventfilter, blkmem, set, SCIP_EVENTTYPE_BESTSOLFOUND, (*conflictstore)->eventhdlr,
            (SCIP_EVENTDATA*)(*conflictstore), -1) );
   }

   if( (*conflictstore)->orderqueue != NULL )
   {
      assert((*conflictstore)->slotqueue != NULL);

      while( !SCIPqueueIsEmpty((*conflictstore)->orderqueue) )
      {
         int idx;

         idx = ((int) (size_t) SCIPqueueRemove((*conflictstore)->orderqueue)) - 1;
         assert(idx >= 0 && idx < (*conflictstore)->conflictsize);

         if( (*conflictstore)->conflicts[idx] == NULL )
            continue;

         /* get the conflict */
         conflict = (*conflictstore)->conflicts[idx];

         SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );
         --(*conflictstore)->nconflicts;
      }

      /* free the queues */
      SCIPqueueFree(&(*conflictstore)->slotqueue);
      SCIPqueueFree(&(*conflictstore)->orderqueue);
   }
   assert((*conflictstore)->nconflicts == 0);

   if( (*conflictstore)->dualrays != NULL )
   {
      int i;
      for( i = 0; i < (*conflictstore)->ndualrays; i++ )
      {
         SCIP_CONS* dualray = (*conflictstore)->dualrays[i];
         SCIP_CALL( SCIPconsRelease(&dualray, blkmem, set) );
      }
      BMSfreeBlockMemoryArray(blkmem, &(*conflictstore)->dualrays, DEFAULT_CONFLICTSTORE_DUALSIZE);
   }

   BMSfreeMemoryArrayNull(&(*conflictstore)->conflicts);
   BMSfreeMemoryArrayNull(&(*conflictstore)->primalbounds);
   BMSfreeMemory(conflictstore);

   return SCIP_OKAY;
}

/** adds a constraint to the pool of dual rays
 *
 *  note: this methods captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddDualray(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_CONS*            dualray,            /**< dual ray to add */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob           /**< transformed problem */
   )
{
   SCIP_CONS* olddualray;

   assert(conflictstore != NULL);

   /* mark the constraint to be a conflict */
   SCIPconsMarkConflict(dualray);

   /* create an array to stre constraints based on dual rays */
   if( conflictstore->dualrays == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->dualrays, DEFAULT_CONFLICTSTORE_DUALSIZE) );
   }

   /* the store is not full */
   if( conflictstore->ndualrays < DEFAULT_CONFLICTSTORE_DUALSIZE )
   {
      conflictstore->dualrays[conflictstore->ndualrays] = dualray;
      SCIPconsCapture(dualray);

      ++conflictstore->ndualrays;
   }

   /* the store is full, we proceed as follows
    *
    * 1. check whether some constraints are marked as deleted and remove those
    * 2. if no constraint is marked as deleted: remove the oldest
    */
   else
   {
      SCIP_Real maxage;
      int idx;
      int i;

      maxage = 0.0;
      idx = -1;

      for( i = 0; i < conflictstore->ndualrays; )
      {
         SCIP_CONS* olddualray = conflictstore->dualrays[i];

         assert(olddualray != NULL);

         /* remove dual rays that are already marked as deleted */
         if( SCIPconsIsDeleted(olddualray) )
         {
            /* remove the deleted dual ray */
            SCIP_CALL( SCIPconsRelease(&olddualray, blkmem, set) );
            --conflictstore->ndualrays;

            SCIPsetDebugMessagePrint(set, "-> remove deleted dual ray: pos=%d\n", i);

            /* replace with the last stored dual ray */
            if( i < conflictstore->ndualrays )
            {
               conflictstore->dualrays[i] = conflictstore->dualrays[conflictstore->ndualrays];
               SCIPsetDebugMessagePrint(set, "->           replaced with: pos=%d\n", conflictstore->ndualrays);
            }

            continue;
         }

         /* remember the position of the oldest dual ray within the array */
         if( SCIPsetIsGT(set, SCIPconsGetAge(olddualray), maxage) )
         {
            maxage = SCIPconsGetAge(olddualray);
            idx = i;
         }

         ++i;
      }
      assert(idx >= 0);

      SCIPsetDebugMessagePrint(set, "-> remove oldest dual ray: pos=%d, age=%g\n", idx, maxage);

      /* remove the oldest dual ray */
      SCIP_CALL( SCIPconsDelete(conflictstore->dualrays[idx], blkmem, set, stat, transprob) );
      SCIP_CALL( SCIPconsRelease(&conflictstore->dualrays[idx], blkmem, set) );

      /* add the new dual ray at the position of the oldest */
      conflictstore->dualrays[idx] = dualray;
      SCIPconsCapture(dualray);
   }

   return SCIP_OKAY;
}

/** adds a conflict to the conflict storage
 *
 *  note: this method captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_PROB*            transprob,          /**< transformed problem */
   SCIP_EVENTFILTER*     eventfilter,        /**< eventfiler */
   SCIP_CONS*            cons,               /**< constraint representing the conflict */
   SCIP_NODE*            node,               /**< node to add conflict (or NULL if global) */
   SCIP_NODE*            validnode,          /**< node at whichaddConf the constraint is valid (or NULL) */
   SCIP_CONFTYPE         conftype,           /**< type of the conflict */
   SCIP_Bool             cutoffinvolved,     /**< is a cutoff bound invaled in this conflict */
   SCIP_Real             primalbound         /**< primal bound the conflict depend on (or -SCIPinfinity) */
   )
{
   int nconflicts;
   int idx;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(transprob != NULL);
   assert(cons != NULL);
   assert(node != NULL);
   assert(validnode != NULL);
   assert(set->conf_allowlocal || SCIPnodeGetDepth(validnode) == 0);
   assert(conftype != SCIP_CONFTYPE_UNKNOWN);
   assert(conftype != SCIP_CONFTYPE_BNDEXCEEDING || cutoffinvolved);
   assert(!cutoffinvolved || (cutoffinvolved && !SCIPsetIsInfinity(set, REALABS(primalbound))));

   nconflicts = conflictstore->nconflicts;

   /* mark the constraint to be a conflict */
   SCIPconsMarkConflict(cons);

   /* initialize the storage */
   if( conflictstore->maxstoresize == -1 )
   {
      SCIP_CALL( initConflictstore(conflictstore, set, transprob, eventfilter, blkmem) );
   }
   assert(conflictstore->initstoresize >= 0);
   assert(conflictstore->initstoresize <= conflictstore->maxstoresize);

   /* return if conflict pool is disabled */
   if( conflictstore->maxstoresize == 0 )
      return SCIP_OKAY;

   SCIP_CALL( conflictstoreEnsureMem(conflictstore, set, nconflicts+1) );

   /* return if the store has size zero */
   if( conflictstore->conflictsize == 0 )
   {
      assert(conflictstore->maxstoresize == 0);
      return SCIP_OKAY;
   }

   /* clean up the storage if we are at a new node or the storage is full */
   if( conflictstore->lastnodenum != SCIPnodeGetNumber(SCIPtreeGetFocusNode(tree)) || SCIPqueueIsEmpty(conflictstore->slotqueue) )
   {
      SCIP_CALL( conflictstoreCleanUpStorage(conflictstore, set, stat, transprob, blkmem) );
   }

   /* update the last seen node */
   conflictstore->lastnodenum = SCIPnodeGetNumber(SCIPtreeGetFocusNode(tree));

   /* get a free slot */
   assert(!SCIPqueueIsEmpty(conflictstore->slotqueue));
   idx = ((int) (size_t) SCIPqueueRemove(conflictstore->slotqueue)-1);
   assert(idx >= 0 && idx < conflictstore->conflictsize);
   assert(conflictstore->conflicts[idx] == NULL);
   assert(conflictstore->primalbounds[idx] == -SCIPsetInfinity(set));

   SCIPconsCapture(cons);
   conflictstore->conflicts[idx] = cons;
   conflictstore->primalbounds[idx] = primalbound;
   conflictstore->ncbconflicts += (SCIPsetIsInfinity(set, REALABS(primalbound)) ? 0 : 1);

   /* add idx shifted by +1 to the ordering queue */
   SCIP_CALL( SCIPqueueInsert(conflictstore->orderqueue, (void*) (size_t) (idx+1)) );

   ++conflictstore->nconflicts;
   ++conflictstore->nconflictsfound;

   SCIPdebugMessage("add conflict <%s> to conflict store at position %d\n", SCIPconsGetName(cons), idx);
   SCIPdebugMessage(" -> conflict type: %d, cutoff involved = %u\n", conftype, cutoffinvolved);
   if( cutoffinvolved )
      SCIPdebugMessage(" -> current primal bound: %g\n", primalbound);
   SCIPdebugMessage(" -> found at node %llu (depth: %d), valid at node %llu (depth: %d)\n", SCIPnodeGetNumber(node),
         SCIPnodeGetDepth(node), SCIPnodeGetNumber(validnode), SCIPnodeGetDepth(validnode));

   return SCIP_OKAY;
}

/** delete all conflicts depending a cutoff bound larger than the given bound */
SCIP_RETCODE SCIPconflictstoreCleanNewIncumbant(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem*/
   SCIP_Real             cutoffbound         /**< current cutoff bound */
   )
{
   SCIP_CONS* conflict;
   SCIP_Real improvement;
   int nseenconfs;
   int ndelconfs;
   int ndelconfs_del;
   int idx;

   assert(conflictstore != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(blkmem != NULL);
   assert(transprob != NULL);

   ndelconfs = 0;
   ndelconfs_del = 0;

   /* return if we do not want to use the storage */
   if( set->conf_maxstoresize == 0 )
      return SCIP_OKAY;

   /* return if we do not want to remove conflicts related to an older cutoff bound */
   if( !set->conf_cleanbnddepend )
      return SCIP_OKAY;

   /* TODO we may want to introduce a paramter */
   if( SCIPsetIsPositive(set, cutoffbound) )
      improvement = (1 - 0.05);
   else
      improvement = (1 + 0.05);

   /* remove all conflicts depending on the cutoffbound */
   nseenconfs = 0;
   while( nseenconfs < conflictstore->nconflicts )
   {
      /* get an index from the queue */
      assert(!SCIPqueueIsEmpty(conflictstore->orderqueue));
      idx = ((int) (size_t) SCIPqueueRemove(conflictstore->orderqueue)) - 1;
      assert(idx >= 0 && idx < conflictstore->conflictsize);

      if( conflictstore->conflicts[idx] == NULL )
      {
         /* add the id shifted by +1 to the queue of empty slots */
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );
         continue;
      }

      /* get the conflict */
      conflict = conflictstore->conflicts[idx];

      ++nseenconfs;

      /* check if the conflict depends on the cutofbound */
      if( SCIPsetIsGT(set, improvement * conflictstore->primalbounds[idx], cutoffbound) )
      {
         SCIP_CALL( SCIPconsDelete(conflict, blkmem, set, stat, transprob) );
         SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );

         conflictstore->conflicts[idx] = NULL;
         conflictstore->primalbounds[idx] = -SCIPsetInfinity(set);

         ++ndelconfs;
         --conflictstore->ncbconflicts;

         /* add the index to queue of free slots */
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );
      }
      /* we remove all as deleted marked constraints too */
      else if( SCIPconsIsDeleted(conflict) )
      {
         /* release the constraint */
         SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );

         conflictstore->ncbconflicts -= (SCIPsetIsInfinity(set, REALABS(conflictstore->primalbounds[idx])) ? 0 : 1);

         /* clean the conflict and primal bound array */
         conflictstore->conflicts[idx] = NULL;
         conflictstore->primalbounds[idx] = -SCIPsetInfinity(set);

         /* add the id shifted by +1 to the queue of empty slots */
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );

         ++ndelconfs_del;
      }
      else
      {
         /* reinsert the id */
         SCIP_CALL( SCIPqueueInsert(conflictstore->orderqueue, (void*) (size_t) (idx+1)) );
      }
   }
   assert(conflictstore->ncbconflicts >= 0);

   SCIPdebugMessage("-> removed %d/%d conflicts, %d depending on cutoff bound\n", ndelconfs+ndelconfs_del, conflictstore->nconflicts, ndelconfs);
   conflictstore->nconflicts -= (ndelconfs+ndelconfs_del);

   return SCIP_OKAY;
}

/* return the maximal size of the conflict pool */
int SCIPconflictstoreGetMaxPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   )
{
   assert(conflictstore != NULL);

   return MAX(conflictstore->storesize, conflictstore->maxstoresize);
}

/* return the initial size of the conflict pool */
int SCIPconflictstoreGetInitPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->initstoresize;
}
