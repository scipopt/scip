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
#define DEFAULT_CONFLICTSTORE_MAXSIZE    60000 /* maximal size of a dynamic conflict storage (multiplied by 3) */
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

/* comparison method for constraints */
static
SCIP_DECL_SORTPTRCOMP(compareConss)
{
   /*lint --e{715}*/
   SCIP_CONS* cons1 = (SCIP_CONS*)elem1;
   SCIP_CONS* cons2 = (SCIP_CONS*)elem2;

   assert(cons1 != NULL);
   assert(cons2 != NULL);

   if( SCIPconsGetAge(cons1) > SCIPconsGetAge(cons2) )
      return -1;
   else if ( SCIPconsGetAge(cons1) < SCIPconsGetAge(cons2) )
      return +1;
   else
      return 0;
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
      SCIP_CALL( SCIPsetGetIntParam(set, "conflict/maxstoresize", &conflictstore->maxstoresize) );

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
         conflictstore->maxstoresize = (int)(MIN(3.0 * conflictstore->initstoresize, DEFAULT_CONFLICTSTORE_MAXSIZE));
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
#ifndef NDEBUG
      int i;
#endif
      /* initialize the complete data structure */
      if( conflictstore->conflictsize == 0 )
      {
         newsize = MIN(conflictstore->storesize, DEFAULT_CONFLICTSTORE_SIZE);
         SCIP_ALLOC( BMSallocMemoryArray(&conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSallocMemoryArray(&conflictstore->primalbounds, newsize) );
      }
      else
      {
         newsize = MIN(conflictstore->storesize, SCIPsetCalcMemGrowSize(set, num));
         SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->primalbounds, newsize) );
      }

#ifndef NDEBUG
      for( i = conflictstore->nconflicts; i < newsize; i++ )
      {
         conflictstore->conflicts[i] = NULL;
         conflictstore->primalbounds[i] = -SCIPsetInfinity(set);
      }
#endif
      conflictstore->conflictsize = newsize;
   }
   assert(num <= conflictstore->conflictsize);

   return SCIP_OKAY;
}

/* increase the dynamic storage if we could not delete enough conflicts
 *
 * we want to have at least set->conf_maxconss free slots in the conflict array, because this is the maximal number
 * of conflict generated at a node. we increase the size by the minimum of set->conf_maxconss and 1% of the current
 * store size. nevertheless, we don't exceed conflictstore->maxstoresize.
 */
static
void adjustStorageSize(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   ndelconfs           /**< number of deleted conflicts */
   )
{
   assert(conflictstore != NULL);

   /* increase storage */
   if( conflictstore->storesize - conflictstore->nconflicts <= set->conf_maxconss
      && conflictstore->storesize < conflictstore->maxstoresize )
   {
      conflictstore->storesize += MIN(set->conf_maxconss, (int)(ceil(0.001 * conflictstore->storesize)));
      conflictstore->storesize = MIN(conflictstore->storesize, conflictstore->maxstoresize);
   }

   return;
}

/* removes conflict at position pos */
static
SCIP_RETCODE delPosConflict(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   pos                 /**< position to remove */
)
{
   SCIP_CONS* conflict;
   int lastpos;

   assert(conflictstore != NULL);

   lastpos = conflictstore->nconflicts-1;
   conflict = conflictstore->conflicts[pos];
   assert(conflict != NULL);

   /* decrease number of conflicts depending an a cutoff bound */
   conflictstore->ncbconflicts -= (SCIPsetIsInfinity(set, REALABS(conflictstore->primalbounds[pos])) ? 0 : 1);

   SCIPdebugMessage("-> remove conflict at pos=%d with age=%g\n", pos, SCIPconsGetAge(conflict));

   /* mark the constraint as deleted */
   if( !SCIPconsIsDeleted(conflict) )
   {
      assert(transprob != NULL);
      SCIP_CALL( SCIPconsDelete(conflict, blkmem, set, stat, transprob) );
   }
   SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );

   /* replace with conflict at the last position */
   if( pos < lastpos )
   {
      conflictstore->conflicts[pos] = conflictstore->conflicts[lastpos];
      conflictstore->primalbounds[pos] = conflictstore->primalbounds[lastpos];

#ifndef NDEBUG
      conflictstore->conflicts[lastpos] = NULL;
      conflictstore->primalbounds[lastpos] = -SCIPsetInfinity(set);
#endif
   }

   /* decrease number of conflicts */
   --conflictstore->nconflicts;

   return SCIP_OKAY;
}

/* removes dual ray at position pos */
static
SCIP_RETCODE delPosDualray(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   pos                 /**< position to remove */
)
{
   SCIP_CONS* dualray;
   int lastpos;

   assert(conflictstore != NULL);

   lastpos = conflictstore->ndualrays-1;
   dualray = conflictstore->dualrayconss[pos];
   assert(dualray != NULL);

   SCIPdebugMessage("-> remove dual ray at pos=%d with age=%g\n", pos, SCIPconsGetAge(dualray));

   /* mark the constraint as deleted */
   if( !SCIPconsIsDeleted(dualray) )
   {
      assert(transprob != NULL);
      SCIP_CALL( SCIPconsDelete(dualray, blkmem, set, stat, transprob) );
   }
   SCIP_CALL( SCIPconsRelease(&dualray, blkmem, set) );

   /* replace with dual ray at the last position */
   if( pos < lastpos )
   {
      conflictstore->dualrayconss[pos] = conflictstore->dualrayconss[lastpos];

#ifndef NDEBUG
      conflictstore->dualrayconss[lastpos] = NULL;
#endif
   }

   /* decrease number of dual rays */
   --conflictstore->ndualrays;

   return SCIP_OKAY;
}

/** removes all deleted conflicts from the storage */
static
SCIP_RETCODE cleanDeletedConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int*                  ndelconfs           /**< pointer to store the number of deleted conflicts */
   )
{
   int i;

   assert(conflictstore != NULL);

   (*ndelconfs) = 0;

   for( i = 0; i < conflictstore->nconflicts; )
   {
      assert(conflictstore->conflicts[i] != NULL);

      /* check whether the constraint is already marked as deleted */
      if( SCIPconsIsDeleted(conflictstore->conflicts[i]) )
      {
         /* remove conflict at current position
          *
          * don't increase i because delPosConflict will swap the last pointer to the i-th position
          */
         SCIP_CALL( delPosConflict(conflictstore, set, stat, NULL, blkmem, i) );

         ++(*ndelconfs);
      }
      else
         /* increase i */
         i++;
   }

   SCIPdebugMessage("removed %d/%d as deleted marked conflicts.\n", *ndelconfs, conflictstore->nconflicts);

   return SCIP_OKAY;
}

/** cleans up the storage */
static
SCIP_RETCODE conflictstoreCleanUpStorage(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,          /**< transformed problem */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   SCIP_Real maxage;
   int i;
   int ndelconfs;
   int ndelconfstmp;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);

   /* the storage is empty  */
   if( conflictstore->nconflicts == 0 )
      return SCIP_OKAY;
   assert(conflictstore->nconflicts >= 1);

   /* increase the number of clean up */
   ++conflictstore->ncleanups;

   ndelconfs = 0;
   ndelconfstmp = 0;

   /* remove all as deleted marked conflicts */
   SCIP_CALL( cleanDeletedConflicts(conflictstore, set, stat, blkmem, &ndelconfstmp) );
   ndelconfs += ndelconfstmp;
   ndelconfstmp = 0;

   /* return if at least one conflict could be deleted */
   if( ndelconfs > 0 )
      goto TERMINATE;

   /* only clean up the storage if it is filled enough */
   if( conflictstore->nconflicts < conflictstore->conflictsize )
      goto TERMINATE;

   /* sort conflict */
   SCIPsortPtrReal((void**)conflictstore->conflicts, conflictstore->primalbounds, compareConss, conflictstore->nconflicts);
   assert(SCIPsetIsGE(set, SCIPconsGetAge(conflictstore->conflicts[0]),
         SCIPconsGetAge(conflictstore->conflicts[conflictstore->nconflicts-1])));

   maxage = SCIPconsGetAge(conflictstore->conflicts[0]);

   /* all conflicts have the same age, remove only the first */
   if( SCIPsetIsEQ(set, maxage, SCIPconsGetAge(conflictstore->conflicts[conflictstore->nconflicts-1])) )
   {
      /* remove conflict at first position */
      SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, 0) );
      ++ndelconfstmp;
   }
   else
   {
      assert(ndelconfstmp == 0);

      /* remove all conflicts with max. age */
      for( i = 0; i < conflictstore->nconflicts; i++ )
      {
         if( SCIPsetIsEQ(set, maxage, SCIPconsGetAge(conflictstore->conflicts[i])) )
         {
            /* remove conflict at first position */
            SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, 0) );
            ++ndelconfstmp;

            /* delPos() decreases number of stored conflicts, we need to increase the number of conflicts again and
             * subtract it later, otherwise it could happen that i > nconflicts but there are conflicts with
             * max. age left.
             */
            ++conflictstore->nconflicts;
         }
         else
         {
            /* subtract number of deleted conflicts */
            conflictstore->nconflicts -= ndelconfstmp;
            break;
         }
      }
   }

   /* update statistics */
   ndelconfs += ndelconfstmp;

   /* adjust size of the storage if we use a dynamic store */
   if( set->conf_maxstoresize == -1 )
      adjustStorageSize(conflictstore, set, ndelconfs);
   assert(conflictstore->initstoresize <= conflictstore->storesize);
   assert(conflictstore->storesize <= conflictstore->maxstoresize);

  TERMINATE:
   SCIPdebugMessage("clean-up #%lld: removed %d/%d conflicts, %d depending on cutoff bound\n",
         conflictstore->ncleanups, ndelconfs, conflictstore->nconflicts+ndelconfs, conflictstore->ncbconflicts);

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
   (*conflictstore)->dualrayconss = NULL;
   (*conflictstore)->primalbounds = NULL;
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
   assert(conflictstore != NULL);
   assert(*conflictstore != NULL);

   if( (*conflictstore)->nconflictsfound > 0 && set->conf_cleanbnddepend )
   {
      /* remove solution event from eventfilter */
      SCIP_CALL( SCIPeventfilterDel(eventfilter, blkmem, set, SCIP_EVENTTYPE_BESTSOLFOUND, (*conflictstore)->eventhdlr,
            (SCIP_EVENTDATA*)(*conflictstore), -1) );
   }

   if( (*conflictstore)->conflicts != NULL )
   {
      int i;
      for( i = 0; i < (*conflictstore)->nconflicts; i++ )
      {
         SCIP_CONS* conflict = (*conflictstore)->conflicts[i];
         SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );
      }
      BMSfreeMemoryArray(&(*conflictstore)->conflicts);
      BMSfreeMemoryArray(&(*conflictstore)->primalbounds);
   }

   if( (*conflictstore)->dualrayconss != NULL )
   {
      int i;
      for( i = 0; i < (*conflictstore)->ndualrays; i++ )
      {
         SCIP_CONS* dualray = (*conflictstore)->dualrayconss[i];
         SCIP_CALL( SCIPconsRelease(&dualray, blkmem, set) );
      }
      BMSfreeBlockMemoryArray(blkmem, &(*conflictstore)->dualrayconss, DEFAULT_CONFLICTSTORE_DUALSIZE);
   }

   BMSfreeMemory(conflictstore);

   return SCIP_OKAY;
}

/** adds a constraint to the pool of dual rays
 *
 *  note: this methods captures the constraint
 */
SCIP_RETCODE SCIPconflictstoreAddDualraycons(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_CONS*            dualraycons,        /**< constraint based on a dual ray */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob           /**< transformed problem */
   )
{
   assert(conflictstore != NULL);

   /* mark the constraint to be a conflict */
   SCIPconsMarkConflict(dualraycons);

   /* create an array to stre constraints based on dual rays */
   if( conflictstore->dualrayconss == NULL )
   {
      SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &conflictstore->dualrayconss, DEFAULT_CONFLICTSTORE_DUALSIZE) );
   }

   /* the store is full, we proceed as follows
    *
    * 1. check whether some constraints are marked as deleted and remove those
    * 2. if no constraint is marked as deleted: remove the oldest
    */
   if( conflictstore->ndualrays >= DEFAULT_CONFLICTSTORE_DUALSIZE )
   {
      int i;

      for( i = 0; i < conflictstore->ndualrays; )
      {
         if( SCIPconsIsDeleted(conflictstore->dualrayconss[i]) )
         {
            /* remove dual ray at current position
             *
             * don't increase i because delPosDualray will swap the last pointer to the i-th position
             */
            SCIP_CALL( delPosDualray(conflictstore, set, stat, transprob, blkmem, i) );
         }
         else
            ++i;
      }

      /* sort dual rays */
      SCIPsortPtr((void**)conflictstore->dualrayconss, compareConss, conflictstore->ndualrays);
      assert(SCIPsetIsGE(set, SCIPconsGetAge(conflictstore->dualrayconss[0]),
            SCIPconsGetAge(conflictstore->dualrayconss[conflictstore->ndualrays-1])));

      SCIP_CALL( delPosDualray(conflictstore, set, stat, transprob, blkmem, 0) );
   }

   /* add the new constraint based on a dual ray at the last position */
   SCIPconsCapture(dualraycons);
   conflictstore->dualrayconss[conflictstore->ndualrays] = dualraycons;
   ++conflictstore->ndualrays;

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
   if( conflictstore->lastnodenum != SCIPnodeGetNumber(SCIPtreeGetFocusNode(tree))
         || conflictstore->nconflicts == conflictstore->conflictsize )
   {
      SCIP_CALL( conflictstoreCleanUpStorage(conflictstore, set, stat, transprob, blkmem) );
   }

   /* update the last seen node */
   conflictstore->lastnodenum = SCIPnodeGetNumber(SCIPtreeGetFocusNode(tree));

   SCIPconsCapture(cons);
   conflictstore->conflicts[conflictstore->nconflicts] = cons;
   conflictstore->primalbounds[conflictstore->nconflicts] = primalbound;
   conflictstore->ncbconflicts += (SCIPsetIsInfinity(set, REALABS(primalbound)) ? 0 : 1);

   ++conflictstore->nconflicts;
   ++conflictstore->nconflictsfound;

   SCIPdebugMessage("add conflict <%s> to conflict store at position %d\n", SCIPconsGetName(cons), conflictstore->nconflicts-1);
   SCIPdebugMessage(" -> conflict type: %d, cutoff involved = %u\n", conftype, cutoffinvolved);
   if( cutoffinvolved )
      SCIPdebugMessage(" -> current primal bound: %g\n", primalbound);
   SCIPdebugMessage(" -> found at node %llu (depth: %d), valid at node %llu (depth: %d)\n", SCIPnodeGetNumber(node),
         SCIPnodeGetDepth(node), SCIPnodeGetNumber(validnode), SCIPnodeGetDepth(validnode));

   return SCIP_OKAY;
}

/** deletes all conflicts depending a cutoff bound larger than the given bound */
SCIP_RETCODE SCIPconflictstoreCleanNewIncumbant(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PROB*            transprob,          /**< transformed problem*/
   SCIP_Real             cutoffbound         /**< current cutoff bound */
   )
{
   SCIP_Real improvement;
   int ndelconfs;
   int ndelconfs_del;
   int i;

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

   /* remove all conflicts depending on a primalbound*improvement > cutoffbound */
   for( i = 0; i < conflictstore->nconflicts; )
   {
      SCIP_CONS* conflict;

      /* get the conflict */
      conflict = conflictstore->conflicts[i];
      assert(conflict != NULL);

      /* check if the conflict depends on the cutoff bound */
      if( SCIPsetIsGT(set, improvement * conflictstore->primalbounds[i], cutoffbound) )
      {
         /* remove conflict at current position
          *
          * don't increase i because delPosConflict will swap the last pointer to the i-th position
          */
         SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, i) );
         ++ndelconfs;
      }
      /* we remove all as deleted marked constraints too */
      else if( SCIPconsIsDeleted(conflict) )
      {
         /* remove conflict at current position
          *
          * don't increase i because delPosConflict will swap the last pointer to the i-th position
          */
         SCIP_CALL( delPosConflict(conflictstore, set, stat, transprob, blkmem, i) );
         ++ndelconfs_del;
      }
      else
         /* increase i */
         ++i;
   }
   assert(conflictstore->ncbconflicts >= 0);
   assert(conflictstore->nconflicts >= 0);

   SCIPdebugMessage("-> removed %d/%d conflicts, %d depending on cutoff bound\n", ndelconfs+ndelconfs_del,
         conflictstore->nconflicts+ndelconfs+ndelconfs_del, ndelconfs);

   return SCIP_OKAY;
}

/** returns the maximal size of the conflict pool */
int SCIPconflictstoreGetMaxPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   )
{
   assert(conflictstore != NULL);

   return MIN(conflictstore->storesize, conflictstore->maxstoresize);
}

/** return the initial size of the conflict pool */
int SCIPconflictstoreGetInitPoolSize(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->initstoresize;
}

/** returns the number of stored conflicts on the conflict pool
 *
 *  note: the number of active conflicts can be less
 */
int SCIPconflictstoreGetNConflictsInStore(
   SCIP_CONFLICTSTORE*   conflictstore       /**< conflict storage */
   )
{
   assert(conflictstore != NULL);

   return conflictstore->nconflicts;
}

/** returns all active conflicts stored in the conflict store */
SCIP_RETCODE SCIPconflictstoreGetConflicts(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_CONS**           conflicts,          /**< array to store conflicts */
   int                   conflictsize,       /**< site of the conflict array */
   int*                  nconflicts          /**< pointer to store the number of conflicts */
   )
{
   int i;

   assert(conflictstore != NULL);

   /* return if the allocated memory is obviously to small */
   if( conflictstore->nconflicts > conflictsize )
   {
      (*nconflicts) = conflictstore->nconflicts;
      return SCIP_OKAY;
   }

   (*nconflicts) = 0;
   for( i = 0; i < conflictstore->nconflicts; i++ )
   {
      SCIP_CONS* conflict;

      conflict = conflictstore->conflicts[i];
      assert(conflict != NULL);

      /* skip deactivated and deleted constraints */
      if( !SCIPconsIsActive(conflict) || SCIPconsIsDeleted(conflict) )
         continue;

      /* count exact number conflicts */
      if( *nconflicts > conflictsize )
         ++(*nconflicts);
      else
      {
         conflicts[*nconflicts] = conflict;
         ++(*nconflicts);
      }
   }

   return SCIP_OKAY;
}
