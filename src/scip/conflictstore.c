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
#define SCIP_DEBUG
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"

#include "scip/struct_conflictstore.h"


#define DEFAULT_CONFLICTSTORE_SIZE       10000 /* maximal size of conflict storage */

/*
 * dynamic memory arrays
 */

/** resizes cuts and score arrays to be able to store at least num entries */
static
SCIP_RETCODE conflictstoreEnsureCutsMem(
   SCIP_CONFLICTSTORE*   conflictstore,      /**< conflict storage */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflictstore != NULL);
   assert(set != NULL);

   /* we do not allocate more memory as allowed */
   if( conflictstore->conflictsize == set->conf_maxstoresize )
      return SCIP_OKAY;

   if( num > conflictstore->conflictsize )
   {
      int newsize;
      int i;

      /* initialize the complete data structure */
      if( conflictstore->conflictsize == 0 )
      {
         newsize = MIN(set->conf_maxstoresize, DEFAULT_CONFLICTSTORE_SIZE);
         SCIP_CALL( SCIPqueueCreate(&conflictstore->slotqueue, newsize, 2) );
         SCIP_CALL( SCIPqueueCreate(&conflictstore->orderqueue, newsize, 2) );
         SCIP_ALLOC( BMSallocMemoryArray(&conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSallocMemoryArray(&conflictstore->primalbounds, newsize) );
      }
      else
      {
         newsize = MIN(set->conf_maxstoresize, conflictstore->conflictsize * 2);
         SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->conflicts, newsize) );
         SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->primalbounds, newsize) );
      }

      /* add all new slots (oldsize,...,newsize-1) withdeclaration a shift of +1 to the slotqueue */
      for( i = conflictstore->conflictsize; i < newsize; i++ )
      {
#ifdef SCIP_DEBUG
         conflictstore->conflicts[i] = NULL;
         conflictstore->primalbounds[i] = SCIP_INVALID;
#endif
         SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (i+1)) );
      }
      conflictstore->conflictsize = newsize;
   }
   assert(num <= conflictstore->conflictsize);

   return SCIP_OKAY;
}


#if 0
/** remove all deleted conflicts from the storage */
static
SCIP_RETCODE cleanDeletedConflicts(
   SCIP_CONFLICTSTORE*   conflictstore
   )
{
   assert(conflictstore != NULL);

   return SCIP_OKAY;
}
#endif

/** clean up the storage */
static
SCIP_RETCODE conflictstoreCleanUpStorage(
   SCIP_CONFLICTSTORE*   conflictstore,
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set,
   SCIP_STAT*            stat,
   SCIP_PROB*            transprob
   )
{
   SCIP_CONS* conflict;
   SCIP_Bool storagecleaned;
   int ndelconfs;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);

   /* the storage is empty  */
   if( conflictstore->nconflicts == 0 )
   {
      assert(SCIPqueueNElems(conflictstore->slotqueue) == conflictstore->conflictsize);
      return SCIP_OKAY;
   }

   assert(conflictstore->nconflicts >= 1);

   /* the storage is not full */
   if( conflictstore->nconflicts < conflictstore->conflictsize )
      return SCIP_OKAY;

   ndelconfs = 0;
   storagecleaned = FALSE;
   do
   {
      int idx;

      idx = ((int) (size_t) SCIPqueueRemove(conflictstore->orderqueue)) - 1;
      assert(idx >= 0 && idx < conflictstore->conflictsize);
      assert(conflictstore->conflicts[idx] != NULL);

      /* get the oldest conflict */
      conflict = conflictstore->conflicts[idx];

      /* check whether the constraint is already marked as deleted */
      if( !SCIPconsIsDeleted(conflict) )
      {
         /* mark the constraint as deleted */
         SCIP_CALL( SCIPconsDelete(conflict, blkmem, set, stat, transprob) );

         storagecleaned = TRUE;
      }

      /* release the constraint */
      SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );

#ifdef SCIP_DEBUG
      /* clean the conflict and primal bound array */
      conflictstore->conflicts[idx] = NULL;
      conflictstore->primalbounds[idx] = SCIP_INVALID;
#endif

      /* add the id shifted by +1 to the queue of empty slots */
      SCIP_CALL( SCIPqueueInsert(conflictstore->slotqueue, (void*) (size_t) (idx+1)) );

      ++ndelconfs;
      --conflictstore->nconflicts;
   }
   while( !storagecleaned );

   SCIPdebugMessage("clean up conflict store: removed %d conflicts, %d currently stored\n",
         ndelconfs, conflictstore->nconflicts);

   return SCIP_OKAY;
}

/** creates conflict storage */
SCIP_RETCODE SCIPconflictstoreCreate(
   SCIP_CONFLICTSTORE**  conflictstore       /**< pointer to store conflict storage */
   )
{
   assert(conflictstore != NULL);

   SCIP_ALLOC( BMSallocMemory(conflictstore) );

   (*conflictstore)->conflicts = NULL;
   (*conflictstore)->primalbounds = NULL;
   (*conflictstore)->slotqueue = NULL;
   (*conflictstore)->orderqueue = NULL;
   (*conflictstore)->conflictsize = 0;
   (*conflictstore)->nconflicts = 0;
   (*conflictstore)->nconflictsfound = 0;
   (*conflictstore)->lastnodenum = -1;

   return SCIP_OKAY;
}

/** frees conflict storage */
SCIP_RETCODE SCIPconflictstoreFree(
   SCIP_CONFLICTSTORE**  conflictstore,      /**< pointer to store conflict storage */
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   SCIP_CONS* conflict;

   assert(conflictstore != NULL);
   assert(*conflictstore != NULL);

   if( (*conflictstore)->orderqueue != NULL )
   {
      assert((*conflictstore)->slotqueue != NULL);

      while( !SCIPqueueIsEmpty((*conflictstore)->orderqueue) )
      {
         int idx;

         idx = ((int) (size_t) SCIPqueueRemove((*conflictstore)->orderqueue)) - 1;
         assert(idx >= 0 && idx < (*conflictstore)->conflictsize);
         assert((*conflictstore)->conflicts[idx] != NULL);

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

   BMSfreeMemoryArrayNull(&(*conflictstore)->conflicts);
   BMSfreeMemoryArrayNull(&(*conflictstore)->primalbounds);
   BMSfreeMemory(conflictstore);

   return SCIP_OKAY;
}

/** adds a conflict to the conflict storage */
SCIP_RETCODE SCIPconflictstoreAddConflict(
   SCIP_CONFLICTSTORE*   conflictstore,
   BMS_BLKMEM*           blkmem,
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic SCIP statistics */
   SCIP_PROB*            transprob,
   SCIP_CONS*            cons,
   SCIP_NODE*            node,
   SCIP_NODE*            validnode,
   SCIP_Bool             global,
   SCIP_CONFTYPE         conftype,
   SCIP_Bool             cutoffinvolved,     /**< is a cutoff bound invaled in this conflict */
   SCIP_Real             primalbound
   )
{
   int nconflicts;
   int idx;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(cons != NULL);
   assert(node != NULL);
   assert(validnode != NULL);
   assert(set->conf_allowlocal || SCIPnodeGetDepth(validnode) == 0);
   assert(conftype != SCIP_CONFTYPE_UNKNOWN);
   assert(conftype != SCIP_CONFTYPE_BNDEXCEEDING || cutoffinvolved);
   assert(!cutoffinvolved || (cutoffinvolved && !SCIPsetIsInfinity(set, REALABS(primalbound))));

   nconflicts = conflictstore->nconflicts;

   SCIP_CALL( conflictstoreEnsureCutsMem(conflictstore, set, nconflicts+1) );

   /* return if the store has size zero */
   if( conflictstore->conflictsize == 0 )
   {
      assert(set->conf_maxstoresize == 0);
      return SCIP_OKAY;
   }

   /* clean up the storage */
   SCIP_CALL( conflictstoreCleanUpStorage(conflictstore, blkmem, set, stat, transprob) );

   /* update the last seen node */
   conflictstore->lastnodenum = SCIPnodeGetNumber(validnode);

   /* get a free slot */
   idx = ((int) (size_t) SCIPqueueRemove(conflictstore->slotqueue)-1);
   assert(idx >= 0 && idx < conflictstore->conflictsize);

   SCIPconsCapture(cons);
   conflictstore->conflicts[idx] = cons;
   conflictstore->primalbounds[idx] = primalbound;

   /* add idx shifted by +1 to the ordering queue */
   SCIP_CALL( SCIPqueueInsert(conflictstore->orderqueue, (void*) (size_t) (idx+1)) );

   ++conflictstore->nconflicts;
   ++conflictstore->nconflictsfound;

   SCIPdebugMessage("add conflict <%s> to conflict store at position %d\n", SCIPconsGetName(cons), idx);
   SCIPdebugMessage(" -> conflict type: %d, cutoff involved = %u\n", conftype, cutoffinvolved);
   SCIPdebugMessage(" -> current primal bound: %g\n", primalbound);
   SCIPdebugMessage(" -> found at node %llu (depth: %d), valid at node %llu (depth: %d)\n", SCIPnodeGetNumber(node),
         SCIPnodeGetDepth(node), SCIPnodeGetNumber(validnode), SCIPnodeGetDepth(validnode));

   return SCIP_OKAY;
}