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

#include "scip/conflictstore.h"
#include "scip/cons.h"
#include "scip/set.h"
#include "scip/tree.h"

#include "scip/struct_conflictstore.h"


#define DEFAULT_CONFLICTSTORE_SIZE        8192 /* maximal size of conflict storage */

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
      int j;

      if( conflictstore->conflictsize == 0 )
         newsize = MIN(set->conf_maxstoresize, DEFAULT_CONFLICTSTORE_SIZE);
      else
         newsize = MIN(set->conf_maxstoresize, conflictstore->conflictsize * 2);

      SCIP_ALLOC( BMSreallocMemoryArray(&conflictstore->conflicts, newsize) );

      /* move all conflicts stored in slots before firstused to the end */
      j = conflictstore->conflictsize;
      conflictstore->firstfree = conflictstore->conflictsize;
      for( i = 0; i < conflictstore->firstused; i++ )
      {
         assert(j < newsize);
         conflictstore->conflicts[j] = conflictstore->conflicts[i];
         conflictstore->conflicts[i] = NULL;
         ++conflictstore->firstfree;
      }
      conflictstore->conflictsize = newsize;
   }
   assert(num <= conflictstore->conflictsize);

   return SCIP_OKAY;
}


/* update pointer to the first used slot */
static
void conflictstoreUpdateFirstused(
   SCIP_CONFLICTSTORE*   conflictstore
   )
{
   assert(conflictstore != NULL);

   if( conflictstore->nconflicts == 0 )
      conflictstore->firstused = -1;
   else
   {
      ++conflictstore->firstused;
      conflictstore->firstused %= conflictstore->conflictsize;
   }
}


/* update pointer to the first free slot */
static
void conflictstoreUpdateFirstfree(
   SCIP_CONFLICTSTORE*   conflictstore
   )
{
   assert(conflictstore != NULL);

   if( conflictstore->nconflicts == conflictstore->conflictsize )
      conflictstore->firstfree = -1;
   else
   {
      ++conflictstore->firstfree;
      conflictstore->firstfree %= conflictstore->conflictsize;
   }
}

/** remove all deleted conflicts from the storage
 *
 * TODO currently, this only counts the number of conflicts marked as deleted
 */
static
SCIP_RETCODE cleanDeletedConflicts(
   SCIP_CONFLICTSTORE*   conflictstore
   )
{
   int ndeletedconflicts;
   int i;
   int j;

   assert(conflictstore != NULL);

   ndeletedconflicts = 0;
   for( i = 0; i < conflictstore->conflictsize; i++ )
   {
      j = (i + conflictstore->firstused) % conflictstore->conflictsize;

      if( j == conflictstore->firstfree )
         break;

      if( SCIPconsIsDeleted(conflictstore->conflicts[j]) )
         ++ndeletedconflicts;
   }

   printf("%d/%d conflicts marked as deleted at node %lld\n", ndeletedconflicts, conflictstore->nconflicts, conflictstore->lastnodenum);

   return SCIP_OKAY;
}

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

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);

   /* the storage is empty  */
   if( conflictstore->firstused == -1 )
      return SCIP_OKAY;

   assert(conflictstore->nconflicts >= 1);

   SCIP_CALL( cleanDeletedConflicts(conflictstore) );

   /* the storage is not full */
   if( conflictstore->nconflicts < conflictstore->conflictsize )
      return SCIP_OKAY;

   conflict = conflictstore->conflicts[conflictstore->firstused];
   assert(conflict != NULL);

   SCIPdebugMessage("clean up conflict store\n");

   SCIP_CALL( SCIPconsDelete(conflict, blkmem, set, stat, transprob) );
   SCIP_CALL( SCIPconsRelease(&conflict, blkmem, set) );
   conflictstore->conflicts[conflictstore->firstused] = NULL;
   --conflictstore->nconflicts;

   /* update pointer to the first free and first used slot */
   conflictstore->firstfree = conflictstore->firstused;
   conflictstoreUpdateFirstused(conflictstore);

   assert(conflictstore->firstused != conflictstore->firstfree);
   assert(conflictstore->firstused >= 0);
   assert(conflictstore->firstfree >= 0);
   assert(conflictstore->firstused < conflictstore->firstfree
       || conflictstore->conflictsize - conflictstore->nconflicts == conflictstore->firstused - conflictstore->firstfree);
   assert(conflictstore->firstused > conflictstore->firstfree
       || conflictstore->nconflicts == conflictstore->firstfree - conflictstore->firstused);

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
   (*conflictstore)->conflictsize = 0;
   (*conflictstore)->nconflicts = 0;
   (*conflictstore)->nconflictsfound = 0;
   (*conflictstore)->firstfree = 0;
   (*conflictstore)->firstused = -1;
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
   int i;

   assert(conflictstore != NULL);
   assert(*conflictstore != NULL);

   for( i = 0; i < (*conflictstore)->conflictsize; i++ )
   {
      int j;

      j = (i + (*conflictstore)->firstused) % (*conflictstore)->conflictsize;

      if( j == (*conflictstore)->firstfree )
         break;

      SCIP_CALL( SCIPconsRelease(&(*conflictstore)->conflicts[j], blkmem, set) );

      --(*conflictstore)->nconflicts;
   }
   assert((*conflictstore)->nconflicts == 0);

   BMSfreeMemoryArrayNull(&(*conflictstore)->conflicts);
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
   SCIP_Bool             global
   )
{
   int nconflicts;

   assert(conflictstore != NULL);
   assert(blkmem != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(cons != NULL);
   assert(node != NULL);

   nconflicts = conflictstore->nconflicts;

   SCIP_CALL( conflictstoreEnsureCutsMem(conflictstore, set, nconflicts+1) );

   /* return if the store has size zero */
   if( conflictstore->conflictsize == 0 )
   {
      assert(set->conf_maxstoresize == 0);

      /* we need to release the constraint */
      SCIPconsRelease(&cons, blkmem, set);

      return SCIP_OKAY;
   }

   /* clean up the storage */
   SCIP_CALL( conflictstoreCleanUpStorage(conflictstore, blkmem, set, stat, transprob) );

   /* update the last seen node */
   conflictstore->lastnodenum = SCIPnodeGetNumber(node);

   SCIPconsCapture(cons);
   conflictstore->conflicts[conflictstore->firstfree] = cons;
   ++conflictstore->nconflicts;
   ++conflictstore->nconflictsfound;

   if( conflictstore->nconflicts == 1 )
      conflictstore->firstused = 0;

   SCIPdebugMessage("add conflict <%s> to conflict store at position %d\n", SCIPconsGetName(cons), conflictstore->firstfree);

   /* update pointer to the first free slot */
   conflictstoreUpdateFirstfree(conflictstore);

   return SCIP_OKAY;
}