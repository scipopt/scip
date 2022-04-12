/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   stpprioqueue.c
 * @brief  priority queue with integer keys
 * @author Daniel Rehfeldt
 *
 * Implements a (minimum) priority queue with integer keys.
 * A list of all interface methods can be found in stpprioqueue.h
 ** todo: if needed for other key type, either use template pattern
 * or intrusive design with macros...
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "stpprioqueue.h"


#define STPPQ_MIN_KEY  INT_MIN


/** PQ entry */
typedef struct stp_priority_queue_entry
{
   void*                 data;
   int                   key;
} PQENTRY;


/** PQ with integer keys */
struct stp_priority_queue
{
   int                   capacity;           /**< maximum size */
   int                   size;               /**< size */
   PQENTRY*              entries;            /**< entries  */
};


/*
 * Local methods
 */


/*
 * Interface methods
 */


/** clean the priority queue */
void stpprioqueue_clean(
   STP_PQ*               prioqueue           /**< the priority queue  */
   )
{
   assert(prioqueue);

   prioqueue->size = 0;
   assert(stpprioqueue_isClean(prioqueue));
}


/** is the priority queue clean? */
SCIP_Bool stpprioqueue_isClean(
   const STP_PQ*          prioqueue          /**< the priority queue  */
   )
{
   assert(prioqueue);

   if( prioqueue->size != 0 )
   {
      SCIPdebugMessage("prioqueue size not 0! (=%d)\n", prioqueue->size);
      return FALSE;
   }

   return TRUE;
}


/** creates new prioqueue. If entries array is provided, it must be of size capacity + 2  */
SCIP_RETCODE stpprioqueue_create(
   SCIP*                 scip,               /**< SCIP */
   int                   capacity,           /**< initial capacity */
   STP_PQ**              prioqueue           /**< the priority queue  */
   )
{
   PQENTRY* entries_heap;

   assert(scip && prioqueue);
   assert(capacity >= 1);

   SCIP_CALL( SCIPallocMemory(scip, prioqueue) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(entries_heap), capacity + 1) );

   (*prioqueue)->capacity = capacity;
   (*prioqueue)->entries = entries_heap;

   /* sentinel */
   entries_heap[0].key = STPPQ_MIN_KEY;

   stpprioqueue_clean(*prioqueue);

   return SCIP_OKAY;
}

/** frees the priority queue */
void stpprioqueue_free(
   SCIP*                 scip,               /**< SCIP */
   STP_PQ**              prioqueue           /**< the priority queue  */
   )
{
   assert(scip && prioqueue);

   SCIPfreeBlockMemoryArray(scip, &((*prioqueue)->entries), (*prioqueue)->capacity + 1);
   SCIPfreeMemory(scip, prioqueue);
}

/** shows top data */
const void* stpprioqueue_peakMinData(
   const STP_PQ*         prioqueue           /**< the priority queue  */
   )
{
   assert(prioqueue);
   assert(prioqueue->size > 0);

   return prioqueue->entries[1].data;
}


/** shows top key */
int stpprioqueue_peakMinKey(
   const STP_PQ*         prioqueue           /**< the priority queue  */
   )
{
   assert(prioqueue);
   assert(prioqueue->size > 0);

   return prioqueue->entries[1].key;
}


/** deletes minimum and gives key and data */
void stpprioqueue_deleteMin(
   void**                data,               /**< pointer to data of minimum */
   int*                  key,                /**< pointer to key of minimum */
   STP_PQ*               prioqueue           /**< the priority queue  */
   )
{
   assert(data && key && prioqueue);

   *key = prioqueue->entries[1].key;
   *data = stpprioqueue_deleteMinReturnData(prioqueue);
}


/** deletes minimum and returns corresponding data */
void* stpprioqueue_deleteMinReturnData(
   STP_PQ*                prioqueue          /**< the priority queue  */
   )
{
   PQENTRY* const RESTRICT entries = prioqueue->entries;
   int fill;
   int parent;
   int hole = 1;
   int child = 2;
   void* data;
   const int lastentry = prioqueue->size--;

   assert(prioqueue && entries);
   assert(prioqueue->size >= 0);

   data = entries[1].data;

   /* move down along min-path */
   while( child < lastentry )
   {
      const SCIP_Real key1 = entries[child].key;
      const SCIP_Real key2 = entries[child + 1].key;
      assert(hole >= 1);

      /* second child with smaller key? */
      if( key1 > key2 )
      {
         entries[hole].key = key2;
         child++;
      }
      else
      {
         entries[hole].key = key1;
      }

      entries[hole].data = entries[child].data;

      hole = child;
      child *= 2;
   }

   /* now hole is at last tree level, fill it with last PQ entry and move it up */

   fill = entries[lastentry].key;
   parent = hole / 2;

   while( entries[parent].key > fill )
   {
      assert(hole >= 1);

      entries[hole] = entries[parent];

      hole = parent;
      parent /= 2;
   }

   /* finally, fill the hole */
   entries[hole].key = fill;
   entries[hole].data = entries[lastentry].data;

   return data;
}


/** inserts data */
SCIP_RETCODE stpprioqueue_insert(
   SCIP*                 scip,               /**< SCIP */
   void*                 data,               /**< the data */
   int                   newkey,             /**< the key */
   STP_PQ*               prioqueue           /**< the priority queue  */
   )
{
   PQENTRY* RESTRICT entries;
   int hole;
   int parent;
   int parentkey;

   assert(scip && prioqueue);
   assert(prioqueue->entries);
   assert(newkey > STPPQ_MIN_KEY);
   assert(prioqueue->size <= prioqueue->capacity);

   if( prioqueue->size == prioqueue->capacity )
   {
      const int oldsize = prioqueue->capacity + 1;
      prioqueue->capacity *= 2;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(prioqueue->entries), oldsize, prioqueue->capacity + 1) );
   }

   entries = prioqueue->entries;
   hole = ++(prioqueue->size);
   parent = hole / 2;
   parentkey = entries[parent].key;

   /* move hole up */
   while( parentkey > newkey )
   {
      assert(hole >= 1);

      entries[hole].key = parentkey;
      entries[hole].data = entries[parent].data;
      hole = parent;
      parent /= 2;
      parentkey = entries[parent].key;
   }

   /* fill the hole */
   entries[hole].key = newkey;
   entries[hole].data = data;

   return SCIP_OKAY;
}
