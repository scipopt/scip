/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sort.h
 * @brief  datastructures and algorithms for sorting and queueing data elements
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SORT_H__
#define __SORT_H__


#if 0
typedef struct PQueue PQUEUE;           /**< priority queue */
#endif

typedef struct HashTable HASHTABLE;     /**< hash table */


#include "def.h"
#include "retcode.h"
#include "memory.h"



/** compares two data element pointers
 *  result:
 *    < 0: elem1 comes before (is better than) elem2
 *    = 0: both elements have the same value
 *    > 0: elem2 comes after (is worse than) elem2
 */
#define DECL_SORTPTRCOMP(x) int x (void* elem1, void* elem2)

/** gets the key of the given element */
#define DECL_HASHGETKEY(x) void* x (void* elem)

/** returns TRUE iff both keys are equal */
#define DECL_HASHKEYEQ(x) Bool x (void* key1, void* key2)

/** returns the hash value of the key */
#define DECL_HASHKEYVAL(x) unsigned int x (void* key)


#if 0 /* PRIORITY QUEUE NOT NEEDED */
extern
RETCODE SCIPpqueueInit(                 /**< initializes priority queue */
   PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int              initsize,           /**< initial number of available element slots */
   Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

extern
void SCIPpqueueFree(                    /**< frees priority queue, but not the data elements themselves */
   PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

extern
RETCODE SCIPpqueueInsert(               /**< inserts element into priority queue */
   PQUEUE*          pqueue,             /**< pointer to a priority queue */
   void*            elem                /**< element to be inserted */
   );

extern
void* SCIPpqueueRemove(                 /**< removes and returns best element from the priority queue */
   PQUEUE*          pqueue              /**< pointer to a priority queue */
   );

extern
void* SCIPpqueueFirst(                  /**< returns the best element of the queue without removing it */
   const PQUEUE*    pqueue              /**< pointer to a priority queue */
   );
#endif


/*
 * Hash Table
 */

extern
RETCODE SCIPhashtableCreate(            /**< creates a hash table */
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   int              tablesize,          /**< size of the hash table */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval))       /**< returns the hash value of the key */
   );

extern
void SCIPhashtableFree(                 /**< frees the hash table */
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPhashtableInsert(            /**< inserts element in hash table (multiple inserts of same element possible) */
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to append to the list */
   );

extern
void* SCIPhashtableRetrieve(            /**< retrieve element with key from hash table, returns NULL if not existing */
   HASHTABLE*       hashtable,          /**< hash table */
   void*            key                 /**< key to retrieve */
   );

extern
RETCODE SCIPhashtableRemove(            /**< removes existing element from the hash table */
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to remove from the table */
   );

extern
void SCIPhashtablePrintStatistics(      /**< prints statistics about hash table usage */
   HASHTABLE*       hashtable           /**< hash table */
   );

extern
DECL_HASHKEYEQ(SCIPhashKeyEqString);    /**< standard hash key comparator for string keys */

extern
DECL_HASHKEYVAL(SCIPhashKeyValString);  /**< standard hashing function for string keys */



/*
 * Sorting algorithms
 */

extern
void SCIPbsortPtr(                      /**< bubble sort of an array of pointers */
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

extern
void SCIPbsortPtrDbl(                   /**< bubble sort of two joint arrays of pointers/Reals, sorted by first array */
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

extern
void SCIPbsortPtrDblInt(                /**< bubble sort of three joint arrays of pointers/Reals/Ints, sorted by first */
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

#endif
