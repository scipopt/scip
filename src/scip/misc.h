/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   misc.h
 * @brief  miscellaneous datastructures and methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MISC_H__
#define __MISC_H__


#if 0
typedef struct PQueue PQUEUE;           /**< priority queue */
#endif

typedef struct HashTable HASHTABLE;     /**< hash table */
typedef struct RealArray REALARRAY;     /**< dynamic array for storing Real values */
typedef struct IntArray INTARRAY;       /**< dynamic array for storing int values */
typedef struct BoolArray BOOLARRAY;     /**< dynamic array for storing Bool values */



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



#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "set.h"



#if 0 /* PRIORITY QUEUE NOT NEEDED */

/** initializes priority queue */
extern
RETCODE SCIPpqueueInit(
   PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int              initsize,           /**< initial number of available element slots */
   Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

/** frees priority queue, but not the data elements themselves */
extern
void SCIPpqueueFree(
   PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** inserts element into priority queue */
extern
RETCODE SCIPpqueueInsert(
   PQUEUE*          pqueue,             /**< pointer to a priority queue */
   void*            elem                /**< element to be inserted */
   );

/** removes and returns best element from the priority queue */
extern
void* SCIPpqueueRemove(
   PQUEUE*          pqueue              /**< pointer to a priority queue */
   );

/** returns the best element of the queue without removing it */
extern
void* SCIPpqueueFirst(
   const PQUEUE*    pqueue              /**< pointer to a priority queue */
   );
#endif


/*
 * Hash Table
 */

/** creates a hash table */
extern
RETCODE SCIPhashtableCreate(
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   int              tablesize,          /**< size of the hash table */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval))       /**< returns the hash value of the key */
   );

/** frees the hash table */
extern
void SCIPhashtableFree(
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   MEMHDR*          memhdr              /**< block memory */
   );

/** inserts element in hash table (multiple inserts of same element possible) */
extern
RETCODE SCIPhashtableInsert(
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to append to the list */
   );

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
extern
RETCODE SCIPhashtableSafeInsert(
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to insert into the table */
   );

/** retrieve element with key from hash table, returns NULL if not existing */
extern
void* SCIPhashtableRetrieve(
   HASHTABLE*       hashtable,          /**< hash table */
   void*            key                 /**< key to retrieve */
   );

/** removes existing element from the hash table */
extern
RETCODE SCIPhashtableRemove(
   HASHTABLE*       hashtable,          /**< hash table */
   MEMHDR*          memhdr,             /**< block memory */
   void*            element             /**< element to remove from the table */
   );

/** prints statistics about hash table usage */
extern
void SCIPhashtablePrintStatistics(
   HASHTABLE*       hashtable           /**< hash table */
   );

/** standard hash key comparator for string keys */
extern
DECL_HASHKEYEQ(SCIPhashKeyEqString);

/** standard hashing function for string keys */
extern
DECL_HASHKEYVAL(SCIPhashKeyValString);



/*
 * Dynamic Arrays
 */

/** creates a dynamic array of real values */
extern
RETCODE SCIPrealarrayCreate(
   REALARRAY**      realarray,          /**< pointer to store the real array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
extern
RETCODE SCIPrealarrayCopy(
   REALARRAY**      realarray,          /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   REALARRAY*       sourcerealarray     /**< dynamic real array to copy */
   );

/** frees a dynamic array of real values */
extern
RETCODE SCIPrealarrayFree(
   REALARRAY**      realarray           /**< pointer to the real array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPrealarrayExtend(
   REALARRAY*       realarray,          /**< dynamic real array */
   const SET*       set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic real array */
extern
void SCIPrealarrayClear(
   REALARRAY*       realarray           /**< dynamic real array */
   );

/** gets value of entry in dynamic array */
extern
Real SCIPrealarrayGetVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPrealarraySetVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Real             val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
RETCODE SCIPrealarrayIncVal(
   REALARRAY*       realarray,          /**< dynamic real array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   Real             incval              /**< value to increase array index */
   );

/** creates a dynamic array of int values */
extern
RETCODE SCIPintarrayCreate(
   INTARRAY**       intarray,           /**< pointer to store the int array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
extern
RETCODE SCIPintarrayCopy(
   INTARRAY**       intarray,           /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   INTARRAY*        sourceintarray      /**< dynamic real array to copy */
   );

/** frees a dynamic array of int values */
extern
RETCODE SCIPintarrayFree(
   INTARRAY**       intarray            /**< pointer to the int array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPintarrayExtend(
   INTARRAY*        intarray,           /**< dynamic int array */
   const SET*       set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic int array */
extern
void SCIPintarrayClear(
   INTARRAY*        intarray            /**< dynamic int array */
   );

/** gets value of entry in dynamic array */
extern
int SCIPintarrayGetVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPintarraySetVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   int              val                 /**< value to set array index to */
   );

/** increases value of entry in dynamic array */
extern
RETCODE SCIPintarrayIncVal(
   INTARRAY*        intarray,           /**< dynamic int array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to increase value for */
   int              incval              /**< value to increase array index */
   );

/** creates a dynamic array of bool values */
extern
RETCODE SCIPboolarrayCreate(
   BOOLARRAY**      boolarray,          /**< pointer to store the bool array */
   MEMHDR*          memhdr              /**< block memory */
   );

/** creates a copy of a dynamic array of real values */
extern
RETCODE SCIPboolarrayCopy(
   BOOLARRAY**      boolarray,          /**< pointer to store the copied real array */
   MEMHDR*          memhdr,             /**< block memory */
   BOOLARRAY*       sourceboolarray     /**< dynamic real array to copy */
   );

/** frees a dynamic array of bool values */
extern
RETCODE SCIPboolarrayFree(
   BOOLARRAY**      boolarray           /**< pointer to the bool array */
   );

/** extends dynamic array to be able to store indices from minidx to maxidx */
extern
RETCODE SCIPboolarrayExtend(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   const SET*       set,                /**< global SCIP settings */
   int              minidx,             /**< smallest index to allocate storage for */
   int              maxidx              /**< largest index to allocate storage for */
   );

/** clears a dynamic bool array */
extern
void SCIPboolarrayClear(
   BOOLARRAY*       boolarray           /**< dynamic bool array */
   );

/** gets value of entry in dynamic array */
extern
Bool SCIPboolarrayGetVal(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int              idx                 /**< array index to get value for */
   );

/** sets value of entry in dynamic array */
extern
RETCODE SCIPboolarraySetVal(
   BOOLARRAY*       boolarray,          /**< dynamic bool array */
   const SET*       set,                /**< global SCIP settings */
   int              idx,                /**< array index to set value for */
   Bool             val                 /**< value to set array index to */
   );





/*
 * Sorting algorithms
 */

/** bubble sort of an array of pointers */
extern
void SCIPbsortPtr(
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

/** bubble sort of two joint arrays of pointers/Reals, sorted by first array */
extern
void SCIPbsortPtrDbl(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

/** bubble sort of three joint arrays of pointers/Reals/Ints, sorted by first */
extern
void SCIPbsortPtrDblInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );

/** bubble sort of four joint arrays of pointers/Reals/Ints/Ints, sorted by first */
extern
void SCIPbsortPtrDblIntInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            dblarray,           /**< Real array to be permuted in the same way */
   int*             intarray1,          /**< first int array to be permuted in the same way */
   int*             intarray2,          /**< second int array to be permuted in the same way */
   int              len,                /**< length of both arrays */
   DECL_SORTPTRCOMP((*ptrcmp))          /**< data element comparator */
   );




/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
Real SCIPcalcMachineEpsilon(
   void
   );

/** calculates the greatest common divisor of the two given values */
extern
Longint SCIPcalcGreComDiv(
   Longint          val1,               /**< first value of greatest common devisor calculation */
   Longint          val2                /**< second value of greatest common devisor calculation */
   );

/** calculates the smallest common multiple of the two given values */
extern
Longint SCIPcalcSmaComMul(
   Longint          val1,               /**< first value of greatest common devisor calculation */
   Longint          val2                /**< second value of greatest common devisor calculation */
   );



/*
 * File methods
 */

/** returns, whether the given file exists */
extern
Bool SCIPfileExists(
   const char*      filename            /**< file name */
   );

#endif
