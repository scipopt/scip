/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_misc.h
 * @ingroup PUBLICMETHODS
 * @brief  public miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MISC_H__
#define __SCIP_PUB_MISC_H__



#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_misc.h"
#include "scip/type_message.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Priority Queue
 */

/** creates priority queue */
extern
SCIP_RETCODE SCIPpqueueCreate(
   SCIP_PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int                   initsize,           /**< initial number of available element slots */
   SCIP_Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** frees priority queue, but not the data elements themselves */
extern
void SCIPpqueueFree(
   SCIP_PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** clears the priority queue, but doesn't free the data elements themselves */
extern
void SCIPpqueueClear(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** inserts element into priority queue */
extern
SCIP_RETCODE SCIPpqueueInsert(
   SCIP_PQUEUE*          pqueue,             /**< priority queue */
   void*                 elem                /**< element to be inserted */
   );

/** removes and returns best element from the priority queue */
extern
void* SCIPpqueueRemove(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the best element of the queue without removing it */
extern
void* SCIPpqueueFirst(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the number of elements in the queue */
extern
int SCIPpqueueNElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
extern
void** SCIPpqueueElems(
   SCIP_PQUEUE*          pqueue              /**< priority queue */
   );




/*
 * Hash Table
 */

/** returns a reasonable hash table size (a prime number) that is at least as large as the specified value */
extern
int SCIPcalcHashtableSize(
   int                   minsize             /**< minimal size of the hash table */
   );

/** creates a hash table */
extern
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval)),      /**< returns the hash value of the key */
   void*                 userptr             /**< user pointer */
   );

/** frees the hash table */
extern
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
   );

/** removes all elements of the hash table
 *
 *  @note From a performance point of view you should not fill and clear a hash table too often since the clearing can
 *        be expensive. Clearing is done by looping over all buckets and removing the hash table lists one-by-one.
 */
extern
void SCIPhashtableClear(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** inserts element in hash table (multiple inserts of same element possible) */
extern
SCIP_RETCODE SCIPhashtableInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
extern
SCIP_RETCODE SCIPhashtableSafeInsert(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to insert into the table */
   );

/** retrieve element with key from hash table, returns NULL if not existing */
extern
void* SCIPhashtableRetrieve(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 key                 /**< key to retrieve */
   );

/** retrieve element with key from hash table, returns NULL if not existing
 * can be used to retrieve all entries with the same key (one-by-one) */
extern
void* SCIPhashtableRetrieveNext(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_HASHTABLELIST**  hashtablelist,      /**< input: entry in hash table list from which to start searching, or NULL; output: entry in hash table list corresponding to element after retrieved one, or NULL */
   void*                 key                 /**< key to retrieve */
   );

/** returns whether the given element exists in the table */
extern
SCIP_Bool SCIPhashtableExists(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to search in the table */
   );

/** removes element from the hash table, if it exists */
extern
SCIP_RETCODE SCIPhashtableRemove(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   void*                 element             /**< element to remove from the table */
   );

/** clears hash table */
extern
void SCIPhashtableRemoveAll(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
   );

/** prints statistics about hash table usage */
extern
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable,          /**< hash table */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** standard hash key comparator for string keys */
extern
SCIP_DECL_HASHKEYEQ(SCIPhashKeyEqString);

/** standard hashing function for string keys */
extern
SCIP_DECL_HASHKEYVAL(SCIPhashKeyValString);




/*
 * Hash Map
 */

/** creates a hash map mapping pointers to pointers */
extern
SCIP_RETCODE SCIPhashmapCreate(
   SCIP_HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash map entries */
   int                   mapsize             /**< size of the hash map */
   );

/** frees the hash map */
extern
void SCIPhashmapFree(
   SCIP_HASHMAP**        hashmap             /**< pointer to the hash map */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
extern
SCIP_RETCODE SCIPhashmapInsert(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** retrieves image of given origin from the hash map, or NULL if no image exists */
extern
void* SCIPhashmapGetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to retrieve image for */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
extern
SCIP_RETCODE SCIPhashmapSetImage(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin,             /**< origin to set image for */
   void*                 image               /**< new image for origin */
   );

/** checks whether an image to the given origin exists in the hash map */
extern
SCIP_Bool SCIPhashmapExists(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to search for */
   );

/** removes origin->image pair from the hash map, if it exists */
extern
SCIP_RETCODE SCIPhashmapRemove(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   void*                 origin              /**< origin to remove from the list */
   );

/** prints statistics about hash map usage */
extern
void SCIPhashmapPrintStatistics(
   SCIP_HASHMAP*         hashmap,            /**< hash map */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   );

/** indicates whether a hash map has no entries */
extern
SCIP_Bool SCIPhashmapIsEmpty(
   SCIP_HASHMAP*      hashmap          /**< hash map */
);

/** gives the number of entries in a hash map */ 
extern
int SCIPhashmapGetNEntries(
   SCIP_HASHMAP*      hashmap          /**< hash map */
);

/** gives the number of lists (buckets) in a hash map */ 
extern
int SCIPhashmapGetNLists(
   SCIP_HASHMAP*      hashmap          /**< hash map */
);

/** gives a specific list (bucket) in a hash map */
extern
SCIP_HASHMAPLIST* SCIPhashmapGetList(
   SCIP_HASHMAP*     hashmap,          /**< hash map */
   int               listindex         /**< index of hash map list */
);

/** gives the number of entries in a list of a hash map */ 
extern
int SCIPhashmapListGetNEntries(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list, can be NULL */
);

/** retrieves origin of given entry in a hash map */ 
extern
void* SCIPhashmapListGetOrigin(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list */
);

/** retrieves image of given entry in a hash map */ 
extern
void* SCIPhashmapListGetImage(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list */
);

/** retrieves next entry from given entry in a hash map list, or NULL if at end of list. */ 
extern
SCIP_HASHMAPLIST* SCIPhashmapListGetNext(
   SCIP_HASHMAPLIST* hashmaplist       /**< hash map list */
);

/** removes all entries in a hash map. */ 
extern
SCIP_RETCODE SCIPhashmapRemoveAll(
   SCIP_HASHMAP*     hashmap           /**< hash map */
);



/*
 * Sorting algorithms
 */

/* first all upwards-sorting methods */

/** sort an indexed element set in non-decreasing order, resulting in a permutation index array */
extern
void SCIPsort(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-decreasing order */
extern
void SCIPsortInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-decreasing order */
extern
void SCIPsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );


/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Reals in non-decreasing order */
extern
void SCIPsortReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,           /**< pointer array to be permuted in the same way */
   void**                ptrarray2,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
extern
void SCIPsortRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-decreasing order */
extern
void SCIPsortInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortIntIntIntPtr(
   int*                  intarray1,           /**< int array to be sorted */
   int*                  intarray2,           /**< int array to be permuted in the same way */
   int*                  intarray3,           /**< int array to be permuted in the same way */
   void**                ptrarray,            /**< pointer array to be permuted in the same way */
   int                   len                  /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortIntPtrIntReal(
   int*                  intarray1,           /**< int array to be sorted */
   void**                ptrarray,            /**< pointer array to be permuted in the same way */
   int*                  intarray2,           /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,           /**< SCIP_Real array to be permuted in the same way */
   int                   len                  /**< length of arrays */
   );

/** sort an array of Longints in non-decreasing order */
extern
void SCIPsortLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
extern
void SCIPsortLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
extern
void SCIPsortLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/* now all downwards-sorting methods */

/** sort an indexed element set in non-increasing order, resulting in a permutation index array */
extern
void SCIPsortDown(
   int*                  perm,               /**< pointer to store the resulting permutation */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< number of elements to be sorted (valid index range) */
   );

/** sort an index array in non-increasing order */
extern
void SCIPsortDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len                 /**< length of array */
   );

/** sort of an array of pointers in non-increasing order */
extern
void SCIPsortDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of array */
   );

/** sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   len                 /**< length of arrays */
   );

/** sort an array of Reals in non-increasing order */
extern
void SCIPsortDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );


/** sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,           /**< pointer array to be permuted in the same way */
   void**                ptrarray2,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort array of ints in non-increasing order */
extern
void SCIPsortDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of ints/reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntIntIntPtr(
   int*                  intarray1,           /**< int array to be sorted */
   int*                  intarray2,           /**< int array to be permuted in the same way */
   int*                  intarray3,           /**< int array to be permuted in the same way */
   void**                ptrarray,            /**< pointer array to be permuted in the same way */
   int                   len                  /**< length of arrays */
   );

/** sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortDownIntPtrIntReal(
   int*                  intarray1,           /**< int array to be sorted */
   void**                ptrarray,            /**< pointer array to be permuted in the same way */
   int*                  intarray2,           /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,           /**< SCIP_Real array to be permuted in the same way */
   int                   len                  /**< length of arrays */
   );

/** sort an array of Longints in non-increasing order */
extern
void SCIPsortDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   len                 /**< length of arrays */
   );

/** sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
extern
void SCIPsortDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of two three arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
extern
void SCIPsortDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len                 /**< length of arrays */
   );

/** sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
extern
void SCIPsortDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   len                 /**< length of arrays */
   );


/*
 * Sorted vectors
 */

/* upwards insertion */

/** insert a new element into an index array in non-decreasing order */
extern
void SCIPsortedvecInsertInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-decreasing order */
extern
void SCIPsortedvecInsertPtr(
   void**                ptrarray,           /**< pointer to the pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,          /**< pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   intval,             /**<  additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an arrays of Reals, sorted in non-decreasing order */
extern
void SCIPsortedvecInsertReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted  */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Longint          field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealRealIntInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< third SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   void*                 field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of ints in non-decreasing order */
extern
void SCIPsortedvecInsertInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntIntPtr(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntIntIntPtr(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int*                  intarray3,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertIntPtrIntReal(
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Longints, sorted in non-decreasing order */
extern
void SCIPsortedvecInsertLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecInsertLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecInsertLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecInsertLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/* downwards insertion */

/** insert a new element into an index array in non-increasing order */
extern
void SCIPsortedvecInsertDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be inserted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of pointers in non-increasing order */
extern
void SCIPsortedvecInsertDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Reals, sorted in non-increasing order */
extern
void SCIPsortedvecInsertDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Bool             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray1,           /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,           /**< second pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array  where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Longint          field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Real             keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   void*                 field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of ints in non-increasing order */
extern
void SCIPsortedvecInsertDownInt(
   int*                  intarray,           /**< int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntReal(
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   SCIP_Real             field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Longint          field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   SCIP_Real             field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   int*                  intarray3,          /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   int                   field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   void*                 field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/** insert a new element into four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray2,          /**< int array where an element is to be inserted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be inserted */
   int                   keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   SCIP_Real             field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into an array of Longints, sorted in non-increasing order */
extern
void SCIPsortedvecInsertDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into two joint arrays of Long/pointer, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecInsertDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into three joint arrays of Long/pointer/ints, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecInsertDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   int                   field2val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   int*                  intarray1,          /**< first int array where an element is to be inserted */
   int*                  intarray2,          /**< second int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   int                   field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );

/** insert a new element into five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecInsertDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be inserted */
   void**                ptrarray1,          /**< first pointer array where an element is to be inserted */
   void**                ptrarray2,          /**< second pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   int*                  intarray,           /**< int array where an element is to be inserted */
   SCIP_Longint          keyval,             /**< key value of new element */
   void*                 field1val,          /**< additional value of new element */
   void*                 field2val,          /**< additional value of new element */
   SCIP_Bool             field3val,          /**< additional value of new element */
   int                   field4val,          /**< additional value of new element */
   int*                  len,                /**< pointer to length of arrays (will be increased by 1) */
   int*                  pos                 /**< pointer to store the insertion position, or NULL */
   );


/* upwards position deletion */

/** delete the element at the given position from an index array in non-decreasing order */
extern
void SCIPsortedvecDelPosInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-decreasing order */
extern
void SCIPsortedvecDelPosPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an arrays of Reals, sorted in non-decreasing order */
extern
void SCIPsortedvecDelPosReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealPtrPtrInt(
   SCIP_Real*            realarray,         /**<  first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,           /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,           /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array  where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of ints in non-decreasing order */
extern
void SCIPsortedvecDelPosInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/pointers/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntPtrReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted by in non-decreasing order */
extern
void SCIPsortedvecDelPosLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Long/pointer, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/* downwards position deletion */

/** delete the element at the given position from an index array in non-increasing order */
extern
void SCIPsortedvecDelPosDownInd(
   int*                  indarray,           /**< pointer to the index array where an element is to be deleted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of pointers in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtr(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrReal(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of pointers/Bools, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrBool(
   void**                ptrarray,           /**< pointer array where an element is to be inserted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be inserted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be increased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrRealInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** deletes the element at the given position from five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Reals, sorted in non-increasing order */
extern
void SCIPsortedvecDelPosDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from three joint arrays of Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealPtrPtr(
   SCIP_Real*            realarray,         /**< first SCIP_Real array where an element is to be deleted */
   void**                ptrarray1,           /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,           /**< second pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Long/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Real*            realarray3,         /**< SCIP_Real array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of ints in non-increasing order */
extern
void SCIPsortedvecDelPosDownInt(
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntInt(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntReal(
   int*                  intarray,           /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntLong(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from two joint arrays of ints/pointers, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownIntPtr(
   int*                  intarray,           /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/** delete the element at the given position from four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosDownIntIntIntPtr(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   int*                  intarray3,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-decreasing order */
extern
void SCIPsortedvecDelPosDownIntPtrIntReal(
   int*                  intarray1,          /**< int array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray2,          /**< int array where an element is to be deleted */
   SCIP_Real*            realarray,          /**< SCIP_Real array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from an array of Longints, sorted in non-increasing order */
extern
void SCIPsortedvecDelPosDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three two arrays of Long/pointer, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from three joint arrays of Long/pointer/int, sorted by the first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray,           /**< pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   int*                  intarray1,          /**< first int array where an element is to be deleted */
   int*                  intarray2,          /**< second int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );

/** delete the element at the given position from five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order */
extern
void SCIPsortedvecDelPosDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array where an element is to be deleted */
   void**                ptrarray1,          /**< first pointer array where an element is to be deleted */
   void**                ptrarray2,          /**< second pointer array where an element is to be deleted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array where an element is to be deleted */
   int*                  intarray,           /**< int array where an element is to be deleted */
   int                   pos,                /**< array position of element to be deleted */
   int*                  len                 /**< pointer to length of arrays (will be decreased by 1) */
   );


/* upwards binary search */

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindInd(
   int*                  indarray,           /**< index array to be searched */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindPtr(
   void**                ptrarray,           /**< pointer array to be searched */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be searched */
   SCIP_Real             val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindInt(
   int*                  intarray,           /**< int array to be searched */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be searched */
   SCIP_Longint          val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );


/* downwards binary search */

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindDownInd(
   int*                  indarray,           /**< index array to be searched */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindDownPtr(
   void**                ptrarray,           /**< pointer array to be searched */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   void*                 val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be searched */
   SCIP_Real             val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindDownInt(
   int*                  intarray,           /**< int array to be searched */
   int                   val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );

/** Finds the position at which 'val' is located in the sorted vector by binary search.
 *  If the element exists, the method returns TRUE and stores the position of the element in '*pos'.
 *  If the element does not exist, the method returns FALSE and stores the position of the element that follows
 *  'val' in the ordering in '*pos', i.e., '*pos' is the position at which 'val' would be inserted.
 *  Note that if the element is not found, '*pos' may be equal to len if all existing elements are smaller than 'val'.
 */
extern
SCIP_Bool SCIPsortedvecFindDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be searched */
   SCIP_Longint          val,                /**< value to search */
   int                   len,                /**< length of array */
   int*                  pos                 /**< pointer to store position of element */
   );


/*
 * Stair map
 */

/** creates stair map */
extern
SCIP_RETCODE SCIPstairmapCreate(
   SCIP_STAIRMAP**       stairmap,           /**< pointer to store the created stair map */
   int                   upperbound,         /**< upper bound of the stairmap */
   int                   ntimepoints         /**< minimum size to ensure */
   );

/** frees given stair map */
extern
void SCIPstairmapFree(
   SCIP_STAIRMAP**       stairmap            /**< pointer to the stair map */
   );

/** resizes the stair map arrays */
extern
SCIP_RETCODE SCIPstairmapResize(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to resize */
   int                   ntimepoints         /**< minimum size to ensure */
   );

/** output of the given stair map */
extern
void SCIPstairmapPrint(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to output */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** insert a core into stair map; if core is non-empty the stair map will be updated otherwise nothing happens */
extern
void SCIPstairmapInsertCore(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer to store if the core does not fit due to capacity */
   );

/** subtracts the height from the stair map during core time */
extern
void SCIPstairmapDeleteCore(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   left,               /**< left side of the core  */
   int                   right,              /**< right side of the core */
   int                   height              /**< height of the core */
   );

/** returns the time point at the given position */   
extern
int SCIPstairmapGetTimepoint(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   pos                 /**< position */
   );

/** returns TRUE if the core  (given by its height and during) can be inserted at the given time point; otherwise FALSE */
extern
SCIP_Bool SCIPstairmapIsFeasibleStart(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   timepoint,          /**< time point to start */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   int*                  pos                 /**< pointer to store the earliest position where the core does not fit */
   );

/** return the earliest possible starting point within the time interval [lb,ub] for a given core (given by its height
 *  and duration)
 */
extern
int SCIPstairmapGetEarliestFeasibleStart(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   lb,                 /**< earliest starting time of the given core */
   int                   ub,                 /**< latest starting time of the given core */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   );

/** return the latest possible starting point within the time interval [lb,ub] for a given core (given by its height and
 *  duration)
 */
extern
int SCIPstairmapGetLatestFeasibleStart(
   SCIP_STAIRMAP*        stairmap,           /**< stair map to use */
   int                   lb,                 /**< earliest possible start point */
   int                   ub,                 /**< latest possible start point */
   int                   duration,           /**< duration of the core */
   int                   height,             /**< height of the core */
   SCIP_Bool*            infeasible          /**< pointer store if the core cannot be inserted */
   );

/*
 * Directed graph
 */

/** creates directed graph structure */
extern
SCIP_RETCODE SCIPdigraphCreate(
   SCIP_DIGRAPH**        digraph,            /**< pointer to store the created directed graph */
   int                   nnodes              /**< number of nodes */
   );

/** sets the sizes of the adjacency lists for the nodes in a directed graph and allocates memory for the lists */
extern
SCIP_RETCODE SCIPdigraphSetSizes(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int*                  sizes               /**< sizes of the adjacency lists */
   );

/** frees given directed graph structure */
extern
void SCIPdigraphFree(
   SCIP_DIGRAPH**        digraph             /**< pointer to the directed graph */
   );

/** add (directed) edge to the directed graph structure
 *  @note: if the edge is already contained, it is added a second time
 */
extern
SCIP_RETCODE SCIPdigraphAddEdge(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the edge */
   int                   endnode             /**< start node of the edge */
   );

/** add (directed) edge to the directed graph structure, if it is not contained, yet */
extern
SCIP_RETCODE SCIPdigraphAddEdgeSafe(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   startnode,          /**< start node of the edge */
   int                   endnode             /**< start node of the edge */
   );

/** returns the number of edges originating at the given node */
extern
int SCIPdigraphGetNOutEdges(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the number of outgoing edges is returned */
   );

/** returns the array of edges originating at the given node; this array must not be changed from outside */
extern
int* SCIPdigraphGetOutEdges(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   node                /**< node for which the array of outgoing edges is returned */
   );

/** Compute undirected connected components on the given graph.
 *
 *  @note For each edge, its reverse is added, so the graph does not need
 *        to be the directed representation of an undirected graph.
 */
extern
SCIP_RETCODE SCIPdigraphComputeUndirectedComponents(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   minsize,            /**< all components with less nodes are ignored */
   int*                  components,         /**< array with as many slots as there are nodes in the directed graph
                                              *   to store for each node the component to which it belongs
                                              *   (components are numbered 0 to ncomponents - 1); or NULL, if components
                                              *   are accessed one-by-one using SCIPdigraphGetComponent() */
   int*                  ncomponents         /**< pointer to store the number of components; or NULL, if the
                                              *   number of components is accessed by SCIPdigraphGetNComponents() */
   );

/** Performes an (almost) topological sort on the undirected components of the directed graph.
 *  The undirected components should be computed before using SCIPdigraphComputeUndirectedComponents().
 *
 *  Note, that in general a topological sort is not unique.
 *  Note, that there might be directed cycles, that are randomly broken,
 *  which is the reason for having only almost topologically sorted arrays.
 */
extern
SCIP_RETCODE SCIPdigraphTopoSortComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** returns the number of previously computed undirected components for the given directed graph */
extern
int SCIPdigraphGetNComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/** Returns the previously computed undirected components of the given number for the given directed graph.
 *  If the components were sorted using SCIPdigraphTopoSortComponents(), the component is (almost) topologically sorted.
 */
extern
void SCIPdigraphGetComponent(
   SCIP_DIGRAPH*         digraph,            /**< directed graph */
   int                   compidx,            /**< number of the component to return */
   int**                 nodes,              /**< pointer to store the nodes in the component */
   int*                  nnodes              /**< pointer to store the number of nodes in the component */
   );

/** frees the component information for the given directed graph */
extern
void SCIPdigraphFreeComponents(
   SCIP_DIGRAPH*         digraph             /**< directed graph */
   );

/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
extern
SCIP_Real SCIPcalcMachineEpsilon(
   void
   );

/** calculates the greatest common divisor of the two given values */
extern
SCIP_Longint SCIPcalcGreComDiv(
   SCIP_Longint          val1,               /**< first value of greatest common devisor calculation */
   SCIP_Longint          val2                /**< second value of greatest common devisor calculation */
   );

/** calculates the smallest common multiple of the two given values */
extern
SCIP_Longint SCIPcalcSmaComMul(
   SCIP_Longint          val1,               /**< first value of smallest common multiple calculation */
   SCIP_Longint          val2                /**< second value of smallest common multiple calculation */
   );

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
extern
SCIP_Bool SCIPrealToRational(
   SCIP_Real             val,                /**< real value r to convert into rational number */
   SCIP_Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   SCIP_Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** tries to find a value, such that all given values, if scaled with this value become integral */
extern
SCIP_RETCODE SCIPcalcIntegralScalar(
   SCIP_Real*            vals,               /**< values to scale */
   int                   nvals,              /**< number of values to scale */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   );

/** given a (usually very small) interval, tries to find a rational number with simple denominator (i.e. a small
 *  number, probably multiplied with powers of 10) out of this interval; returns TRUE iff a valid rational
 *  number inside the interval was found
 */
extern
SCIP_Bool SCIPfindSimpleRational(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed for resulting rational number */
   SCIP_Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   SCIP_Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** given a (usually very small) interval, selects a value inside this interval; it is tried to select a rational number
 *  with simple denominator (i.e. a small number, probably multiplied with powers of 10);
 *  if no valid rational number inside the interval was found, selects the central value of the interval
 */
extern
SCIP_Real SCIPselectSimpleValue(
   SCIP_Real             lb,                 /**< lower bound of the interval */
   SCIP_Real             ub,                 /**< upper bound of the interval */
   SCIP_Longint          maxdnom             /**< maximal denominator allowed for resulting rational number */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
extern
SCIP_Real SCIPrelDiff(
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPrelDiff(val1, val2)         ( ((val1)-(val2))/(MAX3(1.0,REALABS(val1),REALABS(val2))) )

#endif




/*
 * Random Numbers
 */

/** returns a random integer between minrandval and maxrandval */
extern
int SCIPgetRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );

/** returns a random real between minrandval and maxrandval */
extern
SCIP_Real SCIPgetRandomReal(
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   );




/*
 * Additional math functions
 */

/** calculates a binomial coefficient n over m, choose m elements out of n, maximal value will be 33 over 16 (because
 *  the n=33 is the last line in the Pascal's triangle where each entry fits in a 4 byte value), an error occurs due to
 *  big numbers or an negative value m (and m < n) and -1 will be returned
 */
extern
SCIP_Longint SCIPcalcBinomCoef(
   int                   n,                  /**< number of different elements */
   int                   m                   /**< number to choose out of the above */
   );



/*
 * Permutations / Shuffling
 */

/** swaps the addresses of two pointers */
extern
void SCIPswapPointers(
   void**                pointer1,           /**< first pointer */
   void**                pointer2            /**< second pointer */
   );

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
extern
void SCIPpermuteArray(
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end,                /**< last index that should be subject to shuffling (array size for whole array) */
   unsigned int*         randseed            /**< pointer to seed value for the random generator */
   ); 

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
extern
SCIP_RETCODE SCIPgetRandomSubset(
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems,          /**< number of elements that should be drawn and stored */
   unsigned int          randseed            /**< seed value for random generator */
   );


/*
 * Strings
 */

/** copies characters from 'src' to 'dest', copying is stopped when either the 'stop' character is reached or after
 *  'cnt' characters have been copied, whichever comes first.
 *
 *  @note undefined behaviuor on overlapping arrays
 */
extern
int SCIPmemccpy(
   char*                 dest,               /**< destination pointer to copy to */
   const char*           src,                /**< source pointer to copy to */
   char                  stop,               /**< character when found stop copying */
   unsigned int          cnt                 /**< maximal number of characters to copy too */
   );

/** prints an error message containing of the given string followed by a string describing the current system error;
 *  prefers to use the strerror_r method, which is threadsafe; on systems where this method does not exist,
 *  NO_STRERROR_R should be defined (see INSTALL), in this case, srerror is used which is not guaranteed to be
 *  threadsafe (on SUN-systems, it actually is) 
 */
extern
void SCIPprintSysError(
   const char*           message             /**< first part of the error message, e.g. the filename */
   );

/** extracts tokens from strings - wrapper method for strtok_r() */
extern
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
   );

/** translates the given string into a string where symbols ", ', and spaces are escaped with a \ prefix */
extern
void SCIPescapeString(
   char*                 t,                  /**< target buffer to store escaped string */
   int                   bufsize,            /**< size of buffer t */
   const char*           s                   /**< string to transform into escaped string */
   );

/** safe version of snprintf */
extern
int SCIPsnprintf(
   char*                 t,                  /**< target string */
   int                   len,                /**< length of the string to copy */
   const char*           s,                  /**< source string */
   ...                                       /**< further parameters */
   );

/** extract the next token as a double value if it is one; in case a value is parsed the endptr is set to NULL */
extern 
SCIP_Bool SCIPstrToRealValue(
   const char*           str,                /**< string to search */
   SCIP_Real*            value,              /**< pointer to store the parsed value */
   char**                endptr              /**< pointer to store the final string position if successfully parsed
                                              *   otherwise NULL */
   );

/** copies the string between a start and end character */
extern
void SCIPstrCopySection(
   const char*           str,                /**< string to search */
   char                  startchar,          /**< character which defines the beginning */
   char                  endchar,            /**< character which defines the ending */
   char*                 token,              /**< string to store the copy */
   int                   size,               /**< size of the token char array */
   char**                endptr              /**< pointer to store the final string position if successfully parsed
                                              *   otherwise NULL */
   );

/*
 * File methods
 */

/** returns, whether the given file exists */
extern
SCIP_Bool SCIPfileExists(
   const char*           filename            /**< file name */
   );

/** splits filename into path, name, and extension */
extern
void SCIPsplitFilename(
   char*                 filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**                path,               /**< pointer to store path, or NULL if not needed */
   char**                name,               /**< pointer to store name, or NULL if not needed */
   char**                extension,          /**< pointer to store extension, or NULL if not needed */
   char**                compression         /**< pointer to store compression extension, or NULL if not needed */
   );

#ifdef __cplusplus
}
#endif

#endif
