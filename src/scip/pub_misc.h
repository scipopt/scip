/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_misc.h,v 1.10 2005/01/18 09:26:52 bzfpfend Exp $"

/**@file   pub_misc.h
 * @brief  public miscellaneous methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_MISC_H__
#define __PUB_MISC_H__



#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_misc.h"



/*
 * Priority Queue
 */

/** creates priority queue */
extern
RETCODE SCIPpqueueCreate(
   PQUEUE**         pqueue,             /**< pointer to a priority queue */
   int              initsize,           /**< initial number of available element slots */
   Real             sizefac,            /**< memory growing factor applied, if more element slots are needed */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** frees priority queue, but not the data elements themselves */
extern
void SCIPpqueueFree(
   PQUEUE**         pqueue              /**< pointer to a priority queue */
   );

/** clears the priority queue, but doesn't free the data elements themselves */
extern
void SCIPpqueueClear(
   PQUEUE*          pqueue              /**< priority queue */
   );

/** inserts element into priority queue */
extern
RETCODE SCIPpqueueInsert(
   PQUEUE*          pqueue,             /**< priority queue */
   void*            elem                /**< element to be inserted */
   );

/** removes and returns best element from the priority queue */
extern
void* SCIPpqueueRemove(
   PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the best element of the queue without removing it */
extern
void* SCIPpqueueFirst(
   PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the number of elements in the queue */
extern
int SCIPpqueueNElems(
   PQUEUE*          pqueue              /**< priority queue */
   );

/** returns the elements of the queue; changing the returned array may destroy the queue's ordering! */
extern
void** SCIPpqueueElems(
   PQUEUE*          pqueue              /**< priority queue */
   );




/*
 * Hash Table
 */

/** creates a hash table */
extern
RETCODE SCIPhashtableCreate(
   HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   MEMHDR*          memhdr,             /**< block memory used to store hash table entries */
   int              tablesize,          /**< size of the hash table */
   DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   DECL_HASHKEYVAL((*hashkeyval))       /**< returns the hash value of the key */
   );

/** frees the hash table */
extern
void SCIPhashtableFree(
   HASHTABLE**      hashtable           /**< pointer to the hash table */
   );

/** inserts element in hash table (multiple inserts of same element possible) */
extern
RETCODE SCIPhashtableInsert(
   HASHTABLE*       hashtable,          /**< hash table */
   void*            element             /**< element to append to the list */
   );

/** inserts element in hash table (multiple insertion of same element is checked and results in an error) */
extern
RETCODE SCIPhashtableSafeInsert(
   HASHTABLE*       hashtable,          /**< hash table */
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
 * Hash Map
 */

/** creates a hash map mapping pointers to pointers */
extern
RETCODE SCIPhashmapCreate(
   HASHMAP**        hashmap,            /**< pointer to store the created hash map */
   MEMHDR*          memhdr,             /**< block memory used to store hash map entries */
   int              mapsize             /**< size of the hash map */
   );

/** frees the hash map */
extern
void SCIPhashmapFree(
   HASHMAP**        hashmap             /**< pointer to the hash map */
   );

/** inserts new origin->image pair in hash map (must not be called for already existing origins!) */
extern
RETCODE SCIPhashmapInsert(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin,             /**< origin to set image for */
   void*            image               /**< new image for origin */
   );

/** retrieves image of given origin from the hash map, or NULL if no image exists */
extern
void* SCIPhashmapGetImage(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin              /**< origin to retrieve image for */
   );

/** sets image for given origin in the hash map, either by modifying existing origin->image pair or by appending a
 *  new origin->image pair
 */
extern
RETCODE SCIPhashmapSetImage(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin,             /**< origin to set image for */
   void*            image               /**< new image for origin */
   );

/** removes existing origin->image pair from the hash map */
extern
RETCODE SCIPhashmapRemove(
   HASHMAP*         hashmap,            /**< hash map */
   void*            origin              /**< origin to remove from the list */
   );

/** prints statistics about hash map usage */
extern
void SCIPhashmapPrintStatistics(
   HASHMAP*         hashmap             /**< hash map */
   );




/*
 * Sorting algorithms
 */

/** bubble sort an indexed element set, resulting in a permutation index array */
extern
void SCIPbsort(
   void*            dataptr,            /**< pointer to data field that is given to the external compare method */
   int              len,                /**< number of elements to be sorted (valid index range) */
   DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   int*             indarray            /**< pointer to store the sorted index array */
   );

/** bubble sort of an array of pointers */
extern
void SCIPbsortPtr(
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len,                /**< length of array */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of two joint arrays of pointers/Reals, sorted by first array */
extern
void SCIPbsortPtrReal(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            realarray,          /**< Real array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of two joint arrays of pointers/ints, sorted by first array */
extern
void SCIPbsortPtrInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of three joint arrays of pointers/ints/ints, sorted by first array */
extern
void SCIPbsortPtrIntInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   int*             intarray1,          /**< first int array to be permuted in the same way */
   int*             intarray2,          /**< second int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of three joint arrays of pointers/Reals/ints, sorted by first */
extern
void SCIPbsortPtrRealInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            realarray,          /**< Real array to be permuted in the same way */
   int*             intarray,           /**< int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of four joint arrays of pointers/Reals/ints/ints, sorted by first */
extern
void SCIPbsortPtrRealIntInt(
   void**           ptrarray,           /**< pointer array to be sorted */
   Real*            realarray,          /**< Real array to be permuted in the same way */
   int*             intarray1,          /**< first int array to be permuted in the same way */
   int*             intarray2,          /**< second int array to be permuted in the same way */
   int              len,                /**< length of arrays */
   DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of two joint arrays of Reals/pointers, sorted s.t. the Real array is in non-decreasing order */
extern
void SCIPbsortRealPtr(
   Real*            realarray,          /**< Real array to be permuted in the same way */
   void**           ptrarray,           /**< pointer array to be sorted */
   int              len                 /**< length of arrays */
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
   Longint          val1,               /**< first value of smallest common multiple calculation */
   Longint          val2                /**< second value of smallest common multiple calculation */
   );

/** converts a real number into a (approximate) rational representation, and returns TRUE iff the conversion was
 *  successful
 */
extern
Bool SCIPrealToRational(
   Real             val,                /**< real value r to convert into rational number */
   Real             mindelta,           /**< minimal allowed difference r - q of real r and rational q = n/d */
   Real             maxdelta,           /**< maximal allowed difference r - q of real r and rational q = n/d */
   Longint          maxdnom,            /**< maximal denominator allowed */
   Longint*         nominator,          /**< pointer to store the nominator n of the rational number */
   Longint*         denominator         /**< pointer to store the denominator d of the rational number */
   );

/** tries to find a value, such that all given values, if scaled with this value become integral */
extern
RETCODE SCIPcalcIntegralScalar(
   Real*            vals,               /**< values to scale */
   int              nvals,              /**< number of values to scale */
   Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   Real             maxscale,           /**< maximal allowed scalar */
   Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   Bool*            success             /**< stores whether returned value is valid */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
extern
Real SCIPrelDiff(
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPrelDiff(val1, val2)         ( ((val1)-(val2))/(MAX3(1.0,REALABS(val1),REALABS(val2))) )

#endif




/*
 * File methods
 */

/** returns, whether the given file exists */
extern
Bool SCIPfileExists(
   const char*      filename            /**< file name */
   );

/** splits filename into path, name, and extension */
extern
void SCIPsplitFilename(
   char*            filename,           /**< filename to split; is destroyed (but not freed) during process */
   char**           path,               /**< pointer to store path, or NULL if not needed */
   char**           name,               /**< pointer to store name, or NULL if not needed */
   char**           extension           /**< pointer to store extension, or NULL if not needed */
   );

#endif
