/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_misc.h,v 1.24 2006/05/22 15:51:53 bzfheinz Exp $"

/**@file   pub_misc.h
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

/** creates a hash table */
extern
SCIP_RETCODE SCIPhashtableCreate(
   SCIP_HASHTABLE**      hashtable,          /**< pointer to store the created hash table */
   BMS_BLKMEM*           blkmem,             /**< block memory used to store hash table entries */
   int                   tablesize,          /**< size of the hash table */
   SCIP_DECL_HASHGETKEY((*hashgetkey)),      /**< gets the key of the given element */
   SCIP_DECL_HASHKEYEQ ((*hashkeyeq)),       /**< returns TRUE iff both keys are equal */
   SCIP_DECL_HASHKEYVAL((*hashkeyval))       /**< returns the hash value of the key */
   );

/** frees the hash table */
extern
void SCIPhashtableFree(
   SCIP_HASHTABLE**      hashtable           /**< pointer to the hash table */
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

/** prints statistics about hash table usage */
extern
void SCIPhashtablePrintStatistics(
   SCIP_HASHTABLE*       hashtable           /**< hash table */
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
   SCIP_HASHMAP*         hashmap             /**< hash map */
   );




/*
 * Sorting algorithms
 */

/** bubble sort an indexed element set, resulting in a permutation index array */
extern
void SCIPbsort(
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   len,                /**< number of elements to be sorted (valid index range) */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   int*                  indarray            /**< pointer to store the sorted index array */
   );

/** bubble sort of an array of pointers */
extern
void SCIPbsortPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   int                   len,                /**< length of array */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of two joint arrays of pointers/Reals, sorted by first array */
extern
void SCIPbsortPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   len,                /**< length of arrays */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of two joint arrays of pointers/ints, sorted by first array */
extern
void SCIPbsortPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len,                /**< length of arrays */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of three joint arrays of pointers/ints/ints, sorted by first array */
extern
void SCIPbsortPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len,                /**< length of arrays */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of three joint arrays of pointers/Reals/ints, sorted by first */
extern
void SCIPbsortPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   len,                /**< length of arrays */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of four joint arrays of pointers/Reals/ints/ints, sorted by first */
extern
void SCIPbsortPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   len,                /**< length of arrays */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp))         /**< data element comparator */
   );

/** bubble sort of two joint arrays of Reals/pointers, sorted s.t. the SCIP_Real array is in non-decreasing order */
extern
void SCIPbsortRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be sorted */
   int                   len                 /**< length of arrays */
   );




/*
 * Numerical methods
 */

/** returns the machine epsilon: the smallest number eps > 0, for which 1.0 + eps > 1.0 */
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
 * Strings
 */

/** extracts tokens from strings - wrapper method for strtok_r() */
extern
char* SCIPstrtok(
   char*                 s,                  /**< string to parse */
   const char*           delim,              /**< delimiters for parsing */
   char**                ptrptr              /**< pointer to working char pointer - must stay the same while parsing */
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

#endif
