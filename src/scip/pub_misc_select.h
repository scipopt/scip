/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_misc_select.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for selecting (weighted) k-medians
 * @author Gregor Hendel
 *
 * This file contains headers for selecting (weighted) k-medians
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_MISC_SELECT_H__
#define __SCIP_PUB_MISC_SELECT_H__

#include "scip/def.h"
#include "type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Selection and weighted selection algorithms
 */

/**@defgroup SelectionAlgorithms Algorithms for (Weighted) Median Selection
 * @ingroup MiscellaneousMethods
 * @brief public methods for the selection of (weighted) k-median.
 *
 * The methods in this group perform a selection of the (weighted)  \f$ k \f$-median from an unsorted array of elements.
 * The necessary element swaps are performed in-place on the array of keys.
 * The necessary permutations are also performed on up to six associated arrays.
 *
 * For methods that perform complete in place sorting, see \ref SortingAlgorithms.
 *
 * For an array a containing n elements \f$ a[0], ..., a[n-1] \f$ and an integer \f$ 0 \leq k \leq n - 1 \f$ , we call an element
 *  \f$ a[i] \f$   \f$ k \f$-median if
 * there exists a permutation  \f$ \pi \f$  of the array indices such that  \f$ \pi(i) = k  \f$
 * and  \f$ a[\pi^{-1}(j)] \leq a[i]  \f$
 * for  \f$ j = 0, \dots, k-1  \f$ and  \f$ a[\pi^{-1}(j)] > a[i]  \f$ for  \f$ j = k + 1,\dots,n - 1 \f$.
 * The  \f$ k \f$-median is hence an element that would appear at position  \f$ k  \f$ after sorting the input array.
 * Note that there may exist several  \f$ k \f$-medians if the array elements are not unique, only its key value  \f$ a[i] \f$.
 *
 * In order to determine the  \f$ k \f$-median, the algorithm selects a pivot element and determines the array position for
 * this pivot like quicksort. In contrast to quicksort, however, one recursion can be saved during the selection process.
 * After a single iteration that placed the pivot at position  \f$ p \f$ , the algorithm either terminates if  \f$ p = k \f$,
 * or it continues in the left half of the array if  \f$ p > k \f$, or in the right half of the array if  \f$ p < k \f$.
 *
 * After the algorithm terminates, the  \f$ k \f$-median can be accessed by accessing the array element at position  \f$ k \f$.
 *
 * A weighted median denotes the generalization of the  \f$ k \f$-median to arbitrary, nonnegative associated
 * weights \f$ w[0], \dots, w[n-1] \in \mathbb{R}\f$ and a capacity  \f$ 0 \leq C \in \mathbb{R} \f$. An element  \f$ a[i] \f$
 * is called weighted median if there exists a permutation that satisfies the same weak sorting as above and in addition
 * \f$ W:= \sum_{j = 0}^{k - 1}w[\pi^{-1}(j)] < C\f$, but  \f$ W + w[i] \geq C\f$. In other words, the weighted median
 * is the first element in the weak sorting such that its weight together with the sum of all preceding item weights
 * reach or exceed the given capacity \f$ C \f$. If all weights are equal to \f$ 1 \f$ and the capacity is  \f$ C = k + 0.5\f$,
 * the weighted median becomes the  \f$ k \f$-median.
 *
 * @{
 */

/** partial sort an index array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an index array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of an array of pointers in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of an array of pointers in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrRealRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrRealRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Reals in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Reals in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort array of ints in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort array of ints in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedInt(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Longints in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Longints in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an index array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an index array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of an array of pointers in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of an array of pointers in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrPtrLongIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Reals in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Reals in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );

/** partial sort of three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order around the \p k-th element */
SCIP_EXPORT extern
void SCIPselectDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );

/** partial sort of three joint arrays of Reals/Reals/Pointer/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealPtrPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );

/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownRealRealRealBoolBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort array of ints in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort array of ints in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownInt(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort an array of Longints in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Longints in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrRealRealIntBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownPtrIntIntBoolBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the \p k-th element,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity,
 *  see \ref SelectionAlgorithms for more information.
 */
SCIP_EXPORT extern
void SCIPselectWeightedDownIntPtrIntIntBoolBool(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< thrid int array to be permuted in the same way */
   SCIP_Bool*            boolarray1,         /**< first SCIP_Bool array to be permuted in the same way */
   SCIP_Bool*            boolarray2,         /**< second SCIP_Bool array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
