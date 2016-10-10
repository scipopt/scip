/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_misc_select.h
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

/**@defgroup SelectionAlgorithms Selection and weighted selection algorithms
 *
 * @{
 */

/** partial sort an index array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an index array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of an array of pointers in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of an array of pointers in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrRealReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/Bools, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
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


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort an array of Reals in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Reals in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedReal(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second  SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
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


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort array of ints in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort array of ints in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedInt(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Longints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntRealLong(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/pointers/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedIntPtrReal(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/pointers/ints/reals, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort an array of Longints in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Longints in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order around the \p k-th element */
EXTERN
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


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
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


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
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


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort an index array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an index array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownInd(
   int*                  indarray,           /**< pointer to the index array to be sorted */
   SCIP_DECL_SORTINDCOMP((*indcomp)),        /**< data element comparator */
   void*                 dataptr,            /**< pointer to data field that is given to the external compare method */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of an array of pointers in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of an array of pointers in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownPtr(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownPtrPtr(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownPtrReal(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownPtrInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of pointers/Bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownPtrBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrRealInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrRealBool(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/Reals/Bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtrInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtrReal(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of pointers/pointers/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtrIntInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrRealIntInt(
   void**                ptrarray,           /**< pointer array to be sorted */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointers/Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtrRealInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtrRealBool(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Reals/bools, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownPtrPtrLongInt(
   void**                ptrarray1,          /**< first pointer array to be sorted */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_DECL_SORTPTRCOMP((*ptrcomp)),        /**< data element comparator */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of pointer/pointer/Longs/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
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


/** partial sort of five joint arrays of pointer/pointer/Longs/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort an array of Reals in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Reals in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownReal(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealInt(
   SCIP_Real*            realarray,          /**< SCIP_Real  array to be sorted */
   int*                  intarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Bools/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealBoolPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Longs, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealIntLong(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/ints/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealIntPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   int*                  intarray,           /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealRealInt(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< integer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/Reals/Pointer, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealRealPtr(
   SCIP_Real*            realarray1,         /**< first SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealPtrPtrInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/pointers/pointers/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealPtrPtrIntInt(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/pointers/pointers/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealLongRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array  to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Longs/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealRealIntInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray1,          /**< int array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealRealRealInt(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealRealRealPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Reals/Reals/Reals/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of Reals/pointers, sorted by first array in non-decreasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownRealPtrPtr(
   SCIP_Real*            realarray,          /**< SCIP_Real array to be sorted */
   void**                ptrarray1,          /**< pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownRealRealRealBoolPtr(
   SCIP_Real*            realarray1,         /**< SCIP_Real array to be sorted */
   SCIP_Real*            realarray2,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray3,         /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Reals/Reals/Reals/Bools/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
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


/** partial sort of six joint arrays of Reals/Reals/Reals/Bools/Bools/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort array of ints in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownInt(
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort array of ints in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownInt(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntPtr(
   int*                  intarray,           /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of ints/reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntReal(
   int*                  intarray,           /**< int array to be sorted */
   SCIP_Real*            realarray,          /**< real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntIntInt(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int*                  intarray3,          /**< third int  array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/SCIP_Longint, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntIntLong(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three joint arrays of ints/ints/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntIntIntPtr(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/pointers, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntIntIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   int*                  intarray3,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/ints/ints/reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownIntPtrIntReal(
   int*                  intarray1,          /**< int array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray2,          /**< int array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of ints/pointers/ints/Reals, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort an array of Longints in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort an array of Longints in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownLong(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of two joint arrays of Long/pointer, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownLongPtr(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of three arrays of Long/pointer/ints, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
void SCIPselectWeightedDownLongPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   SCIP_Real*            weights,            /**< (optional), nonnegative weights array for weighted median, or NULL (all weights are equal to 1) */
   SCIP_Real             capacity,           /**< the maximum capacity that is exceeded by the median */
   int                   len,                /**< length of arrays */
   int*                  medianpos           /**< pointer to store the index of the weighted median, or NULL, if not needed */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtrRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four arrays of Long/pointer/Real/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtrRealRealBool(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray,           /**< pointer array to be permuted in the same way */
   SCIP_Real*            realarray,          /**< first SCIP_Real array to be permuted in the same way */
   SCIP_Real*            realarray2,         /**< second SCIP_Real array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five arrays of Long/pointer/Real/Real/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order around the \p k-th element */
EXTERN
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


/** partial sort of six arrays of Long/pointer/Real/Real/int/Bool, sorted by the first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtrPtrInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray,           /**< int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of four joint arrays of Long/pointer/pointer/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtrPtrIntInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   int*                  intarray1,          /**< first int array to be permuted in the same way */
   int*                  intarray2,          /**< second int array to be permuted in the same way */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/ints/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
void SCIPselectDownLongPtrPtrBoolInt(
   SCIP_Longint*         longarray,          /**< SCIP_Longint array to be sorted */
   void**                ptrarray1,          /**< first pointer array to be permuted in the same way */
   void**                ptrarray2,          /**< second pointer array to be permuted in the same way */
   SCIP_Bool*            boolarray,          /**< SCIP_Bool array to be permuted in the same way */
   int*                  intarray,           /**< int array to be sorted */
   int                   k,                  /**< the index of the desired element, must be between 0 (search for maximum/minimum) and len - 1 */
   int                   len                 /**< length of arrays */
   );


/** partial sort of five joint arrays of Long/pointer/pointer/Bool/ints, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
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


/** partial sort of five joint arrays of pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the \p k-th element */
EXTERN
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


/** partial sort of six joint arrays of ints/pointer/ints/ints/Bool/Bool, sorted by first array in non-increasing order around the weighted median w.r.t. \p weights and capacity */
EXTERN
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
