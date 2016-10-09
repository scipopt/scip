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

/**@file   random.h
 * @brief  data structures for random number generator
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RANDOM_H__
#define __SCIP_RANDOM_H__

#include "scip/def.h"
#include "scip/struct_random.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates and initialzes a random number generator */
extern
SCIP_RETCODE SCIPrandomCreate(
   SCIP_RANDGEN**        randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   unsigned int          initialseed         /**< initial random seed (> 0) */
   );

/** creates and initialzes a random number generator */
extern
void SCIPrandomFree(
   SCIP_RANDGEN**        randnumgen,         /**< random number generator */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** returns a random integer between minrandval and maxrandval */
extern
int SCIPrandomGetInt(
   SCIP_RANDGEN*         randgen,
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval          /**< maximal value to return */
   );

/** returns a random real between minrandval and maxrandval */
extern
SCIP_Real SCIPrandomGetReal(
   SCIP_RANDGEN*         randgen,
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval          /**< maximal value to return */
   );


/*
 * Permutations / Shuffling
 */

/**@defgroup PermutationsShuffling Permutations Shuffling
 *
 *@{
 */

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm */
extern
void SCIPrandomPermuteIntArray(
   SCIP_RANDGEN*         randgen,            /**< random number generator data */
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end                 /**< last index that should be subject to shuffling (array size for whole
                                               *   array)
                                               */
   );

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
extern
void SCIPrandomPermuteArray(
   SCIP_RANDGEN*         randgen,            /**< random number generator data */
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end                 /**< last index that should be subject to shuffling (array size for whole
                                              *   array)
                                              */
   );

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
extern
SCIP_RETCODE SCIPrandomGetSubset(
   SCIP_RANDGEN*         randgen,            /**< random number generator data */
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems           /**< number of elements that should be drawn and stored */
   );

/**@} */


#ifdef __cplusplus
}
#endif

#endif
