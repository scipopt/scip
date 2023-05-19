/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pattern.h
 * @brief  pattern data for ringpacking problem
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PATTERN__
#define __SCIP_PATTERN__

#include "scip/scip.h"

/*
 * data structures
 */

enum SCIP_Packable
{
   SCIP_PACKABLE_NO      = 0,                  /**< pattern is definitely packable */
   SCIP_PACKABLE_YES     = 1,                  /**< pattern is definitely not packable */
   SCIP_PACKABLE_UNKNOWN = 2                   /**< it is unknown whether pattern is packable */
};
typedef enum SCIP_Packable SCIP_PACKABLE;

enum SCIP_Patterntype
{
   SCIP_PATTERNTYPE_CIRCULAR    = 0,           /**< circular pattern */
   SCIP_PATTERNTYPE_RECTANGULAR = 1            /**< rectangular pattern */
};
typedef enum SCIP_Patterntype SCIP_PATTERNTYPE;

struct SCIP_Pattern
{
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SCIP_PATTERNTYPE      patterntype;        /**< pattern type */
   SCIP_PACKABLE         packable;           /**< packable status */
   SCIP_Real*            xs;                 /**< array containing the x-coordinate of each element */
   SCIP_Real*            ys;                 /**< array containing the y-coordinate of each element */
   int*                  types;              /**< array storing the type of each element */
   int                   size;               /**< size of types, xs, and ys arrays */
   int                   nelems;             /**< number of elements stored */
   int                   nlocks;             /**< number of locks */
   int                   type;               /**< type of the boundary circle */
};
typedef struct SCIP_Pattern SCIP_PATTERN;

/** creates an empty circular pattern */
SCIP_RETCODE SCIPpatternCreateCircular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern,            /**< pointer to store pattern */
   int                   type                /**< circle type (not needed for rectangular patterns) */
   );

/** creates an empty rectangular pattern */
SCIP_RETCODE SCIPpatternCreateRectangular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   );

/** captures a pattern */
void SCIPpatternCapture(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/* frees a pattern */
void SCIPpatternRelease(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to free pattern */
   );

/** copies a pattern */
SCIP_RETCODE SCIPpatternCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN*         pattern,            /**< pattern to copy */
   SCIP_PATTERN**        copy                /**< pointer to store the copy */
   );

/** adds an element of a given type to a pattern; packable status does not change */
SCIP_RETCODE SCIPpatternAddElement(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type,               /**< element of a given type */
   SCIP_Real             x,                  /**< x-coordinate (SCIP_INVALID: unknown) */
   SCIP_Real             y                   /**< y-coordinate (SCIP_INVALID: unknown) */
   );

/** removes the last k elements */
void SCIPpatternRemoveLastElements(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   k                   /**< number of elements to remove */
   );

/** returns the total number of elements of a given type in the pattern */
int SCIPpatternGetNElemens(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/** returns the type of the i-th element */
int SCIPpatternGetElementType(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   i                   /**< i-th element */
   );

/** returns the total number of elements of a given type */
int SCIPpatternCountElements(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< type */
   );

/** returns the x-coordinate of an element */
SCIP_Real SCIPpatternGetElementPosX(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem                /**< index of the element */
   );

/** returns the y-coordinate of an element */
SCIP_Real SCIPpatternGetElementPosY(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem                /**< index of the element */
   );

/** sets the (x,y) position of an element */
void SCIPpatternSetElementPos(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem,               /**< index of the element */
   SCIP_Real             x,                  /**< x-coordinate */
   SCIP_Real             y                   /**< y-coordinate */
   );

/** returns the type of a pattern */
SCIP_PATTERNTYPE SCIPpatternGetPatternType(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/** returns the type of the boundary circle
 *
 * @note this function can only be called for circular patterns
 */
int SCIPpatternGetCircleType(
   SCIP_PATTERN *pattern             /**< pattern */
);

/** sets the type of the boundary circle
 *
 * @note this function can only be called for circular patterns
 */
void SCIPpatternSetType(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< type */
   );

/** returns the packable status of a pattern */
SCIP_PACKABLE SCIPpatternGetPackableStatus(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/** sets the packable status of a pattern */
void SCIPpatternSetPackableStatus(
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_PACKABLE         packable            /**< packable status */
   );

#endif /* __SCIP_PATTERN__ */
