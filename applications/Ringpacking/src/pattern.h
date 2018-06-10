/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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
extern
SCIP_RETCODE SCIPpatternCreateCircular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern,            /**< pointer to store pattern */
   int                   type                /**< circle type (not needed for rectangular patterns) */
   );

/** creates an empty rectangular pattern */
extern
SCIP_RETCODE SCIPpatternCreateRectangular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   );

/** captures a pattern */
extern
void SCIPpatternCapture(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/* frees a pattern */
extern
void SCIPpatternRelease(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to free pattern */
   );

/** copies a pattern */
extern
SCIP_RETCODE SCIPpatternCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN*         pattern,            /**< pattern to copy */
   SCIP_PATTERN**        copy                /**< pointer to store the copy */
   );

/** adds an element of a given type to a pattern; packable status does not change */
extern
SCIP_RETCODE SCIPpatternAddElement(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type,               /**< element of a given type */
   SCIP_Real             x,                  /**< x-coordinate (SCIP_INVALID: unknown) */
   SCIP_Real             y                   /**< y-coordinate (SCIP_INVALID: unknown) */
   );

/** removes the last k elements */
extern
void SCIPpatternRemoveLastElements(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   k                   /**< number of elements to remove */
   );

/** returns the total number of elements of a given type in the pattern */
extern
int SCIPpatternGetNElemens(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/** returns the type of the i-th element */
extern
int SCIPpatternGetElementType(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   i                   /**< i-th element */
   );

/** returns the total number of elements of a given type */
extern
int SCIPpatternCountElements(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< type */
   );

/** returns the x-coordinate of an element */
extern
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
extern
void SCIPpatternSetElementPos(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem,               /**< index of the element */
   SCIP_Real             x,                  /**< x-coordinate */
   SCIP_Real             y                   /**< y-coordinate */
   );

/** returns the type of a pattern */
extern
SCIP_PATTERNTYPE SCIPpatternGetPatternType(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/** returns the type of the boundary circle
 *
 * @note this function can only be called for circular patterns
 */
extern
int SCIPpatternGetCircleType(
   SCIP_PATTERN *pattern             /**< pattern */
);

/** sets the type of the boundary circle
 *
 * @note this function can only be called for circular patterns
 */
extern
void SCIPpatternSetType(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< type */
   );

/** returns the packable status of a pattern */
extern
SCIP_PACKABLE SCIPpatternGetPackableStatus(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

/** sets the packable status of a pattern */
extern
void SCIPpatternSetPackableStatus(
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_PACKABLE         packable            /**< packable status */
   );

#endif /* __SCIP_PATTERN__ */
