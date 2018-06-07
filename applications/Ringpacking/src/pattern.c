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

/**@file   pattern.c
 * @brief  pattern data for Ringpacking Problem
 * @author Benjamin Mueller
 *
 *
 * This file implements the handling of patterns. Each pattern has a <code>SCIP_PATTERNTYPE</code>, accessible by
 * <code>SCIPpatternGetPatternType()</code>, which indicates whether it is a circular or rectangular pattern.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pattern.h"
#include "probdata_rpa.h"

/*
 * local methods
 */

/** ensures that there is enough memory to store elements */
static
SCIP_RETCODE ensureElemSize(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   size                /**< required size */
   )
{
   assert(pattern != NULL);

   if( size > pattern->size )
   {
      int newsize = MAX(4, 2*size);
      assert(newsize > size);

      SCIP_ALLOC( BMSreallocBlockMemoryArray(pattern->blkmem, &pattern->types, pattern->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(pattern->blkmem, &pattern->xs, pattern->size, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(pattern->blkmem, &pattern->ys, pattern->size, newsize) );
      pattern->size = newsize;
   }

   return SCIP_OKAY;
}

/** auxiliary function to create a pattern */
static
SCIP_RETCODE createPattern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern,            /**< pointer to store pattern */
   SCIP_PATTERNTYPE      patterntype,        /**< pattern type */
   int                   type                /**< circle type (not needed for rectangular patterns) */
   )
{
   assert(scip != NULL);
   assert(pattern != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, pattern) );
   BMSclearMemory(*pattern);

   (*pattern)->blkmem = SCIPblkmem(scip);
   (*pattern)->type = type;
   SCIPpatternCapture(*pattern);

   (*pattern)->packable = SCIP_PACKABLE_UNKNOWN;
   (*pattern)->patterntype = patterntype;
   (*pattern)->type = type;

   return SCIP_OKAY;
}

/*
 * interface methods
 */

/** creates an empty circular pattern */
SCIP_RETCODE SCIPpatternCreateCircular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern,            /**< pointer to store pattern */
   int                   type                /**< circle type (not needed for rectangular patterns) */
   )
{
   return createPattern(scip, pattern, SCIP_PATTERNTYPE_CIRCULAR, type);
}

/** creates an empty rectangular pattern */
SCIP_RETCODE SCIPpatternCreateRectangular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   )
{
   return createPattern(scip, pattern, SCIP_PATTERNTYPE_RECTANGULAR, -1);
}

/** captures a pattern */
void SCIPpatternCapture(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);
   assert(pattern->nlocks >= 0);
   ++(pattern->nlocks);
}

/* frees a pattern */
void SCIPpatternRelease(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to free pattern */
   )
{
   assert(scip != NULL);
   assert(pattern != NULL);
   assert(*pattern != NULL);
   assert((*pattern)->nlocks > 0);

   /* reduce locks */
   --((*pattern)->nlocks);

   /* free memory if the pattern is not used anymore */
   if( (*pattern)->nlocks == 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*pattern)->ys, (*pattern)->size);
      SCIPfreeBlockMemoryArrayNull(scip, &(*pattern)->xs, (*pattern)->size);
      SCIPfreeBlockMemoryArrayNull(scip, &(*pattern)->types, (*pattern)->size);
      SCIPfreeBlockMemory(scip, pattern);
   }
   else
      *pattern = NULL;
}

/** copies a pattern */
SCIP_RETCODE SCIPpatternCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN*         pattern,            /**< pattern to copy */
   SCIP_PATTERN**        copy                /**< pointer to store the copy */
   )
{
   int i;

   assert(pattern != NULL);
   assert(copy != NULL);

   SCIP_CALL( createPattern(scip, copy, pattern->patterntype, pattern->type) );
   assert(*copy != NULL);

   /* ensure that we can store all elements */
   SCIP_CALL( ensureElemSize(*copy, pattern->nelems) );

   /* add elements */
   for( i = 0; i < pattern->nelems; ++i )
   {
      SCIP_CALL( SCIPpatternAddElement(*copy, pattern->types[i], pattern->xs[i], pattern->ys[i]) );
   }

   /* copy packable status */
   (*copy)->packable = pattern->packable;

   return SCIP_OKAY;
}

/** adds an element of a given type to a pattern; packable status does not change */
SCIP_RETCODE SCIPpatternAddElement(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type,               /**< element of a given type */
   SCIP_Real             x,                  /**< x-coordinate (SCIP_INVALID: unknown) */
   SCIP_Real             y                   /**< y-coordinate (SCIP_INVALID: unknown) */
   )
{
   assert(pattern != NULL);
   assert(type >= 0);

   SCIP_CALL( ensureElemSize(pattern, pattern->nelems + 1) );
   pattern->types[pattern->nelems] = type;
   pattern->xs[pattern->nelems] = x;
   pattern->ys[pattern->nelems] = y;

   ++(pattern->nelems);

   return SCIP_OKAY;
}

/** removes the last k elements */
void SCIPpatternRemoveLastElements(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   k                   /**< number of elements to remove */
   )
{
   assert(pattern != NULL);
   assert(pattern->nelems >= k);

   pattern->nelems -= k;
}

/** returns the total number of elements */
int SCIPpatternGetNElemens(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);

   return pattern->nelems;
}

/** returns the type of the i-th element */
int SCIPpatternGetElementType(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   i                   /**< index */
   )
{
   assert(pattern != NULL);
   assert(i >= 0 && i < pattern->nelems);

   return pattern->types[i];
}

/** returns the total number of elements of a given type */
int SCIPpatternCountElements(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< type */
   )
{
   int counter = 0;
   int i;

   assert(pattern != NULL);

   for( i = 0; i < pattern->nelems; ++i )
   {
      if( pattern->types[i] == type )
         ++(counter);
   }

   return counter;
}

/** returns the x-coordinate of an element */
SCIP_Real SCIPpatternGetElementPosX(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem                /**< index of the element */
   )
{
   assert(pattern != NULL);
   assert(elem >= 0 && elem < pattern->nelems);

   return pattern->xs[elem];
}

/** returns the y-coordinate of an element */
SCIP_Real SCIPpatternGetElementPosY(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem                /**< index of the element */
   )
{
   assert(pattern != NULL);
   assert(elem >= 0 && elem < pattern->nelems);

   return pattern->ys[elem];
}

/** sets the (x,y) position of an element */
void SCIPpatternSetElementPos(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   elem,               /**< index of the element */
   SCIP_Real             x,                  /**< x-coordinate */
   SCIP_Real             y                   /**< y-coordinate */
   )
{
   assert(pattern != NULL);
   assert(elem >= 0 && elem < pattern->nelems);

   pattern->xs[elem] = x;
   pattern->ys[elem] = y;
}

/** returns the type of a pattern */
SCIP_PATTERNTYPE SCIPpatternGetPatternType(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);

   return pattern->patterntype;
}

/** returns the type of the boundary circle
 *
 * @note this function can only be called for circular patterns
 */
int SCIPpatternGetCircleType(
   SCIP_PATTERN *pattern             /**< pattern */
)
{
   assert(pattern != NULL);
   assert(pattern->patterntype == SCIP_PATTERNTYPE_CIRCULAR);

   return pattern->type;
}

/** sets the type of the boundary circle
 *
 * @note this function can only be called for circular patterns
 */
void SCIPpatternSetType(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< type */
   )
{
   assert(pattern != NULL);
   assert(pattern->patterntype == SCIP_PATTERNTYPE_CIRCULAR);

   pattern->type = type;
}

/** returns the packable status of a pattern */
SCIP_PACKABLE SCIPpatternGetPackableStatus(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);

   return pattern->packable;
}

/** sets the packable status of a pattern */
void SCIPpatternSetPackableStatus(
   SCIP_PATTERN*         pattern,            /**< pattern */
   SCIP_PACKABLE         packable            /**< packable status */
   )
{
   assert(pattern != NULL);

   pattern->packable = packable;
}
