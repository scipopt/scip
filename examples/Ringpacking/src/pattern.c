/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pattern.c
 * @brief  pattern data for ringpacking problem
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "pattern.h"
#include "probdata_rpa.h"

/*
 * local methods
 */

/** auxiliary function to create a pattern */
static
SCIP_RETCODE createPattern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern,            /**< pointer to store pattern */
   SCIP_PATTERNTYPE      patterntype,        /**< pattern type */
   int                   ntypes,             /**< number of different circle types */
   int                   type                /**< circle type (not needed for rectangular patterns) */
   )
{
   assert(scip != NULL);
   assert(pattern != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, pattern) );
   BMSclearMemory(*pattern);

   (*pattern)->type = type;
   SCIPpatternCapture(*pattern);

   (*pattern)->ntypes = ntypes;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*pattern)->nelems, ntypes) );
   BMSclearMemoryArray((*pattern)->nelems, ntypes);

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
   int                   ntypes,             /**< number of different circle types */
   int                   type                /**< circle type (not needed for rectangular patterns) */
   )
{
   return createPattern(scip, pattern, SCIP_PATTERNTYPE_CIRCULAR, ntypes, type);
}

/** creates an empty circular pattern */
SCIP_RETCODE SCIPpatternCreateRectangular(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern,            /**< pointer to store pattern */
   int                   ntypes              /**< number of different circle types */
   )
{
   return createPattern(scip, pattern, SCIP_PATTERNTYPE_RECTANGULAR, ntypes, -1);
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

   /* free memory if the pattern is not used any more */
   if( (*pattern)->nlocks == 0 )
   {
      SCIPfreeBlockMemoryArray(scip, &(*pattern)->nelems, (*pattern)->ntypes);
      SCIPfreeBlockMemory(scip, pattern);
   }
   else
      *pattern = NULL;
}

/** adds an element of a given type to a pattern */
void SCIPpatternAddElement(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< element of a given type */
   )
{
   assert(pattern != NULL);
   assert(type >= 0);

   ++(pattern->nelems[type]);
}

/** removes an element of a given type from a pattern */
void SCIPpatternRemoveElement(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< element of a given type */
   )
{
   assert(pattern != NULL);
   assert(type >= 0);

   --(pattern->nelems[type]);
}

/** returns the number of elements of a given type in the pattern */
int SCIPpatternGetNElemens(
   SCIP_PATTERN*         pattern,            /**< pattern */
   int                   type                /**< element of a given type */
   )
{
   assert(pattern != NULL);
   assert(pattern->nelems != NULL);
   assert(type >= 0 && type < pattern->ntypes);

   return pattern->nelems[type];
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
int SCIPpatternGetType(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);
   assert(pattern->patterntype == SCIP_PATTERNTYPE_CIRCULAR);

   return pattern->type;
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
