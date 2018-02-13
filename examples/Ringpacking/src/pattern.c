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

/*
 * local methods
 */

/*
 * interface methods
 */

/** creates an empty circular pattern */
SCIP_RETCODE SCIPpatternCreateCircularEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   )
{
   assert(scip != NULL);
   assert(pattern != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, pattern) );
   BMSclearMemory(*pattern);

   (*pattern)->type = SCIP_PATTERNTYPE_CIRCULAR;
   SCIPpatternCapture(*pattern);

   return SCIP_OKAY;
}

/** creates an empty rectangular pattern */
SCIP_RETCODE SCIPpatternCreateRectangularEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, pattern) );
   BMSclearMemory(*pattern);

   (*pattern)->type = SCIP_PATTERNTYPE_RECTANGULAR;
   SCIPpatternCapture(*pattern);

   return SCIP_OKAY;
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
      SCIPfreeBlockMemory(scip, pattern);
   }
   else
      *pattern = NULL;
}

/** returns the type of a pattern */
SCIP_PATTERNTYPE SCIPpatternGetType(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);

   return pattern->type;}

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
