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
 * data structures
 */

struct SCIP_Pattern
{
   SCIP_PATTERNTYPE      type;               /**< pattern type */
};

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

   return SCIP_OKAY;
}

/* frees a pattern */
void SCIPpatternFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to free pattern */
   )
{
   assert(scip != NULL);
   assert(pattern != NULL);
   assert(*pattern != NULL);

   SCIPfreeBlockMemory(scip, pattern);
}

/** returns the type of a pattern */
SCIP_PATTERNTYPE SCIPpatternGetType(
   SCIP_PATTERN*         pattern             /**< pattern */
   )
{
   assert(pattern != NULL);

   return pattern->type;
}
