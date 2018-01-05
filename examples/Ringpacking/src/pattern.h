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

/**@file   pattern.h
 * @brief  pattern data for ringpacking problem
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PATTERN__
#define __SCIP_PATTERN__

#include "scip/scip.h"
#include "probdata_rpa.h"

typedef struct SCIP_Pattern SCIP_PATTERN;

enum SCIP_Patterntype
{
   SCIP_PATTERNTYPE_CIRCULAR    = 0,           /**< circular pattern */
   SCIP_PATTERNTYPE_RECTANGULAR = 1            /**< rectangular pattern */
};
typedef enum SCIP_Patterntype SCIP_PATTERNTYPE;

/** creates an empty circular pattern */
extern
SCIP_RETCODE SCIPpatternCreateCircularEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   );

/** creates an empty rectangular pattern */
extern
SCIP_RETCODE SCIPpatternCreateRectangularEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to store pattern */
   );

/* frees a pattern */
extern
void SCIPpatternFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PATTERN**        pattern             /**< pointer to free pattern */
   );

/** returns the type of a pattern */
extern
SCIP_PATTERNTYPE SCIPpatternGetType(
   SCIP_PATTERN*         pattern             /**< pattern */
   );

#endif /* __SCIP_PATTERN__ */
