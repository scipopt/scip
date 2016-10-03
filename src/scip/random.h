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

/**@file   struct_random.h
 * @brief  data structures for random number generator
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RANDOM_H__
#define __SCIP_RANDOM_H__

#include "scip/def.h"
#include "scip/struct_random.h"

/** creates a random number generator */
extern
SCIP_RETCODE SCIPrandomInit(
   SCIP_RANDGEN*         randgen,
   int                   initseed
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

#ifdef __cplusplus
}
#endif

#endif
