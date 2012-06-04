/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_solex.h
 * @brief  public methods for exact primal CIP solutions
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SOLEX_H__
#define __SCIP_PUB_SOLEX_H__

#include "scip/def.h"
#include "scip/type_solex.h"

#ifdef NDEBUG
#include "scip/struct_solex.h"
#endif

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets origin of exact solution */
extern
SCIP_SOLORIGIN SCIPsolexGetOrigin(
   SCIP_SOLEX*           sol                 /**< exact primal CIP solution */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */
#define SCIPsolexGetOrigin(sol)           ((sol)->solorigin)
#endif

#endif
