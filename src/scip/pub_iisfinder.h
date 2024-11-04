/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   pub_iisfinder.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for irreducible infeasible subsytems (IIS) finders
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_IISFINDER_H__
#define __SCIP_PUB_IISFINDER_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_iisfinder.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicIISfinderMethods
 *
 * @{
 */

/** gets name of IIS finder */
SCIP_EXPORT
const char* SCIPiisfinderGetName(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   );

/** gets user data of IIS finder */
SCIP_EXPORT
SCIP_IISFINDERDATA* SCIPiisfinderGetData(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   );

/** gets description of IIS finder */
SCIP_EXPORT
const char* SCIPiisfinderGetDesc(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   );

/** gets priority of IIS finder */
SCIP_EXPORT
int SCIPiisfinderGetPriority(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   );

/** sets user data of IIS finder; user has to free old data in advance! */
SCIP_EXPORT
void SCIPiisfinderSetData(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< new IIS finder user data */
   );

/** gets time in seconds used in this IIS finder */
SCIP_EXPORT
SCIP_Real SCIPiisfinderGetTime(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   );

/** get number of times the IIS finder was called */
SCIP_Longint SCIPiisfinderGetNCalls(
   SCIP_IISFINDER*       iisfinder           /**< IIS finder */
   );

/** compares two IIS finders w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPiisfinderComp);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
