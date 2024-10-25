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

/**@file   pub_iis.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for irreducible infeasible subsytems (IIS)
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_IIS_H__
#define __SCIP_PUB_IIS_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_iis.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicIISMethods
 *
 * @{
 */

/** gets name of IIS */
SCIP_EXPORT
const char* SCIPiisGetName(
   SCIP_IIS*             iis                 /**< IIS */
   );

/** gets user data of IIS */
SCIP_EXPORT
SCIP_IISDATA* SCIPiisGetData(
   SCIP_IIS*             iis                 /**< IIS */
   );

/** gets description of IIS */
SCIP_EXPORT
const char* SCIPiisGetDesc(
   SCIP_IIS*             iis                 /**< IIS */
   );

/** gets priority of IIS */
SCIP_EXPORT
int SCIPiisGetPriority(
   SCIP_IIS*             iis                 /**< IIS */
   );

/** sets user data of IIS; user has to free old data in advance! */
SCIP_EXPORT
void SCIPiisSetData(
   SCIP_IIS*             iis ,               /**< IIS */
   SCIP_IISDATA*         iisdata             /**< new IIS user data */
   );

/** gets time in seconds used in this IIS */
SCIP_EXPORT
SCIP_Real SCIPiisGetTime(
   SCIP_IIS*             iis                 /**< IIS */
   );

/** get number of times the IIS was called */
SCIP_Longint SCIPiisGetNCalls(
   SCIP_IIS*             iis                 /**< IIS */
   );

/** compares two IIS w. r. to their priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPiisComp);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
