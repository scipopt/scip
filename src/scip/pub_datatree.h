/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file    pub_datatree.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for managing data trees
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_DATATREE_H__
#define __SCIP_PUB_DATATREE_H__

#include "scip/def.h"
#include "scip/type_datatree.h"
#include "scip/type_set.h"
#include "scip/type_table.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicDatatreeMethods
 *
 * @{
 */

/** gets a SCIP_Bool value from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetBool(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Bool*            value               /**< buffer to store value */
   );

/** gets a long value from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetLong(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Longint*         value               /**< buffer to store value */
   );

/** gets a SCIP_Real value from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetReal(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Real*            value               /**< buffer to store value */
   );

/** gets a string value from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetString(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   const char**          value               /**< buffer to store pointer to string */
   );

/** gets a SCIP_Bool array from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetBoolArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Bool**           values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   );

/** gets a SCIP_Longint array from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetLongArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Longint**        values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   );

/** gets a SCIP_Real array from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetRealArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_Real**           values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   );

/** gets a string array from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetStringArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   char***               values,             /**< buffer to store pointer to values */
   int*                  nvalues             /**< buffer to store number of values */
   );

/** gets a data tree value from a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPdatatreeGetTree(
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name to look up */
   SCIP_DATATREE**       value               /**< buffer to store pointer to data tree */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
