/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   scip_datatree.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for data tree structure
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_DATATREE_H__
#define __SCIP_SCIP_DATATREE_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_set.h"
#include "scip/type_table.h"
#include "scip/type_datatree.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicDatatreeMethods
 *
 * @{
 */

/** creates a new SCIP_DATATREE */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateDatatree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE**       datatree,           /**< buffer to store created data tree */
   int                   capacity            /**< desired capacity, or -1 for default */
   );

/** creates a new SCIP_DATATREE and inserts it into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateDatatreeInTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree where to insert new data tree */
   SCIP_DATATREE**       newtree,            /**< buffer to store pointer to created store */
   const char*           name,               /**< name of entry to add */
   int                   capacity            /**< capacity of new store, or -1 for default */
   );

/** inserts a bool value into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeBool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_Bool             value               /**< value to add */
   );

/** inserts an int value into a SCIP_DATATREE object
 *
 *  The value will be stored as SCIP_Longint.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   int                   value               /**< value to add */
   );

/** inserts a long value into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeLong(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_Longint          value               /**< value to add */
   );

/** inserts a SCIP_Real value into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeReal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_Real             value               /**< value to add */
   );

/** inserts a string value into a SCIP_DATATREE object
 *
 *  The string value will be copied.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeString(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const char*           value               /**< value to add */
   );

/** inserts a SCIP_Bool array into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeBoolArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const SCIP_Bool*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts an int array into a SCIP_DATATREE object
 *
 *  The value will be stored as array of SCIP_Longint.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeIntArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const int*            values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a SCIP_Longint array into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeLongArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const SCIP_Longint*   values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a SCIP_Real array into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeRealArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const SCIP_Real*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a string array into a SCIP_DATATREE object */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeStringArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   const char* const*    values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a data tree value into a SCIP_DATATREE object
 *
 *  The data tree assumes ownership of value.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPinsertDatatreeTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   const char*           name,               /**< name of entry to add */
   SCIP_DATATREE*        value               /**< value to add */
   );

/** frees a SCIP_DATATREE object */
SCIP_EXPORT
void SCIPfreeDatatree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE**       datatree            /**< pointer to data tree to free */
   );

/** writes a SCIP_DATATREE object as JSON to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPwriteDatatreeJson(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file to write to, or NULL for stdout */
   SCIP_DATATREE*        datatree            /**< data tree to write */
   );

/** prints a generic table from a data store */
SCIP_EXPORT
SCIP_RETCODE SCIPprintDatatreeAsTable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DATATREE*        datatree,           /**< data tree */
   FILE*                 file,               /**< output file */
   const char*           sectionname,        /**< section name to process, e.g., "plugins" */
   const char*           tablename           /**< table name to process, e.g., "heuristics" */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
