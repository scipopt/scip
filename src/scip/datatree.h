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

/**@file    datatree.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for handling data trees
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_DATATREE_H__
#define __SCIP_DATATREE_H__

#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/type_datatree.h"
#include "scip/type_paramset.h"
#include "scip/type_message.h"
#include "scip/type_mem.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a new SCIP_DATATREE with a given capacity for items */
SCIP_RETCODE SCIPdatatreeCreate(
   SCIP_DATATREE**       datatree,           /**< buffer to store pointer to created data tree */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   capacity            /**< initial capacity */
   );

/** frees a SCIP_DATATREE object */
void SCIPdatatreeFree(
   SCIP_DATATREE**       datatree,           /**< pointer to datatree to free */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** inserts a SCIP_Bool value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertBool(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_Bool             value               /**< value of entry */
   );

/** inserts a long value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertLong(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_Longint          value               /**< value of entry */
   );

/** inserts a SCIP_Real value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertReal(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_Real             value               /**< value of entry */
   );

/** inserts a string value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertString(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const char*           value               /**< value of entry */
   );

/** inserts a SCIP_Bool array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertBoolArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const SCIP_Bool*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a SCIP_Real array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertRealArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const SCIP_Real*      values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a SCIP_Longint array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertLongArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const SCIP_Longint*   values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a string array into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertStringArray(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   const char* const*    values,             /**< values of entry */
   int                   nvalues             /**< number of values */
   );

/** inserts a store value into a SCIP_DATATREE object */
SCIP_RETCODE SCIPdatatreeInsertTree(
   SCIP_DATATREE*        datatree,           /**< data tree*/
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   const char*           name,               /**< name of entry */
   SCIP_DATATREE*        value               /**< value of entry */
   );

/** writes a SCIP_DATATREE object as JSON to file */
SCIP_RETCODE SCIPdatatreeWriteJson(
   SCIP_DATATREE*        datatree,           /**< data tree */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   FILE*                 file                /**< file to write to, or NULL for stdout */
   );

#ifdef __cplusplus
}
#endif

#endif
