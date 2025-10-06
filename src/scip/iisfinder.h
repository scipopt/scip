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

/**@file   iisfinder.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for IIS finder
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_IISFINDER_H__
#define __SCIP_IISFINDER_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/pub_iisfinder.h"
#include "scip/lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates an IIS finder */
SCIP_RETCODE SCIPiisfinderCreate(
   SCIP_IISFINDER**      iisfinder,          /**< pointer to store IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of IIS finder */
   const char*           desc,               /**< description of IIS finder */
   int                   priority,           /**< priority of the IIS finder in standard mode */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)), /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)), /**< destructor of IIS finder */
   SCIP_DECL_IISFINDEREXEC ((*iisfinderexec)), /**< IIS finder execution method */
   SCIP_IISFINDERDATA*   iisfinderdata       /**< IIS finder data */
   );

/** enables or disables all clocks of @p iisfinder, depending on the value of the flag */
void SCIPiisfinderEnableOrDisableClocks(
   SCIP_IISFINDER*       iisfinder,          /**< the IIS finder for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the IIS finder be enabled? */
   );

/** calls the IIS finder generation method */
SCIP_RETCODE SCIPiisGenerate(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** copies the given IIS finder to a new scip */
SCIP_RETCODE SCIPiisfinderCopyInclude(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** sets copy method of IIS finder */
void SCIPiisfinderSetCopy(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)) /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** frees memory of IIS finder */
SCIP_RETCODE SCIPiisfinderFree(
   SCIP_IISFINDER**      iisfinder,          /**< pointer to IIS finder data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** sets destructor method of IIS finder */
void SCIPiisfinderSetFree(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)) /**< destructor of IIS finder */
   );

/** sets priority of IIS finder */
void SCIPiisfinderSetPriority(
   SCIP_IISFINDER*       iisfinder,          /**< IIS finder */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the IIS finder */
   );

/** creates and captures a new IIS */
SCIP_RETCODE SCIPiisCreate(
   SCIP_IIS**            iis,                /**< pointer to return the created IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** releases an IIS */
SCIP_RETCODE SCIPiisFree(
   SCIP_IIS**            iis,                /**< pointer to the IIS */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** reset an IIS (in case one exists from a previous solve) */
SCIP_RETCODE SCIPiisReset(
   SCIP_IIS**            iis                 /**< pointer to the IIS */
   );

#ifdef __cplusplus
}
#endif

#endif
