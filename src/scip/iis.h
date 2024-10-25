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

/**@file   iis.h
 * @ingroup INTERNALAPI
 * @brief  internal methods for IIS
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_IIS_H__
#define __SCIP_IIS_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_set.h"
#include "scip/pub_iis.h"
#include "scip/lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates an IIS */
SCIP_RETCODE SCIPiisCreate(
   SCIP_IIS**            iis,                /**< pointer to store IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of IIS */
   const char*           desc,               /**< description of IIS */
   int                   priority,           /**< priority of the IIS in standard mode */
   SCIP_DECL_IISCOPY     ((*iiscopy)),       /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFREE     ((*iisfree)),       /**< destructor of IIS */
   SCIP_DECL_IISGENERATE ((*iisgenerate)),   /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   );

/** enables or disables all clocks of @p iis, depending on the value of the flag */
void SCIPiisEnableOrDisableClocks(
   SCIP_IIS*             iis,                /**< the IIS for which all clocks should be enabled or disabled */
   SCIP_Bool             enable              /**< should the clocks of the IIS be enabled? */
   );

/** calls the IIS generation method */
SCIP_RETCODE SCIPiisGenerate(
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** copies the given IIS to a new scip */
SCIP_RETCODE SCIPiisCopyInclude(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** sets copy method of IIS */
void SCIPiisSetCopy(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISCOPY     ((*iiscopy))        /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** frees memory of IIS */
SCIP_RETCODE SCIPiisFree(
   SCIP_IIS**            iis,                /**< pointer to IIS data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** sets destructor method of IIS */
void SCIPiisSetFree(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISFREE     ((*iisfree))        /**< destructor of IIS */
   );

/** sets priority of IIS */
void SCIPiisSetPriority(
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the IIS */
   );

#ifdef __cplusplus
}
#endif

#endif
