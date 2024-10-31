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

/**@file   scip_iis.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for IIS plugins
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_IIS_H__
#define __SCIP_SCIP_IIS_H__


#include "scip/def.h"
#include "scip/type_iis.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicIISMethods
 *
 * @{
 */

/** creates an IIS and includes it in SCIP
 *
 *  @note this method has all IIS callbacks as arguments and is thus changed every time a new
 *        callback is added in future releases; consider using SCIPincludeIISBasic() and setter functions
 *        if you seek for a method which is less likely to change in future releases
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeIIS(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of cut selector */
   const char*           desc,               /**< description of cut selector */
   int                   priority,           /**< priority of the cut selector */
   SCIP_DECL_IISCOPY     ((*iiscopy)),        /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFREE     ((*iisfree)),       /**< destructor of IIS */
   SCIP_DECL_IISGENERATE ((*iisgenerate)),   /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   );

/** Creates a cut selector and includes it in SCIP with its most fundamental callbacks.
 *
 *  All non-fundamental (or optional) callbacks as, e.g., copy and free callbacks, will be set to NULL. Optional
 *  callbacks can be set via specific setter functions, see SCIPsetIISCopy() and SCIPsetIISFree(),
 *
 *  @note if you want to set all callbacks with a single method call, consider using SCIPincludeIIS() instead
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeIISBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS**            iis,                /**< reference to an IIS, or NULL */
   const char*           name,               /**< name of IIS */
   const char*           desc,               /**< description of IIS */
   int                   priority,           /**< priority of the IIS in standard mode */
   SCIP_DECL_IISGENERATE((*iisgenerate)),    /**< IIS generation method */
   SCIP_IISDATA*         iisdata             /**< IIS data */
   );

/** sets copy method of IIS */
SCIP_EXPORT
SCIP_RETCODE SCIPsetIISCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISCOPY     ((*iiscop))         /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor method of IIS */
SCIP_EXPORT
SCIP_RETCODE SCIPsetIISFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS*             iis,                /**< IIS */
   SCIP_DECL_IISFREE     ((*iisfree))        /**< destructor of IIS */
   );

/** the execution method that iterates over the IIS plugins */
SCIP_EXPORT
SCIP_RETCODE SCIPgenerateIIS(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the IIS of the given name, or NULL if not existing */
SCIP_EXPORT
SCIP_IIS* SCIPfindIIS(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of IIS */
   );

/** returns the array of currently available IISs */
SCIP_EXPORT
SCIP_IIS** SCIPgetIIS(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns the number of currently available IISs */
SCIP_EXPORT
int SCIPgetNIIS(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the priority of an IIS */
SCIP_EXPORT
SCIP_RETCODE SCIPsetIISPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_IIS*             iis,                /**< IIS */
   int                   priority            /**< new priority of the IIS */
   );

/** Gets the IIS storage.
 *
 *  @return the \ref SCIP_IISSTORE iis storage.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_IISSTORE* SCIPgetIISstore(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
