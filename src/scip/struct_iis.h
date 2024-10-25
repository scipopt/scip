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

/**@file   struct_iis.h
 * @ingroup INTERNALAPI
 * @brief  data structures for irreducible infeasible subsystems (IIS)
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_IIS_H__
#define __SCIP_STRUCT_IIS_H__


#include "scip/def.h"
#include "scip/type_iis.h"

#ifdef __cplusplus
extern "C" {
#endif

/** IIS */
struct SCIP_IIS
{
   char*                 name;               /**< name of IIS */
   char*                 desc;               /**< description of IIS */
   SCIP_DECL_IISCOPY     ((*iiscopy));       /**< copy method of IIS or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFREE     ((*iisfree));       /**< destructor of IIS */
   SCIP_DECL_IISGENERATE ((*iisgenerate));   /**< IIS generation method */
   SCIP_CLOCK*           iistime;            /**< IIS execution time */
   SCIP_IISDATA*         iisdata;            /**< IIS data */
   int                   priority;           /**< priority of the IIS */
   SCIP_Longint          ncalls;             /**< number of times, this IIS was called */
};

#ifdef __cplusplus
}
#endif

#endif
