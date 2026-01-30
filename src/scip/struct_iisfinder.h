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

/**@file   struct_iisfinder.h
 * @ingroup INTERNALAPI
 * @brief  data structures for irreducible infeasible subsystems (IIS)
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_IISFINDER_H__
#define __SCIP_STRUCT_IISFINDER_H__


#include "scip/def.h"
#include "scip/type_iisfinder.h"

#ifdef __cplusplus
extern "C" {
#endif

/** IIS */
struct SCIP_IISfinder
{
   char*                 name;               /**< name of IIS finder */
   char*                 desc;               /**< description of IIS finder */
   SCIP_DECL_IISFINDERCOPY ((*iisfindercopy)); /**< copy method of IIS finder or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_IISFINDERFREE ((*iisfinderfree)); /**< destructor of IIS finder */
   SCIP_DECL_IISFINDEREXEC ((*iisfinderexec)); /**< IIS finder generation method */
   SCIP_CLOCK*           iisfindertime;      /**< IIS finder execution time */
   SCIP_IISFINDERDATA*   iisfinderdata;      /**< IIS finder data */
   int                   priority;           /**< priority of the IIS finder */
};

/** IIS */
struct SCIP_IIS
{
   SCIP*                 subscip;            /**< The subscip that stores the IIS */
   SCIP_HASHMAP*         varsmap;            /**< The variable hashmap from the original SCIP to IIS subscip */
   SCIP_HASHMAP*         conssmap;           /**< The constraints hashmap from the original SCIP to IIS subscip */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator */
   SCIP_CLOCK*           iistime;            /**< IIS total execution time */
   int                   niismessagecalls;   /**< The number of times an iis info message has been displayed */
   SCIP_Longint          nnodes;             /**< The number of nodes used over all IIS solves */
   SCIP_Bool             infeasible;         /**< Whether the subscip is currently infeasible, i.e., a valid IS */
   SCIP_Bool             irreducible;        /**< Whether the subscip is an irreducible infeasible subsystem, i.e., an IIS */
};

#ifdef __cplusplus
}
#endif

#endif
