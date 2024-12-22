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

/**@file   type_iisfinder.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for IIS
 * @author Mark Turner
 */

/** @defgroup DEFPLUGINS_IISFINDER Default irreducible infeasible subsystems
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default IISs (irreducible infeasible subsystems) of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_IISFINDER_H__
#define __SCIP_TYPE_IISFINDER_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_lp.h"
#include "scip/type_result.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_IISfinder SCIP_IISFINDER;         /**< IIS finder data structure */
typedef struct SCIP_IISfinderData SCIP_IISFINDERDATA; /**< IIS finder specific data */
typedef struct SCIP_IIS SCIP_IIS;                     /**< IIS storage data structure */


/** copy method for IIS finder plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - iisfinder       : the IIS (irreducible infeasible subsystem) finder itself
 */
#define SCIP_DECL_IISFINDERCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_IISFINDER* iisfinder)

/** destructor of IIS finder plugins to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - iisfinder       : the IIS finder itself
 */
#define SCIP_DECL_IISFINDERFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_IISFINDER* iisfinder)

/** IIS finder execution method to generate an IIS
 *
 *  This method is called to generate an IIS (irreducible infeasible subsystem) for an infeasible problem.
 *  It creates a copy of the SCIP instance and performs different algorithms to create a still infeasible
 *  reduced problem.
 *
 *  input:
 *  - iis             : The IIS data structure. It contains a subscip.
 *  - iisfinder       : the IIS finder itself
 *  - timelim         : The time limit for the total run
 *  - nodelim         : The node limit for the total run
 *  - removebounds    : Whether bounds should also be considered for removal by the IIS finder
 *  - silent          : Whether algorithm specific information should be output during calculations
 *  - result          : pointer to store the result of the IIS finder call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_SUCCESS    : the IIS finder succeeded
 *  - SCIP_DIDNOTFIND : the IIS finder did not find a small enough infeasible subsystem.
 *  - SCIP_DIDNOTRUN  : the IIS finder did not run because some criteria was not satisfied
 */
#define SCIP_DECL_IISFINDEREXEC(x) SCIP_RETCODE x (SCIP_IIS* iis, SCIP_IISFINDER* iisfinder, SCIP_Real timelim, SCIP_Longint nodelim, SCIP_Bool removebounds, SCIP_Bool silent, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
