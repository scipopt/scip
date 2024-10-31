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

/**@file   type_iis.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for IIS
 * @author Mark Turner
 */

/** @defgroup DEFPLUGINS_IIS Default irreducible infeasible subsystems
 *  @ingroup DEFPLUGINS
 *  @brief implementation files (.c files) of the default IISs (irreducible infeasible subsystems) of SCIP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_IIS_H__
#define __SCIP_TYPE_IIS_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_lp.h"
#include "scip/type_result.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_IIS SCIP_IIS;               /**< IIS data structure */
typedef struct SCIP_IISData SCIP_IISDATA;       /**< IIS specific data */
typedef struct SCIP_IISSTORE SCIP_IISSTORE;     /**< IIS storage data structure */


/** copy method for IIS plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - iis             : the IIS (irreducible infeasible subsystem) itself
 */
#define SCIP_DECL_IISCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_IIS* iis)

/** destructor of IIS plugins to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - iis             : the IIS itself
 */
#define SCIP_DECL_IISFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_IIS* iis)

/** IIS generation method of IIS
 *
 *  This method is called to generate an IIS (irreducible infeasible subsystem) for an infeasible problem.
 *  It creates a copy of the SCIP instance and performs different algorithms to create a still infeasible
 *  reduced problem.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - iis             : the IIS itself
 *  - result          : pointer to store the result of the IIS call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_SUCCESS    : the IIS succeeded
 *  - SCIP_DIDNOTFIND : the IIS did not find a small enough infeasible subsystem.
 */
#define SCIP_DECL_IISGENERATE(x) SCIP_RETCODE x (SCIP* scip, SCIP_IIS* iis, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
