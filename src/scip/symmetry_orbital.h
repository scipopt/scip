/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
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

/**@file   symmetry_orbital.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries based on orbits
 * @author Jasper van Doornmalen
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYMMETRY_ORBITAL_H__
#define __SCIP_SYMMETRY_ORBITAL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "scip/type_event.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Data structures
 */

/** data for dynamic orbital fixing propagator */
struct SCIP_OrbitalFixingData;
typedef struct SCIP_OrbitalFixingData SCIP_ORBITALFIXINGDATA;

/*
 * Interface methods
 */

/** propagate orbital fixing */
SCIP_RETCODE SCIPorbitalFixingPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALFIXINGDATA* orbifixdata,      /**< orbitopal fixing data structure */
   SCIP_Bool*            infeasible,         /**< whether infeasibility is found */
   int*                  nred                /**< number of domain reductions */
   );


/** adds component for orbital fixing */
SCIP_RETCODE SCIPorbitalFixingAddComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALFIXINGDATA* orbifixdata,      /**< orbital fixing data structure */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int**                 perms,              /**< permutations in the component */
   int                   nperms              /**< number of permutations in the component */
   );


/** resets orbital fixing data structure (clears all components) */
SCIP_RETCODE SCIPorbitalFixingReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALFIXINGDATA* orbifixdata       /**< orbital fixing data structure */
   );


/** free orbital fixing data */
SCIP_RETCODE SCIPorbitalFixingFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALFIXINGDATA** orbifixdata      /**< orbital fixing data structure */
   );


/** initializes structures needed for orbital fixing
 * This is only done exactly once.
 */
SCIP_RETCODE SCIPorbitalFixingInclude(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALFIXINGDATA** orbifixdata,     /**< pointer to orbital fixing data structure to populate */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   );

#ifdef __cplusplus
}
#endif

#endif
