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

/**@file   symmetry_orbital.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries based on orbits
 * @author Jasper van Doornmalen
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

/** data for orbital reduction propagator */
struct SCIP_OrbitalReductionData;
typedef struct SCIP_OrbitalReductionData SCIP_ORBITALREDDATA; /**< data for orbital reduction propagator */

/*
 * Interface methods
 */

/** prints orbital reduction data */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitalReductionGetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA*  orbireddata,        /**< orbital reduction data structure */
   int*                  nred,               /**< pointer to store the total number of reductions applied */
   int*                  ncutoff             /**< pointer to store the total number of cutoffs applied */
   );


/** prints orbital reduction data */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitalReductionPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA*  orbireddata         /**< orbital reduction data structure */
   );


/** propagates orbital reduction */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitalReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA*  orbireddata,        /**< orbital reduction data structure */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nred,               /**< pointer to store the number of domain reductions */
   SCIP_Bool*            didrun              /**< a global pointer maintaining if any symmetry propagator has run
                                              *   only set this to TRUE when a reduction is found, never set to FALSE */
   );


/** adds component for orbital reduction */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitalReductionAddComponent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA*  orbireddata,        /**< orbital reduction data structure */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int**                 perms,              /**< permutations in the component */
   int                   nperms,             /**< number of permutations in the component */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   );


/** resets orbital reduction data structure (clears all components) */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitalReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA*  orbireddata         /**< orbital reduction data structure */
   );


/** frees orbital reduction data */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitalReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA** orbireddata         /**< orbital reduction data structure */
   );


/** initializes structures needed for orbital reduction
 *
 *  This is only done exactly once.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeOrbitalReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITALREDDATA** orbireddata,        /**< pointer to orbital reduction data structure to populate */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   );

#ifdef __cplusplus
}
#endif

#endif
