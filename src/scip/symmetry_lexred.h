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

/**@file   symmetry_lexred.h
 * @ingroup OTHER_CFILES
 * @brief  methods for handling symmetries by dynamic lexicographic ordering reduction
 * @author Jasper van Doornmalen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYMMETRY_LEXRED_H__
#define __SCIP_SYMMETRY_LEXRED_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "scip/type_event.h"
#include "symmetry/type_symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif


/*
 * Data structures
 */

/** data for dynamic lexicographic reduction propagator */
struct SCIP_LexRedData;
typedef struct SCIP_LexRedData SCIP_LEXREDDATA; /**< data for dynamic lexicographic reduction propagator */

/*
 * Interface methods
 */

/** prints lexicographic reduction propagation data */
SCIP_EXPORT
SCIP_RETCODE SCIPlexicographicReductionGetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   int*                  nred,               /**< total number of reductions applied */
   int*                  ncutoff             /**< total number of cutoffs applied */
   );


/** prints lexicographic reduction propagation data */
SCIP_EXPORT
SCIP_RETCODE SCIPlexicographicReductionPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata          /**< pointer to global data for lexicographic reduction propagator */
   );


/** applies lexicographic reduction propagation */
SCIP_EXPORT
SCIP_RETCODE SCIPlexicographicReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nred,               /**< pointer to store the number of domain reductions */
   SCIP_Bool*            didrun              /**< a global pointer maintaining if any symmetry propagator has run
                                              *   only set this to TRUE when a reduction is found, never set to FALSE */
   );

/** adds permutation for lexicographic reduction propagation */
SCIP_EXPORT
SCIP_RETCODE SCIPlexicographicReductionAddPermutation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_VAR**            permvars,           /**< variable array of the permutation */
   int                   npermvars,          /**< number of variables in that array */
   int*                  perm,               /**< permutation */
   SYM_SYMTYPE           symtype,            /**< type of symmetries in perm */
   SCIP_Real*            permvardomaincenter, /**< array containing center point for each variable domain */
   SCIP_Bool             usedynamicorder,    /**< whether a dynamic variable order shall be used */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   );

/** resets lexicographic reduction propagation (removes all permutations) */
SCIP_EXPORT
SCIP_RETCODE SCIPlexicographicReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA*      masterdata          /**< pointer to global data for lexicographic reduction propagator */
   );


/** frees lexicographic reduction data */
SCIP_EXPORT
SCIP_RETCODE SCIPlexicographicReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA**     masterdata          /**< pointer to global data for lexicographic reduction propagator */
   );


/** initializes structures needed for lexicographic reduction propagation
 *
 *  This is only done exactly once.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeLexicographicReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LEXREDDATA**     masterdata,         /**< pointer to global data for lexicographic reduction propagator */
   SCIP_EVENTHDLR*       shadowtreeeventhdlr /**< pointer to the shadow tree eventhdlr */
   );

#ifdef __cplusplus
}
#endif

#endif
