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

/**@file   symmetry_orbitopal.c
 * @ingroup OTHER_CFILES
 * @brief  methods for handling orbitopal symmetries
 * @author Jasper van Doornmalen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SYMMETRY_ORBITOPAL_H__
#define __SCIP_SYMMETRY_ORBITOPAL_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_var.h"
#include "symmetry/type_symmetry.h"

#ifdef __cplusplus
extern "C" {
#endif


/** variants for orbitope column ordering */
enum SCIP_ColumnOrdering
{
   SCIP_COLUMNORDERING_NONE        = 0,      /**< do not order the columns */
   SCIP_COLUMNORDERING_FIRST       = 1,      /**< choose first possible column */
   SCIP_COLUMNORDERING_LAST        = 2,      /**< choose last possible column */
   SCIP_COLUMNORDERING_CENTRE      = 3,      /**< choose centremost possible column */
   SCIP_COLUMNORDERING_MEDIAN      = 4       /**< choose median column */
};
typedef enum SCIP_ColumnOrdering SCIP_COLUMNORDERING; /**< variants for orbitope column ordering*/

/** variants for orbitope row ordering */
enum SCIP_RowOrdering
{
   SCIP_ROWORDERING_NONE           = 0,      /**< do not order the rows */
   SCIP_ROWORDERING_BRANCHING      = 1       /**< choose rows based on branching variables */
};
typedef enum SCIP_RowOrdering SCIP_ROWORDERING; /**< variants for orbitope row ordering*/


/** data for orbitopal reduction */
struct SCIP_OrbitopalReductionData;
typedef struct SCIP_OrbitopalReductionData SCIP_ORBITOPALREDDATA; /**< data for orbitopal reduction */


/** gets the number of reductions */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitopalReductionGetStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< orbitopal reduction data structure */
   int*                  nred,               /**< total number of reductions applied */
   int*                  ncutoff             /**< total number of cutoffs applied */
   );


/** prints orbitopal reduction data */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitopalReductionPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata        /**< orbitopal reduction data structure */
   );

/** propagates orbitopal reduction */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitopalReductionPropagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< orbitopal reduction data structure */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility is found */
   int*                  nred,               /**< pointer to store the number of domain reductions */
   SCIP_Bool*            didrun              /**< a global pointer maintaining if any symmetry propagator has run
                                              *   only set this to TRUE when a reduction is found, never set to FALSE */
   );


/** adds orbitopal component to orbitopal symmetry handler */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitopalReductionAddOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata,       /**< orbitopal reduction data structure */
   SCIP_ROWORDERING      rowordering,        /**< specifies how rows of orbitope are ordered */
   SCIP_COLUMNORDERING   colordering,        /**< specifies how columnss of orbitope are ordered */
   SCIP_VAR**            vars,               /**< matrix of variables on which the symmetry acts */
   int                   nrows,              /**< number of rows */
   int                   ncols,              /**< number of columns */
   SCIP_Bool*            success             /**< to store whether the component is successfully added */
   );

/** resets orbitopal reduction data structure (clears all orbitopes) */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitopalReductionReset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA* orbireddata        /**< pointer to orbitopal reduction structure to populate */
   );


/** frees orbitopal reduction data */
SCIP_EXPORT
SCIP_RETCODE SCIPorbitopalReductionFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA** orbireddata       /**< pointer to orbitopal reduction structure to populate */
   );


/** initializes structures needed for orbitopal reduction
 *
 *  This is only done exactly once.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeOrbitopalReduction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ORBITOPALREDDATA** orbireddata       /**< pointer to orbitopal reduction structure to populate */
   );


/** returns the default column ordering */
SCIP_EXPORT
SCIP_COLUMNORDERING SCIPorbitopalReductionGetDefaultColumnOrdering(
   SCIP_ORBITOPALREDDATA* orbireddata        /**< pointer to orbitopal reduction structure to populate */
   );

#ifdef __cplusplus
}
#endif

#endif
