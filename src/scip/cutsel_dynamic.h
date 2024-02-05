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

/**@file   cutsel_dynamic.h
 * @ingroup CUTSELECTORS
 * @brief  dynamic cut selector
 * @author Christoph Graczyk
 *
 *  The dynamic cut selector is an extension of the hybrid cut selector to employ a dynamic range for the orthogonality
 *  filtering, dependent on the efficacy ratio between cuts. It allows the user to specify a minimum efficacy gain in
 *  percentage to filter cuts, which corresponds geometrically to the idealized shift from the LP solution adding only
 *  the first cut to the intersection with the second cut.
 *  In this way we can enforce a minimum orthogonality between cuts through the efficacy ratio. This should, in theory,
 *  avoid the selection of cutting plane pairs that do not improve the lp solution, but only increase the number of cuts
 *  in the lp.
 *
 *  The selector allows the user to specify a filtering strategy during cut selection ('d'ynamic- and 'f'ull dynamic),
 *  which determines how the orthogonality filtering is applied.
 *  In both cases, after a cut is selected, the remaining cuts are resorted by computing their relative efficacy to the
 *  selected cut. The efficacy ratio is then used to filter the cuts based on the filtering strategy:
 *  - 'd'ynamic-parallelism: The dynamic parallelism strategy filters cuts based on the efficacy ratio between cut
 *  pairs. It only filters cuts that are not theoretically improving the efficacy of the pair given the first cut, i.e.,
 *  the efficacy ratio is below the minimum efficacy gain in percentage, without calculating or resorting the cut array.
 *
 *  - 'f'ull dynamic parallelism: The full dynamic parallelism strategy filters cuts based on the efficacy ratio between
 *  cut pairs by iteratively recomputing the score given the previously selected cut. This results in a set of cuts
 *  that is a sequence of theoretically most effective almost orthogonal cuts, whereby almost orthogonal means that
 *  the efficacy ratio is above the minimum efficacy gain in percentage.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTSEL_DYNAMIC_H__
#define __SCIP_CUTSEL_DYNAMIC_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dynamic hybrid separator and includes it in SCIP
 *
 * @ingroup CutSelectorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutselDynamic(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CUTSELECTORS
 *
 * @{
 */

/** perform a cut selection algorithm for the given array of cuts
 *
 *  This is an extension of the hybrid cutselector to employ a dynamic range
 *  when applying orthogonality filtering, dependent on the efficacy ratio between cuts.
 *
 *  The input cuts array should be resorted such that the selected cuts come first.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPselectCutsDynamic(
    SCIP*                 scip,               /**< SCIP data structure */
    SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
    SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
    SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
    char                  filtermode,         /**< filtering strategy during cut selection (
                                              *  'd'ynamic- and 'f'ull dynamic parallelism) */
    SCIP_Real             mingain,            /**< minimum efficacy gain in percentage to filter cuts */
    SCIP_Real             maxparall,          /**< maximal parallelism for all cuts that are not good */
    SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
    SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
    SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
    SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
    int                   ncuts,              /**< number of cuts in cuts array */
    int                   nforcedcuts,        /**< number of forced cuts */
    int                   maxselectedcuts,    /**< maximal number of cuts from cuts array to select */
    int*                  nselectedcuts      /**< pointer to return number of selected cuts from cuts array */
  );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
