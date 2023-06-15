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

/**@file   cutsel_ensemble.h
 * @ingroup CUTSELECTORS
 * @brief  ensemble cut selector
 * @author Mark Turner
 *
 * The ensemble cut selector scores cuts by using a weighted sum of normalised efficacy,
 * normalised directed cutoff distance (only at the root node), normalised expected objective improvement,
 * objective parallelism, objective orthogonality, integer support, density, dynamism,
 * and normalised number of locks.
 * It also has a variety of filtering methods. If density filtering is enabled, then it filters all cuts before
 * scoring over some relative density threshold. After scoring, it selects the cuts with the highest score,
 * where after each selection, the remaining cuts are either filtered or penalised via parallelism checks.
 * Whether the cuts are filtered or penalised is a users choice.
 *
 * If a cut is given by \f$ a^T x \leq b \f$, then
 *  - the efficacy is defined as the distance between the LP solution and the hyperplane \f$ a^T x = b \f$.
 *    It is normalised by the largest efficacy from the given array of cuts. ((log(eff(cut) + 1) / log(maxeff + 1))**2 ;
 *  - the directed cutoff distance is defined as the distance between the LP solution and the hyperplane \f$ a^T x = b \f$
 *    restricted to the line segment joining the LP solution to the currently best primal solution; therefore, it is only
 *    defined when a primal solution is available. It is normalised by the largest efficacy from the given array of
 *    cuts. ((log(dcd(cut) + 1) / log(maxdcd + 1))**2;
 *  - the objective parallelism is how parallel the vector \f$ a \f$ is w.r.t. the objective function \f$ c \f$. That
 *    is, the objective parallelism is given by \f$ \frac{a^T c}{\|a\| \|c\|} \f$. Notice that the vectors are parallel
 *    when this formula returns 1;
 *  - the objective orthogonality is how orthogonal the vector \f$ a \f$ is w.r.t. the objective function \f$ c \f$. That
 *    is, the objective orthogonality is given by \f$ 1 - \frac{a^T c}{\|a\| \|c\|} \f$.
 *    Notice that the vectors are orthogonal when this formula returns 1;
 *  - the expected objective improvement is defined by the difference of the objective evaluated at the LP solution
 *    and when evaluated at the orthogonal projection onto the cut. As we normalise the value, we remove the
 *    objective vector multiplication from its calculation. We calculate it as efficacy * objective parallelism.
 *    We also normalise it according to the equation ((log(expimprov(cut) + 1) / log(expimprov + 1))**2;
 *  - the integer support of a cut is the ratio between the number of nonzero integer columns and the number of nonzero
 *    columns.
 *  - the density of a cut is the ratio of non-zero entries divided by the number of columns in the LP;
 *  - the dynamism of a cut is the ratio between the maximum absolute value over all coefficients and the
 *    minimum absolute value over all coefficients. If this ratio is below a threshold, we give the cut a reward for
 *    good numerics;
 *  - the number of locks is the amount
 *
 * These features of a cut can be recovered and/or computed with the functions @ref SCIPgetCutEfficacy(), @ref
 * SCIPgetCutLPSolCutoffDistance(), @ref SCIPgetRowObjParallelism(), and @ref SCIPgetRowNumIntCols(), @ref
 * SCIProwGetNNonz(), @ref SCIProwGetVals(), @ref SCIProwGetCols(), @ref SCIPvarGetNLocksUp(),
 * @ref SCIPvarGetNLocksDown().
 *
 * The filtering (density) works as follows.
 *
 * The non-forced cuts are scanned through, and any cut over a density threshold (0,1) is dropped from
 * consideration.
 *
 * The filtering / penalise (parallelism) step works as follows.
 *
 * First, the forced cuts --- cuts that are going to enter the LP no matter what --- are used to filter / penalise
 * the non-forced cuts. This means that for each forced cut, @p fcut, the parallelism between @p fcut and
 * every non-forced cut, @p cut, is computed (the parallelism between two cuts \f$ a^T x \leq b \f$ and \f$ d^T x \leq e\f$
 * is \f$ \frac{a^T d}{\|a\| \|d\|} \f$).
 * If the parallelism with @p fcut is larger or equal than some maximum threshold then it is either removed from
 * consideration (if filter by parallelism), or its score is decreased (if penalise by parallelism).
 *
 * @note The maximum parallelism is a parameter that can be set, as well as the weights for the score.
 *
 * @note In the case of no primal solution, the weight assigned to the directed cutoff distance is transferred to the
 * efficacy.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTSEL_ENSEMBLE_H__
#define __SCIP_CUTSEL_ENSEMBLE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the ensemble separator and includes it in SCIP
 *
 * @ingroup CutSelectorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutselEnsemble(
   SCIP*                 scip                /**< SCIP data structure */
);

/**@addtogroup CUTSELECTORS
 *
 * @{
 */

/** perform a cut selection algorithm for the given array of cuts
 *
 *  This is the selection method of the ensemble cut selector. It uses a weighted sum of normalised efficacy,
 *  normalised directed cutoff distance, normalised expected improvements, objective parallelism,
 *  objective orthogonality, integer support, sparsity, and dynamism.
 *  As well as the weighted sum scoring there is optional parallelism based filtering, parallelism based penalties,
 *  and density filtering.
 *  There are also additional budget constraints on the amount of cuts that should be added.
 *  The input cuts array gets resorted s.t the selected cuts come first and the remaining ones are the end.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPselectCutsEnsemble(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
   SCIP_CUTSELDATA*      cutseldata,         /**< cut selector data */
   SCIP_Bool             root,               /**< whether we are currently at the root node or not */
   int                   ncuts,              /**< number of cuts in cuts array */
   int                   nforcedcuts,        /**< number of forced cuts */
   int                   maxselectedcuts,    /**< maximal number of cuts from cuts array to select */
   int*                  nselectedcuts       /**< pointer to return number of selected cuts from cuts array */
);

/** @} */

#ifdef __cplusplus
}
#endif

#endif
