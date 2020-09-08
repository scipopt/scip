/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cutsel_default.h
 * @ingroup CUTSELECTORS
 * @brief  default cut selector
 * @author Leona Gottwald
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTSEL_DEFAULT_H__
#define __SCIP_CUTSEL_DEFAULT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the default separator and includes it in SCIP
 *
 * @ingroup CutSelectorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutselDefault(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CUTSELECTOR
 *
 * @{
 */

/** perform a cut selection algorithm for the given array of cuts;
 *  this is the selection method of the default cut selector which does a weighted sum of the
 *  efficacy, parallelism, directed cutoff distance, and the integral support.
 *  The input cuts array gets resorted s.t the selected cuts comes first and the remaining
 *  ones are the end.
 */
SCIP_RETCODE SCIPselectCutsDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
   SCIP_RANDNUMGEN*      randnumgen,         /**< random number generator for tie-breaking, or NULL */
   SCIP_Real             goodscorefac,       /**< factor of best score among the given cuts to consider a cut good
                                              *   and filter with less strict settings of the maximum parallelism */
   SCIP_Real             badscorefac,        /**< factor of best score among the given cuts to consider a cut bad
                                              *   and discard it regardless of its parallelism to other cuts */
   SCIP_Real             goodmaxparall,      /**< maximum parallelism for good cuts */
   SCIP_Real             maxparall,          /**< maximum parallelism for non-good cuts */
   SCIP_Real             dircutoffdistweight,/**< weight of directed cutoff distance in cut score calculation */
   SCIP_Real             efficacyweight,     /**< weight of efficacy in cut score calculation */
   SCIP_Real             objparalweight,     /**< weight of objective parallelism in cut score calculation */
   SCIP_Real             intsupportweight,   /**< weight of integral support in cut score calculation */
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
