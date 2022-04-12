/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cutsel_xyz.h
 * @ingroup CUTSELECTORS
 * @brief  xyz cut selector
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTSEL_XYZ_H__
#define __SCIP_CUTSEL_XYZ_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the xyz separator and includes it in SCIP
 *
 * @ingroup CutSelectorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeCutselXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup CUTSELECTORS
 *
 * @{
 */

/** perform a cut selection algorithm for the given array of cuts
 *
 *  The input cuts array should be resorted such that the selected cuts come first.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPselectCutsXyz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            cuts,               /**< array with cuts to perform selection algorithm */
   SCIP_ROW**            forcedcuts,         /**< array with forced cuts */
   int                   nforcedcuts,        /**< number of forced cuts */
   int                   maxselectedcuts,    /**< maximal number of cuts from cuts array to select */
   int*                  nselectedcuts       /**< pointer to return number of selected cuts from cuts array */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
