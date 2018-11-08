/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_fuzzyround.h
 * @ingroup PRIMALHEURISTICS
 * @brief  primal heuristic that constructs a feasible solution from the lp-relaxation. Round only on the state-variables (binvars)
 * and then reconstruct the rest of the variables accordingly.
 * @author Leon Eifler
 *
 */

#ifndef __SCIP_HEUR_FUZZYROUND_H__
#define __SCIP_HEUR_FUZZYROUND_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the fuzzy rounding primal heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeHeurFuzzyround(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
