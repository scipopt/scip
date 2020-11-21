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

/**@file   benderscut_nogood.h
 * @ingroup BENDERSCUTS
 * @brief  Generates a no-good cut for solutions that are integer infeasible
 * @author Stephen J. Maher
 *
 * The no-good cut is generated for the Benders' decomposition master problem if an integer solution is identified as
 * infeasible in at least one CIP subproblems. The no-good cut is required, because the classical Benders' decomposition
 * feasibility cuts (see benderscut_feas.c) will only cut off the solution \f$\bar{x}\f$ if the LP relaxation of the CIP
 * is infeasible.
 *
 * Consider a Benders' decomposition subproblem that is a CIP and it infeasible. Let \f$S_{r}\f$ be the set of indices
 * for master problem variables that are 1 in \f$\bar{x}\f$. The no-good cut is given by
 *
 * \f[
 * 1 \leq \sum_{i \in S_{r}}(1 - x_{i}) + \sum_{i \notin S_{r}}x_{i}
 * \f]
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_NOGOOD_H__
#define __SCIP_BENDERSCUT_NOGOOD_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the no good Benders' decomposition cut and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBenderscutNogood(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

#ifdef __cplusplus
}
#endif

#endif
