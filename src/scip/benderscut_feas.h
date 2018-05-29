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

/**@file   benderscut_feas.h
 * @ingroup BENDERSCUTS
 * @brief  Standard feasibility cuts for Benders' decomposition
 * @author Stephen J. Maher
 *
 * The classical Benders' decomposition feasibility cuts arise from an infeasible instance of the Benders' decomposition
 * subproblem.
 * Consider the Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : Ty \geq h - H\bar{x}, y \geq 0\}
 * \f]
 * If the subproblem is infeasible as a result of the solution \f$\bar{x}\f$, then the Benders' decomposition
 * feasibility cut can be generated from the dual ray. Let \f$w\f$ be the vector corresponding to the dual ray of the
 * Benders' decomposition subproblem. The resulting cut is:
 * \f[
 * 0 \geq w^{T}(h - Hx)
 * \f]
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_FEAS_H__
#define __SCIP_BENDERSCUT_FEAS_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the Standard Feasibility Benders' decomposition cuts and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBenderscutFeas(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

#ifdef __cplusplus
}
#endif

#endif
