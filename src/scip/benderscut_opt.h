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

/**@file   benderscut_opt.h
 * @ingroup BENDERSCUTS
 * @brief  Generates a standard Benders' decomposition optimality cut
 * @author Stephen J. Maher
 *
 * The classical Benders' decomposition optimality cuts arise from a feasible instance of the Benders' decomposition
 * subproblem. The optimality cuts are an underestimator of the subproblem objective function value. Auxiliary
 * variables, \f$\varphi\f$ are added to the master problem as alower bound on the subproblem objective function value.
 *
 * Consider the Benders' decomposition subproblem that takes the master problem solution \f$\bar{x}\f$ as input:
 * \f[
 * z(\bar{x}) = \min\{d^{T}y : Ty \geq h - H\bar{x}, y \geq 0\}
 * \f]
 * If the subproblem is feasible, and \f$z(\bar{x}) > \varphi\f$ (indicating that the current underestimators are not
 * optimal) then the Benders' decomposition optimality cut can be generated from the optimal dual solution of the
 * subproblem. Let \f$w\f$ be the vector corresponding to the optimal dual solution of the Benders' decomposition
 * subproblem. The resulting cut is:
 * \f[
 * \varphi \geq w^{T}(h - Hx)
 * \f]
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSCUT_OPT_H__
#define __SCIP_BENDERSCUT_OPT_H__


#include "scip/def.h"
#include "scip/type_benders.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the optimality Benders' decomposition cut and includes it in SCIP
 *
 *  @ingroup BenderscutIncludes
 */
EXTERN
SCIP_RETCODE SCIPincludeBenderscutOpt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition */
   );

/* @} */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
