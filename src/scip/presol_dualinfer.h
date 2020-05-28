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

/**@file   presol_dualinfer.h
 * @ingroup PRESOLVERS
 * @brief  dual inference presolver
 * @author Dieter Weninger
 * @author Patrick Gemander
 *
 * This presolver does bound strengthening on continuous variables (columns) for getting bounds on dual variables y.
 * The bounds of the dual variables are then used to fix primal variables or change the side of constraints.
 * For ranged rows one needs to decide which side (rhs or lhs) determines the equality.
 *
 * We distinguish two cases concerning complementary slackness:
 * i)  reduced cost fixing:       c_j - sup_y(y^T A_{.j}) > 0 => x_j = l_j
 *                                c_j - inf_y(y^T A_{.j}) < 0 => x_j = u_j
 * ii) positive dual lower bound: y_i > 0 =>  A_{i.}x = b_i
 *
 * Further information on this presolving approach are given in
 * Achterberg et al. "Presolve reductions in mixed integer programming"
 * and for a two-column extension in
 * Chen et al. "Two-row and two-column mixed-integer presolve using hasing-based pairing methods".
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_DUALINFER_H__
#define __SCIP_PRESOL_DUALINFER_H__

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dual inference presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePresolDualinfer(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
