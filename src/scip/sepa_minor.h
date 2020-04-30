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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_minor.h
 * @ingroup SEPARATORS
 * @brief  principal minor separator
 * @author Benjamin Mueller
 *
 * This separator detects all principal minors of the matrix \f$ xx' \f$ for which all auxiliary variables \f$ X \f$
 * exist, i.e., two indices \f$ i \neq j \f$ such that \f$ X_{ii} \f$, \f$ X_{jj} \f$, and \f$ X_{ij} \f$ exist. Because
 * \f$ xx' - X \f$ is required to be positive semi-definite, it follows that the matrix
 *
 * \f[
 *    A(x,X) = \begin{bmatrix} 1 & x_i & x_j \\ x_i & X_{ii} & X_{ij} \\ x_j & X_{ij} & X_{jj} \end{bmatrix}
 * \f]
 *
 * is also required to be positive semi-definite. Let \f$ v \f$ be a negative eigenvector for \f$ A(x^*,X^*) \f$ in a
 * point \f$ (x^*,X^*) \f$, which implies that \f$ v' A(x^*,X^*) v < 0 \f$. To cut off \f$ (x^*,X^*) \f$, the separator
 * computes the globally valid linear inequality \f$ v' A(x,X) v \ge 0 \f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_MINOR_H__
#define __SCIP_SEPA_MINOR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the minor separator and includes it in SCIP
 *
 * @ingroup SeparatorIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeSepaMinor(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup SEPARATORS
 *
 * @{
 */

/* TODO place other public methods in this group to facilitate navigation through the documentation */

/* @} */

#ifdef __cplusplus
}
#endif

#endif
