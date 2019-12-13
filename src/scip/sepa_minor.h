/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
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
 * This separator detects all principal minors of the matrix xx' for which all auxiliary variables X exist, i.e., two
 * indices i != j such that X_ii, X_jj, and X_ij exist. Because xx' - X is required to be positive semi-definite, it
 * follows that the matrix
 *
 *              [ 1    x_i  x_j  ]
 *    A(x,X) =  [ x_i  X_ii X_ij ]
 *              [ x_j  X_ij X_jj ]
 *
 * is also required to be positive semi-definite. Let v be a negative eigenvector for A(x*,X*) in a point (x*,X*), which
 * implies that v' A(x*,X*) v < 0. To cut off (x*,X*), the separator computes the globally valid linear inequality
 * v' A(x,X) v >= 0.
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
