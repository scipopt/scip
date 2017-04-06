/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_sparsify.h
 * @brief  cancel non-zeros of the constraint matrix
 * @author Dieter Weninger
 *
 * This presolver attempts to cancel non-zero entries of the constraint
 * matrix by adding equalities to other constraints.
 * It supports two cases:
 * a) the variable index set from non-zero entries of the other constraint
 *    is a superset of the variable index set of the non-zero entries of the
 *    equality
 * b) same as case a) except for one variable index
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_SPARSIFY_H__
#define __SCIP_PRESOL_SPARSIFY_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the sparsify presolver and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludePresolSparsify(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
