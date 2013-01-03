/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_dualfix.h
 * @ingroup PRESOLVERS
 * @brief  fixing roundable variables to best bound
 * @author Tobias Achterberg
 *
 * This presolver fixes variables, that have no restrictions in direction of their objective coefficient, to the best
 * possible value. If the objective coefficient of a variable is \f$0\f$ and it may be rounded both up and down, then
 * this variable will be fixed to the closest feasible value to \f$0\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_DUALFIX_H__
#define __SCIP_PRESOL_DUALFIX_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dual fixing presolver and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludePresolDualfix(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
