/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_dualinfer.h
 * @brief  dual inference presolver
 * @author Dieter Weninger
 *
 * This presolver exploits dual information for primal variable fixings:
 * a) The first method is an enhanced dual fixing technique.
 * b) The second method does dual bound strengthening on continuous primal
 *    variables and applies complementary slackness (yA-c)_i > 0 => x_i = 0
 *    for fixing primal variables at their lower bound.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_DUALINFER_H__
#define __SCIP_PRESOL_DUALINFER_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the dual inference presolver and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePresolDualinfer(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
