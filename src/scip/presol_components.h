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

/**@file   presol_components.h
 * @ingroup PRESOLVERS
 * @brief  components presolver
 * @author Dieter Weninger
 * @author Gerald Gamrath
 *
 * This presolver looks for independent components at the end of the presolving.
 * If independent components are found in which a maximum number of discrete variables
 * is not exceeded, the presolver tries to solve them in advance as subproblems.
 * Afterwards, if a subproblem was solved to optimality, the corresponding
 * variables/constraints can be fixed/deleted in the main problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_COMPONENTS_H__
#define __SCIP_PRESOL_COMPONENTS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the components presolver and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludePresolComponents(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
