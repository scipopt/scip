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

/**@file   benders_xyz.h
 * @ingroup BENDERSDECOMPOSITION
 * @brief  xyz Benders' decomposition
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_XYZ_H__
#define __SCIP_BENDERS_XYZ_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the xyz Benders' decomposition and includes it in SCIP
 *
 *  @ingroup BendersIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBendersXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup BENDERSS
 *
 * @{
 */

/** TODO: add public methods to this group for documentation purposes */

/** @} */

#ifdef __cplusplus
}
#endif

#endif
