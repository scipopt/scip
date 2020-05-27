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

/**@file   benders_default.h
 * @ingroup BENDERSDECOMPOSITION
 * @brief  default Benders' decomposition plugin
 * @author Stephen J. Maher
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERS_DEFAULT_H__
#define __SCIP_BENDERS_DEFAULT_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the default Benders' decomposition and includes it in SCIP
 *
 *  @ingroup BendersIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBendersDefault(
   SCIP*                 scip                /**< SCIP data structure */
   );


/**@addtogroup BENDERS
 *
 * @{
 */

/** Creates a default Benders' decomposition algorithm and activates it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateBendersDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   int                   nsubproblems        /**< the number of subproblems in the Benders' decomposition */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
