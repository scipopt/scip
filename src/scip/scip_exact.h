/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_exact.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for exact solving
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_EXACT_H__
#define __SCIP_SCIP_EXACT_H__


#include "scip/def.h"
#include "scip/type_cuts.h"
#include "scip/type_cons.h"
#include "scip/type_heur.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"
#include "scip/type_lpexact.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicExactMethods
 *
 * @{
 */

/** enables or disables exact solving mode
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INIT
 */
SCIP_EXPORT
SCIP_RETCODE SCIPenableExactSolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             enable              /**< enable exact solving (TRUE) or disable it (FALSE) */
   );

/** returns whether the solution process is arithmetically exact, i.e., not subject to roundoff errors
 *
 *  @return Returns TRUE if \SCIP is in exact solving mode, otherwise FALSE
 */
SCIP_EXPORT
SCIP_Bool SCIPisExact(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether aggregation is allowed to use negative slack in exact solving mode
 *
 *  @return Returns TRUE if \SCIP is not in exact solving mode or aggregation is allowed to use negative slack
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Bool SCIPallowNegSlack(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** branches on an LP solution exactly; does not call branching rules, since fractionalities are assumed to small;
 *  if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPbranchLPExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/** adds row to exact separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        rowexact            /**< exact row to add */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
