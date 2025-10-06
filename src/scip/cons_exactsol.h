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

/**@file   cons_exactsol.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for ensuring that primal solution is exact
 * @author Antonia Chmiela
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/lp.h"
#include "scip/lpexact.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
#include "scip/scip_exact.h"
#include "scip/scip_lpexact.h"
#include "scip/scip_sol.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/tree.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for ExactSol constraints and includes it in SCIP
 *
 * @ingroup ConshdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrExactSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif
