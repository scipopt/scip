/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   nlhdlr_signomial.h
 * @ingroup NLHDLRS
 * @brief  signomial nonlinear handler
 * @author Liding Xu
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_SIGNOMIAL_H__
#define __SCIP_NLHDLR_SIGNOMIAL_H__

#include "scip/scip.h"
#include "scip/pub_nlhdlr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes signomial nonlinear handler to nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrSignomial(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLHDLRS
 * @{
 *
 * @name Signomial nonlinear handler
 *
 * This nonlinear handler detects and collects signomial terms and provides specialized propagation and estimation functionality.
 *
 * @{
 */

/** @}
 * @}
 */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_SIGNOMIAL_H__ */
