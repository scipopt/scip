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

/* NLP interface for the CONOPT solver. */

/**@file    nlpi_conopt.h
 * @brief   CONOPT NLP interface
 * @ingroup NLPIS
 * @author  Ksenia Bestuzheva
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_CONOPT_H__
#define __SCIP_NLPI_CONOPT_H__

#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for CONOPT solver and includes it into SCIP
 *
 * @ingroup NLPIIncludes
 */
SCIP_RETCODE SCIPincludeNlpSolverConopt(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLPIS
 *
 * @{
 */

/** sets the license to be passed to CONOPT's COIDEF_License */
SCIP_EXPORT
void SCIPsetLicenseConopt(
   SCIP_NLPI*            nlpi,               /**< CONOPT NLPI */
   int                   integer_1,          /**< CONOPT_LICENSE_INT_1 */
   int                   integer_2,          /**< CONOPT_LICENSE_INT_2 */
   int                   integer_3,          /**< CONOPT_LICENSE_INT_3 */
   const char*           text                /**< CONOPT_LICENSE_TEXT */
   );

/** gets string that identifies CONOPT */
SCIP_EXPORT
const char* SCIPgetSolverNameConopt(
   void
   );

/** gets string that describes CONOPT */
SCIP_EXPORT
const char* SCIPgetSolverDescConopt(
   void
   );

/** returns whether CONOPT is available, i.e., whether it has been linked in */
SCIP_EXPORT
SCIP_Bool SCIPisConoptAvailableConopt(
   void
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_CONOPT_H__ */
