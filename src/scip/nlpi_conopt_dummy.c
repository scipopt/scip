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

/**@file    nlpi_conopt_dummy.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   dummy CONOPT NLP interface
 * @author  Ksenia Bestuzheva
 *
 * This file contains dummy implementations of the interface methods for the CONOPT interface.
 * It is used when SCIP is build without CONOPT.
 */

#include "scip/nlpi_conopt.h"

/** create solver interface for CONOPT solver and includes it into SCIP, if CONOPT is available */  /*lint -e{715}*/
SCIP_RETCODE SCIPincludeNlpSolverConopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{  /*lint --e{715}*/
   assert(scip != NULL);

   return SCIP_OKAY;
}

/** gets string that identifies CONOPT */
const char* SCIPgetSolverNameConopt(void)
{
   return "CONOPT";
}

/** gets string that describes CONOPT */
const char* SCIPgetSolverDescConopt(void)
{
   return "CONOPT not available";
}

/** returns whether CONOPT is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisConoptAvailableConopt(void)
{
   return FALSE;
}

/** sets the license to be passed to CONOPT's COIDEF_License */
void SCIPsetLicenseConopt(
   SCIP_NLPI*            nlpi,               /**< CONOPT NLPI */
   int                   integer_1,          /**< CONOPT_LICENSE_INT_1 */
   int                   integer_2,          /**< CONOPT_LICENSE_INT_2 */
   int                   integer_3,          /**< CONOPT_LICENSE_INT_3 */
   const char*           text                /**< CONOPT_LICENSE_TEXT */
   )
{ } /*lint !e715*/
