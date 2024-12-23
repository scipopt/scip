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

/**@file    nlpi_ipopt_dummy.c
 * @ingroup DEFPLUGINS_NLPI
 * @brief   dummy Ipopt NLP interface for the case that Ipopt is not available
 * @author  Stefan Vigerske
 * @author  Benjamin Müller
 *
 * This code has been separated from nlpi_ipopt.cpp, so the SCIP build system recognizes it as pure C code,
 * thus the linker does not need to be changed to C++.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "scip/nlpi_ipopt.h"
#include "blockmemshell/memory.h"

/** create solver interface for Ipopt solver and includes it into SCIP, if Ipopt is available */  /*lint -e{715}*/
SCIP_RETCODE SCIPincludeNlpSolverIpopt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIP_OKAY;
}  /*lint !e715*/

/** gets string that identifies Ipopt (version number) */
const char* SCIPgetSolverNameIpopt(void)
{
   return "";
}

/** gets string that describes Ipopt */
const char* SCIPgetSolverDescIpopt(void)
{
   return "";
}

/** returns whether Ipopt is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisIpoptAvailableIpopt(void)
{
   return FALSE;
}

/** gives a pointer to the NLPIORACLE object stored in Ipopt-NLPI's NLPI problem data structure */ /*lint -e715*/
void* SCIPgetNlpiOracleIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   )
{
   SCIPerrorMessage("Ipopt not available!\n");
   SCIPABORT();
   return NULL;  /*lint !e527*/
}  /*lint !e715*/

SCIP_RETCODE SCIPcallLapackDsyevIpopt(
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   )
{
   SCIPerrorMessage("Ipopt not available, cannot use it's Lapack link!\n");
   return SCIP_ERROR;
}  /*lint !e715*/

SCIP_RETCODE SCIPsolveLinearEquationsIpopt(
   int                   N,                  /**< dimension */
   SCIP_Real*            A,                  /**< matrix data on input (size N*N); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size N) */
   SCIP_Real*            x,                  /**< buffer to store solution (size N) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   )
{
   SCIPerrorMessage("Ipopt not available, cannot use it's Lapack link!\n");
   return SCIP_ERROR;
}  /*lint !e715*/
