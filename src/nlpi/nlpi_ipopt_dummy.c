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

/**@file    nlpi_ipopt_dummy.c
 * @brief   dummy Ipopt NLP interface for the case that Ipopt is not available
 * @author  Stefan Vigerske
 *
 * This code has been separate from nlpi_ipopt.cpp, so the SCIP build system recognizes it as pure C code,
 * thus the linker does not need to be changed to C++.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "nlpi/nlpi_ipopt.h"

/** create solver interface for Ipopt solver */
SCIP_RETCODE SCIPcreateNlpSolverIpopt(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   )
{
   assert(nlpi != NULL);

   *nlpi = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/

/** gets string that identifies Ipopt (version number) */
const char* SCIPgetSolverNameIpopt(void)
{
   return "";
}

/** gets string that describes Ipopt (version number) */
const char* SCIPgetSolverDescIpopt(void)
{
   return "";
}

/** returns whether Ipopt is available, i.e., whether it has been linked in */
SCIP_Bool SCIPisIpoptAvailableIpopt(void)
{
   return FALSE;
}

/** gives a pointer to the IpoptApplication object stored in Ipopt-NLPI's NLPI problem data structure */
void* SCIPgetIpoptApplicationPointerIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   )
{
   SCIPerrorMessage("Ipopt not available!\n");
   SCIPABORT();
   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** gives a pointer to the NLPIORACLE object stored in Ipopt-NLPI's NLPI problem data structure */
void* SCIPgetNlpiOracleIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   )
{
   SCIPerrorMessage("Ipopt not available!\n");
   SCIPABORT();
   return NULL;  /*lint !e527*/
}  /*lint !e715*/

/** sets modified default settings that are used when setting up an Ipopt problem
 *
 * Do not forget to add a newline after the last option in optionsstring.
 */
void SCIPsetModifiedDefaultSettingsIpopt(
   SCIP_NLPI*            nlpi,               /**< Ipopt NLP interface */
   const char*           optionsstring       /**< string with options as in Ipopt options file */
   )
{
   SCIPerrorMessage("Ipopt not available!\n");
   SCIPABORT();
}  /*lint !e715*/

/** Calls Lapacks Dsyev routine to compute eigenvalues and eigenvectors of a dense matrix. 
 * It's here, because Ipopt is linked against Lapack.
 */
SCIP_RETCODE LapackDsyev(
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   )
{
   SCIPerrorMessage("Ipopt not available, cannot use it's Lapack link!\n");
   return SCIP_ERROR;
}  /*lint !e715*/
