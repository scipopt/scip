/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    nlpi_ipopt.h
 * @brief   Ipopt NLP interface
 * @ingroup NLPIS
 * @author  Stefan Vigerske
 * @author  Benjamin Müller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_IPOPT_H__
#define __SCIP_NLPI_IPOPT_H__

#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup NLPIS
 *
 * @{
 */

/** create solver interface for Ipopt solver
 * sets *nlpi to NULL if Ipopt is not available
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateNlpSolverIpopt(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
   );

/** gets string that identifies Ipopt (version number) */
SCIP_EXPORT
const char* SCIPgetSolverNameIpopt(void);

/** gets string that describes Ipopt */
SCIP_EXPORT
const char* SCIPgetSolverDescIpopt(void);

/** returns whether Ipopt is available, i.e., whether it has been linked in */
SCIP_EXPORT
SCIP_Bool SCIPisIpoptAvailableIpopt(void);

/** gives a pointer to the IpoptApplication object stored in Ipopt-NLPI's NLPI problem data structure */
SCIP_EXPORT
void* SCIPgetIpoptApplicationPointerIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   );

/** gives a pointer to the NLPIORACLE object stored in Ipopt-NLPI's NLPI problem data structure */
SCIP_EXPORT
void* SCIPgetNlpiOracleIpopt(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   );

/** sets modified default settings that are used when setting up an Ipopt problem
 *
 * Do not forget to add a newline after the last option in optionsstring.
 */
SCIP_EXPORT
void SCIPsetModifiedDefaultSettingsIpopt(
   SCIP_NLPI*            nlpi,               /**< Ipopt NLP interface */
   const char*           optionsstring,      /**< string with options as in Ipopt options file */
   SCIP_Bool             append              /**< whether to append to modified default settings or to overwrite */
   );

/** Calls Lapacks Dsyev routine to compute eigenvalues and eigenvectors of a dense matrix. 
 * It's here, because Ipopt is linked against Lapack.
 */
SCIP_EXPORT
SCIP_RETCODE LapackDsyev(
   SCIP_Bool             computeeigenvectors,/**< should also eigenvectors should be computed ? */
   int                   N,                  /**< dimension */
   SCIP_Real*            a,                  /**< matrix data on input (size N*N); eigenvectors on output if computeeigenvectors == TRUE */
   SCIP_Real*            w                   /**< buffer to store eigenvalues (size N) */
   );

/** solves a linear problem of the form Ax = b for a regular matrix A
 *
 *  Calls Lapacks IpLapackDgetrf routine to calculate a LU factorization and uses this factorization to solve
 *  the linear problem Ax = b.
 *  It's here, because Ipopt is linked against Lapack.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsolveLinearProb(
   int                   N,                  /**< dimension */
   SCIP_Real*            A,                  /**< matrix data on input (size N*N); filled column-wise */
   SCIP_Real*            b,                  /**< right hand side vector (size N) */
   SCIP_Real*            x,                  /**< buffer to store solution (size N) */
   SCIP_Bool*            success             /**< pointer to store if the solving routine was successful */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_IPOPT_H__ */
