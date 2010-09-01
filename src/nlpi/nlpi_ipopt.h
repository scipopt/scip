/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: nlpi_ipopt.h,v 1.5 2010/09/01 12:30:23 bzfviger Exp $"

/**@file    nlpi_ipopt.h
 * @brief   Ipopt NLP interface
 * @ingroup NLPINTERFACES
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_IPOPT_H__
#define __SCIP_NLPI_IPOPT_H__

#ifdef WITH_IPOPT

#include "nlpi/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for Ipopt solver */
extern
SCIP_RETCODE SCIPcreateNlpSolverIpopt(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
);

/** gets string that identifies Ipopt (version number) */
extern
const char* SCIPgetSolverNameIpopt(void);

/** gets string that describes Ipopt (version number) */
extern
const char* SCIPgetSolverDescIpopt(void);

/** gives a pointer to the IpoptApplication object stored in Ipopt-NLPI's NLPI problem data structure */
extern
void* SCIPgetIpoptApplication(
   SCIP_NLPIPROBLEM*     nlpiproblem         /**< NLP problem of Ipopt-NLPI */
   );

#ifdef __cplusplus
}
#endif

#else /* WITH_IPOPT not defined */

#define SCIPcreateNlpSolverIpopt(blkmem, nlpi) SCIP_OKAY
#define SCIPgetSolverNameIpopt()               ""
#define SCIPgetSolverDescIpopt()               ""

#endif

#endif /* __SCIP_NLPI_IPOPT_H__ */
