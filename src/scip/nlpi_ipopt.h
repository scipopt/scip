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

/**@file    nlpi_ipopt.h
 * @brief   Ipopt NLP interface
 * @ingroup NLPINTERFACES
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_IPOPT_H__
#define __SCIP_NLPI_IPOPT_H__

#include "scip/type_nlpi.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for Ipopt solver */
extern
SCIP_RETCODE SCIPcreateNlpSolverIpopt( 
   SCIP*                 scip,               /**< central scip datastructure */
   SCIP_NLPI**           nlpi                /**< pointer to buffer for nlpi address */
);

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_IPOPT_H__ */
