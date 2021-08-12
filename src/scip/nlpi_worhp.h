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

/**@file    nlpi_worhp.h
 * @brief   Worhp NLP interface
 * @ingroup NLPIS
 * @author  Benjamin Mueller
 * @author  Renke Kuhlmann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_WORHP_H__
#define __SCIP_NLPI_WORHP_H__

#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for Worhp solver and includes it into SCIP, if Worhp is available
 *
 * @ingroup NLPIIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlpSolverWorhp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             useip               /**< TRUE for using Interior Point, FALSE for SQP */
   );

/**@addtogroup NLPIS
 *
 * @{
 */

/** gets string that identifies Worhp (version number) */
SCIP_EXPORT
const char* SCIPgetSolverNameWorhp(
   void
   );

/** gets string that describes Worhp (version number) */
SCIP_EXPORT
const char* SCIPgetSolverDescWorhp(
   void
   );

/** returns whether Worhp is available, i.e., whether it has been linked in */
SCIP_EXPORT
SCIP_Bool SCIPisWorhpAvailableWorhp(
   void
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_WORHP_H__ */
