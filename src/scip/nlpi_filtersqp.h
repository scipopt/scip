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

/**@file    nlpi_filtersqp.h
 * @brief   filterSQP NLP interface
 * @ingroup NLPIS
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_FILTERSQP_H__
#define __SCIP_NLPI_FILTERSQP_H__

#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for filterSQP solver and include it into SCIP, if filterSQP is available
 *
 * @ingroup NLPIIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlpSolverFilterSQP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLPIS
 *
 * @{
 */

/** gets string that identifies filterSQP */
SCIP_EXPORT
const char* SCIPgetSolverNameFilterSQP(void);

/** gets string that describes filterSQP */
SCIP_EXPORT
const char* SCIPgetSolverDescFilterSQP(void);

/** returns whether filterSQP is available, i.e., whether it has been linked in */
SCIP_EXPORT
SCIP_Bool SCIPisFilterSQPAvailableFilterSQP(void);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_FILTERSQP_H__ */
