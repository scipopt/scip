/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* This is a TEMPLATE for a NLPI. Use this as starting point to implement your own NLPI.
 * Copy the file, rename it, and replace all occurences of XYZ by the name of your NLP solver.
 */

/**@file    nlpi_xyz.h
 * @brief   XYZ NLP interface
 * @ingroup NLPIS
 * @author  you
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLPI_XYZ_H__
#define __SCIP_NLPI_XYZ_H__

#include "scip/type_nlpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create solver interface for Xyz solver and includes it into SCIP
 *
 * @ingroup NLPIIncludes
 */
SCIP_RETCODE SCIPincludeNlpSolverXyz(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLPIS
 *
 * @{
 */

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLPI_XYZ_H__ */
