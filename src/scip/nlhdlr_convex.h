/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_convex.h
 * @ingroup NLHDLRS
 * @brief  nonlinear handlers for convex and concave expressions, respectively
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_CONVEX_H__
#define __SCIP_NLHDLR_CONVEX_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes convex nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrConvex(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** includes concave nonlinear handler in nonlinear constraint handler
 *
 * @ingroup NlhdlrIncludes
 */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrConcave(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@addtogroup NLHDLRS
 * @{
 *
 * @name Convex and concave nonlinear handlers
 *
 * These nonlinear handler detect convex and concave subexpressions and provide specialized estimation functionality.
 *
 * @{
 */

/** checks whether a given expression is convex or concave w.r.t. the original variables
 *
 * This function uses the methods that are used in the detection algorithm of the convex nonlinear handler.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPhasExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPRCURV         curv,               /**< curvature to check for */
   SCIP_Bool*            success,            /**< buffer to store whether expression has curvature curv (w.r.t. original variables) */
   SCIP_HASHMAP*         assumevarfixed      /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   );

/** @}
  * @}
  */

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_CONVEX_H__ */
