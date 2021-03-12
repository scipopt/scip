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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr_soc.h
 * @brief  soc nonlinear handler
 *
 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NLHDLR_SOC_H__
#define __SCIP_NLHDLR_SOC_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes SOC nonlinear handler in nonlinear constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeNlhdlrSoc(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** checks whether a given constraint is SOC-representable
 *
 * This function uses the methods that are used in the detection algorithm of the SOC nonlinear handler.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPisSOCNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool             compeigenvalues,    /**< whether eigenvalues should be computed to detect complex cases */
   SCIP_Bool*            success,            /**< pointer to store whether SOC structure has been detected */
   SCIP_SIDETYPE*        sidetype,           /**< pointer to store which side of cons is SOC representable; only
                                               valid when success is TRUE */
   SCIP_VAR***           vars,               /**< variables that appear on both sides (x) */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int*                  nvars,              /**< total number of variables appearing */
   int*                  nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   );

/** frees arrays created by SCIPisSOCNonlinear() */
SCIP_EXPORT
void SCIPfreeSOCArraysNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variables that appear on both sides (x) */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   );

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_NLHDLR_SOC_H__ */
