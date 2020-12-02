/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nlhdlr.h
 * @brief  private functions of nonlinear handlers of nonlinear constraints
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_NLHDLR_H_
#define SCIP_NLHDLR_H_

#include "scip/pub_nlhdlr.h"

/** creates a nonlinear handler */
SCIP_EXPORT
SCIP_RETCODE SCIPnlhdlrCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr,             /**< buffer to store pointer to created nonlinear handler */
   const char*           name,               /**< name of nonlinear handler (must not be NULL) */
   const char*           desc,               /**< description of nonlinear handler (can be NULL) */
   int                   detectpriority,     /**< detection priority of nonlinear handler */
   int                   enfopriority,       /**< enforcement priority of nonlinear handler */
   SCIP_DECL_NLHDLRDETECT((*detect)),        /**< structure detection callback of nonlinear handler */
   SCIP_DECL_NLHDLREVALAUX((*evalaux)),      /**< auxiliary evaluation callback of nonlinear handler */
   SCIP_NLHDLRDATA*      nlhdlrdata          /**< data of nonlinear handler (can be NULL) */
   );

/** frees a nonlinear handler */
SCIP_EXPORT
SCIP_RETCODE SCIPnlhdlrFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlr              /**< pointer to nonlinear handler to be freed */
   );

/** call the handler copy callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRCOPYHDLR(SCIPnlhdlrCopyhdlr);

/** call the free expression specific data callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRFREEEXPRDATA(SCIPnlhdlrFreeexprdata);

/** call the initialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRINIT(SCIPnlhdlrInit);

/** call the deinitialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLREXIT(SCIPnlhdlrExit);

/** call the detect callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRDETECT(SCIPnlhdlrDetect);

/** call the auxiliary evaluation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLREVALAUX(SCIPnlhdlrEvalaux);

/** calls the interval evaluation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRINTEVAL(SCIPnlhdlrInteval);

/** calls the reverse propagation callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRREVERSEPROP(SCIPnlhdlrReverseprop);

/** calls the separation initialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRINITSEPA(SCIPnlhdlrInitsepa);

/** calls the separation deinitialization callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLREXITSEPA(SCIPnlhdlrExitsepa);

/** calls the enforcement callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRENFO(SCIPnlhdlrEnfo);

/** calls the estimator callback of a nonlinear handler */
SCIP_EXPORT
SCIP_DECL_NLHDLRESTIMATE(SCIPnlhdlrEstimate);

/** reset number of detections counter for last round */
SCIP_EXPORT
void SCIPnlhdlrResetNDetectionslast(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
   );

/** print statistics for nonlinear handlers */
SCIP_EXPORT
void SCIPnlhdlrPrintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLR**         nlhdlrs,            /**< nonlinear handlers */
   int                   nnlhdlrs,           /**< number of nonlinear handlers */
   FILE*                 file                /**< file handle, or NULL for standard out */
   );

#endif /* SCIP_NLHDLR_H_ */
