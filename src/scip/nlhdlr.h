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

#endif /* SCIP_NLHDLR_H_ */
