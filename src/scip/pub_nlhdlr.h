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

/**@file   pub_nlhdlr.h
 * @brief  public functions of nonlinear handlers of nonlinear constraints
 * @ingroup PUBLICCOREAPI
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#ifndef SCIP_PUB_NLHDLR_H_
#define SCIP_PUB_NLHDLR_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_nlhdlr.h"

/**@addtogroup PublicNlhdlrInterfaceMethods
 * @{
 */

/** sets the copy handler callback of a nonlinear handler */
SCIP_EXPORT
void SCIPnlhdlrSetCopyHdlr(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRCOPYHDLR((*copy))         /**< copy callback (can be NULL) */
);

/** sets the nonlinear handler callback to free the nonlinear handler data */
SCIP_EXPORT
void SCIPnlhdlrSetFreeHdlrData(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEHDLRDATA((*freehdlrdata)) /**< handler free callback (can be NULL) */
);

/** sets the nonlinear handler callback to free expression specific data of nonlinear handler */
SCIP_EXPORT
void SCIPnlhdlrSetFreeExprData(
   SCIP_NLHDLR*          nlhdlr,            /**< nonlinear handler */
   SCIP_DECL_NLHDLRFREEEXPRDATA((*freeexprdata)) /**< nonlinear handler expression data free callback (can be NULL if data does not need to be freed) */
);

/** sets the initialization and deinitialization callback of a nonlinear handler */
SCIP_EXPORT
void SCIPnlhdlrSetInitExit(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINIT((*init)),            /**< initialization callback (can be NULL) */
   SCIP_DECL_NLHDLREXIT((*exit))             /**< deinitialization callback (can be NULL) */
);

/** sets the propagation callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPnlhdlrSetProp(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINTEVAL((*inteval)),      /**< interval evaluation callback (can be NULL) */
   SCIP_DECL_NLHDLRREVERSEPROP((*reverseprop)) /**< reverse propagation callback (can be NULL) */
);

/** sets the enforcement callbacks of a nonlinear handler */
SCIP_EXPORT
void SCIPnlhdlrSetSepa(
   SCIP_NLHDLR*          nlhdlr,             /**< nonlinear handler */
   SCIP_DECL_NLHDLRINITSEPA((*initsepa)),    /**< separation initialization callback (can be NULL) */
   SCIP_DECL_NLHDLRENFO((*enfo)),            /**< enforcement callback (can be NULL if estimate is not NULL) */
   SCIP_DECL_NLHDLRESTIMATE((*estimate)),    /**< estimation callback (can be NULL if sepa is not NULL) */
   SCIP_DECL_NLHDLREXITSEPA((*exitsepa))     /**< separation deinitialization callback (can be NULL) */
);

/** gives name of nonlinear handler */
SCIP_EXPORT
const char* SCIPnlhdlrGetName(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** gives description of nonlinear handler, can be NULL */
SCIP_EXPORT
const char* SCIPnlhdlrGetDesc(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** gives detection priority of nonlinear handler */
SCIP_EXPORT
int SCIPnlhdlrGetDetectPriority(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** gives enforcement priority of nonlinear handler */
SCIP_EXPORT
int SCIPnlhdlrGetEnfoPriority(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler is enabled */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrIsEnabled(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** gives handler data of nonlinear handler */
SCIP_EXPORT
SCIP_NLHDLRDATA* SCIPnlhdlrGetData(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the interval evaluation callback */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrHasIntEval(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the reverse propagation callback */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrHasReverseProp(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the separation initialization callback */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrHasInitSepa(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the separation deinitialization callback */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrHasExitSepa(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the enforcement callback */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrHasEnfo(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** returns whether nonlinear handler implements the estimator callback */
SCIP_EXPORT
SCIP_Bool SCIPnlhdlrHasEstimate(
   SCIP_NLHDLR*          nlhdlr              /**< nonlinear handler */
);

/** compares two nonlinear handlers by detection priority
 *
 * if handlers have same detection priority, then compare by name
 */
SCIP_DECL_SORTPTRCOMP(SCIPnlhdlrComp);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* SCIP_PUB_NLHDLR_H_ */
